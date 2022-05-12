# Read and pre-process data from Jonsson et al.
# Available in supplemental data from doi: 10.1038/nature24018

decode_data <- read.table("decode_DNMs.tsv", header = TRUE)
decode_data$mut_class <- paste(decode_data$Ref, decode_data$Alt, sep = "_")
decode_data$ref_triplet <- as.character(getSeq(Hsapiens, decode_data$Chr,
                                               decode_data$Pos_hg38 - 1, decode_data$Pos_hg38 + 1))

# Remove indels
decode_data <- decode_data[!(sapply(decode_data$mut_class, nchar) > 3),]
phased_decode <- subset(decode_data, !is.na(Phase_combined))

# define mutation classes and collapse complementary classes in decode data
mut_classes <- c("A_C", "A_G", "A_T", "C_A", "C_G", "C_T")
complementary_classes <- c("T_G", "T_C", "T_A", "G_T", "G_C", "G_A")
class_map <- setNames(mut_classes, complementary_classes)

collapse_class <- function(mut_class) {
  if(mut_class %in% mut_classes)
    return(mut_class)
  return(class_map[mut_class])
}

decode_data$mut_class <- sapply(decode_data$mut_class, collapse_class)

decode_noCpG <- subset(decode_data, (substr(ref_triplet, 1, 2) != "CG") &
                         (substr(ref_triplet, 2, 3) != "CG"))
decode_noCpG_noEuroTCtriplets <- subset(decode_noCpG,
                                        !(mut_class == "C_T" & (ref_triplet == "TCC" | ref_triplet == "GGA" |
                                                                  ref_triplet == "ACC" | ref_triplet == "GGT" |
                                                                  ref_triplet == "TCT" | ref_triplet == "AGA" |
                                                                  ref_triplet == "CCC" | ref_triplet == "GGG")))

phased_noCpG_noEuroTCtriplets <- subset(decode_noCpG_noEuroTCtriplets, !is.na(Phase_combined))

human_muttable <- table(decode_data$mut_class)
human_muttable_noCpG <- table(decode_noCpG$mut_class)
human_muttable_noCpG_noEuroTCtriplets <- table(decode_noCpG_noEuroTCtriplets$mut_class)

human_spectrum <- human_muttable / sum(human_muttable)
human_spectrum_noCpG <- human_muttable_noCpG / sum(human_muttable_noCpG)
human_spectrum_noCpG_noEuroTCtriplets <- human_muttable_noCpG_noEuroTCtriplets /
  sum(human_muttable_noCpG_noEuroTCtriplets)


###
# aggregate spectra by proband

aggregate_spectra <- function(mutations) {
  mutations_byproband <- aggregate(mutations$mut_class,
                                   by = list(Proband = mutations$Proband_nr), length)  
  ages_byproband <- aggregate(mutations[,c("Fathers_age_at_conception",
                                           "Mothers_age_at_conception")],
                              by = list(Proband = mutations$Proband_nr), unique)
  mutations_byproband[,c("Fathers_age_at_conception", "Mothers_age_at_conception")] <-
    ages_byproband[,c(2,3)]
  
  for(class_i in mut_classes) {
    class_byproband <- aggregate(mutations$mut_class,
                                 by = list(Proband = mutations$Proband_nr),
                                 FUN = function(classes) { return(sum(classes == class_i)) })
    names(class_byproband)[2] <- class_i
    mutations_byproband <- cbind(mutations_byproband, class_byproband[2])
  }
  
  return(mutations_byproband)
}

agg_all <- aggregate_spectra(decode_data)
agg_noCpG_noEuroTCtriplets <- aggregate_spectra(decode_noCpG_noEuroTCtriplets)
agg_phased <- aggregate_spectra(phased_noCpG_noEuroTCtriplets)

# agg_unphased_filtered <- agg_noCpG_noEuroTCtriplets[
#   rowSums(agg_noCpG_noEuroTCtriplets[,mut_classes]) > 20,]
agg_unphased_filtered <- agg_noCpG_noEuroTCtriplets

strictCT_model <- MGLMreg(formula = cbind(A_C, A_G, A_T, C_A, C_G, C_T) ~
                            Fathers_age_at_conception + Mothers_age_at_conception,
                          data = agg_unphased_filtered, dist = "DM")

# filter probands with very few phased mutations to avoid creating indefinite Hessian
agg_phased_filtered <- agg_phased[rowSums(agg_phased[,mut_classes]) > 10,]
phasedstrictCT_model <- MGLMreg(formula = cbind(A_C, A_G, A_T, C_A, C_G, C_T) ~
                                  Fathers_age_at_conception + Mothers_age_at_conception,
                                data = agg_phased_filtered, dist = "DM")

predict_spectrum <- function(model, paternal_age, maternal_age) {
  spectrum <- as.vector(predict(model, c(1, paternal_age, maternal_age)))
  return(setNames(spectrum, mut_classes))
}

predict_parentalages <- function(model, spectrum_vector) {
  return(optim(c(30,30), function(ages) {
    male_age <- ages[1]
    female_age <- ages[2]
    return(aitchison.distance(predict_spectrum(model, male_age, female_age), spectrum_vector))
  }, lower = 0, upper = 100, method = "L-BFGS-B"))
}

aitchison.distance <- function(v1, v2) {
  return(norm(clr(v1) - clr(v2)))
}

LL_DM <- function(obs, coefficients, paternal_age, maternal_age) {
  x <- c(1, paternal_age, maternal_age)
  alpha <- exp(t(x) %*% coefficients)
  
  return(ddirmn(obs, alpha))
}