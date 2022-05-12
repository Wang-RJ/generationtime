# Analyze and make plots for each population separately

AFR_all <- read.table("AFR_allcounts.txt")
names(AFR_all) <- count_header
AFR_all <- cbind(AFR_all[,mut_classes] + AFR_all[,complementary_classes],
                 CpG = AFR_all$CpG, bin_age = AFR_all$bin_age)

load("bootstraps/unphased_allboot_estsAFR.RData")
load("bootstraps/phased_allboot_estsAFR.RData")
AFR_phasedplotframes <- make_plotframes(phasedstrictCT_model, AFR_all, AFRphased_allboot_ests)
AFR_phasedallplots <- make_tgpplot(AFR_phasedplotframes[[1]], AFR_phasedplotframes[[2]])
AFR_phasedallplots[[1]]
AFR_phasedallplots[[2]] + ggtitle("AFR phased")

AFR_unphasedplotframes <- make_plotframes(strictCT_model, AFR_all, AFRunphased_allboot_ests)
AFR_unphasedallplots <- make_tgpplot(AFR_unphasedplotframes[[1]], AFR_unphasedplotframes[[2]])
AFR_unphasedallplots[[1]] + ggtitle("AFR")
AFR_unphasedallplots[[2]] + ggtitle("AFR unphased")


EAS_all <- read.table("tgppops/EAS_allcounts.txt")
names(EAS_all) <- count_header
EAS_all <- cbind(EAS_all[,mut_classes] + EAS_all[,complementary_classes],
                 CpG = EAS_all$CpG, bin_age = EAS_all$bin_age)

load("bootstraps/unphased_allboot_estsEAS.RData")
load("bootstraps/phased_allboot_estsEAS.RData")
EAS_phasedplotframes <- make_plotframes(phasedstrictCT_model, EAS_all, EASphased_allboot_ests)
EAS_phasedallplots <- make_tgpplot(EAS_phasedplotframes[[1]], EAS_phasedplotframes[[2]])
EAS_phasedallplots[[1]]
EAS_phasedallplots[[2]] + ggtitle("EAS phased")

EAS_unphasedplotframes <- make_plotframes(strictCT_model, EAS_all, EASunphased_allboot_ests)
EAS_unphasedallplots <- make_tgpplot(EAS_phasedplotframes[[1]], EAS_phasedplotframes[[2]])
EAS_unphasedallplots[[1]]
EAS_unphasedallplots[[2]] + ggtitle("EAS unphased")


EUR_all <- read.table("tgppops/EUR_allcounts.txt")
names(EUR_all) <- count_header
EUR_all <- cbind(EUR_all[,mut_classes] + EUR_all[,complementary_classes],
                 CpG = EUR_all$CpG, bin_age = EUR_all$bin_age)

load("bootstraps/unphased_allboot_estsEUR.RData")
load("bootstraps/phased_allboot_estsEUR.RData")
EUR_phasedplotframes <- make_plotframes(phasedstrictCT_model, EUR_all, EURphased_allboot_ests)
EUR_phasedallplots <- make_tgpplot(EUR_phasedplotframes[[1]], EUR_phasedplotframes[[2]])
EUR_phasedallplots[[1]]
EUR_phasedallplots[[2]] + ggtitle("EUR phased")

EUR_unphasedplotframes <- make_plotframes(strictCT_model, EUR_all, EURunphased_allboot_ests)
EUR_unphasedallplots <- make_tgpplot(EUR_unphasedplotframes[[1]], EUR_unphasedplotframes[[2]])
EUR_unphasedallplots[[1]]
EUR_unphasedallplots[[2]] + ggtitle("EUR unphased")

SAS_all <- read.table("SAS_allcounts.txt")
names(SAS_all) <- count_header
SAS_all <- cbind(SAS_all[,mut_classes] + SAS_all[,complementary_classes],
                 CpG = SAS_all$CpG, bin_age = SAS_all$bin_age)

load("bootstraps/unphased_allboot_estsSAS.RData")
load("bootstraps/phased_allboot_estsSAS.RData")
SAS_phasedplotframes <- make_plotframes(phasedstrictCT_model, SAS_all, SASphased_allboot_ests)
SAS_phasedallplots <- make_tgpplot(SAS_phasedplotframes[[1]], SAS_phasedplotframes[[2]])
SAS_phasedallplots[[1]]
SAS_phasedallplots[[2]] + ggtitle("SAS phased")

SAS_unphasedplotframes <- make_plotframes(strictCT_model, SAS_all, SASunphased_allboot_ests)
SAS_unphasedallplots <- make_tgpplot(SAS_unphasedplotframes[[1]], SAS_unphasedplotframes[[2]])
SAS_unphasedallplots[[1]]
SAS_unphasedallplots[[2]] + ggtitle("SAS unphased")

# Repeat for only 10k generations

popsframe_all <-
  data.frame(mean_interval = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "phasedplotframes", sep = "_"))[[2]]$sexave })),
             bin_age = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "phasedplotframes", sep = "_"))[[2]]$bin_age })),
             sexdiff = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "phasedplotframes", sep = "_"))[[2]]$sexdiff })),
             sexave_upper = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "phasedplotframes", sep = "_"))[[2]]$sexave_upper })),
             sexave_lower = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "phasedplotframes", sep = "_"))[[2]]$sexave_lower })),
             pop = rep(pops, each = 100))

ggplot(popsframe_all, aes(x = bin_age, y = mean_interval, col = pop)) +
  geom_smooth(method = "loess", span = 0.65, se = FALSE, lwd = 1.1) +
  geom_ribbon(aes(ymin = sexave_lower, ymax = sexave_upper, fill = pop), color = NA, alpha = 0.1) +
  axis_formatting + coord_cartesian(ylim = c(13.5, 38.5), xlim = c(100,80000)) + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "generation interval (years)", x = "allele age (generations ago)")

AFR_phased10kgaplots <- make_tgpplot(AFR_phasedplotframes10kga[[1]], AFR_phasedplotframes10kga[[2]])
AFR_phased10kgaplots[[1]]

EAS_phased10kgaplots <- make_tgpplot(EAS_phasedplotframes10kga[[1]], EAS_phasedplotframes10kga[[2]])
EAS_phased10kgaplots[[1]]

EUR_phased10kgaplots <- make_tgpplot(EUR_phasedplotframes10kga[[1]], EUR_phasedplotframes10kga[[2]])
EUR_phased10kgaplots[[1]]

SAS_phased10kgaplots <- make_tgpplot(SAS_phasedplotframes10kga[[1]], SAS_phasedplotframes10kga[[2]])
SAS_phased10kgaplots[[1]]

# Plots for 1k generation and 1k generation private only

for(pop in pops) {
  pop_private <- read.table(paste(pop, "_1kgaprivatecounts.txt", sep = ""))
  pop_shared <- read.table(paste(pop, "_1kgacounts.txt", sep = ""))
  names(pop_private) <- count_header
  names(pop_shared) <- count_header
  pop_private <- cbind(pop_private[,mut_classes] + pop_private[,complementary_classes],
                       CpG = pop_private$CpG, bin_age = pop_private$bin_age)
  pop_shared <- cbind(pop_shared[,mut_classes] + pop_shared[,complementary_classes],
                      CpG = pop_shared$CpG, bin_age = pop_shared$bin_age)
  assign(paste(pop, "1kga", sep = "_"), pop_shared)
  assign(paste(pop, "private", sep = "_"), pop_private)
}

for(pop in pops) {
  load(paste("bootstraps/phased_10kgaboot_ests", pop, ".RData", sep = ""))
  load(paste("bootstraps/phased_1kgaboot_ests", pop, ".RData", sep = ""))
  assign(paste(pop, "phased_1kgaboot_ests", sep = ""), pop_1kga_est)
  load(paste("bootstraps/phased_priv1kgaboot_ests", pop, ".RData", sep = ""))
}

for(pop in pops) {
  assign(paste(pop, "phasedplotframes10kga", sep = "_"),
         make_plotframes(phasedstrictCT_model, get(paste(pop, "10kga", sep = "_")),
                         get(paste(pop, "phased_10kgaboot_ests", sep = ""))))
  assign(paste(pop, "phasedplotframespriv1kga", sep = "_"),
         make_plotframes(phasedstrictCT_model, get(paste(pop, "private", sep = "_")),
                         get(paste(pop, "phased_priv1kgaboot_ests", sep = ""))))
  assign(paste(pop, "phasedplotframes1kga", sep = "_"),
         make_plotframes(phasedstrictCT_model, get(paste(pop, "1kga", sep = "_")),
                         get(paste(pop, "phased_1kgaboot_ests", sep = ""))))
  
}

popsframe_10kga <-
  data.frame(mean_interval = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "phasedplotframes10kga", sep = "_"))[[2]]$sexave })),
             bin_age = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "phasedplotframes10kga", sep = "_"))[[2]]$bin_age })),
             sexdiff = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "phasedplotframes10kga", sep = "_"))[[2]]$sexdiff })),
             sexave_upper = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "phasedplotframes10kga", sep = "_"))[[2]]$sexave_upper })),
             sexave_lower = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "phasedplotframes10kga", sep = "_"))[[2]]$sexave_lower })),
             sexdiff_upper = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "phasedplotframes10kga", sep = "_"))[[2]]$sexdiff_upper })),
             sexdiff_lower = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "phasedplotframes10kga", sep = "_"))[[2]]$sexdiff_lower })),
             pop = rep(pops, each = 100))

ggplot(popsframe_10kga, aes(x = bin_age, y = mean_interval, col = pop)) +
  geom_smooth(method = "loess", span = 0.65, se = FALSE, lwd = 1.1) +
  geom_ribbon(aes(ymin = sexave_lower, ymax = sexave_upper, fill = pop), color = NA, alpha = 0.1) +
  coord_cartesian(ylim = c(13.5,38.5)) + axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "generation interval (years)", x = "allele age (generations ago)")

ggplot(popsframe_10kga, aes(x = bin_age, y = sexdiff, col = pop)) +
  geom_smooth(method = "loess", span = 0.65, se = FALSE, lwd = 1.1) +
  geom_ribbon(aes(ymin = sexdiff_lower, ymax = sexdiff_upper, fill = pop), color = NA, alpha = 0.1) +
  axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 10, 1:9 * 100, 1:9 * 1000)) +
  labs(y = "male-female difference (years)", x = "allele age (generations ago)")

popsframe_priv1kga <-
  data.frame(mean_interval = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "phasedplotframespriv1kga", sep = "_"))[[2]]$sexave })),
             bin_age = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "phasedplotframespriv1kga", sep = "_"))[[2]]$bin_age })),
             sexdiff = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "phasedplotframespriv1kga", sep = "_"))[[2]]$sexdiff })),
             sexave_upper = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "phasedplotframespriv1kga", sep = "_"))[[2]]$sexave_upper })),
             sexave_lower = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "phasedplotframespriv1kga", sep = "_"))[[2]]$sexave_lower })),
             pop = rep(pops, each = 100))

ggplot(popsframe_priv1kga, aes(x = bin_age, y = mean_interval, col = pop)) +
  geom_smooth(method = "loess", span = .65, se = FALSE, lwd = 1.1) +
  geom_ribbon(aes(ymin = sexave_lower, ymax = sexave_upper, fill = pop), color = NA, alpha = 0.1) +
  axis_formatting + coord_cartesian(ylim = c(13.5,38.5), xlim = c(50,1e3))

popsframe_1kga <-
  data.frame(mean_interval = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "phasedplotframes1kga", sep = "_"))[[2]]$sexave })),
             bin_age = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "phasedplotframes1kga", sep = "_"))[[2]]$bin_age })),
             sexdiff = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "phasedplotframes1kga", sep = "_"))[[2]]$sexdiff })),
             sexave_upper = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "phasedplotframes1kga", sep = "_"))[[2]]$sexave_upper })),
             sexave_lower = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "phasedplotframes1kga", sep = "_"))[[2]]$sexave_lower })),
             pop = rep(pops, each = 100))

ggplot(popsframe_1kga, aes(x = bin_age, y = mean_interval, col = pop)) +
  geom_smooth(method = "loess", span = .65, se = FALSE, lwd = 1.1) + 
  geom_ribbon(aes(ymin = sexave_lower, ymax = sexave_upper, fill = pop), color = NA, alpha = 0.1) +
  axis_formatting + coord_cartesian(ylim = c(13.5,38.5)) +
  labs(y = "generation interval (years)", x = "allele age (generations ago)")

## Calculate averages

t.test(c(EAS_phasedplotframes10kga[[2]]$sexave,
         EUR_phasedplotframes10kga[[2]]$sexave,
         SAS_phasedplotframes10kga[[2]]$sexave),
       AFR_phasedplotframes10kga[[2]]$sexave)

sum(diff(c(0, tgp_phased10kgaplots[[2]]$bin_age)) * tgp_phased10kgaplots[[2]]$sexave) /
  sum(diff(c(0, tgp_phased10kgaplots[[2]]$bin_age)))
sum(diff(c(0, tgp_phased10kgaplots[[2]]$bin_age)) * tgp_phased10kgaplots[[2]]$sexave_sd) /
  sum(diff(c(0, tgp_phased10kgaplots[[2]]$bin_age)))

x <- subset(tgp_phased10kgaplots[[1]], sex == "M")
sum(diff(c(0, x$bin_age)) * x$gen_time) / sum(diff(c(0, x$bin_age)))
sum(diff(c(0, x$bin_age)) * x$gen_sd) / sum(diff(c(0, x$bin_age)))
x <- subset(tgp_phased10kgaplots[[1]], sex == "F")
sum(diff(c(0, x$bin_age)) * x$gen_time) / sum(diff(c(0, x$bin_age)))
sum(diff(c(0, x$bin_age)) * x$gen_sd) / sum(diff(c(0, x$bin_age)))

x <- AFR_phasedplotframes10kga[[2]]
sum(diff(c(0, x$bin_age)) * x$sexave) / sum(diff(c(0, x$bin_age)))
sum(diff(c(0, x$bin_age)) * x$sexave_sd) / sum(diff(c(0, x$bin_age)))
