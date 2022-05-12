# Distribute bootstraps using snow package

library(snow)

clust <- makeCluster(15)
clusterCall(clust, fun = function() { library(MGLM); library(compositions) })
clusterExport(clust, list("predict_parentalages", "TGP_bootsall", "youngspecdiff", "aitchison.distance",
                          "predict_spectrum", "mut_classes"))
phased_boot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(TGP_bootsall[,1:6] / rowSums(TGP_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(phased_boot_ests, file = "phased_boot_estsTGP.RData")
unphased_boot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(TGP_bootsall[,1:6] / rowSums(TGP_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(unphased_boot_ests, file = "unphased_boot_estsTGP.RData")
stopCluster(clust)

clust <- makeCluster(15)
clusterCall(clust, fun = function() { library(MGLM); library(compositions) })
clusterExport(clust, list("predict_parentalages", "TGP_bootsall", "AFR_bootsall", "EAS_bootsall", "EUR_bootsall", "SAS_bootsall",
                          "youngspecdiff", "aitchison.distance", "predict_spectrum", "mut_classes"))
AFRphased_allboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(AFR_bootsall[,1:6] / rowSums(AFR_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(AFRphased_allboot_ests, file = "phased_allboot_estsAFR.RData")
AFRunphased_allboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(AFR_bootsall[,1:6] / rowSums(AFR_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(AFRunphased_allboot_ests, file = "unphased_allboot_estsAFR.RData")

EASphased_allboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EAS_bootsall[,1:6] / rowSums(EAS_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EASphased_allboot_ests, file = "phased_allboot_estsEAS.RData")
EASunphased_allboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EAS_bootsall[,1:6] / rowSums(EAS_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EASunphased_allboot_ests, file = "unphased_allboot_estsEAS.RData")

EURphased_allboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EUR_bootsall[,1:6] / rowSums(EUR_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EURphased_allboot_ests, file = "phased_allboot_estsEUR.RData")
EURunphased_allboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EUR_bootsall[,1:6] / rowSums(EUR_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EURunphased_allboot_ests, file = "unphased_allboot_estsEUR.RData")

SASphased_allboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(SAS_bootsall[,1:6] / rowSums(SAS_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(SASphased_allboot_ests, file = "phased_allboot_estsSAS.RData")
SASunphased_allboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(SAS_bootsall[,1:6] / rowSums(SAS_bootsall[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(SASunphased_allboot_ests, file = "unphased_allboot_estsSAS.RData")
stopCluster(clust)


clust <- makeCluster(15)
clusterCall(clust, fun = function() { library(MGLM); library(compositions) })
clusterExport(clust, list("predict_parentalages", "AFR_boots_10kga", "EAS_boots_10kga", "EUR_boots_10kga", "SAS_boots_10kga",
                          "youngspecdiff", "aitchison.distance", "predict_spectrum", "mut_classes"))
AFRphased_10kgaboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(AFR_boots_10kga[,1:6] / rowSums(AFR_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(AFRphased_10kgaboot_ests, file = "phased_10kgaboot_estsAFR.RData")
AFRunphased_10kgaboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(AFR_boots_10kga[,1:6] / rowSums(AFR_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(AFRunphased_10kgaboot_ests, file = "unphased_10kgaboot_estsAFR.RData")

EASphased_10kgaboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EAS_boots_10kga[,1:6] / rowSums(EAS_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EASphased_10kgaboot_ests, file = "phased_10kgaboot_estsEAS.RData")
EASunphased_10kgaboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EAS_boots_10kga[,1:6] / rowSums(EAS_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EASunphased_10kgaboot_ests, file = "unphased_10kgaboot_estsEAS.RData")

EURphased_10kgaboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EUR_boots_10kga[,1:6] / rowSums(EUR_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EURphased_10kgaboot_ests, file = "phased_10kgaboot_estsEUR.RData")
EURunphased_10kgaboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EUR_boots_10kga[,1:6] / rowSums(EUR_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EURunphased_10kgaboot_ests, file = "unphased_10kgaboot_estsEUR.RData")

SASphased_10kgaboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(SAS_boots_10kga[,1:6] / rowSums(SAS_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(SASphased_10kgaboot_ests, file = "phased_10kgaboot_estsSAS.RData")
SASunphased_10kgaboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(SAS_boots_10kga[,1:6] / rowSums(SAS_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(SASunphased_10kgaboot_ests, file = "unphased_10kgaboot_estsSAS.RData")
stopCluster(clust)


clust <- makeCluster(15)
clusterCall(clust, fun = function() { library(MGLM); library(compositions) })
clusterExport(clust, list("predict_parentalages", "AFR_boots_priv1kga", "EAS_boots_priv1kga", "EUR_boots_priv1kga", "SAS_boots_priv1kga",
                          "youngspecdiff", "aitchison.distance", "predict_spectrum", "mut_classes"))
AFRphased_priv1kgaboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(AFR_boots_priv1kga[,1:6] / rowSums(AFR_boots_priv1kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(AFRphased_priv1kgaboot_ests, file = "phased_priv1kgaboot_estsAFR.RData")
AFRunphased_priv1kgaboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(AFR_boots_priv1kga[,1:6] / rowSums(AFR_boots_priv1kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(AFRunphased_priv1kgaboot_ests, file = "unphased_priv1kgaboot_estsAFR.RData")

EASphased_priv1kgaboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EAS_boots_priv1kga[,1:6] / rowSums(EAS_boots_priv1kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EASphased_priv1kgaboot_ests, file = "phased_priv1kgaboot_estsEAS.RData")
EASunphased_priv1kgaboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EAS_boots_priv1kga[,1:6] / rowSums(EAS_boots_priv1kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EASunphased_priv1kgaboot_ests, file = "unphased_priv1kgaboot_estsEAS.RData")

EURphased_priv1kgaboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EUR_boots_priv1kga[,1:6] / rowSums(EUR_boots_priv1kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EURphased_priv1kgaboot_ests, file = "phased_priv1kgaboot_estsEUR.RData")
EURunphased_priv1kgaboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(EUR_boots_priv1kga[,1:6] / rowSums(EUR_boots_priv1kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(EURunphased_priv1kgaboot_ests, file = "unphased_priv1kgaboot_estsEUR.RData")

SASphased_priv1kgaboot_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(SAS_boots_priv1kga[,1:6] / rowSums(SAS_boots_priv1kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(SASphased_priv1kgaboot_ests, file = "phased_priv1kgaboot_estsSAS.RData")
SASunphased_priv1kgaboot_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(SAS_boots_priv1kga[,1:6] / rowSums(SAS_boots_priv1kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(SASunphased_priv1kgaboot_ests, file = "unphased_priv1kgaboot_estsSAS.RData")
stopCluster(clust)


TGP_boots_10kga <- read.table("Mutations/spectral/bootstraps/TGP_boots_10kga.txt")
names(TGP_boots_10kga) <- c("window", count_header)
TGP_boots_10kga$window <- as.numeric(factor(TGP_boots_10kga$window))
TGP_boots_10kga <- cbind(TGP_boots_10kga[,mut_classes] + TGP_boots_10kga[,complementary_classes],
                         CpG = TGP_boots_10kga$CpG, bin_age = TGP_boots_10kga$bin_age)

clust <- makeCluster(15)
clusterCall(clust, fun = function() { library(MGLM); library(compositions) })
clusterExport(clust, list("predict_parentalages", "TGP_boots_10kga", "youngspecdiff", "aitchison.distance",
                          "predict_spectrum", "mut_classes"))
TGPphased_boot10kga_ests <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(TGP_boots_10kga[,1:6] / rowSums(TGP_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(TGPphased_boot10kga_ests, file = "phased_10kgaboot_estsTGP.RData")
TGPunphased_boot10kga_ests <- clusterApplyLB(clust, unphased_boot_models, fun = function(boot_model) {
  oout <- apply(t(t(TGP_boots_10kga[,1:6] / rowSums(TGP_boots_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = boot_model)
  return(sapply(oout, "[[", 1))
})
save(TGPunphased_boot10kga_ests, file = "unphased_10kgaboot_estsTGP.RData")
stopCluster(clust)

for(pop in pops) {
  pop_boot <- read.table(paste(pop, "_boots_1kga.txt", sep = ""))
  
  names(pop_boot) <- c("window", count_header)
  pop_boot$window <- as.numeric(factor(pop_boot$window))
  pop_boot <- cbind(pop_boot[,mut_classes] + pop_boot[,complementary_classes],
                    CpG = pop_boot$CpG, bin_age = pop_boot$bin_age)
  assign(paste(pop, "boots1kga", sep = "_"), pop_boot)
}

clust <- makeCluster(15)
clusterCall(clust, fun = function() { library(MGLM); library(compositions) })
clusterExport(clust, c("predict_parentalages", "youngspecdiff", "aitchison.distance",
                       "predict_spectrum", "mut_classes"))
for(pop in pops) {
  pop_boots <- get(paste(pop, "boots1kga", sep = "_"))
  clusterExport(clust, "pop_boots")
  pop_1kga_est <- clusterApplyLB(clust, phased_boot_models, fun = function(boot_model) {
    oout <- apply(t(t(pop_boots[,1:6] / rowSums(pop_boots[,1:6])) - youngspecdiff), MARGIN = 1,
                  predict_parentalages, model = boot_model)
    return(sapply(oout, "[[", 1))
  })
  save(pop_1kga_est, file = paste("phased_1kgaboot_ests", pop, ".RData", sep = ""))
}
stopCluster(clust)
