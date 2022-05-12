# Analysis of generation interval with Neanderthal variants masked

tgp_nea <- read.table("neanderthal_masked/TGP_10kgacounts.txt")

names(tgp_nea) <- count_header
tgp_nea <- cbind(tgp_nea[,mut_classes] + tgp_nea[,complementary_classes],
                 CpG = tgp_nea$CpG, bin_age = tgp_nea$bin_age)

tgp_nea10kgaplots <- make_plotframes(phasedstrictCT_model, tgp_nea, TGPphased_boot10kga_ests)
plotframe_noNea <- tgp_nea10kgaplots[[1]]
plotmean_noNea <- tgp_nea10kgaplots[[2]]

tgp_10kga <- read.table("TGP_10kgacounts.txt")
names(tgp_10kga) <- count_header
tgp_10kga <- cbind(tgp_10kga[,mut_classes] + tgp_10kga[,complementary_classes],
                   CpG = tgp_10kga$CpG, bin_age = tgp_10kga$bin_age)

tgp_phased10kgaplots <- make_plotframes(phasedstrictCT_model, tgp_10kga, TGPphased_boot10kga_ests)
plotframe <- tgp_phased10kgaplots[[1]]
plotmean <- tgp_phased10kgaplots[[2]]

ggplot(plotframe_noNea, aes(y = gen_time, x = bin_age, col = sex)) +
  geom_point(size = 2) +
  geom_smooth(data = plotmean_noNea, aes(y = sexave, x = bin_age, col = 'mean'), lwd = 1.2,
              formula = y ~ x, method = "loess", span = 0.65, se = FALSE, color = 'darkgray') +
  geom_smooth(data = plotmean, aes(y = sexave, x = bin_age, col = 'mean'), lwd = 1.2, lty = 2,
              formula = y ~ x, method = "loess", span = 0.65, se = FALSE, color = 'black') +
  geom_ribbon(data = plotmean, aes(x = bin_age, y = sexave, col = 'mean',
                                   ymin = sexave_lower, ymax = sexave_upper),
              alpha = 0.15, color = NA) +
  coord_cartesian(ylim = c(13.5,38.5)) + scale_color_manual(values = rev(pal_aaas("default")(2))) +
  axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "generation interval (years)", x = "allele age (generations ago)", col = "sex")

## 

for(pop in pops) {
  load(paste("bootstraps/phased_10kgaboot_ests", pop, ".RData", sep = ""))
}

tgp_nea <- read.table("neanderthal_masked/TGP_10kgacounts.txt")

names(tgp_nea) <- count_header
tgp_nea <- cbind(tgp_nea[,mut_classes] + tgp_nea[,complementary_classes],
                 CpG = tgp_nea$CpG, bin_age = tgp_nea$bin_age)


for(pop in pops) {
  pop_nea <- read.table(paste("neanderthal_masked/", pop, "_10kgacounts.txt", sep = ""))
  names(pop_nea) <- count_header
  pop_nea <- cbind(pop_nea[,mut_classes] + pop_nea[,complementary_classes],
                   CpG = pop_nea$CpG, bin_age = pop_nea$bin_age)
  
  assign(paste(pop, "noNea", sep = "_"), pop_nea)
  assign(paste(pop, "noNeaphasedplotframes10kga", sep = "_"),
         make_plotframes(phasedstrictCT_model, get(paste(pop, "noNea", sep = "_")),
                         get(paste(pop, "phased_10kgaboot_ests", sep = ""))))
}

popsframe_10kga_noNea <-
  data.frame(mean_interval = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "noNeaphasedplotframes10kga", sep = "_"))[[2]]$sexave })),
             bin_age = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "noNeaphasedplotframes10kga", sep = "_"))[[2]]$bin_age })),
             sexdiff = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "noNeaphasedplotframes10kga", sep = "_"))[[2]]$sexdiff })),
             sexave_upper = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "noNeaphasedplotframes10kga", sep = "_"))[[2]]$sexave_upper })),
             sexave_lower = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "noNeaphasedplotframes10kga", sep = "_"))[[2]]$sexave_lower })),
             sexdiff_upper = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "noNeaphasedplotframes10kga", sep = "_"))[[2]]$sexdiff_upper })),
             sexdiff_lower = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "noNeaphasedplotframes10kga", sep = "_"))[[2]]$sexdiff_lower })),
             pop = rep(pops, each = 100))

loEAS10_noNea <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_10kga_noNea, pop == "EAS"), span = 0.65)
loEUR10_noNea <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_10kga_noNea, pop == "EUR"), span = 0.65)
loSAS10_noNea <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_10kga_noNea, pop == "SAS"), span = 0.65)

predx_noNea <- 10**seq(max(loEAS10_noNea$x[1], loEUR10_noNea$x[1], loSAS10_noNea$x[1]),
                       min(loEAS10_noNea$x[100], loEUR10_noNea$x[100], loSAS10_noNea$x[100]), .01)

predy_noNea <- (predict(loEAS10_noNea, data.frame(bin_age = predx)) + 
                  predict(loEUR10_noNea, data.frame(bin_age = predx)) +
                  predict(loSAS10_noNea, data.frame(bin_age = predx))) / 3

ggplot(subset(popsframe_10kga_noNea, pop == "AFR"), aes(x = bin_age, y = mean_interval, col = pop)) +
  geom_smooth(method = "loess", span = 0.65, se = FALSE, lwd = 1.1) +
  geom_line(data = data.frame(mean_interval = predy_noNea, bin_age = predx_noNea), col = 'black', lwd = 1.1) + 
  coord_cartesian(ylim = c(13.5,38.5)) + axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "generation interval (years)", x = "allele age (generations ago)")

for(pop in pops) {
  pop_10kga <- read.table(paste(pop, "_10kgacounts.txt", sep = ""))
  names(pop_10kga) <- count_header
  pop_10kga <- cbind(pop_10kga[,mut_classes] + pop_10kga[,complementary_classes], CpG = pop_10kga$CpG, bin_age = pop_10kga$bin_age)
  
  oout <- apply(t(t(pop_10kga[,1:6] / rowSums(pop_10kga[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = strictCT_model)
  pop_10kga$mgen_time <- sapply(oout, "[[", 1)[1,]
  pop_10kga$fgen_time <- sapply(oout, "[[", 1)[2,]
  pop_10kga$ss_mean <- (pop_10kga$mgen_time + pop_10kga$fgen_time) / 2
  
  assign(paste(pop, "10kga", sep = "_"), pop_10kga)
  assign(paste(pop, "phasedplotframes10kga", sep = "_"),
         make_plotframes(phasedstrictCT_model, get(paste(pop, "10kga", sep = "_")),
                         get(paste(pop, "phased_10kgaboot_ests", sep = ""))))
  
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

loEAS10 <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_10kga, pop == "EAS"), span = 0.65)
loEUR10 <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_10kga, pop == "EUR"), span = 0.65)
loSAS10 <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_10kga, pop == "SAS"), span = 0.65)

predx <- 10**seq(max(loEAS10$x[1], loEUR10$x[1], loSAS10$x[1]),
                 min(loEAS10$x[100], loEUR10$x[100], loSAS10$x[100]), .01)

predy <- (predict(loEAS10, data.frame(bin_age = predx)) + 
            predict(loEUR10, data.frame(bin_age = predx)) +
            predict(loSAS10, data.frame(bin_age = predx))) / 3

loEAS10upperCI <- loess(sexave_upper ~ log10(bin_age), data = subset(popsframe_10kga, pop == "EAS"), span = 0.65)
loEAS10lowerCI <- loess(sexave_lower ~ log10(bin_age), data = subset(popsframe_10kga, pop == "EAS"), span = 0.65)

loEUR10upperCI <- loess(sexave_upper ~ log10(bin_age), data = subset(popsframe_10kga, pop == "EUR"), span = 0.65)
loEUR10lowerCI <- loess(sexave_lower ~ log10(bin_age), data = subset(popsframe_10kga, pop == "EUR"), span = 0.65)

loSAS10upperCI <- loess(sexave_upper ~ log10(bin_age), data = subset(popsframe_10kga, pop == "SAS"), span = 0.65)
loSAS10lowerCI <- loess(sexave_lower ~ log10(bin_age), data = subset(popsframe_10kga, pop == "SAS"), span = 0.65)

loupperCI <- (predict(loEAS10upperCI, data.frame(bin_age = predx)) + 
                predict(loEUR10upperCI, data.frame(bin_age = predx))+
                predict(loSAS10upperCI, data.frame(bin_age = predx))) / 3

lolowerCI <- (predict(loEAS10lowerCI, data.frame(bin_age = predx)) + 
                predict(loEUR10lowerCI, data.frame(bin_age = predx))+
                predict(loSAS10lowerCI, data.frame(bin_age = predx))) / 3


ggplot(subset(popsframe_10kga, pop == "AFR"), aes(x = bin_age, y = mean_interval, col = pop)) +
  geom_line(data = data.frame(mean_interval = predy, bin_age = predx), col = 'black', lwd = 1.1,
            lty = 2) +
  geom_line(data = data.frame(mean_interval = predy_noNea, bin_age = predx_noNea), col = 'black',
            lwd = 1.1) +
  geom_ribbon(data = subset(popsframe_10kga, pop == "AFR"),
              aes(ymin = sexave_lower, ymax = sexave_upper, fill = pop), color = NA, alpha = 0.1) +
  geom_ribbon(data = data.frame(mean_interval = 0,
                                sexave_lower = lolowerCI, sexave_upper = loupperCI, bin_age = predx),
              aes(ymin = sexave_lower, ymax = sexave_upper), color = NA, alpha = 0.1) +
  geom_smooth(method = "loess", span = 0.65, se = FALSE, lty = 2, lwd = 1.1) +
  geom_smooth(method = "loess", span = 0.65,
              data = subset(popsframe_10kga_noNea, pop == "AFR"), se = FALSE, lwd = 1.1) + 
  coord_cartesian(ylim = c(13.5,38.5)) + axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "generation interval (years)", x = "allele age (generations ago)")

###

for(pop in pops) {
  load(paste("bootstraps/phased_allboot_ests", pop, ".RData", sep = ""))
  
  pop_all <- read.table(paste(pop, "_allcounts.txt", sep = ""))
  names(pop_all) <- count_header
  pop_all <- cbind(pop_all[,mut_classes] + pop_all[,complementary_classes], CpG = pop_all$CpG, bin_age = pop_all$bin_age)
  
  oout <- apply(t(t(pop_all[,1:6] / rowSums(pop_all[,1:6])) - youngspecdiff), MARGIN = 1,
                predict_parentalages, model = strictCT_model)
  pop_all$mgen_time <- sapply(oout, "[[", 1)[1,]
  pop_all$fgen_time <- sapply(oout, "[[", 1)[2,]
  pop_all$ss_mean <- (pop_all$mgen_time + pop_all$fgen_time) / 2
  
  assign(paste(pop, "allcounts", sep = "_"), pop_all)
  assign(paste(pop, "phasedplotframes", sep = "_"),
         make_plotframes(phasedstrictCT_model, get(paste(pop, "allcounts", sep = "_")),
                         get(paste(pop, "phased_allboot_ests", sep = ""))))
}

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

loEASall <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_all, pop == "EAS"), span = 0.65)
loEURall <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_all, pop == "EUR"), span = 0.65)
loSASall <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_all, pop == "SAS"), span = 0.65)

predxall <- 10**seq(max(loEASall$x[1], loEURall$x[1], loSASall$x[1]),
                    min(loEASall$x[100], loEURall$x[100], loSASall$x[100]), .01)

predyall <- (predict(loEASall, data.frame(bin_age = predxall)) + 
               predict(loEURall, data.frame(bin_age = predxall)) +
               predict(loSASall, data.frame(bin_age = predxall))) / 3

for(pop in pops) {
  pop_nea <- read.table(paste("neanderthal_masked/", pop, "_noNea.allcounts.txt", sep = ""))
  names(pop_nea) <- count_header
  pop_nea <- cbind(pop_nea[,mut_classes] + pop_nea[,complementary_classes],
                   CpG = pop_nea$CpG, bin_age = pop_nea$bin_age)
  
  assign(paste(pop, "neaAll", sep = "_"), pop_nea)
  assign(paste(pop, "noNeaphasedplotframesall", sep = "_"),
         make_plotframes(phasedstrictCT_model, get(paste(pop, "neaAll", sep = "_")),
                         get(paste(pop, "phased_allboot_ests", sep = ""))))
}

popsframe_all_noNea <-
  data.frame(mean_interval = unlist(lapply(pops,
                                           function(pop) { get(paste(pop, "noNeaphasedplotframesall", sep = "_"))[[2]]$sexave })),
             bin_age = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "noNeaphasedplotframesall", sep = "_"))[[2]]$bin_age })),
             sexdiff = unlist(lapply(pops,
                                     function(pop) { get(paste(pop, "noNeaphasedplotframesall", sep = "_"))[[2]]$sexdiff })),
             sexave_upper = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "noNeaphasedplotframesall", sep = "_"))[[2]]$sexave_upper })),
             sexave_lower = unlist(lapply(pops,
                                          function(pop) { get(paste(pop, "noNeaphasedplotframesall", sep = "_"))[[2]]$sexave_lower })),
             pop = rep(pops, each = 100))

loEASall_noNea <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_all_noNea, pop == "EAS"), span = 0.65)
loEURall_noNea <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_all_noNea, pop == "EUR"), span = 0.65)
loSASall_noNea <- loess(mean_interval ~ log10(bin_age), data = subset(popsframe_all_noNea, pop == "SAS"), span = 0.65)

predxall_noNea <- 10**seq(max(loEASall_noNea$x[1], loEURall_noNea$x[1], loSASall_noNea$x[1]),
                          min(loEASall_noNea$x[100], loEURall_noNea$x[100], loSASall_noNea$x[100]), .01)

predyall_noNea <- (predict(loEASall_noNea, data.frame(bin_age = predxall_noNea)) + 
                     predict(loEURall_noNea, data.frame(bin_age = predxall_noNea)) +
                     predict(loSASall_noNea, data.frame(bin_age = predxall_noNea))) / 3

ggplot(subset(popsframe_all, pop == "AFR"), aes(x = bin_age, y = mean_interval, col = pop)) +
  geom_smooth(method = "loess", span = 0.65, se = FALSE, lwd = 1.1) +
  geom_smooth(method = "loess", span = 0.65, se = FALSE, lwd = 1.1, lty = 2,
              data = subset(popsframe_all_noNea, pop == "AFR")) +
  geom_line(data = data.frame(bin_age = predxall_noNea, mean_interval = predyall_noNea),
            col = 'black', lwd = 1.1, lty = 2) +
  geom_line(data = data.frame(bin_age = predxall, mean_interval = predyall),
            col = 'black', lwd = 1.1) +
  axis_formatting + coord_cartesian(ylim = c(13.5, 38.5), xlim = c(100,80000)) + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "generation interval (years)", x = "allele age (generations ago)")