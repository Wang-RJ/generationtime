source("age_modeling.R")
source("plot_helper.R")

options(stringsAsFactors = FALSE)

count_header <- c("A_C", "ACC_ATC", "A_G", "A_T", "C_A", "CCC_CTC", "C_G", "CpG", "C_T", "G_A",
                  "G_C", "G_T", "T_A", "T_C", "TCC_TTC", "TCT_TTT", "T_G", "bin_age")

pops <- c("AFR", "EAS", "EUR", "SAS")

## Read in dated TGP dataset, bootstraps, and plot

load("bootstraps/unphased_10kgaboot_estsTGP.RData")
load("bootstraps/phased_10kgaboot_estsTGP.RData")

tgp_all <- read.table("TGP_allcounts.txt")
names(tgp_all) <- count_header
tgp_all <- cbind(tgp_all[,mut_classes] + tgp_all[,complementary_classes],
                 CpG = tgp_all$CpG, bin_age = tgp_all$bin_age)

tgpall_phasedplotframes <- make_plotframes(phasedstrictCT_model, tgp_all, phased_boot_ests)
tgpall_phasedallplots <- make_tgpplot(tgpall_phasedplotframes[[1]], tgpall_phasedplotframes[[2]])
tgp_phasedallplots[[1]]
tgp_phasedallplots[[2]] + ggtitle("TGP phased")

tgp_unphasedallplotframes <- make_plotframes(strictCT_model, tgp_all, unphased_boot_ests)
tgp_unphasedallplots <- make_tgpplot(tgp_unphasedallplotframes[[1]], tgp_unphasedallplotframes[[2]])
tgp_unphasedallplots[[1]]
tgp_unphasedallplots[[2]] + ggtitle("TGP unphased")

tgp_all <- read.table("TGP_allcounts.txt")
names(tgp_all) <- count_header
tgp_all <- cbind(tgp_all[,mut_classes] + tgp_all[,complementary_classes],
                 CpG = tgp_all$CpG, bin_age = tgp_all$bin_age)

tgp_phasedallplots <- make_plotframes(phasedstrictCT_model, tgp_all, phased_boot_ests)
tgp_phasedallplots[[1]]
tgp_phasedallplots[[2]]

tgp_unphasedallplots <- make_plotframes(strictCT_model, tgp_all, unphased_boot_ests)
tgp_unphasedallplots[[1]]
tgp_unphasedallplots[[2]]

tgp_10kga <- read.table("TGP_10kgacounts.txt")
names(tgp_10kga) <- count_header
tgp_10kga <- cbind(tgp_10kga[,mut_classes] + tgp_10kga[,complementary_classes],
                   CpG = tgp_10kga$CpG, bin_age = tgp_10kga$bin_age)

tgp_phased10kgaplots <- make_plotframes(phasedstrictCT_model, tgp_10kga, TGPphased_boot10kga_ests)
plotframe <- tgp_phased10kgaplots[[1]]
plotmean <- tgp_phased10kgaplots[[2]]

ggplot(plotframe, aes(y = gen_time, x = bin_age, col = sex)) +
  geom_point(size = 2) +
  geom_smooth(data = plotmean, aes(y = sexave, x = bin_age, col = 'mean'), lwd = 1.2,
              formula = y ~ x, method = "loess", span = 0.65, se = FALSE, color = 'darkgray') +
  geom_ribbon(data = plotmean, aes(x = bin_age, y = sexave, col = 'mean',
                                   ymin = sexave_lower, ymax = sexave_upper),
              alpha = 0.15, color = NA) +
  geom_errorbar(aes(ymin = gen_time - gen_sd, ymax = gen_time + gen_sd), lwd = .5,
                alpha = 0.2, width = 0) +
  coord_cartesian(ylim = c(13.5,38.5)) + scale_color_manual(values = rev(pal_aaas("default")(2))) +
  axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "generation interval (years)", x = "allele age (generations ago)", col = "sex")

tgp_unphased10kgaplots <- make_plotframes(strictCT_model, tgp_10kga, TGPunphased_boot10kga_ests)
unp_plotframe <- tgp_unphased10kgaplots[[1]]
unp_plotmean <- tgp_unphased10kgaplots[[2]]

ggplot(unp_plotframe, aes(y = gen_time, x = bin_age, col = sex)) +
  geom_point(size = 2) +
  geom_smooth(data = unp_plotmean, aes(y = sexave, x = bin_age, col = 'mean'), lwd = 1.2,
              formula = y ~ x, method = "loess", span = 0.65, se = FALSE, color = 'darkgray') +
  geom_ribbon(data = unp_plotmean, aes(x = bin_age, y = sexave, col = 'mean',
                                       ymin = sexave_lower, ymax = sexave_upper),
              alpha = 0.15, color = NA) +
  geom_errorbar(aes(ymin = gen_time - gen_sd, ymax = gen_time + gen_sd), lwd = .5,
                alpha = 0.2, width = 0) +
  coord_cartesian(ylim = c(13.5,38.5)) + scale_color_manual(values = rev(pal_aaas("default")(2))) +
  axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "generation interval (years)", x = "allele age (generations ago)", col = "sex")

smoothed <- loess(formula = sexave ~ log(bin_age), data = tgp_phased10kgaplots[[2]], span = 0.65)
agebygen <- cumsum(c(rep(predict(smoothed, data.frame(bin_age = 78)), 77),
                     predict(smoothed, data.frame(bin_age = 78:9330)),
                     rep(predict(smoothed, data.frame(bin_age = 9330)), 10e3-9330)))

plotframe$year_age <- agebygen[round(plotframe$bin_age)]
ggplot(plotframe, aes(y = gen_time, x = year_age, col = sex)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = gen_time - gen_sd, ymax = gen_time + gen_sd), lwd = .5,
                alpha = 0.2, width = 0) +
  coord_cartesian(ylim = c(13.5,38.5)) + scale_color_manual(values = rev(pal_aaas("default")(2))) +
  axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "generation interval (years)", x = "years ago", col = "sex")

ggplot(plotmean, aes(y = sexdiff, x = bin_age)) +
  geom_smooth(lwd = 1.2, formula = y ~ x, method = "loess", se = FALSE, span = 0.7) +
  #  geom_ribbon(aes(ymin = sexdiff_upper, ymax = sexdiff_lower), alpha = 0.2) +
  axis_formatting + annotation_logticks(sides = "b") + scale_color_aaas() + scale_fill_aaas() +
  coord_cartesian(ylim = c(0, 10)) + scale_y_continuous(breaks = c(5,10)) +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000))

ggplot(unp_plotmean, aes(y = sexdiff, x = bin_age)) +
  geom_smooth(lwd = 1.2, formula = y ~ x, method = "loess", se = FALSE, span = 0.7) +
  #  geom_ribbon(aes(ymin = sexdiff_upper, ymax = sexdiff_lower), alpha = 0.2) +
  axis_formatting + annotation_logticks(sides = "b") + scale_color_aaas() + scale_fill_aaas() +
  coord_cartesian(ylim = c(0, 12)) + scale_y_continuous(breaks = c(5,10)) +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000))

tgp_unphased10kgaplots <- make_plotframes(strictCT_model, tgp_10kga, TGPunphased_boot10kga_ests)
tgp_unphased10kgaplots[[1]]
tgp_unphased10kgaplots[[2]]

## Functions to plot spectral across 10k generations for each population

freq_plot_abs <- function(pop) {
  popincl <- read.table(paste(pop, "_10kgacounts.txt", sep = ""))
  names(popincl) <- count_header
  popincl <- cbind(popincl[,mut_classes] + popincl[,complementary_classes], CpG = popincl$CpG,
                   ACC_ATC = popincl$ACC_ATC, CCC_CTC = popincl$CCC_CTC,
                   TCC_TTC = popincl$TCC_TTC, TCT_TTT = popincl$TCT_TTT,
                   bin_age = popincl$bin_age)

  hspec_ggframe <- data.frame(delta_frequency =
                                as.vector(
                                  apply(popincl[,1:7]/rowSums(popincl[,1:7]), MARGIN = 1,
                                        function(row) {
                                          return(clr(row))
                                        })),

                              class = rep(c(names(popincl)[1:6], "CpG"), nrow(popincl)),
                              age = rep(popincl$bin_age, each = 7),
                              idx = rep(1:100, each = 7))
  
  hspec_ggframe$fit_frequency = as.vector(t(sapply(mut_classes, FUN = function(c_idx) {
    loess(delta_frequency ~ log10(age), data = subset(hspec_ggframe, class == c_idx))$fitted })))
  
  ggplot(hspec_ggframe, aes(x = age, y = delta_frequency, color = class)) +
    geom_line(lwd = 1.2) +
    scale_color_manual(values = brewer.pal(6, "Dark2")[c(2,3,6,5,4,1,7)]) +
    ylab("frequency") + xlab("generations ago") + annotation_logticks(sides = "b") +
    scale_x_continuous(trans = "log10") +
    coord_cartesian(ylim = c(-.75,1))
}

freq_plot_rel <- function(pop) {
  popincl <- read.table(paste(pop, "_10kgacounts.txt", sep = ""))
  names(popincl) <- count_header
  popincl <- cbind(popincl[,mut_classes] + popincl[,complementary_classes], CpG = popincl$CpG,
                   ACC_ATC = popincl$ACC_ATC, CCC_CTC = popincl$CCC_CTC,
                   TCC_TTC = popincl$TCC_TTC, TCC_TTC = popincl$TCC_TTC,
                   bin_age = popincl$bin_age)
  
  hspec_ggframe <- data.frame(delta_frequency =
                                as.vector(
                                  apply(popincl[,1:6]/rowSums(popincl[,1:6]), MARGIN = 1,
                                        function(row) {
                                          return(row - pd_noCpGnoTCtrip - youngspecdiff)
                                        })),
                              class = rep(names(popincl)[1:6], nrow(popincl)),
                              age = rep(popincl$bin_age, each = 6),
                              idx = rep(1:100, each = 6))
  
  hspec_ggframe$fit_frequency = as.vector(t(sapply(mut_classes, FUN = function(c_idx) {
    loess(delta_frequency ~ log10(age), data = subset(hspec_ggframe, class == c_idx))$fitted })))
  
  ggplot(hspec_ggframe, aes(x = age, y = delta_frequency, col = class)) +
    geom_line(size = 0.75, alpha = 0.3) +
    axis_formatting + geom_smooth(se = FALSE, size = 1.75, span = 0.65) +
    scale_color_manual(values = brewer.pal(6, "Dark2")[c(2,3,6,5,4,1)]) +
    ylab("delta-frequency") + xlab("allele age") + annotation_logticks(sides = "b") +
    coord_cartesian(ylim = c(-.008,0.0075)) + scale_x_continuous(trans = "log10")
}