axis_formatting <- theme_bw() + theme(axis.ticks.length=unit(0.2, "cm")) +
  theme(axis.text.x = element_text(size = 12, color = 'black')) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 12, color = 'black')) +
  theme(axis.title.y = element_text(size = 14)) + theme(panel.border = element_blank(), axis.line = element_line())

make_plotframes <- function(mu_model, snp_data, boot_data) {
  predict_out <- apply(t(t(snp_data[,1:6] / rowSums(snp_data[,1:6])) - youngspecdiff), MARGIN = 1,
                       predict_parentalages, model = mu_model)
  snp_data$mgen_time <- sapply(predict_out, "[[", 1)[1,]
  snp_data$fgen_time <- sapply(predict_out, "[[", 1)[2,]
  snp_data$sexavg_time <- (snp_data$mgen_time + snp_data$fgen_time) / 2
  
  bin_apply <- function(func) { sapply(1:100, FUN = func) }
  
  bootextract <- function(boot_data, sex_flag) {
    bin_apply(function(bin) {
      unlist(lapply(boot_data, function(model_entry) {
        model_entry[sex_flag,((bin - 1)*100+1):(bin*100)]
      }))
    })
  }
  
  male_boots <- bootextract(boot_data, 1)
  female_boots <- bootextract(boot_data, 2)
  
  snp_data$mgen_sd <- bin_apply(function(bin) { sd(male_boots[,bin]) })
  snp_data$fgen_sd <- bin_apply(function(bin) { sd(female_boots[,bin]) })
  
  boot_sexave <- bin_apply(function(bin) { (male_boots[,bin] + female_boots[,bin]) / 2 })
  boot_sexdiff <- bin_apply(function(bin) { (male_boots[,bin] - female_boots[,bin]) })
  
  snp_data$sexave_sd <- apply(boot_sexave, MARGIN = 2, FUN = sd)
  snp_data$sexdiff_sd <- apply(boot_sexdiff, MARGIN = 2, FUN = sd)
  
  plotframe <- data.frame(sex = c(rep("M", nrow(snp_data)), rep("F", nrow(snp_data))),
                          bin_age = c(snp_data$bin_age, snp_data$bin_age),
                          gen_time = c(snp_data$mgen_time, snp_data$fgen_time),
                          gen_sd = c(snp_data$mgen_sd, snp_data$fgen_sd))
  plotmean <- data.frame(bin_age = snp_data$bin_age,
                         sexave = (snp_data$mgen_time + snp_data$fgen_time) / 2,
                         sexdiff = (snp_data$mgen_time - snp_data$fgen_time),
                         sexave_sd = snp_data$sexave_sd,
                         sexdiff_sd = snp_data$sexdiff_sd)
  
  plotframe$gen_upper <-
    c(loess((mgen_time + mgen_sd) ~ log(bin_age), data = snp_data, span = 0.75)$fitted,
      loess((fgen_time + fgen_sd) ~ log(bin_age), data = snp_data, span = 0.75)$fitted)
  plotframe$gen_lower <-
    c(loess((mgen_time - mgen_sd) ~ log(bin_age), data = snp_data, span = 0.75)$fitted,
      loess((fgen_time - fgen_sd) ~ log(bin_age), data = snp_data, span = 0.75)$fitted)
  
  plotmean$sexave_upper <-
    loess((sexave + sexave_sd) ~ log(bin_age), data = plotmean, span = 0.75)$fitted
  plotmean$sexave_lower <-
    loess((sexave - sexave_sd) ~ log(bin_age), data = plotmean, span = 0.75)$fitted
  
  plotmean$sexdiff_upper <-
    loess((sexdiff + sexdiff_sd) ~ log(bin_age), data = plotmean, span = 0.75)$fitted
  plotmean$sexdiff_lower <-
    loess((sexdiff - sexdiff_sd) ~ log(bin_age), data = plotmean, span = 0.75)$fitted  
  
  return(list(plotframe, plotmean))
}

make_tgpplot <- function(plotframe, plotmean) {
  p1 <- ggplot(plotframe, aes(y = gen_time, x = bin_age, col = sex)) +
    geom_point(size = 1.8) +
    geom_smooth(data = plotmean, aes(y = sexave, x = bin_age, col = 'mean'), lwd = 1.2,
                formula = y ~ x, method = "loess", span = 0.65, se = FALSE, color = 'darkgray') +
    geom_ribbon(data = plotmean, aes(x = bin_age, y = sexave, col = 'mean',
                                     ymin = sexave_lower, ymax = sexave_upper),
                alpha = 0.15, color = NA) +
    geom_errorbar(aes(ymin = gen_time - gen_sd, ymax = gen_time + gen_sd), alpha = 0.2, width = 0) +
    scale_color_manual(values = rev(pal_aaas("default")(2))) + coord_cartesian(ylim = c(10,42)) + 
    axis_formatting + annotation_logticks(sides = "b") +
    scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
    labs(y = "generation interval (years)", x = "allele age (generations ago)", col = "sex") 
  
  p2 <- ggplot(plotframe, aes(y = gen_time, x = bin_age, col = sex)) +
    geom_point() +
    geom_smooth(data = plotmean, aes(y = sexave, x = bin_age, col = 'mean'), lwd = 1.2,
                formula = y ~ x, method = "loess", span = 0.65, se = FALSE, color = 'darkgray') +
    geom_ribbon(data = plotmean, aes(x = bin_age, y = sexave, col = 'mean',
                                     ymin = sexave_lower, ymax = sexave_upper),
                alpha = 0.15, color = NA) +
    coord_cartesian(ylim = c(10,42)) + scale_color_manual(values = rev(pal_aaas("default")(2))) +
    axis_formatting + annotation_logticks(sides = "b") +
    scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
    labs(y = "generation interval (years)", x = "allele age (generations ago)", col = "sex") 
  
  return(list(p1,p2))
}