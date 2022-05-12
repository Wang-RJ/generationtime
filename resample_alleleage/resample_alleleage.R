# analyze redrawn allele ages
redraw_frames <- list()
for(i in 25:34) {
  redraw <- read.table(paste("Mutations/spectral/manuscript/submit/revisions/family_redraw/TGP_10kgacounts.redraw", i, sep = ""))
  names(redraw) <- count_header
  redraw <- cbind(redraw[,mut_classes] + redraw[,complementary_classes],
                  CpG = redraw$CpG, bin_age = redraw$bin_age)
  redraw_frames[[i-24]] <- make_plotframes(phasedstrictCT_model, redraw, TGPphased_boot10kga_ests)[[2]]
}

redraw_sexaves <- lapply(redraw_frames, "[", c("sexave"))

redraw_sexaves <- mapply(cbind, redraw_sexaves,
                         data.frame(bin_age = tgp_phased10kgaplots[[2]]$bin_age),
                         SIMPLIFY = FALSE)

redraw_ggframe <- do.call(rbind, mapply(cbind, redraw_sexaves,
                                        lapply(1:10, function(i) {
                                          data.frame(redraw_idx = rep(i, 100)) }),
                                        SIMPLIFY = FALSE))
names(redraw_ggframe)[2] <- "bin_age"

ggplot(redraw_ggframe, aes(y = sexave, x = bin_age, col = factor(redraw_idx))) +
  geom_smooth(method = "loess", se = FALSE, span = 0.65) +
  geom_smooth(data = tgp_phased10kgaplots[[2]], aes(y = sexave, x = bin_age, col = 'mean'), lwd = 1, lty = 3,
              formula = y ~ x, method = "loess", span = 0.65, se = FALSE, color = 'black') +
  geom_ribbon(data = tgp_phased10kgaplots[[2]], aes(x = bin_age, y = sexave, col = 'mean',
                                                    ymin = sexave_lower, ymax = sexave_upper),
              alpha = 0.15, color = NA) +
  coord_cartesian(ylim = c(13.5,38.5)) + #scale_color_manual(values = rev(pal_aaas("default")(2))) +
  axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) + 
  labs(y = "generation interval (years)", x = "allele age (generations ago)") +
  theme(legend.position = "none")
