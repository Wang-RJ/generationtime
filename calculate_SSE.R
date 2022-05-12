# Calculate SSE in bootstraps

bootextract <- function(boot_data, sex_flag) {
  bin_apply <- function(func) { sapply(1:100, FUN = func) }
  bin_apply(function(bin) {
    unlist(lapply(boot_data, function(model_entry) {
      model_entry[sex_flag,((bin - 1)*100+1):(bin*100)]
    }))
  })
}

TGP10_maleboots <- bootextract(TGPphased_boot10kga_ests, 1)
TGP10_femaleboots <- bootextract(TGPphased_boot10kga_ests, 2)

boot_predspecs <- sapply(1:100, function(idx) {
  patageslice <- TGP10_maleboots[,idx]
  matageslice <- TGP10_femaleboots[,idx]
  return(mapply(function(p_age, m_age) { predict_spectrum(phasedstrictCT_model, p_age, m_age) },
                patageslice, matageslice, SIMPLIFY = FALSE))
})

boot_obsdiff <- sapply(1:100, function(idx) {
  bootslice <- do.call(rbind, boot_predspecs[,idx])
  bootSSE <- colSums((t(bootslice) - unlist(obs_specs[idx,]))**2)
  return(bootSSE)
}, simplify = FALSE)

boxplot(boot_obsdiff, names = c(1,rep(NA,18),20,rep(NA,19),40,rep(NA,19),60,rep(NA,19),80,rep(NA,19),100),
        outline = FALSE, ylab = "SSE", xlab = "bin", ylim = c(0,1e-3))

bootobsd_ggframe <- data.frame(bin_age = rep(tgp_10kga$bin_age, each = 10e3),
                               bin = as.factor(rep(tgp_10kga$bin_age, each = 10e3)),
                               SSE = unlist(boot_obsdiff))

ggplot(bootobsd_ggframe, aes(x = bin_age, y = SSE, fill = bin)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.2) + ylim(0, 1e-3) +
  axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "SSE", x = "allele age (generations ago)") + theme(legend.position = "none") +
  scale_fill_manual(values = rep("gray", 100))


# SSE for each population

pop_SSE <- list()
for(pop in pops) {
  pop_10kga <- get(paste(pop, "10kga", sep = "_"))
  pred_spec <- mapply(function(p_age, m_age) { predict_spectrum(phasedstrictCT_model, p_age, m_age) },
                      pop_10kga$mgen_time, pop_10kga$fgen_time)
  obs_spec <- t(t(pop_10kga[,1:6] / rowSums(pop_10kga[,1:6])) - youngspecdiff)
  
  pop_SSE[[pop]] <- data.frame(SSE = rowSums((obs_spec - t(pred_spec))**2),
                               bin_age = pop_10kga$bin_age,
                               pop = pop)
}

popSSE_ggframe <- do.call(rbind, pop_SSE)

ggplot(popSSE_ggframe, aes(x = bin_age, y = SSE, col = pop)) + geom_line(lwd = 1) +
  geom_line(data = data.frame(pop = rep("TGP", 100), bin_age = TGP_10kga$bin_age,
                              SSE = obs_diff), col = "black", lwd = 1.2) +
  axis_formatting + annotation_logticks(sides = "b") +
  scale_x_log10(minor_breaks = c(1:9 * 100, 1:9 * 1000)) +
  labs(y = "SSE", x = "allele age (generations ago)")