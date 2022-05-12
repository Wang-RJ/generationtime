### Simulate family variance
decode_pat_mean <- mean(agg_phased_filtered$Fathers_age_at_conception)
decode_mat_mean <- mean(agg_phased_filtered$Mothers_age_at_conception)

decode_pat_sd <- sd(agg_phased_filtered$Fathers_age_at_conception)
decode_mat_sd <- sd(agg_phased_filtered$Mothers_age_at_conception)

pat_spread <- rnorm(1e3, decode_pat_mean, decode_pat_sd)
mat_spread <- rnorm(1e3, decode_mat_mean, decode_mat_sd)

decode_cov <- cov(agg_phased_filtered$Fathers_age_at_conception,
                  agg_phased_filtered$Mothers_age_at_conception)

decode_sigma <- matrix(c(decode_pat_sd**2, decode_cov,
                         decode_cov, decode_mat_sd**2), ncol = 2)

total_mumodel <- lm(x ~ Fathers_age_at_conception + Mothers_age_at_conception,
                    data = agg_all)

sim_family <- function(nsims, sigma_scalar) {
  sim_ages <- mvrnorm(n = nsims, mu = c(decode_pat_mean,decode_mat_mean),
                      Sigma = sigma_scalar * decode_sigma)
  drop_idx <- which(!apply(sim_ages > 0, MARGIN = 1, all))
  if(length(drop_idx))
    sim_ages <- sim_ages[-drop_idx,]
  nsims <- nrow(sim_ages)
  
  sim_alpha <- t(apply(sim_ages, MARGIN = 1, function(ages) {
    exp(c(1, ages) %*% phasedstrictCT_model@coefficients)
  }))
  sim_total <- apply(sim_ages, MARGIN = 1, function(ages) {
    c(1, ages) %*% total_mumodel$coefficients
  })
  
  sim_muts <- rdirmn(n = nsims, size = sim_total, alpha = sim_alpha)
  return(colSums(sim_muts) / sum(sim_muts))
}

sim_famspecs <- t(sapply(seq(0.5,2,0.1), function(scalar) { sim_family(1e5, scalar) }))

sim_diff <- apply(t(t(sim_famspecs) - as.vector(human_spectrum_noCpG_noEuroTCtriplets)),
                  MARGIN = 1, function(x) {
                    sum(x**2)
                  })

plot(sim_diff ~ seq(0.5,2,0.1), ylim = c(0,6e-5), xlab = "sigma scalar", ylab = "SSE", pch = 21,
     bg = 'grey', cex = 1.2)

pred_specs <- mapply(function(p_age, m_age) { predict_spectrum(phasedstrictCT_model, p_age, m_age) },
                     plotframe$gen_time[1:100], plotframe$gen_time[101:200])
obs_specs <- t(t(tgp_10kga[,1:6] / rowSums(tgp_10kga[,1:6])) - youngspecdiff)

obs_diff <- mapply(function(diff, pred) { sum(diff**2)},
                   split(t(t(obs_specs) - pred_specs), 1:100),
                   split(pred_specs, 1:100))

plot(obs_diff, pch = 21, bg = 'grey', ylim = c(0, 4e-4), ylab = "SSE", xlab = 'bin')