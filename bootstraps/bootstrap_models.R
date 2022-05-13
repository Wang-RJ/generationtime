# Resample mutation data to create 100 bootstrap models
# Bootstrap estimates calculated distributed cluster call in bootstrap_call.R

phased_boot_models <- lapply(replicate(100, list(sample.int(nrow(agg_phased_filtered),
                                                            size = nrow(agg_phased_filtered),
                                                            replace = TRUE))),
                      FUN = function(sampled_row_idx) {
                        dup_idx <- which(duplicated(sampled_row_idx))
                        boot_data <- agg_phased_filtered[sampled_row_idx,]
                        out <-
                          purrr::quietly(MGLMreg)(formula = cbind(A_C, A_G, A_T, C_A, C_G, C_T) ~
                                  Fathers_age_at_conception + Mothers_age_at_conception,
                                dist = "DM",
                                data = boot_data)
                        if(length(grep("hessian", out$warnings)) > 0)
                          return(NA)
                        else
                          return(out$result)
       })

unphased_boot_models <- lapply(replicate(100, list(sample.int(nrow(agg_unphased_filtered),
                                                     size = nrow(agg_unphased_filtered),
                                                     replace = TRUE))),
                      FUN = function(sampled_row_idx) {
                        boot_data <- agg_unphased_filtered[sampled_row_idx,]
                        out <-
                          purrr::quietly(MGLMreg)(formula = cbind(A_C, A_G, A_T, C_A, C_G, C_T) ~
                                  Fathers_age_at_conception + Mothers_age_at_conception,
                                dist = "DM",
                                data = boot_data)
                        if(length(grep("hessian", out$warnings)) > 0)
                          return(NA)
                        else
                          return(out$result)
                      })
