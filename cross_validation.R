# Sample 20% for cross validation

resample_SSE20 <- function() {
  cross_validation_size <- round(nrow(agg_phased_filtered) * 0.2)
  cross_validation_idx1 <- sample.int(nrow(agg_phased_filtered), size = cross_validation_size)
  cross_validation_idx2 <- sample.int(nrow(agg_phased_filtered), size = cross_validation_size)
  
  rows1 <- agg_phased_filtered[cross_validation_idx1,]
  rows2 <- agg_phased_filtered[cross_validation_idx2,]
  
  return(sum((colSums(rows1[,5:10]) / sum(sum(rows1[,5:10])) -
                colSums(rows2[,5:10]) / sum(sum(rows2[,5:10])))**2))
}

resampled_SSE20 <- replicate(1e4, resample_SSE20())
resampled_SSE <- replicate(1e2, resample_SSE())
boxplot(resampled_SSE, outline = FALSE, ylab = "SSE", ylim = c(0,4e-4))