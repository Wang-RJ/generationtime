# h38 map, maps_concat.txt, extracted from plink.GRCh38.map.zip with all .map files concatenated together

h38map <- read.table("maps_concat.txt")
names(h38map) <- c("chr", "snp", "cm", "pos", "rate")
h38map$chr <- paste("chr", h38map$chr, sep = "")

splitmap <- split(h38map, factor(h38map$chr))

phased_noCpG_noEuroTCtriplets$rec_pos <- NA
phased_noCpG_noEuroTCtriplets$rec_rate <- NA

for(map_chr in names(splitmap)) {
  decode_idx <- which(phased_noCpG_noEuroTCtriplets$Chr == map_chr)
  minmap_idx <- sapply(phased_noCpG_noEuroTCtriplets[decode_idx,]$Pos_hg38, function(decode_pos) {
    which.min(abs(splitmap[[map_chr]]$pos - decode_pos))
  })
  phased_noCpG_noEuroTCtriplets[decode_idx,]$rec_pos <- splitmap[[map_chr]][["pos"]][minmap_idx]
  phased_noCpG_noEuroTCtriplets[decode_idx,]$rec_rate <- splitmap[[map_chr]][["rate"]][minmap_idx]
}

recsorted <- phased_noCpG_noEuroTCtriplets[order(phased_noCpG_noEuroTCtriplets$rec_rate),]
recsorted <- subset(recsorted, Chr != "chrX")

begs <- floor(nrow(recsorted) / 5) * 0:4 + 1
ends <- floor(nrow(recsorted) / 5) * 1:5

quint_specs <- mapply(function(start,end) { table(recsorted$mut_class[start:end]) }, begs, ends)
quint_mat <- quint_specs / floor(nrow(recsorted) / 5)

barplot(quint_mat, col = brewer.pal(6, "Dark2"))
legend("topright", rownames(quint_specs), fill = brewer.pal(6, "Dark2"))

sapply(1:5, function(column) { predict_parentalages(phasedstrictCT_model, quint_mat[,column])$par })

recquint_ggframe <- data.frame(proportion = as.vector(quint_mat),
                               mutation = rep(rownames(quint_mat), 5),
                               quintile = rep(1:5, each = 6))

ggplot(recquint_ggframe, aes(x = quintile, y = log10(proportion), color = mutation)) +
  geom_line(lwd = 1.2) + geom_point(size = 3) + ylim(0,0.4) + axis_formatting +
  scale_color_brewer(palette = "Dark2")

ggplot(recquint_ggframe, aes(x = quintile, y = proportion, color = mutation)) +
  geom_line(lwd = 1.2) + geom_point(size = 3) + ylim(0,0.4) + axis_formatting +
  scale_color_brewer(palette = "Dark2")

ggplot(recquint_ggframe, aes(x = mutation, y = proportion, fill = factor(quintile))) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = brewer.pal(7, 'Blues')[2:6]) + axis_formatting

quint_ecdf <- ecdf(resampled_SSE20)
1 - quint_ecdf(sum((quint_mat[,2] - quint_mat[,5])**2))
1 - quint_ecdf(sum((quint_mat[,1] - quint_mat[,4])**2))
1 - quint_ecdf(sum((quint_mat[,2] - quint_mat[,4])**2))
