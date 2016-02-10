##Plot heatmap with dtw distance

setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/')
load('analysis/warp_distance.rda')
source('../GitHub/Proteomics/normalise_script.R')
library('gplots')

rpkm <- read.delim('analysis/me_rpkm.txt', header = T, strings = F)

norm_rpkm <- normalise_rna(rpkm)
avg_norm <- median_rna(norm_rpkm)

col = colorRampPalette(c('blue','grey','yellow'))
hc = hclust(warp_dis)

heatmap.2(as.matrix(norm_rpkm), col = colorRampPalette(c('blue','grey','yellow')) , Colv = F, trace = 'none',
					Rowv = as.dendrogram(hc))

