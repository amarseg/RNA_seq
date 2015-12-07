library('gplots')

load('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/CDC2_RNA_051015.rda')

sam_rpkm <- list()

sam_rpkm[[1]] <- CDC2
sam_rpkm[[2]] <- CDC2_2
sam_rpkm[[3]] <- CDC2_3
sam_rpkm[[4]] <- CDC2_4

for(i in 2:4)
{
	sam_rpkm[[i]][,1:12] <- log2(sam_rpkm[[i]][,1:12]/sam_rpkm[[i]][,1])
	is.na(sam_rpkm[[i]]) <- sapply(sam_rpkm[[i]], is.nan)
	is.na(sam_rpkm[[i]]) <- sapply(sam_rpkm[[i]], is.infinite)
	sam_rpkm[[i]][is.na(sam_rpkm[[i]])] <- 0
	
}

norm_sam <- cbind(sam_rpkm[[2]][,1:12],sam_rpkm[[3]][,1:12],sam_rpkm[[4]][,1:12])

nc_rna <- norm_sam[grep(row.names(norm_sam), pattern = 'NCRNA'),]

heatmap.2(as.matrix(nc_rna), col = colorRampPalette(c('blue','gray','yellow')), trace = 'none', Colv = F)

