setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/')
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
library('gtools')

metabolic <- read.delim('aa_carbohidrate_metabolism', header = T, strings = F)
metabolic_id <- metabolic$ensembl_id
ribosomes <- read.delim('ribosome&biogenesisGenes', header = T, strings = F)
ribosomes_id <- ribosomes$ensembl_id

load('P:/CDC2_RNA_051015.rda')

sam_rpkm <- list()

sam_rpkm[[1]] <- CDC2
sam_rpkm[[2]] <- CDC2_2
sam_rpkm[[3]] <- CDC2_3
sam_rpkm[[4]] <- CDC2_4

norm_sam <- cbind(sam_rpkm[[2]][,1:12],sam_rpkm[[3]][,1:12],sam_rpkm[[4]][,1:12])
norm_sam <- norm_sam[,mixedorder(colnames(norm_sam))]


fractions <- as.data.frame(matrix(nrow = 3, ncol = ncol(norm_sam), NA))
row.names(fractions) <- c('ribosome', 'metabolism', 'rest')
colnames(fractions) <- colnames(norm_sam)

met <- colSums(norm_sam[which(row.names(norm_sam) %in% metabolic_id),])
ribo <- colSums(norm_sam[which(row.names(norm_sam) %in% ribosomes_id), ])
total <- colSums(norm_sam)

fractions[1,] <- ribo/total
fractions[2,] <- met/total
fractions[3,] <- (total - ribo - met)/total

barplot(as.matrix(fractions), col = terrain.colors(3))
