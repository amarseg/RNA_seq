###Annotation generator
###Generates a lists of genes with their pertenence to the different fractions in the 
###modelling based on Andrea's paper
setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/')
source('../GitHub/Proteomics/normalise_script.R')
protein = F
if(protein)
{
	abs_data <- read.delim('analysis/rep_cpc.txt', header = T, strings = F)
	abs_data <- 2^abs_data
	
	data <- read.delim('analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)
	norm_data <- normalise_ProtDataset(data, what = 'heatmap')
	norm_data <- norm_data[,7:42]
	avg_data <- median_rna(norm_data)
}else{
	abs_data <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/rna_cpc.txt',
										 header = T, strings = F)
	data <- read.delim('analysis/me_rpkm.txt', header = T)
	norm_data <- normalise_rna(data)
	avg_data <- median_rna(data)
}


#Get gene_list 



gene_list <- data.frame(gene_name = row.names(norm_data), annotation = NA)

##Annotate ribosomes

ribo_id <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/rproteins.txt', 
											header = T, strings = F)

gene_list[which(gene_list$gene_name %in% ribo_id$Sp_name), 2] <- 'ribosome'

##Annotate polymerases

poly1 <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/polymerase1', 
											header = T, strings = F)
pol2 <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/polymerase2', 
									 header = T, strings = F)
pol3 <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/polymerase3', 
									 header = T, strings = F)

all_pol <- unique(c(poly1$ensembl_id, pol2$ensembl_id, pol3$ensembl_id))
gene_list[which(gene_list$gene_name %in% all_pol), 2] <- 'polymerases'

##Annotate transporters

transporters <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/transporters.txt', header = T, strings = F)
transporters <- transporters$ensembl_id

gene_list[which(gene_list$gene_name %in% transporters), 2] <- 'transporters'


##Chromatin structure
chromatin <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/chromatin_remodelling_histones.txt', header = T, strings = F)
chromatin = chromatin$ensembl_id

gene_list[which(gene_list$gene_name %in% chromatin), 2] <- 'chromatin structure'

##Annotate metabolism

metabolism <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/aa_carbohidrate_metabolism', header = T, strings = F)
gene_list[which(gene_list$gene_name %in% metabolism$ensembl_id), 2] <- 'metabolism'


##Annotate stress, we are gonna get genes that are up or downregulated two fold in more than 5 time points
stress <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/stress', header = T, strings = F)
stress <- stress$ensembl_id

fold_cutoff <- 2
n_samples <- 5

induced_genes <- avg_data[rowSums( avg_data <= -fold_cutoff | avg_data >= fold_cutoff) > n_samples, ]

stress_list <- intersect(stress, row.names(induced_genes))

gene_list[which(gene_list$gene_name %in% stress_list), 2] <- 'stress'


##Constant fraction
sd_cutoff <- 2

gene_list$StandardDeviation <- apply(norm_data, 1, sd)

for(i in 1:nrow(gene_list))
{
	if(gene_list[i,]$StandardDeviation <= sd_cutoff & is.na(gene_list[i,]$annotation))
	{
		gene_list[i,]$annotation <- c('Constant')
	}
}
##Not included in the analysis


gene_list[which(is.na(gene_list$annotation)),]$annotation <- 'NotIncluded'

##Summarize data based on the annotation and plot (fingers crossed)

if(protein)
{
	abs_data <- median_rna(abs_data)
}

agg <- aggregate(abs_data, by = list(gene_list$annotation), median, na.rm = T)
agg[,2:13] <- apply(agg[,2:13],2, as.numeric)

col = rainbow(7)

par(mfrow = c(1,1))
plot(y = agg[1,2:13], x = rep(0:11), ylim = c(0,max(agg[,2:12])),type = 'l',  col = col[1])
for(i in 2:8)
{
	lines(y = agg[i,2:13], x = rep(0:11), col = col[i])
	
}
lines(y = apply((abs_data[,1:12]),2,median), x = rep(0:11), col = 'black')
legend('topleft',legend = agg[,1], col = col, fill = col, cex = 0.5)

