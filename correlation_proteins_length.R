##
library('dtw')
library('lattice')
library('TSdist')
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')

absolute = F

rna_cell <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/rna_percell.txt', header = T, strings = F)
isx_data <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/isx_data_summary.txt')

if(absolute)
{
	rep_cpc <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/rep_cpc.txt', header = T, strings = F)
	rep_cpc <- 2^rep_cpc
}else{
	rep_cpc <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)

	rep_cpc <- normalise_ProtDataset(rep_cpc, what = 'heatmap')
	rep_cpc <- rep_cpc[,7:42]
	#rep_cpc <- reorder_proteomics(rep_cpc)
	rep_cpc 
}


sample_info <- cbind(colnames(rep_cpc), rep(1:3,3))
sample_info <- as.data.frame(sample_info)
colnames(sample_info) <- c('SampleName','Replicate')


pearson_coeff <- rep_cpc[,1:3]
colnames(pearson_coeff) <- c('Replicate1','Replicate2', 'Replicate3')

coeff <- rep_cpc[,1:3]

par(mfrow = c(3,2))
for(i in 1:3)
{
	lengths <- isx_data[isx_data$rep == i,]
	samples <- rep_cpc[, colnames(rep_cpc) %in% sample_info[sample_info$Replicate == i,]$SampleName]
	pearson_coeff[,i] <- cor(t(samples), lengths$Length_Erode.M03..4., method = 'pearson')
	#pearson_coeff[,i] <- DTWDistance(t(samples), lengths$Length_Erode.M03..4.)
	coeff [,i] <- coefficients(lm(t(samples)~lengths$Length_Erode.M03..4.))[2,]
	hist(pearson_coeff[,i], breaks = 50)
	hist(coeff[,i], breaks = 100)
}
par(mfrow = c(1,1))

par(mfrow = c(1,3))
for(i in 1:3)
{
	plot(pearson_coeff[,i], coeff[,i])
}

cutoff_pearson = -0.3
cutoff_slope = -0.025
n_rep = 2

all_coeff <- cbind(pearson_coeff, coeff)
colnames(all_coeff) <- c('Pearson1','Pearson2','Pearson3','Slope1','Slope2','Slope3')
row.names(all_coeff) <- row.names(rep_cpc)

gene_list <- all_coeff[which(rowSums(all_coeff[,1:3] < cutoff_pearson) >= 2),]
gene_list <- gene_list[which(rowSums(gene_list[,4:6] < cutoff_slope) >= 2),]

go_terms <- GOanalysis(row.names(gene_list), go = GOtable, all = 5123)
write.table(gene_list, 'pearson_slope_list.txt', sep = '\t')
