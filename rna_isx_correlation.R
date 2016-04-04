##
rm(list = ls())

setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/rna_isx_correlation/')
library('lattice')
library('TSdist')
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
library('clusterProfiler')
source('C:/Users/am4613/Documents/GitHub/RNA_seq/cell_cycle_plotter.R')

isx_data <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/isx_data_summary.txt')


rep_cpc <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/me_rpkm.txt', header = T, strings = F)
	
rep_cpc <- normalise_rna(rep_cpc)


sample_info <- cbind(colnames(rep_cpc), rep(1:3,3))
sample_info <- as.data.frame(sample_info)
colnames(sample_info) <- c('SampleName','Replicate')


pearson_coeff <- rep_cpc[,1:3]
colnames(pearson_coeff) <- c('Replicate1','Replicate2', 'Replicate3')

coeff <- rep_cpc[,1:3]
##Correlation with length


par(mfrow = c(3,2))
for(i in 1:3)
{
	lengths <- isx_data[isx_data$rep == i,]
	samples <- rep_cpc[, colnames(rep_cpc) %in% sample_info[sample_info$Replicate == i,]$SampleName]
	pearson_coeff[,i] <- cor(t(samples), lengths$Length_Erode.M03..4., method = 'pearson')
	coeff [,i] <- coefficients(lm(t(samples)~lengths$Length_Erode.M03..4.))[2,]
	hist(pearson_coeff[,i], breaks = 50)
	hist(coeff[,i], breaks = 100)
}

cutoff_pearson = -0.5
cutoff_slope = -0.05
n_rep = 2


all_coeff <- cbind(pearson_coeff, coeff)
colnames(all_coeff) <- c('Pearson1','Pearson2','Pearson3','Slope1','Slope2','Slope3')
row.names(all_coeff) <- row.names(rep_cpc)

gene_list_length <- all_coeff[which(rowSums(all_coeff[,1:3] < cutoff_pearson) >= n_rep),]
gene_list_length <- gene_list_length[which(rowSums(gene_list_length[,4:6] < cutoff_slope) >= n_rep),]

write.table(gene_list_length, 'negative_correlation_length.txt', sep = '\t')

par(mfrow = c(1,3))
for(i in 1:3)
{
	plot(pearson_coeff[,i], coeff[,i])
	points(gene_list_length[,i], gene_list_length[,i+3], col = 'red')
}


##Correlation with Area

par(mfrow = c(3,2))
for(i in 1:3)
{
	lengths <- isx_data[isx_data$rep == i,]
	samples <- rep_cpc[, colnames(rep_cpc) %in% sample_info[sample_info$Replicate == i,]$SampleName]
	pearson_coeff[,i] <- cor(t(samples), lengths$Area_Erode.M03..4., method = 'pearson')
	coeff [,i] <- coefficients(lm(t(samples)~lengths$Area_Erode.M03..4.))[2,]
	hist(pearson_coeff[,i], breaks = 50)
	hist(coeff[,i], breaks = 100)
}
par(mfrow = c(1,1))


cutoff_pearson = -0.5
cutoff_slope = -0.01
n_rep = 2

all_coeff <- cbind(pearson_coeff, coeff)
colnames(all_coeff) <- c('Pearson1','Pearson2','Pearson3','Slope1','Slope2','Slope3')
row.names(all_coeff) <- row.names(rep_cpc)

gene_list_area <- all_coeff[which(rowSums(all_coeff[,1:3] < cutoff_pearson) >= n_rep),]
gene_list_area <- gene_list_area[which(rowSums(gene_list_area[,4:6] < cutoff_slope) >= n_rep),]

par(mfrow = c(1,3))
for(i in 1:3)
{
	plot(pearson_coeff[,i], coeff[,i])
	points(gene_list_area[,i], gene_list_area[,i+3], col = 'red')
}

write.table(gene_list_area, 'negative_correlation_area.txt', sep = '\t')

###Correlation with volume

isx_data$Volume <- (pi*(isx_data$Width_Erode.M03..4.^2)*(isx_data$Length_Erode.M03..4. - isx_data$Width_Erode.M03..4./3))/4

par(mfrow = c(3,2))
for(i in 1:3)
{
	lengths <- isx_data[isx_data$rep == i,]
	samples <- rep_cpc[, colnames(rep_cpc) %in% sample_info[sample_info$Replicate == i,]$SampleName]
	pearson_coeff[,i] <- cor(t(samples), lengths$Volume, method = 'pearson')
	coeff [,i] <- coefficients(lm(t(samples)~lengths$Volume))[2,]
	hist(pearson_coeff[,i], breaks = 50)
	hist(coeff[,i], breaks = 100)
}
par(mfrow = c(1,1))

par(mfrow = c(1,3))
for(i in 1:3)
{
	plot(pearson_coeff[,i], coeff[,i])
}

cutoff_pearson = -0.5
cutoff_slope = -0.001
n_rep = 2

all_coeff <- cbind(pearson_coeff, coeff)
colnames(all_coeff) <- c('Pearson1','Pearson2','Pearson3','Slope1','Slope2','Slope3')
row.names(all_coeff) <- row.names(rep_cpc)

gene_list_volume <- all_coeff[which(rowSums(all_coeff[,1:3] < cutoff_pearson) >= n_rep),]
gene_list_volume <- gene_list_volume[which(rowSums(gene_list_volume[,4:6] < cutoff_slope) >= n_rep),]

write.table(gene_list_volume, 'negative_correlation_volume.txt', sep = '\t')

negative_correlation <- list(gene_list_length, gene_list_area, gene_list_volume)
names(negative_correlation) <- c('Length','Area','Volume')
save(negative_correlation, file = 'Negative_correlation.rda')
