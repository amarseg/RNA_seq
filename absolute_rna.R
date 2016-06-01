##Absolute RNA counts
##Assuming rRNA percentage is constant!!! (this might change after experiments)
rm(list = ls())
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
source('C:/Users/am4613/Documents/GitHub/Misc/oneDplot.v4.r')
library('DESeq2')
setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/')
library('gtools')
#Load FISH data
fq_data <- read.delim('C:/Users/am4613/Documents/FISH_QUANT/Results_mature/Collated_results/results_FQ.txt', header = T, strings = F)
agg <- aggregate(Spots ~ Sample_Name*Acession, fq_data, mean)
rrna <- 0.80

#Load rna per cell 
rna_per_cell <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/rna_starting_material.txt', strings = F)
rna_per_cell$RNA_per_cell <- rna_per_cell$RNA_per_cell
rna_per_cell$Starting_RNA <- rna_per_cell$Starting_material*rna_per_cell$volume
rna_per_cell$cells <- rna_per_cell$Starting_RNA/rna_per_cell$RNA_per_cell



#Load Sam rpkms
# load('P:/CDC2_RNA_051015.rda')
# 
# sam_rpkm <- cbind(CDC2_2,CDC2_3,CDC2_4)
# 
# rpkm <- sam_rpkm[,grep(colnames(sam_rpkm), pattern = 'Annot', invert = T)]
# rpkm <- rpkm[,grep(colnames(rpkm), pattern = 'MM', invert = T)]
# rpkm <- rpkm[,grep(colnames(rpkm), pattern = 'common', invert = T)]
# rpkm <- rpkm[which(rowSums(rpkm) > 0),]
# # 
# rpkm <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/me_rpkm.txt', header = T, strings = F)
# 
#rpkm <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/rpk_me.txt', header = T, strings = F)
#rpkm <- rpkm*1000
#rpkm <- reorder_proteomics(rpkm)

all_counts <- read.delim('analysis/all_rev_counts.txt', header = T, strings = F)
row.names(all_counts) <- all_counts$ID
all_counts <- all_counts[,2:ncol(all_counts)]
all_counts <- all_counts[,mixedorder(colnames(all_counts))]

exp_design = data.frame(Time = factor(rep(0:11)), 
												Replicate = rep(1:3,12), 
												group = rep(1,ncol(all_counts)),
												row.names = colnames(all_counts))

cds <- DESeqDataSetFromMatrix(countData = all_counts, colData = exp_design, design = ~Time)
cds <- estimateSizeFactors(cds)
size_fact <- sizeFactors(cds)

rpkm <- all_counts

for(i in 1:36)
{
	rpkm[,i] <- all_counts[,i]/size_fact[i]
}

rpkm <- reorder_proteomics(rpkm)

##This bit of the code does the calibration for the replicates separately

probes_ids <- unique(fq_data$Acession)
rpkm_probes <- rpkm[row.names(rpkm) %in% probes_ids,c(1,13,25)]

parameters <- data.frame(sample = paste0('sho',c(1:3)), intercept = NA, slope = NA)
corrs <- c(1:3)

for(i in 1:3)
{
	sample_fish <- agg[agg$Sample_Name == paste0('sho',i),]
	toDo <- merge(rpkm_probes, sample_fish, by.x = 'row.names', by.y = 'Acession', all.y = T)
	fit <- lm(log2(toDo[,i+1]) ~ log2(toDo$Spots))
	corrs[i] <- summary(fit)$r.squared
	parameters[i,2] <- coefficients(fit)[[1]]
	parameters[i,3] <- coefficients(fit)[[2]]
}

cpc <- rpkm
j = 1
for(i in c(1,13,25))
{
	n_cell <- rna_per_cell[rna_per_cell$Replicate == paste0('sho',j),]
	cpc[,i:(i+11)] <- log2(rpkm[,i:(i+11)])/parameters$slope[j] - parameters$intercept[j]
	j = j+1
}
cpc <- 2^cpc
cpc <- cpc/(1 - rrna)
cpc <- cpc/rna_per_cell$cells
boxplot(log2(cpc), outline = F)

##This bit of the code does the calibration for the median of the samples

avg_cell <- aggregate(RNA_per_cell ~ Time, data = rna_per_cell, median)
avg_cell$fold_change <- avg_cell$RNA_per_cell/avg_cell$RNA_per_cell[1]

avg_fq <- aggregate(Spots ~ Acession, fq_data, mean)
var_fq <- aggregate(Spots ~ Acession, fq_data, var)

avg_rpkm <- rpkm[,1:12]

for(i in 1:12)
{
	avg_rpkm[,i] <- apply(rpkm[,seq(i,i+33,12)], 1, mean)
}
	
avg_rpkm_probes <- merge(avg_rpkm, avg_fq, by.x = 'row.names', by.y = 'Acession', all.y = T)

fit <- lm(log2(rev_tp0_sho1.bam) ~ log2(Spots) + 0 , data= avg_rpkm_probes)
intercept <- coefficients(fit)[[1]]
slope <- coefficients(fit)[[2]]
avg_cpc <- avg_rpkm

avg_norm_rpkm <- avg_rpkm/avg_rpkm[,1]

avg_cpc[,1] <- log2(avg_rpkm[,1])/intercept
avg_cpc[,1] <- 2^avg_cpc[,1]

avg_cpc[,2:12] <- avg_norm_rpkm[,2:12]*avg_cpc[,1]*avg_cell$fold_change[2:12]

boxplot(log2(avg_cpc), outline = F)
log_cpc <- avg_cpc[which(rowSums(avg_cpc) > 0),]

oneDplot(df2list(log_cpc), log  = T, spread = 200, breaks = 50, ylim = c(-5,10), col = 'darkgreen')
write.table(avg_cpc,'C:/Users/am4613/Documents/Summaries_as_timecourses/rna_cpc.txt', sep = '\t')

plot(log2(avg_cpc[,1]), log2(avg_cpc[,12]))
nc_rna <- avg_cpc[grep(row.names(avg_cpc), pattern = 'SPNC'),]
points(log2(nc_rna[,1]), log2(nc_rna[,12]), col = 'green')


##Fit all the points at the same time

meh <- c(1:36)

toDo <- merge(rpkm_probes, agg, by.x = 'row.names', by.y = 'Acession', all.y = T)
for(i in seq(1,34,3))
{
	meh[seq(i,i+2)] <- as.numeric(toDo[i,2:4])
}

agg$rpkm <- meh

fit <- lm(agg$rpkm ~ 0 +agg$Spots)
plot(x = agg$Spots, y = agg$rpkm, xlab = 'Average number of FISH counts', ylab = 'DESeq normalised expression data')
abline(fit)
corrs <- round(summary(fit)$r.squared,3)
legend('topleft', legend = paste0('R squared: ', corrs))
corrs[i] <- summary(fit)$r.squared

new_slope <- coefficients(fit)

new_cpc <- avg_rpkm
new_cpc[,1] <- avg_rpkm[,1]/new_slope
#new_cpc[,1] <- 2^new_cpc[,1]

new_cpc[,2:12] <- new_cpc[,1]*avg_norm_rpkm[,2:12]

for(i in 2:12)
{
	new_cpc[,i] <- new_cpc[,i]*avg_cell$fold_change[i]
}

col = colorRampPalette(c('white', 'darkgray'))
boxplot(new_cpc, outline = F, ylab = 'RNA copies per cell', col = col(12))


write.table(new_cpc,'rna_copies_per_cell.txt', sep ='\t')


