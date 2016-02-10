###Masigpro 

library(maSigPro)
library(gtools)
library(DESeq2)
library(edgeR)
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
library(gplots)
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')

low_exp_cutoff <- 10 

all_counts <- read.delim("C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/all_rev_counts.txt", header= T, strings = F)
row.names(all_counts) <- all_counts$ID
all_counts <- all_counts[,2:37]
all_counts <- all_counts[,mixedorder(colnames(all_counts))]


##Remove genes with little reads in all samples
all_counts <- all_counts[rowSums(all_counts) > low_exp_cutoff,]

exp_design = data.frame(Time = factor(rep(0:11,each = 3)), 
												Replicate = rep(1:3,12), 
												group = rep(1,ncol(all_counts)),
												row.names = colnames(all_counts))



##Data needs to be normalised before going to maSigPro
##DESeq normalisation (divide by size factors)
cds <- DESeqDataSetFromMatrix(countData = all_counts, colData = exp_design, design = ~Time)
cds <- estimateSizeFactors(cds)
size_fact <- sizeFactors(cds)

norm_counts <- all_counts

for(i in 1:36)
{
	norm_counts[,i] <- all_counts[,i]/size_fact[i]
}

exp_design$Time <- rep(0:11, each = 3)

design <- make.design.matrix(exp_design, degree = 3, time.col = 1, repl.col = 2, group.col = 3)

fit <- p.vector(norm_counts, design, counts = T)
tstep <- T.fit(fit)
sigs <- get.siggenes(tstep, vars = 'all')
see.genes(sigs, edesign = exp_design , 
					group.cols = 3, show.fit = T, summary.mode = 'median',
					color.mode = 'gray')


###Plotting of DE genes
rpkms <- read.delim("C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/me_rpkm.txt", header= T, strings = F)
norm_rpkm <- normalise_rna(rpkms)

de_rpkm <- norm_rpkm[row.names(norm_rpkm) %in% sigs$summary,]

cols <- colorRampPalette(c('blue','gray','yellow'))

heatmap.2(as.matrix(de_rpkm), col = cols, Colv = F, trace = 'none')

masigpro_sig <- data.frame(ID = sigs$summary, pval = sigs$sig.genes$sig.pvalues$'p-value', beta0 =sigs$sig.genes$coefficients$beta0, beta1 = sigs$sig.genes$coefficients$betaTime, beta2 = sigs$sig.genes$coefficients$betaTime2,
													 beta3 = sigs$sig.genes$coefficients$betaTime3)

write.table(masigpro_sig,'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/masigpro_sig.txt', sep = '\t')

#Introduce Fold cutoff

fold_cutoff <- 1.5
n_tp <- 1

cut_de_rpkm <- de_rpkm[rowSums(de_rpkm > fold_cutoff | de_rpkm < -fold_cutoff) > n_tp,]

heatmap.2(as.matrix(cut_de_rpkm), col = cols, Colv = F, trace = 'none')

up_masigpro <- de_rpkm[rowSums(de_rpkm > fold_cutoff) > n_tp,]
down_masigpro <- de_rpkm[rowSums(de_rpkm < -fold_cutoff) > n_tp,]

write.table(up_masigpro,'upregulated_masigpro.txt', sep = '\t')
write.table(down_masigpro,'downregulated_masigpro.txt', sep = '\t')



