rm(list = ls())
library('betr')
library(maSigPro)
library(gtools)
library(DESeq2)
library(edgeR)
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
library(gplots)
library(Biobase)
setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/betr_clustering/DESeq_norm/')
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')
source('C:/Users/am4613/Documents/GitHub/Misc/gene_plotter.R')

rpkm = F

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

es_norm <- ExpressionSet(as.matrix(norm_counts))

prob <- betr(eset = es_norm, timepoint = exp_design$Time, replicate = exp_design$Replicate, alpha = 0.05, twoCondition = FALSE)
##Warnings meaning -> This just means that for at least one gene the log ratio is 
##identical for all samples. Since this will give a zero variance (which 
##will end up in the denominator of your statistic and could possibly 
##result in an infinite value for your test statistic) it has been offset
##to a small value to prevent that possibility.

sig_level <- 0.99

de_betr <- prob[prob > sig_level]

write.table(de_betr, 'C:/Users/am4613/Documents/Summaries_as_timecourses/DE_genes/result_betr.txt', sep = '\t')

##See how DE genes look in the data
if(rpkm == T)
{
	rpkms <- read.delim("C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/me_rpkm.txt", header= T, strings = F)
	norm_rpkm <- normalise_rna(rpkms)
}else{
	norm_rpkm <- reorder_proteomics(norm_counts)
	norm_rpkm <- normalise_rna(norm_rpkm)
}

de_rpkm <- norm_rpkm[row.names(norm_rpkm) %in% names(de_betr),]

cols <- colorRampPalette(c('blue','gray','yellow'))

heatmap.2(as.matrix(de_rpkm), col = cols, Colv = F, trace = 'none')

hc.rows <- hclust(dist(de_rpkm))
plot(hc.rows)

k = 2

rect.hclust(hc.rows,k)

clusters<- cutree(hc.rows, k = k)

pval = 0.05

pdf(file  = 'results_betr_deseq.pdf')
for(i in 1:k)
{
	cl <- names(clusters[clusters == i])
	write.table(cl, file = paste0('Cluster_',i,'.txt'), sep = '\t')
	
	go_results <- GOanalysis(cl, GOtable, all = 7102)
	go_results <- go_results[go_results[,2] > pval, ]
	write.table(go_results, file = paste0('GO_cluster',i,'.txt'), sep = '\t')
	
	gene_plotter(cl, what = 'RNA', norm = 'DESeq')
}
dev.off()

