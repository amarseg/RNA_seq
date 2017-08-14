rm(list = ls())
library('betr')
library(maSigPro)
library(gtools)
library(DESeq2)
library(edgeR)
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
library(gplots)
library(Biobase)
setwd('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/early_late_de/')
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')
source('C:/Users/am4613/Documents/GitHub/Misc/gene_plotter.R')

rpkmi= F

low_exp_cutoff <- 10 

all_counts <- read.delim("C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/analysis/all_rev_counts.txt", header= T, strings = F)
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

early <- norm_counts[,1:21]
late <- norm_counts[,c(22:36,1:3)]

es_early <- ExpressionSet(as.matrix(early))
es_late <- ExpressionSet(as.matrix(late))

exp_design_early <- exp_design[1:21,]
exp_design_late <- exp_design[22:36,]
prob_early <- betr(eset = es_early, timepoint = as.numeric(exp_design_early$Time), replicate = as.numeric(exp_design_early$Replicate), alpha = 0.05, twoCondition = FALSE)
prob_late <- betr(eset = es_late, timepoint = as.numeric(exp_design_late$Time), replicate = as.numeric(exp_design_late$Replicate), alpha = 0.05, twoCondition = FALSE)

##Warnings meaning -> This just means that for at least one gene the log ratio is 
##identical for all samples. Since this will give a zero variance (which 
##will end up in the denominator of your statistic and could possibly 
##result in an infinite value for your test statistic) it has been offset
##to a small value to prevent that possibility.
 

sig_level <- 0.99

de_early <- prob_early[prob_early > sig_level]
de_late <- prob_late[prob_late > sig_level]

gene_plotter(names(de_early), what = 'RNA')
gene_plotter(names(de_late), what = 'RNA')

write.table(de_early, sep = '\t', 'betr_early.txt')
write.table(de_late, sep = '\t', 'betr_late.txt')
