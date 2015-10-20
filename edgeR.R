###edgeR analysis

library(edgeR)
library(gtools)

low_exp_cutoff <- 10

all_counts <- read.delim("C:/Users/am4613/Desktop/all_rev_counts.txt", header= T, strings = F)
row.names(all_counts) <- all_counts$ID
all_counts <- all_counts[,2:37]
all_counts <- all_counts[,mixedorder(colnames(all_counts))]

##Remove genes with little reads in all samples
all_counts <- all_counts[rowSums(all_counts) > low_exp_cutoff,]

exp_design = data.frame(Time = rep(0:11,each = 3), 
												Replicate = rep(1:3,12), 
												row.names = colnames(all_counts))


dge <- DGEList(counts = as.matrix(all_counts))
design <- model.matrix(~exp_design$Time)
dge <- estimateGLMCommonDisp(dge, design, verbose = T)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)

lrt <- glmLRT(fit)

topTags(lrt)

go <- goana(lrt)

