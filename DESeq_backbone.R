##DESeq analysis comparing each time point with the previous one
library("DESeq2")

###Read and prepare data
se.path <- ("~/rna_seq0206/seq_0206/Transcripts_1.rda")

counts <- assays(se)$counts[,seq(1,24,2)]
row.names(counts) <- elementMetadata(se)$ID

colData <- as.data.frame(factor(rep(c(0,1,10,11,2,3,4,5,6,7,8,9))))
rownames(colData) <- colnames(counts)
colnames(colData) <- "timepoint"

test <- DESeqDataSetFromMatrix(countData = counts,
															 colData = colData ,
															 design = ~ timepoint)

##Test for differential expression
test <- DESeq(test)
res <- results(test)

#Get significative hits using the adjusted p-value
alpha = 0.05
resSig <- subset(res, padj < alpha)

#Save diagnostic plot as pdf
pdf("Diagnostics.pdf")
plotDispEsts(test)
plotMA(res)
dev.off()

