#DESEq for comparing cdc2-33 and cdc2-as
library("DESeq2")

cdc2.33 <- c("")
cdc2.as <- c("")

load(cdc2.as)
counts.as<- assays(se)$counts[,seq(1,24,2)]
row.names(counts.as) <- elementMetadata(se)$ID
row.names(counts.as) <- substr(elementMetadata(se)$ID, start = 12, stop = nchar(elementMetadata(se)$ID))

counts.33 <- read.csv(cdc2.33, header = T, strings = F)

all.counts <- merge(counts.33, counts.as, by = row.names, all.x = T)

colData <- as.data.frame(factor(rep(c(0,1,10,11,2,3,4,5,6,7,8,9)), c(0:11)))
rownames(colData) <- colnames(counts.as)
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

