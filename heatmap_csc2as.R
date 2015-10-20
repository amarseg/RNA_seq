library("rtracklayer")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library('gplots')
library('edgeR')
library('gtools')

gff.path <- 'C:/Users/am4613/Desktop/am4613/Schizosaccharomyces_pombe.ASM294v2.26.gff3'
all_counts <- read.delim("C:/Users/am4613/Desktop/all_rev_counts.txt", header= T, strings = F)
row.names(all_counts) <- all_counts[,1]
all_counts <- all_counts[,2:37]
all_counts <- all_counts[,mixedorder(colnames(all_counts))]

gff <- import.gff3(gff.path)
gff.exons <- gff[which(elementMetadata(gff)[,"type"] == "exon" | elementMetadata(gff)[,"type.1"] == "LTR") ]
gff.exons[which(elementMetadata(gff.exons)[,"type.1"] == "LTR")]$ID <- paste0("LTR:", gff.exons[which(elementMetadata(gff.exons)[,"type.1"] == "LTR")]$Name)
exons <- gff.exons[which(elementMetadata(gff.exons)[,"type"] == "exon")]$Parent
exons <- substr(exons, 1, nchar(exons) - 2)
gff.exons[which(elementMetadata(gff.exons)[,"type"] == "exon")]$ID <- exons
a <- strsplit(elementMetadata(gff.exons)$ID, "\\:")
vector_with_IDs <- sapply(a,"[[",2)
geneLengths <- width(gff.exons)# Length of exon union per gene in kbp
geneLengths <- data.frame(ID = vector_with_IDs, length = geneLengths)

geneLengths <- aggregate(geneLengths[,2], FUN = sum, na.rm = T, by = list(geneLengths[,1]))

rpkm <- rpkm(all_counts, gene.length = geneLengths[,2])

norm_rpkm <- rpkm

for(i in 1:3)
{
	norm_rpkm[,seq(i,33+i,3)] <- rpkm[,seq(i,33+i,3)]/rpkm[,i]
}

log_rpkm <- log2(norm_rpkm + 0.00001)
is.na(log_rpkm) <- sapply(log_rpkm, is.nan)

heatmap.2(na.omit(log_rpkm), Colv = NULL, trace = 'none', col = colorRampPalette(c('blue','grey','yellow')), symkey = T)
