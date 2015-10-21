##Split the data into three differente tables
library('gtools')
library('rtracklayer')
library('edgeR')

gff.path <- 'C:/Users/am4613/Documents/am4613/Schizosaccharomyces_pombe.ASM294v2.26.gff3'

load('P:/CDC2_RNA_051015.rda')

sam_rpkm <- list()

sam_rpkm[[1]] <- CDC2
sam_rpkm[[2]] <- CDC2_2
sam_rpkm[[3]] <- CDC2_3
sam_rpkm[[4]] <- CDC2_4

all_counts <- read.delim("C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/all_rev_counts.txt", header= T, strings = F)
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

separate_rpkm <- list()

for(i in 1:3)
{
	separate_rpkm[[i]] <- rpkm[,seq(i,33+i,3)]
}

separate_order <- cbind(separate_rpkm[[1]],separate_rpkm[[2]],separate_rpkm[[3]])

write.table(separate_order, sep = '\t', 'C:/Users/am4613/Desktop/me_rpkm.txt')

merged_rpkm <- list()
for(i in 1:3)
{
	merged_rpkm[[i]] <- merge(separate_rpkm[[i]], sam_rpkm[[i+1]], by.x = 'row.names', by.y = 'row.names')
}


plot(merged_rpkm[[1]][,2], merged_rpkm[[1]][,14], log = 'xy')

