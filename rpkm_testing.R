##Split the data into three differente tables
library('gtools')
library('rtracklayer')
library('edgeR')
library('DESeq2')

gff.path <- 'C:/Users/am4613/Documents/am4613/Schizosaccharomyces_pombe.ASM294v2.26.gff3'

load('P:/CDC2_RNA_051015.rda')

sam_rpkm <- cbind(CDC2_2,CDC2_3,CDC2_4)

all_counts <- read.delim("C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/new_all_counts.txt", header= T, strings = F)
row.names(all_counts) <- all_counts[,1]
all_counts <- all_counts[,2:37]
all_counts <- all_counts[,mixedorder(colnames(all_counts))]

gff <- import.gff3(gff.path)
gff.exons <- gff[which(elementMetadata(gff)[,"type"] == "CDS" | elementMetadata(gff)[,"type.1"] == "LTR" | elementMetadata(gff)[,'type'] == 'ncRNA_gene' | elementMetadata(gff)[,'type'] == 'tRNA_gene' | elementMetadata(gff)[,'type'] == 'rRNA_gene' | elementMetadata(gff)[,'type'] == 'snoRNA_gene' | elementMetadata(gff)[,'type'] == 'snRNA_gene' | elementMetadata(gff)[,'type'] == 'pseudogene')]
gff.exons[which(elementMetadata(gff.exons)[,"type.1"] == "LTR")]$ID <- paste0("LTR:", gff.exons[which(elementMetadata(gff.exons)[,"type.1"] == "LTR")]$Name) 
exons <- gff.exons[which(elementMetadata(gff.exons)[,"type"] == "CDS")]$ID
exons <- substr(exons, 1, nchar(exons) - 6)
gff.exons[which(elementMetadata(gff.exons)[,"type"] == "CDS")]$ID <- exons
meh <- gff.exons[which(elementMetadata(gff.exons)[,"type"] == "rRNA" | elementMetadata(gff.exons)[,"type"] == "tRNA" | elementMetadata(gff.exons)[,"type"] == "snoRNA")]$ID
meh <- substr(meh, 1, nchar(meh) - 2)
gff.exons[which(elementMetadata(gff.exons)[,"type"] == "rRNA" | elementMetadata(gff.exons)[,"type"] == "tRNA" | elementMetadata(gff.exons)[,"type"] == "snoRNA")]$ID <- meh
#unique_exonID <- unique(as.character(elementMetadata(gff.exons)$ID))
a <- strsplit(elementMetadata(gff.exons)$ID, "\\:")
vector_with_IDs <- sapply(a,"[[",2)
geneLengths <- width(gff.exons)# Length of exon union per gene in kbp
geneLengths <- data.frame(ID = vector_with_IDs, length = geneLengths)

geneLengths <- aggregate(geneLengths[,2], FUN = sum, na.rm = T, by = list(geneLengths[,1]))

rpkm <- all_counts
sumas <- colSums(all_counts)

rpkm <- rpkm * 10^9

for(i in 1:ncol(rpkm))
{
	rpkm[,i] <- rpkm[,i]/geneLengths[,2]
	rpkm[,i] <- rpkm[,i]/sumas[i]
}


separate_rpkm <- list()

for(i in 1:3)
{
	separate_rpkm[[i]] <- rpkm[,seq(i,33+i,3)]
}

separate_order <- cbind(separate_rpkm[[1]],separate_rpkm[[2]],separate_rpkm[[3]])

par(mfrow = c(1,3))


exp_design = data.frame(Time = factor(rep(0:11,each = 3)), 
												Replicate = rep(1:3,12), 
												group = rep(1,ncol = all_counts),
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

meh <- sam_rpkm[,grep(colnames(sam_rpkm), pattern = 'Annot', invert = T)]
meh <- meh[,grep(colnames(meh), pattern = 'MM', invert = T)]
meh <- meh[,grep(colnames(meh), pattern = 'common', invert = T)]

separate_counts <- list()
for(i in 1:3)
{
	separate_counts[[i]] <- all_counts[,seq(i,33+i,3)]
}

counts_order <- cbind(separate_counts[[1]],separate_counts[[2]],separate_counts[[3]])