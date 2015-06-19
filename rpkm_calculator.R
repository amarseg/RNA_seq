##RPKM calculator
library("rtracklayer")
library("gplots")

gff.path <- c("Z:/GFF/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
count.path <- c("C:/Users/am4613/Desktop/am4613/count_table.tsv")



gff <- import.gff3(gff.path)
gff.genes <- gff[which(elementMetadata(gff)[,"type"] == "gene" | elementMetadata(gff)[,"type"] == "rRNA_gene" | elementMetadata(gff)[,"type"] == "snRNA_gene" | elementMetadata(gff)[,"type"] == "ncRNA_gene" | elementMetadata(gff)[,"type"] == "snoRNA_gene" | elementMetadata(gff)[,"type"] == "tRNA_gene" | elementMetadata(gff)[,"type.1"] == "LTR")]
gff.cds <- gff[which(elementMetadata(gff)[,"type"] == "CDS")]
gffsub = gff.genes

gffsub[which(elementMetadata(gffsub)[,"type.1"] == "LTR")]$ID <-paste0("LTR:",gffsub[which(elementMetadata(gffsub)[,"type.1"] == "LTR")]$Name) 

counts <- read.delim(count.path, header = T, strings = F)

#returnRPKM <- function(counts, gffsub) {
	geneLengthsInKB <- width(gffsub)/1000# Length of exon union per gene in kbp
	IDs <- elementMetadata(gffsub)$ID
	a <- strsplit(IDs, "\\:")
	vector_with_ids <- sapply(a,"[[",2)
	geneLengthsInKB <- as.data.frame(cbind(vector_with_ids, geneLengthsInKB))
	geneLengthsInKB[,1] <- as.character(geneLengthsInKB[,1])
	geneLengthsInKB[,1] <- substr(geneLengthsInKB[,1], start = 1, stop = nchar(geneLengthsInKB[,1]) - 2)
	all <- merge(x = counts, y = geneLengthsInKB, by.x = row.names(counts), by.y = as.vector(geneLengthsInKB[,1]))
	millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
	rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
	rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
	#return(rpkm)
#}

rpkms <- returnRPKM(counts, gff.genes)
norm.rpkms <- log2(rpkms/rpkms[,1])
is.na(norm.rpkms) <- sapply(norm.rpkms, is.nan)
is.na(norm.rpkms) <- sapply(norm.rpkms, is.infinite)



heatmap.2(as.matrix(rpkms))
