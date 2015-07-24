library("rtracklayer")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")

#gff.path <- c("C:/Users/am4613/Desktop/am4613/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
gff.path <- c("~/annot/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
bam.path <- c("~/rna_seq0206/seq_0206/")
count_yes <- FALSE #Do you want a count table or an rda file
output.name <- "Transcript_1"
setwd("~/rna_seq0206/withLTR")

#Import gff and select certain features
gff <- import.gff3(gff.path)
gff.genes <- gff[which(elementMetadata(gff)[,"type"] == "gene" | elementMetadata(gff)[,"type"] == "rRNA_gene" | elementMetadata(gff)[,"type"] == "snRNA_gene" | elementMetadata(gff)[,"type"] == "ncRNA_gene" | elementMetadata(gff)[,"type"] == "snoRNA_gene" | elementMetadata(gff)[,"type"] == "tRNA_gene" | elementMetadata(gff)[,"type.1"] == "LTR")]

gff.genes[which(elementMetadata(gff.genes)[,"type.1"] == "LTR")]$ID <- paste0("LTR:", gff.genes[which(elementMetadata(gff.genes)[,"type.1"] == "LTR")]$Name)

IDs <- elementMetadata(gff.genes)$ID
a <- strsplit(IDs, "\\:")
vector_with_ids <- sapply(a,"[[",2)

#Import list of bam files
bam.list <- list.files(bam.path, pattern = "bam$", full.names = T)
bam <- BamFileList(bam.list, index = character(), obeyQname = T)

#Obtain counts
se <- summarizeOverlaps(gff.genes, bam, mode = "Union")

if(count_yes){
	#Construct table with IDs and counts for all the samples 
	count_table <- cbind(vector_with_ids, assays(se)$counts)
	count_table <- as.data.frame(count_table)
	write.table(count_table, file = paste0(output.name, ".txt"), sep = "\t")  

}else{
	save(se, file = paste0(bam.path, output.name,".rda"))
}

counts <- assays(se)$counts
counts <- assays(se)$counts[,seq(1,24,2)]

##Calculate RPKMS (HOPEFULLY)
geneLengthsInKB <- width(gff.genes)/1000# Length of exon union per gene in kbp
millionsMapped <- colSums(counts)/1e+06 # Factor for converting to million of mapped reads.
rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
rpkm_table <- cbind(vector_with_ids, rpkm)
rpkm_table <- as.data.frame(rpkm_table)

if(count_yes){
	write.table(rpkm_table, file = paste0(output.name, "rpkm.txt"), sep = "\t")
}else{
	save(rpkm_table, file = paste0(output.name, "rpkm.rda"))
}
	 
