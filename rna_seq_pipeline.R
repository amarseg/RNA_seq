library("rtracklayer")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")

##I PRODUCE DATA THAT DOES NOT CORRELATE VERY WELL WITH SAMS COUNTS

gff.path <- c("C:/Users/am4613/Desktop/am4613/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
#gff.path <- c("~/annot/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
bam.path <- c("~/rna_seq0206/seq_0206/")
count_yes <- FALSE #Do you want a count table or an rda file
output.name <- "Transcript_1"
setwd("~/rna_seq0206/withLTR")

#Import gff and select certain features
gff <- import.gff3(gff.path)
gff.exons <- gff[which(elementMetadata(gff)[,"type"] == "exon" | elementMetadata(gff)[,"type.1"] == "LTR") ]
gff.exons[which(elementMetadata(gff.exons)[,"type.1"] == "LTR")]$ID <- paste0("LTR:", gff.exons[which(elementMetadata(gff.exons)[,"type.1"] == "LTR")]$Name)
exons <- gff.exons[which(elementMetadata(gff.exons)[,"type"] == "exon")]$Parent
exons <- substr(exons, 1, nchar(exons) - 2)
gff.exons[which(elementMetadata(gff.exons)[,"type"] == "exon")]$ID <- exons
#unique_exonID <- unique(as.character(elementMetadata(gff.exons)$ID))
a <- strsplit(elementMetadata(gff.exons)$ID, "\\:")
vector_with_IDs <- sapply(a,"[[",2)

#Import list of bam files
bam.list <- list.files(bam.path, pattern = "bam$", full.names = T)
bam <- BamFileList(bam.list, index = character(), obeyQname = T)

#Obtain counts. if inter.feature is false, reads that are included in several features are counted for 
#each one of them. Otherwise they are discarded
se <- summarizeOverlaps(gff.exons, bam, mode = "Union", inter.feature = F, ignore.strand = F) 

#Find ovelapping features and divide the number of reads by the number of features that span the same region
n_overlap <- countOverlaps(gff.exons, gff.exons, type = "within")
counts <- assays(se)$counts/n_overlap 

#counts <- assays(se)$counts[,seq(1,24,2)]
row.names(counts) <- vector_with_IDs
ag <- aggregate(counts, FUN = sum, by = list(row.names(counts)))

if(count_yes){
	#Construct table with IDs and counts for all the samples 
	write_table(counts, file = paste0(output.name, ".txt"), sep = "\t")  
	
}else{
	save(se, file = paste0(bam.path, output.name,".rda"))
}



##FIX ME
geneLengthsInKB <- width(gff.genes)/1000# Length of exon union per gene in kbp
millionsMapped <- colSums(counts)/1e+06 # Factor for converting to million of mapped reads.
rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
rpkm_table <- cbind(vector_with_ids, rpkm)
rpkm_table <- as.data.frame(rpkm_table)

if(count_yes){
	write_table(rpkm_table, file = paste0(output.name, "rpkm.txt"), sep = "\t")
}else{
	save(rpkm_table, file = paste0(output.name, "rpkm.rda"))
}


