library("rtracklayer")
library("GenomicRanges")
library("Rsamtools")


#gff.path <- c("C:/Users/am4613/Desktop/am4613/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
gff.path <- c("Z:/GFF/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
bam.path <- c("C:/Users/am4613/Desktop/simpleRNAtest/merged_bam/")

#Import gff and select certain features
gff <- import.gff3(gff.path)
gff.genes <- gff[which(elementMetadata(gff)[,"type"] == "gene"),]

#Import list of bam files
bam.list <- list.files(bam.path, pattern = "bam$", full.names = T)
bam <- BamFileList(bam.list, index = character(), obeyQname = T)

#Obtain counts
se <- summarizeOverlaps(gff.genes, bam, mode = "Union")

#Construct table with IDs and counts for all the samples 
count_table <- cbind(substr(rowRanges(se)$ID, start = 6, stop = nchar(rowRanges(se)$ID)), assays(se)$counts)
count_table <- as.data.frame(count_table)
count_table[,2:8] <- as.numeric(count_table[,2:8])
count_table[,1] <- as.character(count_table[,1])
