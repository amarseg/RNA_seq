library("rtracklayer")
library("GenomicRanges")
library("GenomicAlignments")
library("GenomicFeatures")
library("easyRNASeq")


#gff.path <- c("C:/Users/am4613/Desktop/am4613/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
gff.path <- c("Z:/GFF/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
bam.path <- c("C:/Users/am4613/Desktop/simpleRNAtest/")

#Remove repeat_regions
gff <- import.gff3(gff.path)
gff.genes <- gff[which(elementMetadata(gff)[,"type"] == "gene"),]
annotParam <- AnnotParam(datasource = GRangesList(gff.genes))

#annotParam <- AnnotParam(datasource = gff.path, type = "gff3")
rnaSeqParam <- RnaSeqParam(annotParam = annotParam, countBy = "genes")
bamFiles <- getBamFileList(list.files(bam.path, full.names = T))

count.table <- simpleRNASeq(bamFiles = bamFiles, param = rnaSeqParam, verbose = T)
			
