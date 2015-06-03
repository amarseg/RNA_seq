library("rtracklayer")
library("GenomicRanges")
library("GenomicAlignments")
library("GenomicFeatures")
library("easyRNASeq")


#gff.path <- c("C:/Users/am4613/Desktop/am4613/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
gff.path <- c("Z:/GFF/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
bam.path <- c("C:/Users/am4613/Desktop/simpleRNAtest/TP11_dna_sorted.bam")

#gff <- import.gff3(gff.path)
tbx <- makeTxDbFromGFF(gff.path, format = "gff3")
exbygene <- exonsBy(tbx, "gene")

bam <- readGAlignments(bam.path)
indexBam(bam.path)

se <- summarizeOverlaps(exbygene, bam, mode = "Union")


