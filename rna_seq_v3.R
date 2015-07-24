library("rtracklayer")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")

##Copied from Rsamtools manual, turns the lists of lists into a single list
.unlist <- function(x)
{
	x1 <- x[[1L]]
	if(is.factor(x1)){
		structure(unlist(x), class = "factor", levels = levels(x1))
	}else{
		do.call(c,x)
	}
}

gff.path <- c("C:/Users/am4613/Desktop/am4613/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
#gff.path <- c("~/annot/Schizosaccharomyces_pombe.ASM294v2.26.gff3")
bam.path <- c("C:/Users/am4613/Desktop/run_1802/")
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
what <- c("rname","strand","pos","qwidth")
param <- ScanBamParam(what = what)

count_table <- matrix(ncol = length(bam.list) + 1, nrow = length(vector_with_IDs), NA)
count_table <- as.data.frame(count_table)
count_table[,1] <- vector_with_IDs

for(i in seq(1,length(bam.list)))
{
	bam <- scanBam(bam.list[i], param = param)
	bam <- unname(bam)
	elts <- setNames(bamWhat(param), bamWhat(param))
	lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))
	rm(bam)
	rm(elts)
	df.bam <- do.call("data.frame", lst)
	df.bam[,5] <- df.bam[,3] + df.bam[,4] 
	df.bam[,1] <- as.character(df.bam[,1])
	df.bam[which(df.bam$rname == "chr1"),1] <- "I"
	df.bam[which(df.bam$rname == "chr2"),1] <- "II"
	df.bam[which(df.bam$rname == "chr3"),1] <- "III"
	df.bam[which(df.bam$rname == "chr4"),1] <- "MT"
	df.bam[which(df.bam$rname == "chr5"),1] <- "MTR"
	df.bam[which(df.bam$rname == "chr6"),1] <- "AB325691"
	gr.bam <- makeGRangesFromDataFrame(df.bam, start.field = "pos", end.field = "V5", strand.field = "strand", seqnames.field = "rname")
	overlaps <- findOverlaps(gr.bam, gff.exons, select = "arbitrary", type = "any")
	count <- table(overlaps)
	count <- as.data.frame(count)
	count <- apply(count, MARGIN = c(1,2), as.numeric)
	count_table[count[,1],i+1] <- count[,2]
	rm(bam)
}


#counts <- assays(se)$counts[,seq(1,24,2)]
agg <- aggregate(count_table[,2:13], FUN = sum, na.rm = T, by = list(count_table[,1]) ) 


if(count_yes){
	#Construct table with IDs and counts for all the samples 0
	write_table(counts, file = paste0(output.name, ".txt"), sep = "\t")  
	
}else{
	save(se, file = paste0(bam.path, output.name,".rda"))
}

x1

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




