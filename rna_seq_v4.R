

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

rna.seq.quant <- function(gff.path, bam.path, count_yes = T)
{
	library("rtracklayer")
	library("GenomicRanges")
	library("Rsamtools")
	library("GenomicAlignments")
	
	#Import gff and select certain features
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
		df.bam[,5] <- df.bam[,3] + df.bam[,4] #I need the end of the sequence (start + length of the sequence) to do the overlaps
		df.bam[,1] <- as.character(df.bam[,1])
		if(df.bam$rname[1] == "chr1") #For Sam bam files 
		{
			df.bam[which(df.bam$rname == "chr1"),1] <- "I"
			df.bam[which(df.bam$rname == "chr2"),1] <- "II"
			df.bam[which(df.bam$rname == "chr3"),1] <- "III"
			df.bam[which(df.bam$rname == "chr4"),1] <- "MT"
			df.bam[which(df.bam$rname == "chr5"),1] <- "MTR"
			df.bam[which(df.bam$rname == "chr6"),1] <- "AB325691"
		}
		gr.bam <- makeGRangesFromDataFrame(df.bam, start.field = "pos", end.field = "V5", strand.field = "strand", seqnames.field = "rname")
		overlaps <- findOverlaps(gr.bam, gff.exons, select = "arbitrary", type = "any")
		count <- table(overlaps)
		count <- as.data.frame(count)
		count <- apply(count, MARGIN = c(1,2), as.numeric)
		count_table[count[,1],i+1] <- count[,2]
		
	}
	
	agg <- aggregate(count_table[,2:ncol(count_table)], FUN = sum, na.rm = T, by = list(count_table[,1]) ) 
	
	
	if(count_yes){
		#Construct table with IDs and counts for all the samples 0
		colnames(agg) <- c("ID",basename(bam.list))
		return(agg)  
		
	}else{
		
		##FIX ME
		geneLengths <- width(gff.exons)# Length of exon union per gene in kbp
		geneLengths <- data.frame(ID = vector_with_IDs, length = geneLengths)
		geneLengths <- aggregate(geneLengths[,2], FUN = sum, na.rm = T, by = list(geneLengths[,1]))
		millionsMapped <- sum(agg[,2:ncol(agg)])/1e+06 # Factor for converting to million of mapped reads.
		rpm <- agg[,2:ncol(agg)]/millionsMapped # RPK: reads per kilobase of exon model.
		rpkm <- rpm/geneLengths[,2] # RPKM: reads per kilobase of exon model per million mapped reads.
		rpkm <- rpkm*10^9
		rpkm <- cbind(geneLengths[,1], rpkm)
		return(rpkm)
		
	}
	
}


