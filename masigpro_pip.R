###Masigpro 

library(maSigPro)
library(gtools)
library(DESeq2)
library(edgeR)

low_exp_cutoff <- 10 

all_counts <- read.delim("C:/Users/am4613/Desktop/all_rev_counts.txt", header= T, strings = F)
row.names(all_counts) <- all_counts$ID
all_counts <- all_counts[,2:37]
all_counts <- all_counts[,mixedorder(colnames(all_counts))]


##Remove genes with little reads in all samples
all_counts <- all_counts[rowSums(all_counts) > low_exp_cutoff,]

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


exp_design$Time <- rep(0:11,each = 3)
design <- make.design.matrix(exp_design, degree = 3, time.col = 1, repl.col = 2, group.col = 3)

fit <- p.vector(norm_counts, design, counts = T)
tstep <- T.fit(fit)
sigs <- get.siggenes(tstep, vars = 'each')
see.genes(sigs$sig.genes[[1]], edesign = exp_design , 
					group.cols = 3, show.fit = T, summary.mode = 'median',
					color.mode = 'gray')
see.genes(norm_counts, edesign = exp_design, group.cols = 3)

top_norm <- norm_counts
######All normalised to time point zero
for(i in 1:3)
{
	top_norm[,seq(i,33+i, 3)] = norm_counts[,seq(i,i + 11, 3)]/norm_counts[,i]
}

is.na(top_norm) <- sapply(top_norm, is.nan)
is.na(top_norm) <- sapply(top_norm, is.infinite)
top_norm[is.na(top_norm)] <- 0

fit <- p.vector(top_norm, design, counts = T)
tstep_tops <- T.fit(fit_tops)
sigs_tops <- get.siggenes(tstep_tops, vars = 'each')
see.genes(sigs_tops$sig.genes[[1]], edesign = exp_design , 
					group.cols = 3, show.fit = T, summary.mode = 'median',
					color.mode = 'gray')
