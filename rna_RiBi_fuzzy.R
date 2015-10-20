##

library('gplots')
library('gtools')
library('Mfuzz')
library('matrixStats')

biogenesis <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/GO_ribosome_biogenesis', header = T, strings = F)
rps <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/rproteins.txt', header = T, strings = F)

load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')
load('P:/CDC2_RNA_051015.rda')

sam_rpkm <- list()

sam_rpkm[[1]] <- CDC2
sam_rpkm[[2]] <- CDC2_2
sam_rpkm[[3]] <- CDC2_3
sam_rpkm[[4]] <- CDC2_4

for(i in 2:4)
{
	sam_rpkm[[i]][,1:12] <- log2(sam_rpkm[[i]][,1:12]/sam_rpkm[[i]][,1])
	is.na(sam_rpkm[[i]]) <- sapply(sam_rpkm[[i]], is.nan)
	is.na(sam_rpkm[[i]]) <- sapply(sam_rpkm[[i]], is.infinite)
	sam_rpkm[[i]][is.na(sam_rpkm[[i]])] <- 0
	
}

norm_sam <- cbind(sam_rpkm[[2]][,1:12],sam_rpkm[[3]][,1:12],sam_rpkm[[4]][,1:12])
median_sam <- matrix(ncol = 12, nrow = nrow(norm_sam), NA)
median_sam <- as.data.frame(median_sam)
colnames(median_sam) <- as.character(c(0:11))
row.names(median_sam) <- row.names(norm_sam)

for(i in 1:12)
{
	median_sam[,i] <- rowMedians(cbind(norm_sam[,i],norm_sam[,i+12],norm_sam[,i+24]), na.rm = T)
}


prot_bio <- subset(median_sam,row.names(median_sam) %in% biogenesis[,1])
prot_rp <- subset(median_sam,row.names(median_sam) %in% rps[,3])

bioset <- ExpressionSet(as.matrix(prot_bio))
rpset <- ExpressionSet(as.matrix(prot_rp))

bioset.r <- filter.NA(bioset,0.5) #Threshold for discarding missing data 
bioset.f <- fill.NA(bioset.r, mode = 'mean') #Fill the gaps with the mean of the other values
rpset.r <- filter.NA(rpset,0.5) #Threshold for discarding missing data 
rpset.f <- fill.NA(rpset.r, mode = 'mean') #Fill the gaps with the mean of the other values

bio_tmp <- filter.std(bioset.f, min.std = 0)
rp_tmp <- filter.std(rpset.f, min.std = 0)

bioset.s <- standardise(bio_tmp)
rpset.s <- standardise(rp_tmp)
n_cl = 3
m_fuzzy = 1.25
biocl <- mfuzz(bioset.s, c = n_cl, m = m_fuzzy)
mfuzz.plot(bioset.s, cl = biocl, mfrow = c(3,3), min.mem = 0.7, colo = rainbow(n=10))

rpcl <- mfuzz(rpset.s, c = n_cl, m = m_fuzzy)
mfuzz.plot(rpset.s, cl = rpcl, mfrow = c(3,3), min.mem = 0.7, colo = rainbow(n=10))

##Code to extract the members of the different clusters and do GO analysis on them 

bio_clusters <- list()

for(i in 1:n_cl)
{
	bio_clusters[[i]] <- names(biocl$cluster[biocl$cluster == i])
}


rp_clusters <- list()

for(i in 1:n_cl)
{
	rp_clusters[[i]] <- names(rpcl$cluster[rpcl$cluster == i])
}


##GO analysis of the different clusters
p = 0.05

results_go <- list()

for(i in 1:n_cl)
{
	temp <- GOanalysis(li = bio_clusters[[i]], go = GOtable, all = 5123, sort = T)
	temp[,6] <- i 
	results_go[[i]] <- temp[temp[,2] <= p,]
}

output_final <- data.frame(Reduce(rbind,results_go))
write.table(output_final, 'C:/Users/am4613/Desktop/tableGO_RiBI_RNA.txt', sep = '\t')
