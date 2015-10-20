library('Mfuzz')
library('gtools')
library('matrixStats')
library('gplots')
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')


transcriptomics <- read.delim('C:/Users/am4613/Desktop/me_rpkm.txt', header = T, strings = F)
transcriptomics[,1:12] <- transcriptomics[,1:12]/transcriptomics[,1]
transcriptomics[,13:24] <- transcriptomics[,13:24]/transcriptomics[,13]
transcriptomics[,25:36] <- transcriptomics[,25:36]/transcriptomics[,25]

transcriptomics <- log2(transcriptomics)
is.na(transcriptomics) <- sapply(transcriptomics, is.nan)
is.na(transcriptomics) <- sapply(transcriptomics, is.infinite)
transcriptomics <- na.omit(transcriptomics)

heatmap.2(as.matrix(transcriptomics), Colv = F, trace = 'none', col = colorRampPalette(c('blue','gray','yellow')))

median_trans <- matrix(ncol = 12, nrow = nrow(transcriptomics), NA) 
median_trans <- as.data.frame(median_trans)
row.names(median_trans) <- row.names(transcriptomics)

for(i in 1:12)
{
	median_trans[,i] <- rowMedians(cbind(transcriptomics[,i], transcriptomics[,i+12], transcriptomics[,i+24]))
}

eset <- ExpressionSet(as.matrix(median_trans))

eset.r <- filter.NA(eset,0.5) #Threshold for discarding missing data 
eset.f <- fill.NA(eset.r, mode = 'mean') #Fill the gaps with the mean of the other values

tmp <- filter.std(eset.f, min.std = 0)

eset.s <- standardise(tmp)
n_cl = 9
m_fuzzy = 1.25
cl <- mfuzz(eset.s, c = n_cl, m = m_fuzzy)
mfuzz.plot(eset.s, cl = cl, mfrow = c(4,4), min.mem = 0.7, colo = rainbow(n=10))

##Code to extract the members of the different clusters and do GO analysis on them 

f_clusters <- list()

for(i in 1:n_cl)
{
	f_clusters[[i]] <- names(cl$cluster[cl$cluster == i])
}

##GO analysis of the different clusters
p = 0.05

results_go <- list()

for(i in 1:n_cl)
{
	temp <- GOanalysis(li = f_clusters[[i]], go = GOtable, all = 5123, sort = T)
	temp[,6] <- i 
	results_go[[i]] <- temp[temp[,2] <= p,]
}

output_final <- data.frame(Reduce(rbind,results_go))
write.table(output_final, 'C:/Users/am4613/Desktop/tableGO_fuzzy_rna.txt', sep = '\t')


