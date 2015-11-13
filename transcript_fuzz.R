library('Mfuzz')
library('gtools')
library('matrixStats')
library('gplots')
library('ggplot2')
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')


load('P:/CDC2_RNA_051015.rda')

sam_rpkm <- cbind(CDC2_2,CDC2_3,CDC2_4)

transcriptomics <- sam_rpkm[,grep(colnames(sam_rpkm), pattern = 'Annot', invert = T)]
transcriptomics <- transcriptomics[,grep(colnames(transcriptomics), pattern = 'MM', invert = T)]
transcriptomics <- transcriptomics[,grep(colnames(transcriptomics), pattern = 'common', invert = T)]

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
n_cl =9
m_fuzzy = mestimate(eset.s)
cl <- mfuzz(eset.s, c = n_cl, m = m_fuzzy)
mfuzz.plot(eset.s, cl = cl, mfrow = c(4,3), min.mem = 0.7, colo = rainbow(n=10))

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
write.table(output_final, 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/tableGO_fuzzy_rna.txt', sep = '\t')
save(f_clusters, file = 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/transcript_clusters.rda')


median_trajectories <- as.data.frame(matrix(nrow = n_cl, ncol = 12, NA))
par(mfrow = c(3,3))
for(i in 1:n_cl)
{
	todo <- exprs(eset.s)[row.names(exprs(eset.s)) %in% f_clusters[[i]],]
	median_trajectories[i,] <- colMedians(as.matrix(todo), na.rm = T)
	plot(y = median_trajectories[i,], x = c(0:11), type = 'l', lwd = 3, xlab = 'Time', ylab = paste0('Cluster ',i))
}


