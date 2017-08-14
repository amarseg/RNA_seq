rm(list = ls())
library('gplots')
library('pheatmap')
load('C:/Users/am4613/OneDrive/Summaries_as_timecourses/analysis/CDC2_RNA_051015.rda')
source('C:/Users/am4613/Documents/GitHub/Misc/gene_plotter.R')
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
norm_counts <- read.delim('C:/Users/am4613/OneDrive/Summaries_as_timecourses/DESeq_norm_counts.txt', header = T, strings = F)
norm_counts <- reorder_proteomics(norm_counts)
norm_counts <- normalise_rna(norm_counts)

norm_prot <- read.delim('C:/Users/am4613/OneDrive/Summaries_as_timecourses/analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)
norm_prot <- normalise_ProtDataset(norm_prot, what = 'heatmap')
norm_prot <- norm_prot[,7:ncol(norm_prot)]
norm_prot <- reorder_proteomics(norm_prot)

###CDC2_as strain plotting non-conding RNAs
sam_rpkm <- list()

sam_rpkm[[1]] <- CDC2
sam_rpkm[[2]] <- CDC2_2
sam_rpkm[[3]] <- CDC2_3
sam_rpkm[[4]] <- CDC2_4

blue_yellow_palette <- colorRampPalette(c('Blue','gray','yellow'))

for(i in 2:4)
{
	sam_rpkm[[i]][,1:12] <- log2(sam_rpkm[[i]][,1:12]/sam_rpkm[[i]][,1])
	is.na(sam_rpkm[[i]]) <- sapply(sam_rpkm[[i]], is.nan)
	is.na(sam_rpkm[[i]]) <- sapply(sam_rpkm[[i]], is.infinite)
	sam_rpkm[[i]][is.na(sam_rpkm[[i]])] <- 0
	
}

norm_sam <- cbind(sam_rpkm[[2]][,1:12],sam_rpkm[[3]][,1:12],sam_rpkm[[4]][,1:12])

nc_rna <- norm_sam[grep(row.names(norm_sam), pattern = 'NCRNA'),]

col = colorRampPalette(c('blue','gray','yellow'))
pheatmap(as.matrix(nc_rna), cluster_cols = F, clustering_method = 'ward.D2')

nc_rna_names <- row.names(nc_rna)

gene_plotter(nc_rna_names, what = 'RNA')

###Thermosensitive strain - plotting non-coding RNAs

cdc2_ts <- read.delim('C:/Users/am4613/Documents/am4613/rna_seq0804/RPKM.tsv', header = T, strings = F)

norm_ts <- log2(cdc2_ts/cdc2_ts[,1])
is.na(norm_ts) <- sapply(norm_ts, is.nan)
is.na(norm_ts) <- sapply(norm_ts, is.infinite)
norm_ts[is.na(norm_ts)] <- 0

nc_ts <- norm_ts[grep(row.names(norm_ts), pattern = 'NCRNA'),]

pheatmap(nc_ts, cluster_cols = F, clustering_method = 'ward.D2')

##Do different types of RNAs cluster together?
nc_types <- read.delim('C:/Users/am4613/OneDrive/Summaries_as_timecourses/analysis/nc_rna_type.txt', header = T, strings = F)
##Remove test data (new unnanotated nc RNAs)
clean_nc_types <- nc_types[grep(nc_types$gene.ID, pattern = 'test', invert = T),]
clean_nc_types <- data.frame(clean_nc_types$gene.ID, clean_nc_types$XUT, clean_nc_types$CUT, clean_nc_types$DUT)
clean_nc_types$name <- NA
clean_nc_types[which(clean_nc_types$clean_nc_types.XUT == 1),]$name <- 'XUT' 
clean_nc_types[which(clean_nc_types$clean_nc_types.CUT == 1),]$name <- 'CUT' 
clean_nc_types[which(clean_nc_types$clean_nc_types.DUT == 1),]$name <- 'RUT' 

##How many of each type? 
barplot(table(clean_nc_types$name))
##Create data frame for labelling the rows
nc_labels <- data.frame(row.names = clean_nc_types$clean_nc_types.gene.ID, names = clean_nc_types$name)

gene_plotter(nc_rna_names, what = 'RNA', rowlabs = nc_labels)
pheatmap(nc_ts, cluster_cols = F, clustering_method = 'ward.D2', annotation_row = nc_labels, color = blue_yellow_palette(100))

gene_plotter(nc_rna_names, what = 'RNA', rowlabs = nc_labels)


#They don't seem to cluster together in any of the datasets...(FUCK ME)
##See if nc_RNAs and their gene pairs anticorrelate at the transcript level
pairs <- read.delim('C:/Users/am4613/OneDrive/Summaries_as_timecourses/analysis/nc_rna_pairs.txt', header = T, strings = F, skip = 1)

t <- merge(pairs, norm_counts, by.x = 'g1', by.y = 'row.names', all.x = T)
t2 <- merge(t, norm_counts, by.x = 'g2', by.y = 'row.names', all.x = T)

##only non_codingi
nc_t2 <- t2[which(t2$biotype.g1 == 'ncRNA' | t2$biotype.g2 == 'ncRNA'),]
nc_t2$orientation <- as.factor(nc_t2$orientation)
nc_t2$overlap.at.5..end <- as.factor(nc_t2$overlap.at.5..end)
nc_t2$overlap.at.3..end <- as.factor(nc_t2$overlap.at.3..end)

plot(nc_t2[,52], nc_t2[,88], log = 'xy', xlab = 'Protein_coding exp', ylab = 'non_coding partner', 
		 col = nc_t2$overlap.at.5..end)

pheatmap(nc_t2[,17:ncol(nc_t2)], cluster_cols = F, clustering_method = 'ward.D2' )

div_t2 <- nc_t2[which(nc_t2$orientation == 'divergent'),]
pheatmap(div_t2[,17:ncol(div_t2)], cluster_cols = F, clustering_method = 'ward.D2')

con_t2 <- nc_t2[which(nc_t2$orientation == 'convergent'),]
pheatmap(con_t2[,17:ncol(con_t2)], cluster_cols = F, clustering_method = 'ward.D2')

##Same for protein level

pt <- merge(pairs, norm_prot, by.x = 'g1', by.y = 'row.names', all.x = T)
pt2 <- merge(pt, norm_counts, by.x = 'g2', by.y = 'row.names', all.x = T)

##only non_coding
pt_nc_t2 <- pt2[which(t2$biotype.g1 == 'ncRNA' | pt2$biotype.g2 == 'ncRNA'),]
pt_nc_t2$orientation <- as.factor(pt_nc_t2$orientation)
pt_nc_t2$overlap.at.5..end <- as.factor(pt_nc_t2$overlap.at.5..end)
pt_nc_t2$overlap.at.3..end <- as.factor(pt_nc_t2$overlap.at.3..end)

plot(pt_nc_t2[,52], pt_nc_t2[,88], log = 'xy', xlab = 'Protein_coding exp', ylab = 'non_coding partner', 
		 col = pt_nc_t2$overlap.at.5..end)

pheatmap(na.omit(pt_nc_t2[,17:ncol(pt_nc_t2)]), cluster_cols = F, clustering_method = 'ward.D2' )

pt_div_t2 <- pt_nc_t2[which(pt_nc_t2$orientation == 'divergent'),]
pheatmap(na.omit(pt_div_t2[,17:ncol(pt_div_t2)]), cluster_cols = F, clustering_method = 'ward.D2')

pt_con_t2 <- pt_nc_t2[which(pt_nc_t2$orientation == 'convergent'),]
pheatmap(na.omit(pt_con_t2[,17:ncol(pt_con_t2)]), cluster_cols = F, clustering_method = 'ward.D2')


##Looking only at differentially expressed non-coding RNAs
de1 <- read.delim('C:/Users/am4613/OneDrive/Summaries_as_timecourses/both_clustering/Cluster_1.txt', header =T, strings = F)
de2 <- read.delim('C:/Users/am4613/OneDrive/Summaries_as_timecourses/both_clustering/Cluster_2.txt', header =T, strings = F)

de_ncrna <- rbind(de1, de2)
test <- nc_t2[which(nc_t2$g2 %in% de_ncrna[,1]),]
pheatmap(na.omit(test[,17:ncol(test)]), cluster_cols = F, clustering_method = 'ward.D2',
				 color = colorRampPalette(c('blue','gray','yellow'))(100))
plot(test[,52], test[,88], xlab = 'Protein_coding exp', ylab = 'non_coding partner', 
		 col = test$orientation)

cor_prot_nc <- c(1:36)
for(i in 17:ncol(test))
{
	cor_prot_nc[i-16] <- cor(test[,i], test[,i+36], use = 'complete.obs')
	
}
plot(cor_prot_nc, type = 'l', lwd = 3)
cor_back <- c(1:36)
for(i in 17:ncol(test))
{
	cor_back[i-16] <- cor(pt2[,i], pt2[,i+36], use = 'complete.obs')
	
}

lines(cor_back, col = 'red', lwd = 3)
##Lets do the same analyisis but on a gene to gene basis
test <- na.omit(test)
test2 <- t(test)
cor_genes <- c(1:ncol(test2))
for(i in 1:ncol(test2))
{
	cor_genes[i] <- cor(as.numeric(test2[17:52,i]), as.numeric(test2[53:88,i]), use = 'complete.obs')
}
pt2_rev <- na.omit(nc_t2)
pt2_rev <- t(pt2_rev)
cor_genes_back <- c(1:ncol(pt2_rev))
for(i in 1:ncol(pt2_rev))
{
	cor_genes_back[i] <- cor(as.numeric(pt2_rev [17:52,i]), as.numeric(pt2_rev [53:88,i]), use = 'complete.obs')
}

hist(cor_genes, lwd = 3, freq = F)
lines(density(cor_genes_back), col = 'red', lwd = 3)

cual <- which(cor_genes < 0)
m <- test[cual,]
write.table(m, 'cool_non_coding_genes.txt', sep = '\t')

###I have a very nice table with types of non-coding RNAs

nc_type <- read.delim('Y:/Pers - Amalia/chip_ncRNA/sophie_lnc_rna.txt', header = T,strings = F) ##Only genes with a pair
danny_table <- read.delim('Y:/Pers - Amalia/chip_ncRNA/ncRNA_Overlap_current_annotation.txt', header = T,strings = F)
intergenic <- danny_table[which(rowSums(danny_table[,9:13] == '')==5),]

intersect(intergenic, nc_type$non_coding)

slim_nc_type <- nc_type[,2:4]
##remove test nc_rnas
slim_nc_type <- slim_nc_type[grep(slim_nc_type$non_coding, pattern = 'test', invert = T),]
type_annot <- data.frame(row.names = nc_rna_names)
type_annot[row.names(type_annot) %in% intergenic$Gene.ID,1] <- 'Intergenic'
type_annot[is.na(type_annot[,1]),] <- 'Genic'

gene_plotter(nc_rna_names, what = 'RNA', rowlabs = type_annot)


##I'm going to do the same analysis of nc_pairs and coding pairs but using sophies table
soph_t <- merge(slim_nc_type, norm_counts, by.x = 'non_coding', by.y = 'row.names', all.x = T)
soph_t2 <- merge(soph_t, norm_counts, by.x = 'coding', by.y = 'row.names', all.x = T)

pheatmap(na.omit(soph_t2[,4:75]), cluster_cols = F, clustering_method = 'ward.D2')

soph_pt <- merge(slim_nc_type, norm_prot, by.x = 'coding', by.y = 'row.names')
soph_pt2 <- merge(soph_pt, norm_counts, by.x = 'non_coding', by.y = 'row.names')

pheatmap(soph_pt2[,4:75], cluster_cols = F, clustering_method = 'ward.D2')

de_soph_pt2 <- soph_pt2[which(de_ncrna[,1] %in% soph_pt2$non_coding),]
pheatmap(na.omit(de_soph_pt2[,4:75]), cluster_cols = F, clustering_method = 'ward.D2', 
				 color = colorRampPalette(c('blue','gray','yellow'))(100))

de_soph_t2 <- soph_t2[which(de_ncrna[,1] %in% soph_t2$non_coding),]
pheatmap(na.omit(de_soph_t2[,4:75]), cluster_cols = F, clustering_method = 'ward.D2', 
				 color = colorRampPalette(c('blue','gray','yellow'))(100))

