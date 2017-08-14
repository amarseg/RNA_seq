###Masigpro 
##I'm adding a comparison between early and late time points?
library(maSigPro)
library(gtools)
library(DESeq2)
library(edgeR)
library(clusterProfiler)
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
source('C:/Users/am4613/Documents/GitHub/Misc/gene_plotter.R')
library(gplots)
library(dendextend)
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')
source('C:/Users/am4613/Documents/GitHub/cluster_plotter.R')
setwd('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/early_late_de/')

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

low_exp_cutoff <- 10 

rpkm = F

all_counts <- read.delim("C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/analysis/all_rev_counts.txt", header= T, strings = F)
row.names(all_counts) <- all_counts$ID
all_counts <- all_counts[,2:37]
all_counts <- all_counts[,mixedorder(colnames(all_counts))]
all_counts <- all_counts[rowSums(all_counts == 0) !=  36 ,]


##Remove genes with little reads in all samples
all_counts <- all_counts[rowSums(all_counts) > low_exp_cutoff,]

exp_design = data.frame(Time = factor(rep(0:11,each = 3)), 
												Replicate = rep(1:3,12), 
												group = rep(1,36),
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
write.table(norm_counts,'deseq_norm_counts.txt', sep = '\t')
early <- norm_counts[,1:21]
early <- early + 0.001

late <- norm_counts[,c(22:36,1:3)]

exp_design$Time <- rep(0:11, each = 3)

design_early <- make.design.matrix(exp_design[1:21,], degree = 3, time.col = 1, repl.col = 2, group.col = 3)
design_late <- make.design.matrix(exp_design[c(22:36,1:3),], degree = 3, time.col = 1, repl.col = 2, group.col = 3)

fit_early <- p.vector(early, design_early, counts = T)
fit_late <- p.vector(late, design_late, counts = T)

tstep_early <- T.fit(fit_early)
tstep_late <- T.fit(fit_late)

sigs_early <- get.siggenes(tstep_early, vars = 'all')
sigs_late <- get.siggenes(tstep_late, vars = 'all')


###Plotting of DE genes
masigpro_sig_early <- data.frame(ID = sigs_early$summary, pval = sigs_early$sig.genes$sig.pvalues$'p-value', beta0 =sigs_early$sig.genes$coefficients$beta0, beta1 = sigs_early$sig.genes$coefficients$betaTime, beta2 = sigs_early$sig.genes$coefficients$betaTime2,
													 beta3 = sigs_early$sig.genes$coefficients$betaTime3)

masigpro_sig_late <- data.frame(ID = sigs_late$summary, pval = sigs_late$sig.genes$sig.pvalues$'p-value', beta0 =sigs_late$sig.genes$coefficients$beta0, beta1 = sigs_late$sig.genes$coefficients$betaTime, beta2 = sigs_late$sig.genes$coefficients$betaTime2,
													 beta3 = sigs_late$sig.genes$coefficients$betaTime3)

write.table(masigpro_sig_early, sep = '\t', 'masig_early.txt')
write.table(masigpro_sig_late, sep = '\t', 'masig_late.txt')


betr_early <- read.delim('betr_early.txt', header = T, strings = F)
betr_late <- read.delim('betr_late.txt', header = T, strings = F)
masigpro_sig_early <- read.delim('masig_early.txt', header = T, strings = F)
masigpro_sig_late <- read.delim('masig_late.txt', header = T, strings = F)


de_early <- intersect(row.names(betr_early), masigpro_sig_early$ID)
de_late <- intersect(row.names(betr_early), masigpro_sig_late$ID)

write.table(de_early, sep = '\t', 'both_early.txt')
write.table(de_late, sep = '\t', 'both_late.txt')

pdf('de_early.pdf')
gene_plotter(de_early, what = 'RNA')
dev.off()

pdf('de_late.pdf')
obj<- gene_plotter(de_late, what = 'RNA')
dev.off()

dend <- obj$tree_row
clusters <- cutree(dend, k = 2)


re_cl <- list()
pdf(file  = 'both_reclustering.pdf')
for(i in 1:max(clusters))
{
	cl <- names(clusters[clusters == i])
	re_cl[[i]] <- cl
	write.table(cl, file = paste0('redown_Cluster_',i,'.txt'), sep = '\t')
	
	#go_results <- GOanalysis(cl, GOtable, all = 5123)
	#go_results <- go_results[go_results[,2] > pval, ]
	#write.table(go_results, file = paste0('both_clustering/DESeq_both/down_GO_cluster',i,'.txt'), sep = '\t')
	
	gene_plotter(cl, what = 'RNA', norm = 'DESeq')
}
dev.off()

up <- enrichKEGG(as.data.frame(re_cl[1])[,1], organism = 'spo')
enrichMap(up)

down <- enrichKEGG(as.data.frame(re_cl[2])[,1], organism = 'spo')
enrichMap(down)

meh$variable <- factor(meh$variable, levels = c('polymerase','ribosome'))

p <- ggplot(data = meh, aes(y = value, x=Time, group = variable, color = variable))
p +geom_line(size = 2, alpha = 0.60) + theme_bw() 

#########################
##Here I'm doing the analysis dividing the time course in two groups, will that work?
##############################

exp_design = data.frame(Time = rep(0:11,each = 3), 
												Replicate = rep(1:3,12), 
												group = c(rep(1,21),rep(2,15)),
												row.names = colnames(all_counts))

design_all <- make.design.matrix(exp_design, degree = 3, time.col = 1, repl.col = 2, group.col = 3)
fit_all <- p.vector(norm_counts, design_all, counts = T)
tstep_all <- T.fit(fit_all)
sigs_all <- get.siggenes(tstep_all, vars = 'groups')

masigpro_sig_all <- data.frame(ID = sigs_all$summary, pval = sigs_all$sig.genes$group$sig.pvalues$'p-value', beta0 =sigs_all$sig.genes$group$coefficients$beta0, beta1 = sigs_all$sig.genes$group$coefficients$betaTime, beta2 = sigs_all$sig.genes$group$coefficients$betaTime2,
																beta3 = sigs_all$sig.genes$group$coefficients$betaTime3)

ids_sigs <- masigpro_sig_all[,1]
t <- enrichKEGG(ids_sigs, organism = 'spo')
enrichMap(t)

all_tree <- gene_plotter(ids = ids_sigs, what = 'RNA')
dend <- all_tree$tree_row
clusters <- cutree(dend, k = 9)
gene_clusters <- split(clusters, clusters)

###venn diagram of late + tp0 and late vs. early

n <- length(intersect(ids_sigs, sigs_late$summary))
grid.newpage()
draw.pairwise.venn(area1 = length(ids_sigs), area2 = length(sigs_late$summary), cross.area = n,
									 category = c('Late vs. Early','Late + tp0'))

int_gene_list <- intersect(ids_sigs, sigs_late$summary)
t2 <- enrichKEGG(int_gene_list, organism = 'spo')
enrichMap(t2)

#############################
#I do all the things, but with the lists in one of the two choices (which one is better?)
#########################

all_data <- read.delim('../mmc1_from_cell_paper.txt', header = T, strings = F, skip = 6)

up <- names(gene_clusters[[1]])
kegg_up <- enrichKEGG(up, organism = 'spo')
enrichMap(kegg_up)

down <- names(gene_clusters[[2]])
kegg_down <- enrichKEGG(down, organism = 'spo')
enrichMap(kegg_down)

all_data$scale <- NA

all_data[all_data$Systematic.name %in% up,]$scale <- 'Gene Up'
all_data[all_data$Systematic.name %in% down,]$scale <- 'Gene down'
all_data[is.na(all_data$scale),]$scale <- 'Scaling'

p <- ggplot(all_data, aes(x = scale, y= log2(MM.mRNA.cpc),color = scale, fill = scale))
print(p + geom_boxplot(notch = T, alpha = 0.5)+ theme_bw() +scale_fill_manual(values=cbPalette)+ scale_colour_manual(values=cbPalette)) 

p <- ggplot(all_data, aes(x = scale, y= log2(as.numeric(MM.protein.cpc)),color = scale, fill = scale))
print(p + geom_boxplot(notch = T, alpha = 0.5)+ theme_bw() +scale_fill_manual(values=cbPalette)+ scale_colour_manual(values=cbPalette)) 

p <- ggplot(all_data, aes(x = scale, y= exonic.length,color = scale, fill = scale))
print(p + geom_boxplot(notch = T, alpha = 0.5)+ theme_bw() +scale_fill_manual(values=cbPalette)+ scale_colour_manual(values=cbPalette)) 

###Merge with half-life data

mrna_half <- read.delim('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/analysis/rna_half_lives.txt', header = T,strings = F, skip = 6)
all_data_half <- merge(all_data, mrna_half, by.x = 'Systematic.name', by.y = 'Gene.name', all.x =T)

p <- ggplot(all_data_half, aes(x = scale, y= Average.half.life,color = scale, fill = scale))
print(p + geom_boxplot(notch = T, alpha = 0.5) + theme_bw() +scale_fill_manual(values=cbPalette)+ scale_colour_manual(values=cbPalette)) 

prot_half <- read.delim('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/analysis/proteins_half_lives.txt', header = T, strings = F)
prot_half <- prot_half[-which(prot_half$ENSG == ''),]

all_data_prot <- merge(all_data, prot_half, by.x = 'Systematic.name', by.y = 'ENSG', all.x = T)
p <- ggplot(all_data_prot, aes(x = scale, y= as.numeric(t1.2..hours.), color= scale, fill = scale),)
print(p + geom_boxplot(notch = T, alpha = 0.5) + theme_bw() +scale_fill_manual(values=cbPalette)+ scale_colour_manual(values=cbPalette) + coord_cartesian(ylim = c(0, 25)))


#############
#Let's do some clusters bitch
############
k = 6

t <- gene_plotter(de_late, what = 'RNA')
cutree(t$tree_row, k = k) %>% cluster_plotter(what = 'RNA')

cutree(t$tree_row, k = k) -> x 
t_cl <- split(x, x)
names_cl <- lapply(t_cl, names)
ck <- compareCluster(names_cl, fun = 'enrichKEGG', organism = 'spo')
dotplot(ck)

t2 <- gene_plotter(ids_sigs, what = 'RNA')
cutree(t2$tree_row, k = k) %>% cluster_plotter(what = 'RNA')


cutree(t2$tree_row, k = k) -> x2 
t_cl2 <- split(x2, x2)
names_cl2 <- lapply(t_cl2, names)
ck2 <- compareCluster(names_cl2, fun = 'enrichKEGG', organism = 'spo')
dotplot(ck2)

#############
#cdc2 targets
#############

cdc2_targets <- read.delim('C:/Users/am4613/Desktop/cdc2_targets.txt', header = F, strings = F)
gene_plotter(cdc2_targets[,1], what = 'RNA')

targets_in_cl <- c(1:k)
target_names <- list()
for(i in 1:k)
{
	targets_in_cl[i] <-  length(intersect(cdc2_targets[,1], names_cl[[i]]))
	target_names[[i]] <- intersect(cdc2_targets[,1], names_cl[[i]])
}

t <- plyr::ldply(target_names, data.frame)
gene_plotter(t[,1], what = 'RNA')

source('C:/Users/am4613/Documents/GitHub/RNA_seq/cell_cycle_plotter.R')
cell_cycle_heatmap(ids_sigs)
cell_cycle_plotter(ids_sigs)

cell_cycle_heatmap(de_late)
cell_cycle_plotter(de_late)

####################
#Get promoters for all the clusters (i'm going to use early versus late but that can be changed)
######################

save(names_cl2, file = 'Y:/Pers - Amalia/Promoter analysis/clusters_RNA_2.rda')

###################
#Cluster 6 seems to be enriched in H3K9 methilation o.o look at overlaps
####################

cluster_6 <- read.delim('Y:/Pers - Amalia/chip_ncRNA/CHIP_SEQ/clusters/cluster_6.txt', header = T, strings = F)
methil_genes <- read.delim('Y:/Pers - Amalia/chip_ncRNA/methil_genes.txt', header = T, strings = F)
grid.newpage()
int <- intersect(cluster_6[,1], methil_genes[,1])

cluster_7 <- read.delim('Y:/Pers - Amalia/chip_ncRNA/CHIP_SEQ/clusters/cluster_7.txt', header = T, strings = F)
int <- intersect(cluster_7[,1], methil_genes[,1])

cluster_8 <- read.delim('Y:/Pers - Amalia/chip_ncRNA/CHIP_SEQ/clusters/cluster_8.txt', header = T, strings = F)
int <- intersect(cluster_8[,1], methil_genes[,1])

########################
#Separate data in clustes, check sam's data in all of them to talk about chromatin modifications
#######################
cl_df <- plyr::ldply(names_cl2, data.frame)
colnames(cl_df) <- c('Clusters', 'ID')

non_regulated <- row.names(all_counts)[!(row.names(all_counts) %in% cl_df$ID)]

non_regulated_df <- data.frame(Clusters = 'Non regulated',ID = non_regulated)

all_genes_df <- as.data.frame(rbind(cl_df, non_regulated_df)) ##Create data frame with all genes

##merge with all_gene_data
merg_cl <- inner_join(all_genes_df, all_data, by = c('ID' = 'Systematic.name'))

p <- ggplot(merg_cl, aes(x = Clusters, y= as.numeric(MM.mRNA.cpc), color= Clusters, fill = Clusters))
p + geom_boxplot(alpha = 0.5) + theme_bw() +scale_fill_brewer(palette = 'Set1') + scale_color_brewer(palette = 'Set1')+ scale_y_log10()

p <- ggplot(merg_cl, aes(x = Clusters, y= as.numeric(MM.protein.cpc), color= Clusters, fill = Clusters))
p + geom_boxplot(alpha = 0.5) + theme_bw() +scale_fill_brewer(palette = 'Set1') + scale_color_brewer(palette = 'Set1')+ scale_y_log10()

p <- ggplot(merg_cl, aes(x = Clusters, y= as.numeric(exonic.length), color= Clusters, fill = Clusters))
p + geom_boxplot(alpha = 0.5) + theme_bw() +scale_fill_brewer(palette = 'Set1') + scale_color_brewer(palette = 'Set1')+ scale_y_log10()

##merge with half_life data
merg_half_rna <- inner_join(all_genes_df, mrna_half, by = c('ID' = 'Gene.name'))

p <- ggplot(merg_half_rna, aes(x = Clusters, y= as.numeric(Average.half.life), color= Clusters, fill = Clusters))
p + geom_boxplot(alpha = 0.5) + theme_bw() +scale_fill_brewer(palette = 'Set1') + scale_color_brewer(palette = 'Set1')

merg_half_prot <- inner_join(all_genes_df, prot_half, by = c('ID' = 'ENSG'))

p <- ggplot(merg_half_prot, aes(x = Clusters, y= as.numeric(t1.2..hours.), color= Clusters, fill = Clusters))
p + geom_boxplot(alpha = 0.5) + theme_bw() +scale_fill_brewer(palette = 'Set1') + scale_color_brewer(palette = 'Set1')
