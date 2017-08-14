###Annotating genes in cell cycle phases (using Sam's GO final table)
rm(list = ls())
setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/')
peak_times <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/fission_timecourses/Peaktimes all genes.txt', header= T, strings = F)
load('../GitHub/Misc/GO.analysis.110914.rda')
source('../GitHub/RNA_seq/cell_cycle_plotter.R')

go_extractor <- function(term_name)
{
	toDo <- GOfinal[,colnames(GOfinal) == term_name]
	genes <- row.names(GOfinal[which(toDo==1),])
	return(genes)
}


colnames_cell_cycle <- grep(colnames(GOfinal), pattern = 'Rustici', value = T)
cell_cycle_genes <- lapply(colnames_cell_cycle, FUN = go_extractor)
names(cell_cycle_genes) <- colnames_cell_cycle

rna_list <- read.delim('both_clustering/Cluster_2.txt', header = T, strings = F)
prot_list <- read.delim('protein_isx_correlation/negative_correlation_length.txt', header = T, strings = F)

gene_list <- c(rna_list[,1], row.names(prot_list))
gene_list <- as.data.frame(unique(gene_list))
gene_list$cell_cycle <- NA

for(i in 1:5)
{
	x <- cell_cycle_genes[[i]]
	gene_list[which(gene_list[,1] %in% x),]$cell_cycle <- names(cell_cycle_genes)[i]
}

only_cell_cycle <- na.omit(gene_list)
only_cell_cycle$cell_cycle <- substr(only_cell_cycle$cell_cycle,start = 1, stop= 2)

barplot(table(only_cell_cycle[,2]), las = 2)
