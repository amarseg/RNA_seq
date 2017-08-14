#analysing sam's clusters
rm(list = ls())
library(ddplyr)
setwd('C:/Users/am4613/OneDrive/gene_list_analysis/')
clusters <- read.delim('clusters_090916.txt', header = T, strings = F)
source('C:/Users/am4613/Documents/GitHub/Misc/half_life_plotter.R')
source('C:/Users/am4613/Documents/GitHub/Misc/gene_plotter.R')


boxplot(MM_prot ~ cluster, data = clusters)

cl_list <- split(clusters, f = clusters$cluster)


half_life_plotter(cl_list[[1]]$names)
half_life_plotter(cl_list[[2]]$names)
half_life_plotter(cl_list[[3]]$names)
half_life_plotter(cl_list[[4]]$names)

half_life_plotter(clusters$names, data = 'rna')
half_life_plotter(clusters$names, data = 'protein')

gene_plotter(cl_list[[1]]$names, what = 'RNA')
gene_plotter(cl_list[[2]]$names, what = 'RNA')
gene_plotter(cl_list[[3]]$names, what = 'RNA')
gene_plotter(cl_list[[4]]$names, what = 'RNA')
gene_plotter(clusters$name, what = 'RNA')

all_data <- read.delim('../Summaries_as_timecourses/mmc1_from_cell_paper.txt', header = T, strings = F, skip = 6)
all_data$scale <- NA
all_data[which(all_data$Systematic.name %in% clusters$names),]$scale <- 'Not Scaling'
all_data[which(is.na(all_data$scale)),]$scale <- 'Scaling'
p <- ggplot(all_data, aes(x = scale, y= log2(MM.mRNA.cpc), color = scale))
p + geom_boxplot(notch = T)
ggsave('copy_number.wmf')

all_data$cluster <- NA
all_data[which(all_data$Systematic.name %in% cl_list[[1]]$names),]$cluster <- '1'
all_data[which(all_data$Systematic.name %in% cl_list[[2]]$names),]$cluster <- '2'
all_data[which(all_data$Systematic.name %in% cl_list[[3]]$names),]$cluster <- '3'
all_data[which(all_data$Systematic.name %in% cl_list[[4]]$names),]$cluster <- '4'
all_data[which(is.na(all_data$cluster)),]$cluster <- 'Not Scaling'

p <- ggplot(all_data, aes(x = cluster, y= log2(MM.mRNA.cpc), color = scale))
p + geom_boxplot()
ggsave('copy_number_clusters.wmf', device = 'wmf')
