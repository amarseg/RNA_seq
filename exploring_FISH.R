##Processing FISH QUANT data 
library('plyr')
library('lattice')
par(mfrow = c(1,1))

probe_list <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/Probe_list.txt', header = T, strings = F)

setwd('C:/Users/am4613/Documents/FISH_QUANT/Results_mature/Collated_results/')

file_list <- dir()
file_list <- file_list[grep(file_list, pattern = 'Results')]
data <- list()

for(i in 1:length(file_list))
{
	data[[i]] <- read.delim(file_list[i], header = T, sep = '\t')
}

data_df <- ldply(data, data.frame)
data_df$Channel <- as.factor(data_df$Channel)

data_df$Gene <- NA
data_df$Acession <- NA

for(i in 1:nrow(probe_list))
{
	gene <- probe_list[i,]
	data_df[data_df$Mix == gene$Mix & data_df$Channel == gene$Channel,]$Gene <- gene$Gene
	data_df[data_df$Mix == gene$Mix & data_df$Channel == gene$Channel,]$Acession <- gene$Accession
}

medians <- aggregate(data_df$Spots ~ data_df$Channel*data_df$Mix, data_df, median)
medians_v <- medians[,3]

mypanel <- function(x,y,...)
{
	panel.xyplot(x,y,...)
	panel.text(30,33, labels = medians_v)
}

densityplot(~data_df$Spots|data_df$Gene)

write.table(data_df,'results_FQ.txt', sep = '\t')
write.table(medians, 'median_result_FQ.txt', sep = '\t')
