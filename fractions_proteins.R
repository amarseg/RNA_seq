setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/')
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')

avg_isx <- read.delim('isx_data_summary.txt')

metabolic <- read.delim('aa_carbohidrate_metabolism', header = T, strings = F)
metabolic_id <- metabolic$ensembl_id
ribosomes <- read.delim('ribosome&biogenesisGenes', header = T, strings = F)
ribosomes_id <- ribosomes$ensembl_id

prot_data <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)

##Preprocessing, and averaging between technical replicates.Reordering

norm_prot <- normalise_ProtDataset(prot_data, what = 'nada')
norm_prot <- norm_prot[,7:42]

separate_prot <- list()

for(i in 1:3)
{
	separate_prot[[i]] <- norm_prot[,seq(i,33+i,3)]
}

fractions <- as.data.frame(matrix(nrow = 3, ncol = ncol(norm_prot), NA))
row.names(fractions) <- c('ribosome', 'metabolism', 'rest')
colnames(fractions) <- colnames(norm_prot)

met <- colSums(norm_prot[which(row.names(norm_prot) %in% metabolic_id),])
ribo <- colSums(norm_prot[which(row.names(norm_prot) %in% ribosomes_id), ])
total <- colSums(norm_prot)

fractions[1,] <- ribo/total
fractions[2,] <- met/total
fractions[3,] <- (total - ribo - met)/total

separate_frac <- list()

for(i in 1:3)
{
	separate_frac[[i]] <- fractions[,seq(i,33+i,3)]
}

par(mfrow = c(1,3))
for(i in 1:3)
{
	dfbar <- barplot(as.matrix(separate_frac[[i]]), col = terrain.colors(3))
	par(new= T)
	length = avg_isx[avg_isx$rep == i, ]
	plot(x = dfbar, y = length$Length_Erode.M03..4., type = 'o', col = 'blue',xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	axis(4)
}

