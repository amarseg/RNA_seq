library("DESeq")


path <- "~"
filename <-"count_table.tsv"

count_table <- read.delim("count_table.tsv", header = T, strings = F)



dataset <- newCountDataSet(count_table, condition = factor(c(0:11)))

dataset <- estimateSizeFactors(dataset)

##NO BIOLOGICAL REPLICATES!!
dataset <- estimateDispersions(dataset, method = "blind", sharingMode = "fit-only")

results <- list()

for (i in 0:10)
{
	results[[i+1]] <- nbinomTest(dataset, as.character(i), as.character(i+1))
}

save(results, file = "results.rda")
