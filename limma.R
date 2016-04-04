library('limma')
library('splines')
setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/')

count_table <- read.delim('all_rev_counts.txt', strings = F, header = T)
row.names(count_table) <- count_table[,1]
count_table <- count_table[,2:37]

group <- factor(rep(0:11,each = 3))
names(group) <- colnames(count_table)

y <- DGEList(counts = count_table, group = group)
y = calcNormFactors(y)

x <- ns(rep(0:11,each = 3), df = 5)

design <- model.matrix(~group)


v <- voom(y, design, plot = T)


fit <- lmFit(v, design)
fit <- eBayes(fit)
