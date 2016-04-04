library('edgeR')
setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/')
rm(list = ls())

count_table <- read.delim('all_rev_counts.txt', strings = F, header = T)
row.names(count_table) <- count_table[,1]
count_table <- count_table[,2:37]

group <- factor(rep(0:11,each = 3))

y <- DGEList(counts = count_table, group = group)
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, ,keep.lib.sizes = F]
y = calcNormFactors(y)
#x <- ns(rep(0:11,each = 3), df = 2)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

pval_cutoff = 0.05

result <- topTags(qlf, p.value = pval_cutoff, n = 300)

write.table(result, file = 'C:/Users/am4613/Documents/Summaries_as_timecourses/DE_genes/result_edgeR.txt', sep = '\t')
