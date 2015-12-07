library('locfit')
library('rootSolve')
library('plyr')

library('polynom')

load('P:/CDC2_RNA_051015.rda')

sam_rpkm <- cbind(CDC2_2,CDC2_3,CDC2_4)

rpkm <- sam_rpkm[,grep(colnames(sam_rpkm), pattern = 'Annot', invert = T)]
rpkm <- rpkm[,grep(colnames(rpkm), pattern = 'MM', invert = T)]
rpkm <- rpkm[,grep(colnames(rpkm), pattern = 'common', invert = T)]

rpkm <- rpkm[which(rowSums(rpkm) > 0),]

time <- rep(0:11,each = 3)

poly1 <- list()
poly2 <- list()
poly3 <- list()

best_model <- list()

time <- rep(0:11, each = 3)

for(i in 1:nrow(rpkm))
{
	test <- cbind(t(rpkm[i,]),time)
	colnames(test) <- c('gene','time')
	test <- as.data.frame(test)
	poly1 <- lm(test$gene ~ poly(test$time,1,raw=TRUE))
	poly2 <- lm(test$gene ~ poly(test$time,2,raw=TRUE))
	poly3 <- lm(test$gene ~ poly(test$time,3,raw=TRUE))
	
	adj_r_squared <- c(summary(poly1)$adj.r.squared,summary(poly2)$adj.r.squared,summary(poly3)$adj.r.squared)
	pos <- which(adj_r_squared == max(adj_r_squared))
	if(pos == 1)
	{
		best_model[[i]] <- poly1
	}else if(pos == 2)
	{
		best_model[[i]] <- poly2
	}else if(pos == 3)
	{
		best_model[[i]] <- poly3
	}
}

model_degree <- c(1:nrow(rpkm))
for(i in 1:nrow(rpkm))
{
	model_degree[i] <-  length(best_model[[i]]$coefficients) -1
}

hist(model_degree)

first_derivative <- list()
second_derivative <- list()
sol_first <- list()
sol_second <- list()
for(i in 1:nrow(rpkm))
{
	first_derivative[[i]] <- deriv(polynomial(best_model[[i]]$coefficients))
	sol_first[[i]] <- Re(polyroot(first_derivative[[i]]))
	second_derivative[[i]] <- deriv(polynomial(first_derivative[[i]]))
	sol_second[[i]] <- Re(polyroot(second_derivative[[i]]))
}

sol_firstdf <- ldply(sol_first, rbind)
sol_firstdf <- as.data.frame(sol_firstdf)
row.names(sol_firstdf) <- row.names(rpkm)

hist(na.omit(sol_firstdf)[,1], breaks = 100, freq = F)
lines(density(sol_firstdf[,1], na.rm = T, to = max(sol_firstdf[,1], na.rm = T)), col = 'red')

sol_seconddf <- ldply(sol_second, rbind)
sol_seconddf <- as.data.frame(sol_seconddf)
row.names(sol_seconddf) <- row.names(rpkm)

hist(sol_seconddf[,1], breaks = 100, freq =F)
lines(density(sol_seconddf[,1], na.rm = T, to = max(sol_seconddf[,1], na.rm = T)), col = 'red')

save(sol_seconddf, file = 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/second_der_rna.rda')
save(sol_firstdf, file = 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/first_der_rna.rda')
