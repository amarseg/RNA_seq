library('locfit')
library('rootSolve')
library('plyr')


load('P:/CDC2_RNA_051015.rda')

sam_rpkm <- cbind(CDC2_2,CDC2_3,CDC2_4)

meh <- sam_rpkm[,grep(colnames(sam_rpkm), pattern = 'Annot', invert = T)]
meh <- meh[,grep(colnames(meh), pattern = 'MM', invert = T)]
meh <- meh[,grep(colnames(meh), pattern = 'common', invert = T)]

time <- rep(0:11,each = 3)

result_loess <- list()
result_spline <- list()
result_locfit <- list()
result_locfit2 <- list()
result_locfit1 <- list()
rev <- list()

for(i in 1:nrow(meh))
{
	test <- cbind(t(meh[i,]),time)
	colnames(test) <- c('gene','time')
	test <- as.data.frame(test)
	result_loess[[i]] <- loess(gene~time, data = test)
	result_locfit[[i]] <- locfit(gene~time, data = test)
	result_locfit2[[i]] <- locfit(gene~time, data = test, deriv = c(1,1)) #That's how you express the second derivative
	result_locfit1[[i]] <- locfit(gene~time, data = test, deriv = 1)
	result_spline[[i]] <- smooth.spline(x = test$time, y= test$gene, all.knots = T)
}

###Code to plot the fit of any gene
plot_loess_fit <- function(gene_name)
{
	j = which(row.names(norm_sq) == gene_name)
	plot(y = norm_sq[j,7:42], x = time)
	lines(time,result_loess[[j]]$fitted, col = 'red', lwd = 3)
	lines(result_spline[[j]], col = 'blue', lwd = 3)
	derivada <- predict.smooth.spline(result[[j]])
	lines(derivada$y, col = 'green', lwd = 3)
	lines(result_locfit[[j]], col = 'purple', lwd = 3)
	plot(result_locfit2[[j]])
}

roots_locfit2 = list()
for(i in 1:length(result_locfit2))
{
	roots_locfit2[[i]] <- uniroot.all(function(t) predict(result_locfit2[[i]],t), interval = c(0,11))
}

roots_locfit2 <- ldply(roots_locfit2, rbind)

raices_2 <- rbind(roots_locfit2[,1],roots_locfit2[,2],roots_locfit2[,3])
hist(raices_2, breaks = 100, main = 'Histogram of roots', freq = F)
lines(density(raices_2, na.rm = T), col = 'red')

roots_locfit1 = list()
for(i in 1:length(result_locfit1))
{
	roots_locfit1[[i]] <- uniroot.all(function(t) predict(result_locfit1[[i]],t), interval = c(0,11))
}

roots_locfit1 <- ldply(roots_locfit1, rbind)

raices_1 <- rbind(roots_locfit1[,1],roots_locfit1[,2],roots_locfit1[,3], roots_locfit1[,4])
hist(raices_1, breaks = 100, main = 'Histogram of roots', freq = F)
lines(density(raices_1, na.rm = T), col = 'red')


myInterval <- rowSums(roots_locfit2 > 3 & roots_locfit2 < 5 , na.rm = T)
genesInterval <- row.names(meh[myInterval > 0,])
write.table(genesInterval, sep = '\t', 'R_3&5_2ndstDerivative.txt')
