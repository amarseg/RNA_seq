##This function gets a list of genes and plots the data from cycle base


cell_cycle_plotter <- function(ids, average = T, ylim = c(-2,2))
{
	load('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/fission_timecourses/cell_cycle_data.rda')
	load('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/fission_timecourses/exp_time.rda')
	
	
	
	if(average){
		
		col <- rainbow(10)
		
		plot(0,type='l',xlim=c(0,350),ylim=ylim,xlab='Time', ylab='Expression Level', xaxt = 'n')
		axis(side = 1, at = c(0,15,30,90,100,115,130,190,200,215,230,290,300,315,330,390))
		
		exp_name <- 1:10
		for(i in 1:10)
		{
			exp_levels <- exp_list[[i]]
			time_levels <- time_list[which(names(time_list) == exp_levels$Experiment[1])]
			exp_name[i] <- exp_levels$Experiment[1]
			time_levels <- as.numeric(time_levels[[1]])
			ids_levels <- exp_levels[which(exp_levels$identifier %in% ids),]
			ids_levels[,4:ncol(ids_levels)] <- apply(ids_levels[,4:ncol(ids_levels)],c(1,2),as.numeric)
			avg_levels <- colMeans(ids_levels[,4:ncol(ids_levels)], na.rm = T)
			lines(y = avg_levels, x = time_levels, col = col[i], fill = T, lty = i, lwd = 3)
		}
		
		legend('topleft', legend = exp_name, fill = col, cex = 0.75, ncol = 4)
	}else{
		par(mfrow = c(5,2))
		for(i in 1:10)
		{
			exp_levels <- exp_list[[i]]
			time_levels <- time_list[which(names(time_list) == exp_levels$Experiment[1])]
			time_levels <- as.numeric(time_levels[[1]])
			ids_levels <- exp_levels[which(exp_levels$identifier %in% ids),]
			ids_levels[,4:ncol(ids_levels)] <- apply(ids_levels[,4:ncol(ids_levels)],c(1,2),as.numeric)
			matplot(y = t(ids_levels[,4:ncol(ids_levels)]), x = time_levels, col = 'gray', type = 'l',
							xlab = 'Cell Cycle Time', ylab = 'Expression Level', main = exp_levels$Experiment[1])
			avg_levels <- colMeans(ids_levels[,4:ncol(ids_levels)], na.rm = T)
			lines(y = avg_levels, x = time_levels, col = 'black')
		}
		par(mfrow = c(1,1))
	}
}

cell_cycle_heatmap <- function(ids)
{
	load('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/fission_timecourses/cell_cycle_data.rda')
	load('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/fission_timecourses/exp_time.rda')
	peak_times <- read.delim('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/fission_timecourses/Peaktimes all genes.txt', header= T, strings = F)
	library(gplots)
	
	rustici_elu <- merge(exp_list[[8]][,c(-1,-2)], exp_list[[9]][,c(-1,-2)], by = 'identifier', all = T)
	rustici_elu <- merge(rustici_elu, exp_list[[10]][,c(-1,-2)], by = 'identifier', all = T)
	rustici_elu <- merge(rustici_elu, peak_times, by.x = 'identifier', by.y = 'Gene.ID')
	
	
	rustici_elu[,-1] <- apply(rustici_elu[,-1], c(1,2), as.numeric)
	
	col = colorRampPalette(c('blue','gray','yellow'))
	
	test <- rustici_elu[which(ids %in% rustici_elu$identifier),]
	ordered_test <- test[order(test$Peak),]
	ordered_test[is.na(ordered_test)] <- 0
	heatmap.2(as.matrix(ordered_test[,2:61]), Colv = F, trace = 'none', col = col, labRow = ordered_test$identifier, Rowv = F)
}


