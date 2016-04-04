##This function gets a list of genes and plots the data from cycle base


cell_cycle_plotter <- function(ids, average = T, ylim = c(-2,2))
{
	load('C:/Users/am4613/Documents/Summaries_as_timecourses/fission_timecourses/cell_cycle_data.rda')
	load('C:/Users/am4613/Documents/Summaries_as_timecourses/fission_timecourses/exp_time.rda')
	
	
	
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
			lines(y = avg_levels, x = time_levels, col = col[i], fill = T, lty = i)
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