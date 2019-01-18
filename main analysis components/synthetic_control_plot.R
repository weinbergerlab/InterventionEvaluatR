#This is the plotting file. Run it after the analysis file to visualize results if not using the rmarkdown report.

plot_list <- setNames(vector(mode = 'list', length = length(syncon$groups)), syncon$groups)
for (group in syncon$groups) {
	plots <- list()
	
	#View scaled covariates
	covars.sub<-syncon$covars$full[[group]][,-c(1:(syncon$n_seasons-1))]
	alpha1=rep(impact_results$full$inclusion_probs[[group]][-c(1:(syncon$n_seasons-1)),'inclusion_probs'], each=nrow(syncon$covars$full[[group]]))
	covar_plot <- ggplot(melt(covars.sub, id.vars = NULL), mapping = aes_string(x = rep(syncon$time_points, ncol(covars.sub)), y = 'value', group = 'variable', alpha = alpha1 )) + 
		geom_line() + 
		labs(x = 'Time', y = 'Scaled Covariates') + 
		ggtitle(paste(group, 'Scaled Covariates Weighted by Inclusion Probability')) +
		theme_bw() +
	  theme(axis.line = element_line(colour = "black"),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.border = element_blank(),
	        panel.background = element_blank()) +
		theme(legend.position = 'none', plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

		#Plot predictions
	min_max <- c(min(c(impact_results$full$pred_quantiles[, , group], syncon$outcome[, group])), max(c(impact_results$full$pred_quantiles[, , group], syncon$outcome[, group])))
	pred_best_plot <-        plotPred(impact_results$best$pred_quantiles[, , group], syncon$time_points, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'Best estimate'))
	pred_full_plot <-        plotPred(impact_results$full$pred_quantiles[, , group], syncon$time_points, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'Synthetic controls estimate'))
	pred_time_plot <-        plotPred(impact_results$time$pred_quantiles[, , group], syncon$time_points, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'Interupted time series estimate'))
	pred_pca_plot <-        plotPred(impact_results$pca$pred_quantiles[, , group], syncon$time_points, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'STL+PCA estimate'))
	if(params$crossval){
	  pred_stack_plot <-        plotPred(pred_quantiles_stack[, , group], syncon$time_points, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'Stacked estimate'))
	  pred_stack_plot_agg <-        plotPredAgg(ann_pred_quantiles_stack[[group]], syncon$time_points, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'Stacked estimate'))
	  
	}
	pred_best_plot_agg <-        plotPredAgg(impact_results$best$ann_pred_quantiles[[group]], syncon$time_points, syncon$year_def, syncon$intervention_date, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'Best estimate'))
	pred_full_plot_agg <-        plotPredAgg(impact_results$full$ann_pred_quantiles[[group]], syncon$time_points, syncon$year_def, syncon$intervention_date, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'Synthetic controls estimate'))
	pred_time_plot_agg <-        plotPredAgg(impact_results$time$ann_pred_quantiles[[group]], syncon$time_points, syncon$year_def, syncon$intervention_date, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'Interupted time series estimate'))
	pred_pca_plot_agg <-        plotPredAgg(impact_results$pca$ann_pred_quantiles[[group]], syncon$time_points, syncon$year_def, syncon$intervention_date, syncon$post_period, min_max, syncon$outcome[, group], title = paste(group, 'STL+PCA estimate'))
	
	if (exists('sensitivity_table_intervals')) {
		pred_sensitivity_plot <- plotPred(impact_results$full$pred_quantiles[, , group], syncon$time_points, syncon$post_period, min_max, syncon$outcome[, group], sensitivity_pred_quantiles = syncon$sensitivity_pred_quantiles[[group]], sensitivity_title = paste(group, 'Sensitivity Plots'), plot_sensitivity = TRUE)
	} else {
	  pred_sensitivity_plot <- NA
	}
	#matplot(pred_quantiles_full[, , 10], ylim=c(0,22000), type='l')	 ##Check
	#points(prelog_data[[10]]$J12_18)
	
	#Plot rolling rate ratio
	min_max <- c(min(impact_results$full$rr_roll[, , group], impact_results$time$rr_roll[, , group]), max(impact_results$full$rr_roll[, , group], impact_results$time$rr_roll[, , group]))
	rr_roll_best_plot <- ggplot(melt(as.data.frame(impact_results$best$rr_roll[, , group]), id.vars = NULL), mapping = aes_string(x = rep(syncon$time_points[(length(syncon$time_points) - nrow(impact_results$best$rr_roll[, , group]) + 1):length(syncon$time_points)], ncol(impact_results$best$rr_roll[, , group])), y = 'value', linetype = 'variable')) + 
	  geom_line() + geom_hline(yintercept = 1, linetype = 4) +
	  labs(x = 'Time', y = 'Rolling Rate Ratio') + 
	  ggtitle(paste(group, 'Synthetic Control Rolling Rate Ratio')) +
	  coord_cartesian(ylim = min_max) +
	  theme_bw() +
	  theme(axis.line = element_line(colour = "black"),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.border = element_blank(),
	        panel.background = element_blank()) +
	  theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	rr_roll_full_plot <- ggplot(melt(as.data.frame(impact_results$full$rr_roll[, , group]), id.vars = NULL), mapping = aes_string(x = rep(syncon$time_points[(length(syncon$time_points) - nrow(impact_results$full$rr_roll[, , group]) + 1):length(syncon$time_points)], ncol(impact_results$full$rr_roll[, , group])), y = 'value', linetype = 'variable')) + 
		geom_line() + geom_hline(yintercept = 1, linetype = 4) +
		labs(x = 'Time', y = 'Rolling Rate Ratio') + 
		ggtitle(paste(group, 'Synthetic Control Rolling Rate Ratio')) +
		coord_cartesian(ylim = min_max) +
		theme_bw() +
	  theme(axis.line = element_line(colour = "black"),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.border = element_blank(),
	        panel.background = element_blank()) +
		theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	rr_roll_time_plot <- ggplot(melt(as.data.frame(impact_results$time$rr_roll[, , group]), id.vars = NULL), mapping = aes_string(x = rep(syncon$time_points[(length(syncon$time_points) - nrow(impact_results$time$rr_roll[, , group]) + 1):length(syncon$time_points)], ncol(impact_results$time$rr_roll[, , group])), y = 'value', linetype = 'variable')) + 
		geom_line() + geom_hline(yintercept = 1, linetype = 4) +
		labs(x = 'Time', y = 'Rolling Rate Ratio') + 
		ggtitle(paste(group, 'TT Rolling Rate Ratio')) +
		coord_cartesian(ylim = min_max) +
		theme_bw() +
	  theme(axis.line = element_line(colour = "black"),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.border = element_blank(),
	        panel.background = element_blank()) +
		theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	rr_roll_pca_plot <- ggplot(melt(as.data.frame(impact_results$pca$rr_roll[, , group]), id.vars = NULL), mapping = aes_string(x = rep(syncon$time_points[(length(syncon$time_points) - nrow(impact_results$pca$rr_roll[, , group]) + 1):length(syncon$time_points)], ncol(impact_results$pca$rr_roll[, , group])), y = 'value', linetype = 'variable')) + 
	  geom_line() + geom_hline(yintercept = 1, linetype = 4) +
	  labs(x = 'Time', y = 'Rolling Rate Ratio') + 
	  ggtitle(paste(group, 'STL+PCA Rolling Rate Ratio')) +
	  coord_cartesian(ylim = min_max) +
	  theme_bw() +
	  theme(axis.line = element_line(colour = "black"),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.border = element_blank(),
	        panel.background = element_blank()) +
	  theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	if(params$crossval){
	  rr_roll_stack_plot <- ggplot(melt(as.data.frame(rr_roll_stack[, , group]), id.vars = NULL), mapping = aes_string(x = rep(syncon$time_points[(length(syncon$time_points) - nrow(rr_roll_stack[, , group]) + 1):length(syncon$time_points)], ncol(rr_roll_stack[, , group])), y = 'value', linetype = 'variable')) + 
	    geom_line() + geom_hline(yintercept = 1, linetype = 4) +
	    labs(x = 'Time', y = 'Rolling Rate Ratio') + 
	    ggtitle(paste(group, 'Stacked Rolling Rate Ratio')) +
	    coord_cartesian(ylim = min_max) +
	    theme_bw() +
	    theme(axis.line = element_line(colour = "black"),
	          panel.grid.major = element_blank(),
	          panel.grid.minor = element_blank(),
	          panel.border = element_blank(),
	          panel.background = element_blank()) +
	    theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	  
	}
	#Plot cumulative sums
	cumsum_prevented_plot <- ggplot(melt(as.data.frame(impact_results$best$cumsum_prevented[, , group]), id.vars = NULL), mapping = aes_string(x = rep(syncon$time_points, ncol(impact_results$best$cumsum_prevented[, , group])), y = 'value', linetype = 'variable')) + 
	  geom_line() + geom_hline(yintercept = 1, linetype = 4) +
		labs(x = 'Time', y = 'Cumulative Sum Prevented') + 
		ggtitle(paste(group, 'Cumulative Number of Cases Prevented')) + 
		theme_bw() +
	  theme(axis.line = element_line(colour = "black"),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.border = element_blank(),
	        panel.background = element_blank()) +
		theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	if( params$crossval){
	  cumsum_prevented_stack_plot <- ggplot(melt(as.data.frame(cumsum_prevented_stack[, , group]), id.vars = NULL), mapping = aes_string(x = rep(syncon$time_points, ncol(cumsum_prevented_stack[, , group])), y = 'value', linetype = 'variable')) + 
	    geom_line() + geom_hline(yintercept = 1, linetype = 4) +
	    labs(x = 'Time', y = 'Cumulative Sum Prevented') + 
	    ggtitle(paste(group, 'Cumulative Number of Cases Prevented (Stacked model)')) + 
	    theme_bw() +
	    theme(axis.line = element_line(colour = "black"),
	          panel.grid.major = element_blank(),
	          panel.grid.minor = element_blank(),
	          panel.border = element_blank(),
	          panel.background = element_blank()) +
	    theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	  
	}
	
	
if(params$crossval){	plot_list[[group]] <- list(covar_plot = covar_plot, 
                                         
                                         pred_stack_plot=pred_stack_plot,
	                           	pred_full_plot = pred_full_plot, 
															pred_time_plot = pred_time_plot, 
															pred_pca_plot = pred_pca_plot, 
															pred_stack_plot_agg=pred_stack_plot_agg,
															pred_full_plot_agg = pred_full_plot_agg, 
															pred_time_plot_agg = pred_time_plot_agg, 
															pred_pca_plot_agg = pred_pca_plot_agg, 
															#pred_sensitivity_plot = pred_sensitivity_plot, 
															rr_roll_stack_plot = rr_roll_stack_plot, 
															rr_roll_full_plot = rr_roll_full_plot, 
															rr_roll_time_plot = rr_roll_time_plot, 
															rr_roll_pca_plot = rr_roll_pca_plot,
															cumsum_prevented_stack_plot=cumsum_prevented_stack_plot,
															cumsum_prevented_plot = cumsum_prevented_plot)
}else {
  plot_list[[group]] <- list(covar_plot = covar_plot, 
                             pred_best_plot = pred_best_plot, 
                             pred_full_plot = pred_full_plot, 
                             pred_time_plot = pred_time_plot, 
                             pred_pca_plot = pred_pca_plot, 
                             pred_best_plot_agg = pred_best_plot_agg, 
                             pred_full_plot_agg = pred_full_plot_agg, 
                             pred_time_plot_agg = pred_time_plot_agg, 
                             pred_pca_plot_agg = pred_pca_plot_agg, 
                           #  pred_sensitivity_plot = pred_sensitivity_plot, 
                            rr_roll_best_plot = rr_roll_best_plot, 
                             rr_roll_full_plot = rr_roll_full_plot, 
                             rr_roll_time_plot = rr_roll_time_plot, 
                             rr_roll_pca_plot = rr_roll_pca_plot, 
                             cumsum_prevented_plot = cumsum_prevented_plot)
}
}

for (group in syncon$groups) {
	cat('###', group, '\n', sep = '')
	for (group_plot in plot_list[[group]]) {
		print(group_plot)
	}
	cat('\n\n')
}
