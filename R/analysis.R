library(R6)

syncon_factory <- R6Class(
	"SyntheticControl",
	private = list(
		# Inputs
		data_full = NA,
		data_time = NA,
		data_time_no_offset = NA,
		data_pca = NA,
		intervention_date = NA,
		time_points = NA,
		n_seasons = NA,
		year_def = NA,
		
		# Produced by cross-validation
		cv_impact_full = NA,
		cv_impact_time = NA,
		cv_impact_time_no_offset = NA,
		cv_impact_pca = NA
	),
	public = list(
		# Produced by cross-validation
		quantiles_stack = NA,
		ann_pred_quantiles_stack = NA, 
		pred_quantiles_stack = NA,
		rr_roll_stack = NA,
		rr_mean_stack = NA,
		rr_mean_stack_intervals = NA,
		cumsum_prevented_stack = NA,
		stacking_weights.all.m = NA,
		
		# Produced by main analysis
		impact_full = NA,
		impact_time = NA,
		impact_time_no_offset = NA,
		impact_pca = NA,
		
		# Produced by sensitivy analysis
		sensitivity_pred_quantiles = NA,
    sensitivity_table = NA,
    sensitivity_table_intervals = NA,
    rr_table = NA,
    rr_table_intervals = NA,
		
		initialize = function(data_full, data_time, data_time_no_offset, data_pca, intervention_date, time_points, n_seasons, year_def) {
			private$data_full <- data_full
			private$data_time <- data_time
			private$data_time_no_offset <- data_time_no_offset
			private$data_pca <- data_pca
			private$intervention_date <- intervention_date
			private$time_points <- time_points
			private$n_seasons <- n_seasons
			private$year_def <- year_def
		},
		impact = function() {
			#Start Cluster for CausalImpact (the main analysis function).
			cl <- makeCluster(n_cores)
			clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE)})
			clusterExport(cl, c('doCausalImpact',  'intervention_date', 'time_points', 'n_seasons'), environment())
			self$impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, var.select.on=TRUE, time_points = time_points), groups)
			self$impact_time <- setNames(parLapply(cl, data_time, doCausalImpact, intervention_date = intervention_date,  var.select.on=FALSE,time_points = time_points, trend = TRUE), groups)
			self$impact_time_no_offset <- setNames(parLapply(cl, data_time_no_offset, doCausalImpact, intervention_date = intervention_date,  var.select.on=FALSE,time_points = time_points,  trend = FALSE), groups)
			self$impact_pca <- setNames(parLapply(cl, data_pca, doCausalImpact, intervention_date = intervention_date, var.select.on=FALSE, time_points = time_points), groups)
			#No covariates, but with random intercept
			#impact_seas_only <- setNames(parLapply(cl, data_null, doCausalImpact, intervention_date = intervention_date, var.select.on=FALSE, time_points = time_points, n_seasons = n_seasons), groups)
			#No covariates except seasonal, no random intercept
			#impact_seas_only_no_re <- setNames(parLapply(cl, data_null, doCausalImpact, intervention_date = intervention_date, var.select.on=FALSE,ri.select=FALSE, time_points = time_points, n_seasons = n_seasons), groups)
			stopCluster(cl)
		},
		crossval = function() {
			#Creates List of lists: 1 entry for each stratum; within this, there are CV datasets for each year left out, and within this, there are 2 lists, one with full dataset, and one with the CV dataset
			cv.data_full<-lapply(data_full, makeCV)
			cv.data_time<-lapply(data_time, makeCV)
			cv.data_time_no_offset<-lapply(data_time_no_offset, makeCV)
			cv.data_pca<-lapply(data_pca, makeCV)
			#zoo_data<-cv.data_time[[1]][[2]]
			#Run the models on each of these datasets
			cl <- makeCluster(n_cores)
			clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE)})
			clusterExport(cl, c('doCausalImpact',  'intervention_date', 'time_points', 'n_seasons','crossval'), environment())
			private$cv_impact_full <-setNames(parLapply(cl, cv.data_full, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, intervention_date = intervention_date,  var.select.on=TRUE, time_points = time_points)), groups)
			private$cv_impact_time_no_offset <-setNames(parLapply(cl, cv.data_time_no_offset, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, trend=FALSE, intervention_date = intervention_date,  var.select.on=FALSE, time_points = time_points)), groups)
			private$cv_impact_time <-setNames(parLapply(cl, cv.data_time, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, trend=TRUE, intervention_date = intervention_date,  var.select.on=FALSE, time_points = time_points)), groups)
			private$cv_impact_pca <-setNames(parLapply(cl, cv.data_pca, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, intervention_date = intervention_date,  var.select.on=FALSE, time_points = time_points)), groups)
			stopCluster(cl)
			
			#Calculate pointwise log likelihood for cross-val prediction sample vs observed
			#These are N_iter*N_obs*N_cross_val array
			ll.cv.full<-lapply(private$cv_impact_full, function(x) lapply(x,crossval.log.lik))
			ll.cv.full2<-lapply(ll.cv.full, reshape.arr)
			#
			ll.cv.time_no_offset<-lapply(private$cv_impact_time_no_offset, function(x) lapply(x,crossval.log.lik))
			ll.cv.time_no_offset2<-lapply(ll.cv.time_no_offset, reshape.arr)
			#
			ll.cv.time<-lapply(private$cv_impact_time, function(x) lapply(x,crossval.log.lik))
			ll.cv.time2<-lapply(ll.cv.time, reshape.arr)
			#
			ll.cv.pca<-lapply(private$cv_impact_pca, function(x) lapply(x,crossval.log.lik))
			ll.cv.pca2<-lapply(ll.cv.pca, reshape.arr)
			#Create list that has model result for each stratum
			ll.compare<- vector("list", length(ll.cv.pca2)) 
			stacking_weights.all<-matrix(NA, nrow=length(ll.cv.pca2), ncol=4)
			
			for(i in 1:length(ll.compare)){
				ll.compare[[i]]<-cbind(ll.cv.full2[[i]],ll.cv.time_no_offset2[[i]],ll.cv.time2[[i]],ll.cv.pca2[[i]])#will get NAs if one of covariates is constant in fitting period (ie pandemic flu dummy)...shoud=ld fix this above
				keep<-complete.cases(ll.compare[[i]])
				ll.compare[[i]]<-ll.compare[[i]][keep,]
				#occasionally if there is a very poor fit, likelihood is very very small, which leads to underflow issue and log(0)...delete these rows to avoid this as a dirty solution. Better would be to fix underflow
				row.min<-apply(exp(ll.compare[[i]]),1,min)
				ll.compare[[i]]<-ll.compare[[i]][!(row.min==0),]
				#if(min(exp(ll.compare[[i]]))>0){
				stacking_weights.all[i,]<-stacking_weights(ll.compare[[i]])
				#}
			}
			stacking_weights.all<-as.data.frame(round(stacking_weights.all,3))
			names(stacking_weights.all)<-c('Synthetic Controls', 'Time trend', 'Time trend (no offset)', 'STL+PCA')
			stacking_weights.all<-cbind.data.frame(groups,stacking_weights.all)
			self$stacking_weights.all.m<-melt(stacking_weights.all, id.vars='groups')
			# stacking_weights.all.m<-stacking_weights.all.m[order(stacking_weights.all.m$groups),]
			
			stacked.ests<-mapply(  FUN=stack.mean,group=groups,impact_full=self$impact_full,impact_time=self$impact_time,impact_time_no_offset=self$impact_time_no_offset,impact_pca=self$impact_pca, MoreArgs=list(stacking_weights.all=stacking_weights.all), SIMPLIFY=FALSE )
			# plot.stacked.ests<-lapply(stacked.ests,plot.stack.est)
			self$quantiles_stack <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = stacked.ests[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period, year_def = private$year_def)}), groups)
			self$pred_quantiles_stack <- sapply(self$quantiles_stack, getPred, simplify = 'array')
			self$rr_roll_stack <- sapply(self$quantiles_stack, FUN = function(quantiles_stack) {quantiles_stack$roll_rr}, simplify = 'array')
			self$rr_mean_stack <- round(t(sapply(self$quantiles_stack, getRR)),2)
			self$rr_mean_stack_intervals <- data.frame('Stacking Estimate (95% CI)'     = makeInterval(self$rr_mean_stack[, 2], self$rr_mean_stack[, 3], self$rr_mean_stack[, 1]), check.names = FALSE, row.names = groups)
			self$cumsum_prevented_stack <- sapply(groups, FUN = cumsum_func, quantiles = self$quantiles_stack, simplify = 'array')
			self$ann_pred_quantiles_stack <- sapply(self$quantiles_stack, getAnnPred, simplify = FALSE)
			#Preds: Compare observed and expected
			pred.cv.full<-lapply(private$cv_impact_full, function(x) sapply(x,pred.cv,simplify='array'))
			pred.cv.pca<-lapply(private$cv_impact_pca, function(x) sapply(x,pred.cv,simplify='array'))
		},
		sensitivity = function() {
			bad_sensitivity_groups <- sapply(covars_full, function (covar) {ncol(covar) <= n_seasons-1+3})
			sensitivity_covars_full <- covars_full[!bad_sensitivity_groups]
			sensitivity_ds <- ds[!bad_sensitivity_groups]
			sensitivity_impact_full <- self$impact_full[!bad_sensitivity_groups]
			sensitivity_groups <- groups[!bad_sensitivity_groups]
			
			if (length(sensitivity_groups)!=0) {
				#Weight Sensitivity Analysis - top weighted variables are excluded and analysis is re-run.
				cl <- makeCluster(n_cores)
				clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
				clusterExport(cl, c('sensitivity_ds', 'weightSensitivityAnalysis', 'sensitivity_groups', 'intervention_date', 'outcome', 'time_points', 'n_seasons',  'eval_period', 'post_period', 'rrPredQuantiles'), environment())
				sensitivity_analysis_full <- setNames(parLapply(cl, sensitivity_groups, weightSensitivityAnalysis, covars = sensitivity_covars_full, ds = sensitivity_ds, impact = sensitivity_impact_full, time_points = time_points, intervention_date = intervention_date, n_seasons = n_seasons, outcome = outcome,  eval_period = eval_period, post_period = post_period, year_def=private$year_def), sensitivity_groups)
				stopCluster(cl)
				
				self$sensitivity_pred_quantiles  <- lapply(sensitivity_analysis_full, FUN = function(sensitivity_analysis) {
					pred_list <- vector(mode = 'list', length = length(sensitivity_analysis))
					for (sensitivity_index in 1:length(sensitivity_analysis)) {
						pred_list[[sensitivity_index]] <- getPred(sensitivity_analysis[[sensitivity_index]])
					}
					return(pred_list)
				})
				
				#Table of rate ratios for each sensitivity analysis level
				self$sensitivity_table <- t(sapply(sensitivity_groups, sensitivityTable, sensitivity_analysis = sensitivity_analysis_full, original_rr = rr_mean_full))
				self$sensitivity_table_intervals <- data.frame('Estimate (95% CI)' = makeInterval(self$sensitivity_table[, 2],  self$sensitivity_table[, 3],  self$sensitivity_table[, 1]),
																									'Top Control 1' = self$sensitivity_table[, 'Top Control 1'],
																									'Inclusion Probability of Control 1' = self$sensitivity_table[, 'Inclusion Probability of Control 1'],
																									'Control 1 Estimate (95% CI)' = makeInterval(self$sensitivity_table[, 7],  self$sensitivity_table[, 8],  self$sensitivity_table[, 6]),
																									'Top Control 2' = self$sensitivity_table[, 'Top Control 2'],
																									'Inclusion Probability of Control 2' = self$sensitivity_table[, 'Inclusion Probability of Control 2'],
																									'Control 2 Estimate (95% CI)' = makeInterval(self$sensitivity_table[, 12],  self$sensitivity_table[, 13],  self$sensitivity_table[, 11]),
																									'Top Control 3' = self$sensitivity_table[, 'Top Control 3'],
																									'Inclusion Probability of Control 3' = self$sensitivity_table[, 'Inclusion Probability of Control 3'],
																									'Control 3 Estimate (95% CI)' = makeInterval(self$sensitivity_table[, 17],  self$sensitivity_table[, 18],  self$sensitivity_table[, 16]), check.names = FALSE)
				self$rr_table <- cbind.data.frame(round(rr_mean_time[!bad_sensitivity_groups, ],2), self$sensitivity_table)
				self$rr_table_intervals <- cbind('Trend Estimate (95% CI)' = rr_mean_time_intervals[!bad_sensitivity_groups, ], self$sensitivity_table_intervals)
			} else {
				self$sensitivity_table_intervals <- NA
			}
		}
	)
)
