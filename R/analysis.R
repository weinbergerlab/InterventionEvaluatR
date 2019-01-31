syncon.init <- function(
  country, data,
  pre_period_start, pre_period_end,
  post_period_start, post_period_end,
  eval_period_start, eval_period_end,
  n_seasons, year_def, 
  group_name, date_name, outcome_name, denom_name
) {
  
  syncon = listenv(
    time_points = NA,

    groups = NA,
    sparse_groups = NA,
    model_size = NA,
    covars = NA,
    outcome = NA,
    
    results = list(
      impact=NA,
      crossval=NA,
      sensitivity=NA
    ),
    
    stacking_weights.all.m = NA, # TODO make private

    .private=listenv(
      # Variants
      variants = list(
        full = list(var.select.on=TRUE, trend=FALSE, name="Synthetic controls"),
        time = list(var.select.on=FALSE, trend=TRUE, name="Time trend"),
        time_no_offset = list(var.select.on=FALSE, trend=FALSE, name="Time trend (no offset)"),
        pca = list(var.select.on=FALSE, trend=FALSE, name="STL+PCA")
      ),

      exclude_covar = NA,
      exclude_group = NA,

      # Computation state
      data = list(),
      data.cv = list(),
      n_cores = NA,
      ds = NA
    )
  )

  syncon$country <- country #Country or region name.
  syncon$n_seasons <- n_seasons #Number of months (seasons) per year. 12 for monthly, 4 for quarterly, 3 for trimester data.
  syncon$year_def <- params$year_def #Can be cal_year to aggregate results by Jan-Dec; 'epi_year' to aggregate July-June

  #MOST DATES MUST BE IN FORMAT "YYYY-MM-01", exception is end of pre period, which is 1 day before end of post period
  syncon$pre_period <- as.Date(c(pre_period_start, pre_period_end)) #Range over which the data is trained for the CausalImpact model.
  syncon$start_date <- syncon$pre_period[1]
  syncon$post_period <- as.Date(c(post_period_start, post_period_end )) #Range from the intervention date to the end date.
  syncon$intervention_date <- syncon$post_period[1]-1
  syncon$end_date <- syncon$post_period[2]

  syncon$eval_period <- as.Date(c(eval_period_start, eval_period_end)) #Range over which rate ratio calculation will be performed.

  syncon$group_name <- group_name #Name of column containing group labels.
  syncon$date_name <- date_name #Name of column containing dates.
  syncon$outcome_name <- outcome_name #Name of column containing outcome.
  syncon$denom_name <- denom_name #Name of column containing denominator to be used in offset.

  syncon$.private$exclude_covar <- c() #User-defined list of covariate columns to exclude from all analyses.
  syncon$.private$exclude_group <- c() #User-defined list of groups to exclude from analyses.

  #Assign variable values
  syncon$input_data <- data
  return(syncon)
}

syncon.impact = function(syncon) {
  syncon.impact.pre(syncon)
  results = list()
  
  #Start Cluster for CausalImpact (the main analysis function).
  cl <- makeCluster(syncon$.private$n_cores)
  clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE)})
  clusterExport(cl, c('doCausalImpact'), environment())

  for (variant in names(syncon$.private$variants)) {
    results[[variant]]$groups <- setNames(pblapply(
      cl=cl, syncon$.private$data[[variant]], FUN=doCausalImpact, 
      syncon$intervention_date, 
      syncon$n_seasons,
      var.select.on=syncon$.private$variants[[variant]]$var.select.on, 
      time_points=syncon$time_points,
      trend=syncon$.private$variants[[variant]]$trend,
      crossval.stage=FALSE
    ), syncon$groups)
  }
  stopCluster(cl)

  for (variant in c('full', 'time')) {
    #Save the inclusion probabilities from each of the models
    results[[variant]]$inclusion_prob <- setNames(lapply(results[[variant]]$groups, inclusionProb), syncon$groups)
  }

  for (variant in names(syncon$.private$variants)) {
    #All model results combined
    results[[variant]]$quantiles <- setNames(lapply(syncon$groups, FUN=function(group) {
      rrPredQuantiles(impact = results[[variant]]$groups[[group]], denom_data = syncon$.private$ds[[group]][, syncon$denom_name], eval_period=syncon$eval_period, post_period=syncon$post_period, year_def=syncon$year_def, time_points=syncon$time_points, n_seasons=syncon$n_seasons)
    }), syncon$groups)
  }

  # Calculate best model
  syncon$model_size <- sapply(results$full$groups, modelsize_func, n_seasons=syncon$n_seasons)
  results$best$quantiles <- vector("list", length(results$full$quantiles)) 
  results$best$quantiles[syncon$model_size>=1] <- results$full$quantiles[syncon$model_size>=1]
  results$best$quantiles[syncon$model_size<1] <- results$pca$quantiles[syncon$model_size<1]
  results$best$quantiles <- setNames(results$best$quantiles, syncon$groups)

  for (variant in c("best", names(syncon$.private$variants))) {
    # Predictions, aggregated by year
    results[[variant]]$pred_quantiles <- sapply(results[[variant]]$quantiles, getPred, simplify = 'array')
    results[[variant]]$ann_pred_quantiles <- sapply(results[[variant]]$quantiles, getAnnPred, simplify = FALSE)
  }

  for (variant in c('full', 'best')) {
    # Pointwise RR and uncertainty for second stage meta variant
    results[[variant]]$log_rr_quantiles <- sapply(results[[variant]]$quantiles, FUN = function(quantiles) {quantiles$log_rr_full_t_quantiles}, simplify = 'array')
    dimnames(results[[variant]]$log_rr_quantiles)[[1]] <- syncon$time_points
    results[[variant]]$log_rr_sd <- sapply(results[[variant]]$quantiles, FUN = function(quantiles) {quantiles$log_rr_full_t_sd}, simplify = 'array')
    results[[variant]]$log_rr_full_t_samples.prec <- sapply(results[[variant]]$quantiles, FUN = function(quantiles) {quantiles$log_rr_full_t_samples.prec}, simplify = 'array')
  }

  for (variant in c("best", names(syncon$.private$variants))) {
  	# Rolling rate ratios
    results[[variant]]$rr_roll <- sapply(results[[variant]]$quantiles, FUN = function(quantiles) {quantiles$roll_rr}, simplify = 'array')
    # Rate ratios for evaluation period.
    results[[variant]]$rr_mean <- t(sapply(results[[variant]]$quantiles, getRR))
  }

  results$best$log_rr <- t(sapply(results$best$quantiles, getsdRR))

  for (variant in c("best", names(syncon$.private$variants))) {
    results[[variant]]$rr_mean_intervals <- setNames(data.frame(makeInterval(
      results[[variant]]$rr_mean[, 2], results[[variant]]$rr_mean[, 3], results[[variant]]$rr_mean[, 1]
    ), check.names = FALSE, row.names = syncon$groups), c(paste(syncon$.private$variants[[variant]]$name, 'Estimate (95% CI)')))
  }

  colnames(results$time$rr_mean) <- paste('Time_trend', colnames(results$time$rr_mean))

  for (variant in c("best", names(syncon$.private$variants))) {
    results[[variant]]$cumsum_prevented <- sapply(syncon$groups, FUN=cumsum_func, quantiles = results[[variant]]$quantiles, outcome=syncon$outcome, syncon$time_points, syncon$post_period, simplify = 'array')
  }

  #Run a classic ITS analysis
  rr.its1 <- lapply(syncon$.private$data$time, its_func, post_period=syncon$post_period, eval_period=syncon$eval_period, time_points=syncon$time_points)
  rr.t <- sapply(rr.its1, `[[`, "rr.q.t", simplify='array')
  results$its = list()
  results$its$rr_end <- t(sapply(rr.its1, `[[`, "rr.q.post", simplify='array')) 
  results$its$rr_mean_intervals <- data.frame('Classic ITS (95% CI)' = makeInterval(results$its$rr_end[, 2], results$its$rr_end[, 3], results$its$rr_end[, 1]), check.names = FALSE, row.names = syncon$groups)

  syncon$results$impact <- results
  return(results)
}

syncon.crossval = function(syncon) {
  results = list()

  #Creates List of lists: 1 entry for each stratum; within this, there are CV datasets for each year left out, and within this, there are 2 lists, one with full dataset, and one with the CV dataset
  for (variant in names(syncon$.private$variants)) {
    syncon$.private$data.cv[[variant]]<-lapply(syncon$.private$data[[variant]], makeCV, syncon$time_points, syncon$intervention_date)
  }

  #Run the models on each of these datasets
  cl <- makeCluster(syncon$.private$n_cores)
  clusterEvalQ(cl, {
    library(pogit, quietly = TRUE); 
    library(lubridate, quietly = TRUE)
  })
  clusterExport(cl, c('doCausalImpact'), environment())
  for (variant in names(syncon$.private$variants)) {
    results[[variant]]$groups <-setNames(pblapply(
      cl=cl, syncon$.private$data.cv[[variant]], FUN=function(x) lapply(
        x, doCausalImpact, 
        syncon$intervention_date, 
        syncon$n_seasons,
        time_points=syncon$time_points,
        crossval.stage=TRUE,
        var.select.on=syncon$.private$variants[[variant]]$var.select.on,
      )), syncon$groups)
  }
  stopCluster(cl)

  ll.cv = list()

  #Calculate pointwise log likelihood for cross-val prediction sample vs observed
  #These are N_iter*N_obs*N_cross_val array
  for (variant in names(syncon$.private$variants)) {
    ll.cv[[variant]] <- lapply(results[[variant]]$groups, function(x) lapply(x, crossval.log.lik))
    ll.cv[[variant]] <- lapply(ll.cv[[variant]], reshape.arr)
  }

  #Create list that has model result for each stratum
  ll.compare<- vector("list", length(ll.cv$pca)) 
  stacking_weights.all<-matrix(NA, nrow=length(ll.cv$pca), ncol=4)

  for(i in 1:length(ll.compare)){
    #will get NAs if one of covariates is constant in fitting period (ie pandemic flu dummy)...should fix this above
    ll.compare[[i]] <- cbind(ll.cv$full[[i]], ll.cv$time_no_offset[[i]], ll.cv$time[[i]], ll.cv$pca[[i]]) 
    keep <- complete.cases(ll.compare[[i]])
    ll.compare[[i]] <- ll.compare[[i]][keep,]
    #occasionally if there is a very poor fit, likelihood is very very small, which leads to underflow issue and log(0)...delete these rows to avoid this as a dirty solution. Better would be to fix underflow
    row.min <- apply(exp(ll.compare[[i]]), 1, min)
    ll.compare[[i]] <- ll.compare[[i]][!(row.min==0),]
    #if(min(exp(ll.compare[[i]]))>0){
    stacking_weights.all[i,] <- stacking_weights(ll.compare[[i]])
    #}
  }
  stacking_weights.all <- as.data.frame(round(stacking_weights.all,3))
  names(stacking_weights.all) <- lapply(syncon$.private$variants, function(v) { v$name })
  stacking_weights.all <- cbind.data.frame(data.frame(groups=syncon$groups), stacking_weights.all)
  syncon$stacking_weights.all.m <- melt(stacking_weights.all, id.vars='groups')
  # stacking_weights.all.m<-stacking_weights.all.m[order(stacking_weights.all.m$groups),]
  
  stacked.ests <- mapply(
    FUN=stack.mean,
    group=syncon$groups,
    impact_full=syncon$results$impact$full$groups,
    impact_time=syncon$results$impact$time$groups,
    impact_time_no_offset=syncon$results$impact$time_no_offset$groups,
    impact_pca=syncon$results$impact$pca$groups,
    MoreArgs=list(stacking_weights.all=stacking_weights.all, outcome=syncon$outcome), 
    SIMPLIFY=FALSE
  )
  # plot.stacked.ests<-lapply(stacked.ests,plot.stack.est)
  results$quantiles_stack <- setNames(lapply(syncon$groups, FUN = function(group) {
    rrPredQuantiles(impact = stacked.ests[[group]], denom_data = syncon$.private$ds[[group]][, syncon$denom_name], syncon$eval_period, syncon$post_period, syncon$n_seasons, syncon$year_def, syncon$time_points)
  }), syncon$groups)
  results$pred_quantiles_stack <- sapply(results$quantiles_stack, getPred, simplify = 'array')
  results$rr_roll_stack <- sapply(results$quantiles_stack, FUN = function(quantiles_stack) {quantiles_stack$roll_rr}, simplify = 'array')
  results$rr_mean_stack <- round(t(sapply(results$quantiles_stack, getRR)), 2)
  results$rr_mean_stack_intervals <- data.frame('Stacking Estimate (95% CI)' = makeInterval(results$rr_mean_stack[, 2], results$rr_mean_stack[, 3], results$rr_mean_stack[, 1]), check.names = FALSE, row.names = syncon$groups)
  results$cumsum_prevented_stack <- sapply(syncon$groups, FUN = cumsum_func, quantiles = results$quantiles_stack, outcome=syncon$outcome, syncon$time_points, syncon$post_period, simplify = 'array')
  results$ann_pred_quantiles_stack <- sapply(results$quantiles_stack, getAnnPred, simplify = FALSE)
  #Preds: Compare observed and expected
  results$full$pred <- lapply(results$impact$full, function(x) sapply(x,pred.cv,simplify='array'))
  results$pca$pred <- lapply(results$impact$pca, function(x) sapply(x,pred.cv,simplify='array'))
  
  syncon$results$crossval <- results
  return(results)
}

syncon.sensitivity = function(syncon) {
  results = list()
  bad_sensitivity_groups <- sapply(syncon$covars$full, function (covar) {ncol(covar) <= syncon$n_seasons-1+3})
  sensitivity_covars_full <- syncon$covars$full[!bad_sensitivity_groups]
  sensitivity_ds <- syncon$.private$ds[!bad_sensitivity_groups]
  sensitivity_impact_full <- syncon$results$impact$full$groups[!bad_sensitivity_groups]
  sensitivity_groups <- syncon$groups[!bad_sensitivity_groups]

  if (length(sensitivity_groups)!=0) {
    #Weight Sensitivity Analysis - top weighted variables are excluded and analysis is re-run.
    cl <- makeCluster(syncon$.private$n_cores)
    clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
    clusterExport(cl, c('sensitivity_ds', 'weightSensitivityAnalysis', 'sensitivity_groups'), environment())
    sensitivity_analysis_full <- setNames(pblapply(cl=cl, sensitivity_groups, FUN=weightSensitivityAnalysis, covars = sensitivity_covars_full, ds = sensitivity_ds, impact = sensitivity_impact_full, time_points = syncon$time_points, intervention_date = syncon$intervention_date, n_seasons = syncon$n_seasons, outcome = syncon$outcome, eval_period = syncon$eval_period, post_period = syncon$post_period, year_def=syncon$year_def), sensitivity_groups)
    stopCluster(cl)
  
    results$sensitivity_pred_quantiles <- lapply(sensitivity_analysis_full, FUN = function(sensitivity_analysis) {
      pred_list <- vector(mode = 'list', length = length(sensitivity_analysis))
      for (sensitivity_index in 1:length(sensitivity_analysis)) {
        pred_list[[sensitivity_index]] <- getPred(sensitivity_analysis[[sensitivity_index]])
      }
      return(pred_list)
    })
  
    #Table of rate ratios for each sensitivity analysis level
    results$sensitivity_table <- t(sapply(sensitivity_groups, sensitivityTable, sensitivity_analysis = sensitivity_analysis_full, original_rr = syncon$results$impact$full$rr_mean))
    results$sensitivity_table_intervals <- data.frame(
      'Estimate (95% CI)' = makeInterval(results$sensitivity_table[, 2], results$sensitivity_table[, 3], results$sensitivity_table[, 1]),
      'Top Control 1' = results$sensitivity_table[, 'Top Control 1'],
      'Inclusion Probability of Control 1' = results$sensitivity_table[, 'Inclusion Probability of Control 1'],
      'Control 1 Estimate (95% CI)' = makeInterval(results$sensitivity_table[, 7], results$sensitivity_table[, 8], results$sensitivity_table[, 6]),
      'Top Control 2' = results$sensitivity_table[, 'Top Control 2'],
      'Inclusion Probability of Control 2' = results$sensitivity_table[, 'Inclusion Probability of Control 2'],
      'Control 2 Estimate (95% CI)' = makeInterval(results$sensitivity_table[, 12], results$sensitivity_table[, 13], results$sensitivity_table[, 11]),
      'Top Control 3' = results$sensitivity_table[, 'Top Control 3'],
      'Inclusion Probability of Control 3' = results$sensitivity_table[, 'Inclusion Probability of Control 3'],
      'Control 3 Estimate (95% CI)' = makeInterval(results$sensitivity_table[, 17], results$sensitivity_table[, 18], results$sensitivity_table[, 16]), check.names = FALSE
    )
    results$rr_table <- cbind.data.frame(round(syncon$results$impact$time$rr_mean[!bad_sensitivity_groups, ],2), results$sensitivity_table)
    results$rr_table_intervals <- cbind('Trend Estimate (95% CI)' = syncon$results$impact$time$rr_mean_intervals[!bad_sensitivity_groups, ], results$sensitivity_table_intervals)
  } else {
    results$sensitivity_table_intervals <- NA
  }
  
  syncon$results$sensitivity <- results
  return(results)
}

syncon.impact.pre = function(syncon) {
  # Setup data
  prelog_data <- syncon$input_data[!is.na(syncon$input_data[, syncon$outcome_name]),]#If outcome is missing, delete 
  prelog_data[, syncon$group_name] = prelog_data[, syncon$group_name] %% 2
  syncon$groups <- as.character(unique(unlist(prelog_data[, syncon$group_name], use.names = FALSE)))
  if (exists('exclude_group')) {
    syncon$groups <- syncon$groups[!(syncon$groups %in% exclude_group)]
  }

  # Format covars
  prelog_data[,syncon$date_name]<-as.Date(as.character(prelog_data[,syncon$date_name]), tryFormats=c("%m/%d/%Y",'%Y-%m-%d' ))

  #test<-split(prelog_data, factor(prelog_data[,syncon$group_name]))
  #outcome.na<-sapply(test, function(x) sum(is.na(x[,syncon$outcome_name])))
  prelog_data[, syncon$date_name] <- formatDate(prelog_data[, syncon$date_name])
  prelog_data <- setNames(lapply(
    syncon$groups, FUN=splitGroup, 
    ungrouped_data=prelog_data, group_name=syncon$group_name, date_name=syncon$date_name, 
    start_date=syncon$start_date, end_date=syncon$end_date, 
    no_filter=c(syncon$group_name, syncon$date_name, syncon$outcome_name, syncon$denom_name)
  ), syncon$groups)
  #if (exists('exclude_group')) {prelog_data <- prelog_data[!(names(prelog_data) %in% exclude_group)]}

  #Log-transform all variables, adding 0.5 to counts of 0.
  syncon$.private$ds <- setNames(lapply(prelog_data, FUN = logTransform, no_log = c(syncon$group_name, syncon$date_name, syncon$outcome_name)), syncon$groups)
  syncon$time_points <- unique(syncon$.private$ds[[1]][, syncon$date_name])

  #Monthly dummies
  if(syncon$n_seasons==4){
    dt <- quarter(as.Date(syncon$time_points))
  }
  if(syncon$n_seasons==12){
    dt <- month(as.Date(syncon$time_points))
  }
  if(syncon$n_seasons==3){
    dt.m <- month(as.Date(syncon$time_points))
    dt <- dt.m
    dt[dt.m %in% c(1,2,3,4)] <- 1
    dt[dt.m %in% c(5,6,7,8)] <- 2
    dt[dt.m %in% c(9,10,11,12)] <- 3
  }
  season.dummies <- dummies::dummy(dt)
  season.dummies <- as.data.frame(season.dummies)
  names(season.dummies) <- paste0('s', 1:syncon$n_seasons)
  season.dummies <- season.dummies[,-syncon$n_seasons]

  syncon$.private$ds <- lapply(syncon$.private$ds, function(ds) {
    if (!(syncon$denom_name %in% colnames(ds))) {
      ds[syncon$denom_name] <- 0
    }
    return(ds)
  })

  syncon$sparse_groups <- sapply(syncon$.private$ds, function(ds) {
    return(ncol(ds[!(colnames(ds) %in% c(syncon$date_name, syncon$group_name, syncon$denom_name, syncon$outcome_name, syncon$.private$exclude_covar))]) == 0)
  })
  syncon$.private$ds <- syncon$.private$ds[!syncon$sparse_groups]
  syncon$groups <- syncon$groups[!syncon$sparse_groups]

  #Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
  syncon$covars = list()
  syncon$covars$full <- setNames(lapply(syncon$.private$ds, function(group) { makeCovars(syncon$country, syncon$time_points, syncon$intervention_date, season.dummies, group) }), syncon$groups)
  syncon$covars$full <- lapply(syncon$covars$full, FUN = function(covars) {covars[, !(colnames(covars) %in% syncon$.private$exclude_covar), drop = FALSE]})
  syncon$covars$time <- setNames(lapply(syncon$covars$full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies,time_index = 1:nrow(covars))))}), syncon$groups)
  syncon$covars$null <- setNames(lapply(syncon$covars$full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies)))}), syncon$groups)

  #Standardize the outcome variable and save the original mean and SD for later analysis.
  syncon$outcome <- sapply(syncon$.private$ds, FUN = function(data) {data[, syncon$outcome_name]})
  offset <- sapply(syncon$.private$ds, FUN=function(data) exp(data[, syncon$denom_name]) ) #offset term on original scale; 1 column per age group

  ##SECTION 1: CREATING SMOOTHED VERSIONS OF CONTROL TIME SERIES AND APPENDING THEM ONTO ORIGINAL DATAFRAME OF CONTROLS
  #EXTRACT LONG TERM TREND WITH DIFFERENT LEVELS OF SMOOTHNESS USING STL
  # Set a list of parameters for STL
  stl.covars<-mapply(smooth_func,ds.list=syncon$.private$ds,covar.list=syncon$covars$full, SIMPLIFY=FALSE, MoreArgs=list(n_seasons=syncon$n_seasons)) 
  post.start.index<-which(syncon$time_points==syncon$post_period[1])

  if (length(syncon$groups)>1){ 
    stl.data.setup<-mapply(stl_data_fun,covars=stl.covars, ds.sub=syncon$.private$ds ,SIMPLIFY=FALSE, MoreArgs=list(n_seasons=syncon$n_seasons, outcome_name=syncon$outcome_name, post.start.index=post.start.index)) #list of lists that has covariates for each regression for each strata
  }else{
    stl.data.setup <- list(mapply(stl_data_fun,covars=stl.covars, ds.sub=syncon$.private$ds, MoreArgs=list(n_seasons=syncon$n_seasons, outcome_name=syncon$outcome_name, post.start.index=post.start.index)))
  }

  ##SECTION 2: run first stage models
  syncon$.private$n_cores <- detectCores() - 1
  glm.results <- vector("list", length=length(stl.data.setup)) #combine models into a list
  cl <- makeCluster(syncon$.private$n_cores)
  clusterEvalQ(cl, {library(lme4, quietly = TRUE)})
  clusterExport(cl, c('stl.data.setup', 'glm.fun', 'post.start.index'), environment())
  for(i in 1:length(stl.data.setup)){
    glm.results[[i]] <- pblapply(cl=cl, stl.data.setup[[i]], FUN=function(d) {glm.fun(d, post.start.index)})
  }
  stopCluster(cl)
  ######################

  # Combine data
  #Combine the outcome, covariates, and time point information.
  syncon$.private$data$full <- setNames(lapply(syncon$groups, makeTimeSeries, outcome=syncon$outcome, covars=syncon$covars$full), syncon$groups)
  syncon$.private$data$time <- setNames(lapply(syncon$groups, makeTimeSeries, outcome=syncon$outcome, covars=syncon$covars$time, trend=TRUE, offset=offset), syncon$groups)
  syncon$.private$data$pca <- mapply(FUN=pca_top_var, glm.results.in=glm.results, covars=stl.covars, ds.in=syncon$.private$ds, SIMPLIFY=FALSE, MoreArgs=list(outcome_name=syncon$outcome_name, season.dummies=season.dummies))
  names(syncon$.private$data$pca) <- syncon$groups
  #Time trend model but without a denominator
  syncon$.private$data$time_no_offset <- setNames(lapply(syncon$groups, makeTimeSeries, outcome=syncon$outcome, covars=syncon$covars$time, trend=FALSE), syncon$groups)
}
