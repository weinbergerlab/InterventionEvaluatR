library(R6)

syncon_factory <- R6Class(
  "SyntheticControl",
  private = list(
    # Variants
    variants = list(
      full = list(var.select.on=TRUE, trend=FALSE, name="Synthetic controls"),
      time = list(var.select.on=TRUE, trend=FALSE, name="Time trend"),
      time_no_offset = list(var.select.on=TRUE, trend=FALSE, name="Time trend (no offset)"),
      pca = list(var.select.on=TRUE, trend=FALSE, name="STL+PCA")
    ),

    # Inputs
    input_data = NA,

    group_name = NA,
    date_name = NA,
    outcome_name = NA,
    denom_name = NA,

    exclude_covar = NA,
    exclude_group = NA,

    # Computation state
    data = list(),
    n_cores = NA,

    impact.pre = function() {
      # Setup data
      prelog_data <- private$input_data[!is.na(private$input_data[, private$outcome_name]),]#If outcome is missing, delete 
      prelog_data[, private$group_name] = prelog_data[, private$group_name] %% 2
      self$groups <- as.character(unique(unlist(prelog_data[, private$group_name], use.names = FALSE)))
      if (exists('exclude_group')) {
        self$groups <- self$groups[!(self$groups %in% exclude_group)]
      }

      # Format covars
      prelog_data[,private$date_name]<-as.Date(as.character(prelog_data[,private$date_name]), tryFormats=c("%m/%d/%Y",'%Y-%m-%d' ))

      #test<-split(prelog_data, factor(prelog_data[,private$group_name]))
      #outcome.na<-sapply(test, function(x) sum(is.na(x[,private$outcome_name])))
      prelog_data[, private$date_name] <- formatDate(prelog_data[, private$date_name])
      prelog_data <- setNames(lapply(
        self$groups, FUN=splitGroup, 
        ungrouped_data=prelog_data, group_name=private$group_name, date_name=private$date_name, 
        start_date=self$start_date, end_date=self$end_date, 
        no_filter=c(private$group_name, private$date_name, private$outcome_name, private$denom_name)
      ), self$groups)
      #if (exists('exclude_group')) {prelog_data <- prelog_data[!(names(prelog_data) %in% exclude_group)]}

      #Log-transform all variables, adding 0.5 to counts of 0.
      ds <- setNames(lapply(prelog_data, FUN = logTransform, no_log = c(private$group_name, private$date_name, private$outcome_name)), self$groups)
      self$time_points <- unique(ds[[1]][, private$date_name])
      print(3)

      #Monthly dummies
      if(self$n_seasons==4){
        dt <- quarter(as.Date(self$time_points))
      }
      if(self$n_seasons==12){
        dt <- month(as.Date(self$time_points))
      }
      if(self$n_seasons==3){
        dt.m <- month(as.Date(self$time_points))
        dt <- dt.m
        dt[dt.m %in% c(1,2,3,4)] <- 1
        dt[dt.m %in% c(5,6,7,8)] <- 2
        dt[dt.m %in% c(9,10,11,12)] <- 3
      }
      season.dummies <- dummies::dummy(dt)
      season.dummies <- as.data.frame(season.dummies)
      names(season.dummies) <- paste0('s', 1:self$n_seasons)
      season.dummies <- season.dummies[,-self$n_seasons]
      print(4)

      ds <- lapply(ds, function(ds) {
        if (!(private$denom_name %in% colnames(ds))) {
          ds[private$denom_name] <- 0
        }
        return(ds)
      })
      print(5)

      self$sparse_groups <- sapply(ds, function(ds) {
        return(ncol(ds[!(colnames(ds) %in% c(private$date_name, private$group_name, private$denom_name, private$outcome_name, private$exclude_covar))]) == 0)
      })
      ds <- ds[!self$sparse_groups]
      self$groups <- self$groups[!self$sparse_groups]

      #Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
      self$covars = list()
      self$covars$full <- setNames(lapply(ds, function(group) { makeCovars(self$country, self$time_points, self$intervention_date, season.dummies, group) }), self$groups)
      self$covars$full <- lapply(self$covars$full, FUN = function(covars) {covars[, !(colnames(covars) %in% private$exclude_covar), drop = FALSE]})
      self$covars$time <- setNames(lapply(self$covars$full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies,time_index = 1:nrow(covars))))}), self$groups)
      self$covars$null <- setNames(lapply(self$covars$full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies)))}), self$groups)
      print(6)

      #Standardize the outcome variable and save the original mean and SD for later analysis.
      self$outcome <- sapply(ds, FUN = function(data) {data[, private$outcome_name]})
      offset <- sapply(ds, FUN=function(data) exp(data[, private$denom_name]) ) #offset term on original scale; 1 column per age group

      ##SECTION 1: CREATING SMOOTHED VERSIONS OF CONTROL TIME SERIES AND APPENDING THEM ONTO ORIGINAL DATAFRAME OF CONTROLS
      #EXTRACT LONG TERM TREND WITH DIFFERENT LEVELS OF SMOOTHNESS USING STL
      # Set a list of parameters for STL
      stl.covars<-mapply(smooth_func,ds.list=ds,covar.list=self$covars$full, SIMPLIFY=FALSE, MoreArgs=list(n_seasons=self$n_seasons)) 
      post.start.index<-which(self$time_points==self$post_period[1])

      if (length(self$groups)>1){ 
        stl.data.setup<-mapply(stl_data_fun,covars=stl.covars, ds.sub=ds ,SIMPLIFY=FALSE, MoreArgs=list(n_seasons=self$n_seasons, outcome_name=private$outcome_name, post.start.index=post.start.index)) #list of lists that has covariates for each regression for each strata
      }else{
        stl.data.setup <- list(mapply(stl_data_fun,covars=stl.covars, ds.sub=ds, MoreArgs=list(n_seasons=self$n_seasons, outcome_name=private$outcome_name, post.start.index=post.start.index)))
      }

      ##SECTION 2: run first stage models
      private$n_cores <- detectCores() - 1
      glm.results <- vector("list", length=length(stl.data.setup)) #combine models into a list
      cl <- makeCluster(private$n_cores)
      clusterEvalQ(cl, {library(lme4, quietly = TRUE)})
      clusterExport(cl, c('stl.data.setup', 'glm.fun', 'post.start.index'), environment())
      for(i in 1:length(stl.data.setup)){
        print(i)
        glm.results[[i]] <- parLapply(cl=cl, stl.data.setup[[i]], fun=function(d) {glm.fun(d, post.start.index)})
      }
      stopCluster(cl)
      ######################

      # Combine data
      #Combine the outcome, covariates, and time point information.
      private$data$full <- setNames(lapply(self$groups, makeTimeSeries, outcome=self$outcome, covars=self$covars$full), self$groups)
      private$data$time <- setNames(lapply(self$groups, makeTimeSeries, outcome=self$outcome, covars=self$covars$time, trend=TRUE, offset=offset), self$groups)
      private$data$pca <- mapply(FUN=pca_top_var, glm.results.in=glm.results, covars=stl.covars, ds.in=ds, SIMPLIFY=FALSE, MoreArgs=list(outcome_name=private$outcome_name, season.dummies=season.dummies))
      names(private$data$pca) <- self$groups
      #Time trend model but without a denominator
      private$data$time_no_offset <- setNames(lapply(self$groups, makeTimeSeries, outcome=self$outcome, covars=self$covars$time, trend=FALSE), self$groups)
    }
    
  ),
  public = list(
    country = NA,

    pre_period = NA,
    intervention_date = NA,
    post_period = NA,
    eval_period = NA,

    start_date = NA,
    end_date = NA,

    n_seasons = NA,
    time_points = NA,
    year_def = NA,

    groups = NA,
    sparse_groups = NA,
    model_size = NA,
    covars = NA,
    outcome = NA,

    initialize = function(
      country, data,
      pre_period_start, pre_period_end,
      post_period_start, post_period_end,
      eval_period_start, eval_period_end,
      n_seasons, year_def, 
      group_name, date_name, outcome_name, denom_name
    ) {

      self$country <- country #Country or region name.
      self$n_seasons <- n_seasons #Number of months (seasons) per year. 12 for monthly, 4 for quarterly, 3 for trimester data.
      self$year_def <- params$year_def #Can be cal_year to aggregate results by Jan-Dec; 'epi_year' to aggregate July-June

      #MOST DATES MUST BE IN FORMAT "YYYY-MM-01", exception is end of pre period, which is 1 day before end of post period
      self$pre_period <- as.Date(c(pre_period_start, pre_period_end)) #Range over which the data is trained for the CausalImpact model.
      self$start_date <- self$pre_period[1]
      self$post_period <- as.Date(c(post_period_start, post_period_end )) #Range from the intervention date to the end date.
      self$intervention_date <- self$post_period[1]-1
      self$end_date <- self$post_period[2]

      self$eval_period <- as.Date(c(eval_period_start, eval_period_end)) #Range over which rate ratio calculation will be performed.

      private$group_name <- group_name #Name of column containing group labels.
      private$date_name <- date_name #Name of column containing dates.
      private$outcome_name <- outcome_name #Name of column containing outcome.
      private$denom_name <- denom_name #Name of column containing denominator to be used in offset.

      private$exclude_covar <- c() #User-defined list of covariate columns to exclude from all analyses.
      private$exclude_group <- c() #User-defined list of groups to exclude from analyses.

      #Assign variable values
      private$input_data <- data
    },

    impact = function() {
      private$impact.pre()
      results = list()
      #Start Cluster for CausalImpact (the main analysis function).
      cl <- makeCluster(private$n_cores)
      clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE)})
      clusterExport(cl, c('doCausalImpact'), environment())

      for (variant in names(private$variants)) {
        results[[variant]]$groups <- setNames(parLapply(
          cl, private$data[[variant]], doCausalImpact, 
          self$intervention_date, 
          self$n_seasons,
          var.select.on=private$variants[[variant]]$var.select.on, 
          time_points=self$time_points,
          trend=private$variants[[variant]]$trend
        ), self$groups)
      }
      stopCluster(cl)
      
      for (variant in c('full', 'time')) {
        #Save the inclusion probabilities from each of the models
        results[[variant]]$inclusion_prob <- setNames(lapply(results[[variant]]$groups, inclusionProb), self$groups)
      }

      for (variant in names(private$variants)) {
        #All model results combined
        results[[variant]]$quantiles <- setNames(lapply(self$groups, FUN=function(group) {
          rrPredQuantiles(impact = results[[variant]]$groups[[group]], denom_data = ds[[group]][, denom_name], eval_period=self$eval_period, post_period=self$post_period, year_def=self$year_def, time_points=self$time_points, n_seasons=self$n_seasons)
        }), self$groups)
      }

      # Calculate best model
      self$model_size <- sapply(results$full$groups, modelsize_func, n_seasons=self$n_seasons)
      results$best$quantiles <- vector("list", length(results$full$quantiles)) 
      results$best$quantiles[self$model_size>=1] <- results$full$quantiles[self$model_size>=1]
      results$best$quantiles[self$model_size<1] <- results$pca$quantiles[self$model_size<1]
      results$best$quantiles <- setNames(results$best$quantiles, self$groups)

      for (variant in c("best", names(private$variants))) {
        # Predictions, aggregated by year
        results[[variant]]$pred_quantiles <- sapply(results[[variant]]$quantiles, getPred, simplify = 'array')
        results[[variant]]$ann_pred_quantiles <- sapply(results[[variant]]$quantiles, getAnnPred, simplify = FALSE)
      }

      for (variant in c('full', 'best')) {
        # Pointwise RR and uncertainty for second stage meta variant
        results[[variant]]$log_rr_quantiles <- sapply(results[[variant]]$quantiles, FUN = function(quantiles) {quantiles$log_rr_full_t_quantiles}, simplify = 'array')
        dimnames(results[[variant]]$log_rr_quantiles)[[1]] <- self$time_points
        results[[variant]]$log_rr_sd <- sapply(results[[variant]]$quantiles, FUN = function(quantiles) {quantiles$log_rr_full_t_sd}, simplify = 'array')
        results[[variant]]$log_rr_full_t_samples.prec <- sapply(results[[variant]]$quantiles, FUN = function(quantiles) {quantiles$log_rr_full_t_samples.prec}, simplify = 'array')
      }

      for (variant in c("best", names(private$variants))) {
      	# Rolling rate ratios
        results[[variant]]$rr_roll <- sapply(results[[variant]]$quantiles, FUN = function(quantiles) {quantiles$roll_rr}, simplify = 'array')
        # Rate ratios for evaluation period.
        results[[variant]]$rr_mean <- t(sapply(results[[variant]]$quantiles, getRR))
      }

      results$best$log_rr <- t(sapply(results$best$quantiles, getsdRR))
      
      for (variant in c("best", names(private$variants))) {
        results[[variant]]$rr_mean_intervals <- setNames(data.frame(makeInterval(
          results[[variant]]$rr_mean[, 2], results[[variant]]$rr_mean[, 3], results[[variant]]$rr_mean[, 1]
        ), check.names = FALSE, row.names = self$groups), c(paste(private$variants[[variant]]$name, 'Estimate (95% CI)')))
      }

      colnames(results$time$rr_mean) <- paste('Time_trend', colnames(results$time$rr_mean))

      for (variant in c("best", names(private$variants))) {
        results[[variant]]$cumsum_prevented <- sapply(self$groups, FUN=cumsum_func, quantiles = results[[variant]]$quantiles, outcome=self$outcome, self$time_points, self$post_period, simplify = 'array')
      }

      #Run a classic ITS analysis
      rr.its1 <- lapply(private$data$time, its_func, post_period=self$post_period, eval_period=self$eval_period, time_points=self$time_points)
      rr.t <- sapply(rr.its1, `[[`, "rr.q.t", simplify='array')
      results$its = list()
      results$its$rr_end <- t(sapply(rr.its1, `[[`, "rr.q.post", simplify='array')) 
      results$its$rr_mean_intervals <- data.frame('Classic ITS (95% CI)' = makeInterval(results$its$rr_end[, 2], results$its$rr_end[, 3], results$its$rr_end[, 1]), check.names = FALSE, row.names = self$groups)

      results
    },
    crossval = function() {
      results = list()
      #Creates List of lists: 1 entry for each stratum; within this, there are CV datasets for each year left out, and within this, there are 2 lists, one with full dataset, and one with the CV dataset
      for (variant in names(private$variants)) {
        private$data.cv[[variant]]<-lapply(private$data[[variant]], makeCV, self$time_points, self$intervention_date)
      }

      #Run the models on each of these datasets
      cl <- makeCluster(private$n_cores)
      clusterEvalQ(cl, {
        library(pogit, quietly = TRUE); 
        library(lubridate, quietly = TRUE)
      })
      clusterExport(cl, c('doCausalImpact'), environment())
      for (variant in variants) {
        results$impact[[variant]] <-setNames(parLapply(
          cl, private$data.cv[[variant]], function(x) lapply(
            x, doCausalImpact, 
            self$intervention_date, 
            self$n_seasons,
            time_points=self$time_points,
            crossval.stage=TRUE,
            var.select.on=private$variants[[variant]]$var.select.on,
          )), self$groups)
      }
      stopCluster(cl)
      
      ll.cv = list()

      #Calculate pointwise log likelihood for cross-val prediction sample vs observed
      #These are N_iter*N_obs*N_cross_val array
      for (variant in variants) {
        ll <- lapply(results$impact[[variant]], function(x) lapply(x, crossval.log.lik))
        ll.cv[[variant]] <- lapply(ll.cv, reshape.arr)
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
      names(stacking_weights.all) <- lapply(private$variants, function(v) { v$name })
      stacking_weights.all <- cbind.data.frame(self$groups, stacking_weights.all)
      self$stacking_weights.all.m <- melt(stacking_weights.all, id.vars='self$groups')
      # stacking_weights.all.m<-stacking_weights.all.m[order(stacking_weights.all.m$groups),]
      
      stacked.ests <- mapply(
        FUN=stack.mean,
        outcome=self$outcome,
        group=self$groups,
        impact_full=self$results$impact$full,
        impact_time=self$results$impact$time,
        impact_time_no_offset=self$results$impact$time_no_offset,
        impact_pca=self$results$impact$pca,
        MoreArgs=list(stacking_weights.all=stacking_weights.all), 
        SIMPLIFY=FALSE
      )
      # plot.stacked.ests<-lapply(stacked.ests,plot.stack.est)
      results$quantiles_stack <- setNames(lapply(self$groups, FUN = function(group) {
        rrPredQuantiles(impact = stacked.ests[[group]], denom_data = ds[[group]][, denom_name], eval_period, post_period, self$year_def, self$time_points)
      }), self$groups)
      results$pred_quantiles_stack <- sapply(results$quantiles_stack, getPred, simplify = 'array')
      results$rr_roll_stack <- sapply(results$quantiles_stack, FUN = function(quantiles_stack) {quantiles_stack$roll_rr}, simplify = 'array')
      results$rr_mean_stack <- round(t(sapply(results$quantiles_stack, getRR)), 2)
      results$rr_mean_stack_intervals <- data.frame('Stacking Estimate (95% CI)' = makeInterval(results$rr_mean_stack[, 2], results$rr_mean_stack[, 3], results$rr_mean_stack[, 1]), check.names = FALSE, row.names = self$groups)
      results$cumsum_prevented_stack <- sapply(self$groups, FUN = cumsum_func, quantiles = results$quantiles_stack, outcome=self$outcome, self$time_points, self$post_period, simplify = 'array')
      results$ann_pred_quantiles_stack <- sapply(results$quantiles_stack, getAnnPred, simplify = FALSE)
      #Preds: Compare observed and expected
      results$pred$full <- lapply(results$impact$full, function(x) sapply(x,pred.cv,simplify='array'))
      results$pred$pca <- lapply(results$impact$pca, function(x) sapply(x,pred.cv,simplify='array'))
      
      results
    },
    sensitivity = function() {
      results = list()
      bad_sensitivity_groups <- sapply(self$covars$full, function (covar) {ncol(covar) <= n_seasons-1+3})
      sensitivity_covars_full <- self$covars$full[!bad_sensitivity_groups]
      sensitivity_ds <- ds[!bad_sensitivity_groups]
      sensitivity_impact_full <- self$impact_full[!bad_sensitivity_groups]
      sensitivity_groups <- self$groups[!bad_sensitivity_groups]
      
      if (length(sensitivity_groups)!=0) {
        #Weight Sensitivity Analysis - top weighted variables are excluded and analysis is re-run.
        cl <- makeCluster(private$n_cores)
        clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
        clusterExport(cl, c('sensitivity_ds', 'weightSensitivityAnalysis', 'sensitivity_groups', 'outcome', 'time_points', 'n_seasons', 'eval_period', 'post_period', 'rrPredQuantiles'), environment())
        sensitivity_analysis_full <- setNames(parLapply(cl, sensitivity_groups, weightSensitivityAnalysis, covars = sensitivity_covars_full, ds = sensitivity_ds, impact = sensitivity_impact_full, time_points = self$time_points, intervention_date = self$intervention_date, n_seasons = n_seasons, outcome = outcome, eval_period = eval_period, post_period = post_period, year_def=self$year_def), sensitivity_groups)
        stopCluster(cl)
        
        results$sensitivity_pred_quantiles <- lapply(sensitivity_analysis_full, FUN = function(sensitivity_analysis) {
          pred_list <- vector(mode = 'list', length = length(sensitivity_analysis))
          for (sensitivity_index in 1:length(sensitivity_analysis)) {
            pred_list[[sensitivity_index]] <- getPred(sensitivity_analysis[[sensitivity_index]])
          }
          return(pred_list)
        })
        
        #Table of rate ratios for each sensitivity analysis level
        results$sensitivity_table <- t(sapply(sensitivity_groups, sensitivityTable, sensitivity_analysis = sensitivity_analysis_full, original_rr = rr_mean_full))
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
        reults$rr_table <- cbind.data.frame(round(rr_mean_time[!bad_sensitivity_groups, ],2), results$sensitivity_table)
        results$rr_table_intervals <- cbind('Trend Estimate (95% CI)' = rr_mean_time_intervals[!bad_sensitivity_groups, ], results$sensitivity_table_intervals)
      } else {
        results$sensitivity_table_intervals <- NA
      }
    },
    generate_plots = function() {
      sc_plots(self)
    }
  )
)
