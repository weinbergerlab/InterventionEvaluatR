#' Initialize analysis
#'
#' @param country A one-word label for output (eg country or state name).
#' @param data Input dataframe. There should be a row for each time point (e.g., month) and category (e.g., age group). There should be variables for the date, for the category (or a column of 1s if only 1 category), for the outcome variable (a count), for the denominator (or a column of 1s if no denominator), and columns for each control variable. 
#' @param pre_period_start Date when analysis starts, YYYY-MM-01. defaults to first date in dataset.
#' @param post_period_start Month when intervention introduced. YYY-MM-01
#' @param set.burnN Number of MCMC iterations for burn in (default 5000),
#' @param set.sampleN Number of MCMC iterations post-burn-in to use for inference (default 10000),
#' @param post_period_end Date when analysis ends, YYYY-MM-01. defaults to first date in dataset. Defaults to end of evaluation period.
#' @param eval_period_start First month of the period when the effect of intervention is evaluated. YYYY-MM-01. typically 12-24 months after post_period_start.
#' @param eval_period_end Last month of the period when the effect of intervention is evaluated. YYYY-MM-01. 
#' @param n_seasons How many observations per year? Defaults to 12 (monthly data) Change to 4 for quarterly
#' @param year_def Should results be aggregated by calendar year ('cal_year': the default) or epidemiological year ('epi_year'; July-June)
#' @param group_name Name of the stratification variable (e.g., age group). If only one age group present, add a column of 1s to the dataset
#' @param log.covars Should the covariate be log transformed? (default: TRUE)
#' @param date_name Name of the variable with the date for the time series
#' @param outcome_name Name of the outcome (y) variable in the 'data' dataframe. Should be a count
#' @param denom_name Name of the denominator variable in the 'data' dataframe. if there is no denominator, include a column of 1s.
#' @param ridge Run ridge regression with AR(1) random intercepts (faster) or spike and slab (with iid random intercept) for variable selection. Logical, Default TRUE. 
#' @param error_dist For the INLA models: use an 'iid' or 'ar1' error on the random intercept. Defaults to iid. Use caution with AR(1) as it can introduce bias in some situations 
#' @param sparse_threshold Threshold for filtering out control variables based on sparsity (mean number of cases per time period). Defaults to 5. 
#' @return Initialized analysis object, `analysis` as described below
#'
#' `analysis$country` as passed to `country`
#' 
#' `analysis$input_data` as passed to `data`
#' 
#' `analysis$n_seasons` as passed to `n_seasons`
#' 
#' `analysis$year_def` as passed to `year_def`
#' 
#' `analysis$pre_period` Range of dates in the pre-intervention period
#' 
#' `analysis$post_period` Range of dates in the post-intervention period
#' 
#' `analysis$eval_period` Range of dates in the evaluation period
#' 
#' `analysis$start_date` First date of the pre-intervention period
#' 
#' `analysis$intervention_date` Last time point before the start of the post-period
#' 
#' `analysis$end_date` Last date in the evaluation period
#' 
#' `analysis$group_name` as passed to `group_name`
#' 
#' `analysis$date_name` as passed in in `date_name`
#' 
#' `analysis$outcome_name` as passed in in `outcome_name`
#' 
#' `analysis$denom_name` as passed in in `denom_name`
#' 
#' `analysis$time_points` Vector of time points in the dataset
#' 
#' `analysis$set.burnN` as passed in in `set.burnN`
#'
#' `analysis$set.sampleN` as passed in in `set.sampleN`
#' 
#' `analysis$log.covars` as passed in in `log.covars`
#'   
#' `analysis$groups` Vector of groups analyzed
#' 
#' `analysis$sparse_groups` Vector indicating which groups were too sparse to analyze
#' 
#' `analysis$model_size` Average number of covariates included in the synthetic control model
#' 
#' `analysis$covars` Matrix of covariates used for analysis
#' 
#' `analysis$outcome` as passeed to `outcome_name`
#' 
#' `analysis$ridge` as passeed to `ridge`
#' `analysis$error_dist` as passeed to `error_dist`
#' `analysis$sparse_threshold` as passed to `sparse_threshold`
#'
#' @importFrom listenv listenv
#' @importFrom lubridate is.Date
#' @export

evaluatr.init <- function(country,
                        data,
                        pre_period_start='start',
                        post_period_start,
                        post_period_end=eval_period_end,
                        eval_period_start,
                        eval_period_end,
                        n_seasons=12,
                        year_def='cal_year',
                        group_name,
                        date_name,
                        outcome_name,
                        set.burnN=5000,
                        set.sampleN=10000,
                        denom_name,
                        log.covars=TRUE,
                        ridge=F,
                        error_dist='iid',
                        sparse_threshold = 5) {
  analysis = listenv(
    time_points = NA,
    
    groups = NA,
    sparse_groups = NA,
    model_size = NA,
    covars = NA,
    outcome = NA,
    sparse_threshold = sparse_threshold,
    
    results = list(
      impact = NA,
      crossval = NA,
      sensitivity = NA,
      univariate = NA
    ),
    
    stacking_weights.all.m = NA,
    
    .private = listenv(
      # Variants
      variants = list(
        full = list(
          var.select.on = TRUE,
          trend = FALSE,
          name = "Synthetic controls"
        ),
        time = list(
          var.select.on = FALSE,
          trend = TRUE,
          name = "Time trend"
        ),
        time_no_offset = list(
          var.select.on = FALSE,
          trend = FALSE,
          name = "Time trend (no offset)"
        ),
        pca = list(
          var.select.on = FALSE,
          trend = FALSE,
          name = "STL+PCA"
        )
      ),
      
      exclude_covar = NA,
      exclude_group = NA,
      
      # Computation state
      data = list(),
      data.cv = list(),
      ds = NA,
      cluster = NULL
    )
  )
  
  analysis$country <- country #Country or region name.
  analysis$n_seasons <-
    n_seasons #Number of months (seasons) per year. 12 for monthly, 4 for quarterly, 3 for trimester data.
  analysis$year_def <-
    year_def #Can be cal_year to aggregate results by Jan-Dec; 'epi_year' to aggregate July-June
  analysis$set.burnN <-set.burnN
  analysis$set.sampleN <-set.sampleN
  analysis$log.covars <- log.covars
  analysis$ridge <- ridge
  analysis$error_dist<-error_dist
    normalizeDate <- function(d) {
    if (is.Date(d)) {
      d
    } else {
      as.Date(as.character(d), tryFormats = c("%m/%d/%Y", '%Y-%m-%d','%Y/%m/%d'))
    }
  }
  

  # Parse date formats of various inputs  
  data[, date_name] <- normalizeDate(data[, date_name])
  post_period_start <- normalizeDate(post_period_start)
  post_period_end <- normalizeDate(post_period_end)
  eval_period_start <- normalizeDate(eval_period_start)
  eval_period_end <- normalizeDate(eval_period_end)
  

  # Identify various significant dates  
  if(!is.Date(pre_period_start) && pre_period_start=='start'){
    first.date.data<-min(data[,date_name])
  } else {
    first.date.data<-normalizeDate(pre_period_start)
  }
  
  #MOST DATES MUST BE IN FORMAT "YYYY-MM-01", exception is end of pre period, which is 1 day before end of post period
  analysis$post_period <-
    as.Date(c(post_period_start, post_period_end)) #Range from the intervention date to the end date.
  analysis$pre_period_end<- 
    analysis$post_period[1] - 1
  analysis$pre_period <-
    as.Date(c(first.date.data, analysis$pre_period_end)) #Range over which the data is trained for the CausalImpact model.
  analysis$start_date <- analysis$pre_period[1]
  analysis$intervention_date <- analysis$post_period[1] - 1
  analysis$end_date <- analysis$post_period[2]
  
  analysis$eval_period <-
    as.Date(c(eval_period_start, eval_period_end)) #Range over which rate ratio calculation will be performed.
  
  analysis$group_name <-
    group_name #Name of column containing group labels.
  analysis$date_name <- date_name #Name of column containing dates.
  analysis$outcome_name <-
    outcome_name #Name of column containing outcome.
  analysis$denom_name <-
    denom_name #Name of column containing denominator to be used in offset.
  
  analysis$.private$exclude_covar <-
    c() #User-defined list of covariate columns to exclude from all analyses.
  analysis$.private$exclude_group <-
    c() #User-defined list of groups to exclude from analyses.
  
  #Assign variable values
  analysis$input_data <- data
  
  # Setup default progress reporting
  analysis$.private$progress_parts = list()
  analysis$.private$progress_done = 0
  analysis$.private$progress = showProgress
  
  # Setup default cluster
  analysis$.private$startCluster = defaultStartCluster
  analysis$.private$stopCluster = defaultStopCluster
  return(analysis)
}

# This is used by the web UI to set up parallel computation 
evaluatr.initParallel = function(analysis, startCluster, stopCluster, progress) {
  if (!is.null(startCluster)) {
    analysis$.private$startCluster = startCluster
    analysis$.private$stopCluster = stopCluster
  } else {
    analysis$.private$startCluster = defaultStartCluster
    analysis$.private$stopCluster = defaultStopCluster
  }
  analysis$.private$progress = progress
}

#' Perform impact analysis
#'
#' @param analysis Analysis object, initialized by evaluatr.init.
#' @param variants List of one or more of the following analysis types to perform: 'time', 'time_no_offset', 'pca', and 'full'. Defaults to all 4.
#' @return Analysis results, `results`, as described below
#'
#' `results$full` Results from synthetic controls model
#' 
#' `results$time` Results from time trend model (with offset term)
#' 
#' `results$time_no_offset` Results from time trend model (no offset term)
#' 
#' `results$pca` Results from STL+PCA model
#' 
#' `results$its` Results from classic interrupted time series model
#' 
#' `results$best` Results from 'best model (either synthetic controls or STL+PCA, depending on )
#' 
#' `results$point.weights` TODO
#' 
#' `results$rr_mean_combo` Rate ratios during the evaluation period, as calculated in each of the models.
#'
#' @importFrom stats AIC as.formula cov dpois glm median poisson prcomp predict quantile rmultinom rnorm rpois sd setNames stl var vcov complete.cases
#' @importFrom loo stacking_weights
#' @importFrom lme4 glmer glmerControl fixef
#' @importFrom lubridate as_date %m+% year month quarter time_length interval yday
#' @importFrom splines bs
#' @importFrom pomp logmeanexp
#' @importFrom reshape melt
#' @importFrom MASS mvrnorm
#' @importFrom HDInterval hdi
#' @importFrom RcppRoll roll_sum
#' @importFrom plyr rbind.fill
#' @importFrom pogit poissonBvs
#' @importFrom parallel makeCluster clusterEvalQ clusterExport stopCluster parLapply
#' @importFrom future availableCores
#' @importFrom coda geweke.diag mcmc
#' @importFrom stats pnorm
#' @export

evaluatr.impact = function(analysis, variants=names(analysis$.private$variants)) {
  addProgress(analysis, sprintf("Impact analysis (%s)", lapply(analysis$.private$variants, function(variant) variant$name)))
  evaluatr.impact.pre(analysis, run.stl= ('pca' %in% variants))
  results1 = list()
  
  #Start Cluster for CausalImpact (the main analysis function).
  clusterEvalQ(cluster(analysis), {
    library(pogit, quietly = TRUE)
    library(lubridate, quietly = TRUE)
    library(RcppRoll, quietly = TRUE)
    library(HDInterval, quietly = TRUE)
    library(plyr, quietly = TRUE)
    library(INLA, quietly = TRUE)
  })
  clusterExport(cluster(analysis), c('doCausalImpact','inla_mods','rrPredQuantiles','cumsum_func','analysis'), environment())
  
  analysis$.private$variants = analysis$.private$variants[variants]
  if(analysis$ridge==F){
    for (variant in variants) {
      progressStartPart(analysis)
      results1[[variant]]$groups <- setNames(
        parLapply(
          cl = cluster(analysis),
          analysis$.private$data[[variant]],
          fun = doCausalImpact,
          intervention_date=analysis$intervention_date,
          n_seasons=analysis$n_seasons,
          var.select.on = analysis$.private$variants[[variant]]$var.select.on,
          time_points = analysis$time_points,
          trend = analysis$.private$variants[[variant]]$trend,
          burnN=analysis$set.burnN,
          sampleN=analysis$set.sampleN,
          crossval.stage = FALSE,
          analysis=analysis
        ),
        analysis$groups
      )
      progressEndPart(analysis)
    }
  }else{
    for (variant in variants) {
      progressStartPart(analysis)
      results1[[variant]]$groups <- setNames(
        pblapply(
          cl = cl,
          analysis$.private$data[['full']],
          FUN = inla_mods,
          model.variant=variant
        ),
        analysis$groups
      )
      progressEndPart(analysis)
    }
  }  

  stopCluster(analysis)

#    return(results1)
#}
#part2<-function(ds){
  results <- vector("list", length(variants))
  names(results)<-variants
  quantiles<-vector("list", length(variants))
  names(quantiles)<-variants
  cumsum1<-vector("list", length(variants))
  names(cumsum1)<-variants  
 # cumsum1.hdi<-vector("list", length(variants))
  #names(cumsum1.hdi)<-variants
  
  for (variant in variants) {
    results[[variant]]$groups<-sapply(results1[[variant]]$groups,function(x) x[['impact']], simplify=F)  
    quantiles[[variant]]<-sapply(results1[[variant]]$groups,function(x) x$quantiles, simplify=F)  
    results[[variant]]$cumsum_prevented <-sapply(results1[[variant]]$groups,function(x) x$cumsum_prevented, simplify='array') 
    results[[variant]]$cumsum_prevented_hdi<-sapply(results1[[variant]]$groups,function(x) x$cumsum_prevented_hdi, simplify='array') 
    results[[variant]]$quantiles <- quantiles[[variant]]
  }
  
  clusterUpdateAnalysis(analysis, function(analysis) {
    for (variant in intersect(c('full', 'time'), variants)) {
      #Save the inclusion probabilities from each of the models
      results[[variant]]$inclusion_prob <-
        setNames(lapply(results[[variant]]$groups, inclusionProb), analysis$groups)
    }
  
    # Calculate best model
    if(analysis$ridge==F){
      if ("full" %in% variants) {
        analysis$model_size <-
          sapply(results$full$groups, modelsize_func, n_seasons = analysis$n_seasons)
      }
    }else{
      analysis$model_size <- NA
    }
  
    if (all(c("full", "pca") %in% variants)) {
      results$best$quantiles <-
        vector("list", length(results$full$quantiles))
      results$best$quantiles[analysis$model_size >= 1] <-
        results$full$quantiles[analysis$model_size >= 1]
      results$best$quantiles[analysis$model_size < 1] <-
        results$pca$quantiles[analysis$model_size < 1]
      results$best$quantiles <-
        setNames(results$best$quantiles, analysis$groups)

      results$best$variant <-
        vector("list", length(results$full$quantiles))
      results$best$variant[analysis$model_size >= 1] <- "full"
      results$best$variant[analysis$model_size < 1] <- "pca"
      results$best$variant <-
        setNames(results$best$variant, analysis$groups)

      variants = c("best", variants)
    }
  
    # if(analysis$model_size >= 1){
      results$best$cumsum_prevented<- results$full$cumsum_prevented
    # }else{
    #   results$best$cumsum_prevented<- results$pca$cumsum_prevented
    # }

    for (variant in variants) {
      # Predictions, aggregated by year
      results[[variant]]$pred_quantiles <-
        sapply(results[[variant]]$quantiles, getPred, simplify = 'array')
      results[[variant]]$pred_quantiles_HDI <-
        sapply(results[[variant]]$quantiles, getPredHDI, simplify = 'array')
      results[[variant]]$ann_pred_quantiles <-
        sapply(results[[variant]]$quantiles, getAnnPred, simplify = FALSE)
      results[[variant]]$ann_pred_HDI <-
        sapply(results[[variant]]$quantiles, getAnnPredHDI, simplify = FALSE)
    }

    for (variant in intersect(c('full', 'best'), variants)) {
      # Pointwise RR and uncertainty for second stage meta variant
      results[[variant]]$log_rr_quantiles <-
        sapply(
          results[[variant]]$quantiles,
          FUN = function(quantiles) {
            quantiles$log_rr_full_t_quantiles
          },
          simplify = 'array'
        )
      dimnames(results[[variant]]$log_rr_quantiles)[[1]] <-
        analysis$time_points
      results[[variant]]$log_rr_sd <-
        sapply(
          results[[variant]]$quantiles,
          FUN = function(quantiles) {
            quantiles$log_rr_full_t_sd
          },
          simplify = 'array'
        )
    
      results[[variant]]$log_rr_hdi <-
        sapply(
          results[[variant]]$quantiles,
          FUN = function(quantiles) {
            quantiles$log_rr_full_t_hdi
          },
          simplify = 'array'
        )
      dimnames(results[[variant]]$log_rr_hdi)[[1]] <-
        analysis$time_points
    
      results[[variant]]$log_rr_full_t_samples.prec <-
        sapply(
          results[[variant]]$quantiles,
          FUN = function(quantiles) {
            quantiles$log_rr_full_t_samples.prec
          },
          simplify = 'array'
        )
    }
  
    for (variant in variants) {
      # Rolling rate ratios
      results[[variant]]$rr_roll <-
        sapply(
          results[[variant]]$quantiles,
          FUN = function(quantiles) {
            quantiles$roll_rr
          },
          simplify = 'array'
        )
      # Rate ratios for evaluation period.
      results[[variant]]$rr_mean <-
        t(sapply(results[[variant]]$quantiles, getRR))
      results[[variant]]$rr_iter <-
        t(sapply(results[[variant]]$quantiles, getRRiter))
      results[[variant]]$rr_mean_hdi <-
        t(sapply(results[[variant]]$quantiles, getRRHDI))
    
      #Convergence status
      trace1<-results[[variant]]$rr_iter
      con.stat<-matrix(NA, nrow=nrow(trace1), ncol=2)
      colnames(con.stat)<- c('geweke.p','status')
      for(i in 1: nrow(trace1)){
        geweke.p<- pnorm(abs(geweke.diag(mcmc(trace1[i,]))$z),lower.tail=FALSE)*2
        con.stat[i,1]<-geweke.p
        if(geweke.p>0.05){
          con.stat[i,2]<-'Model converged'
        }else{
          con.stat[i,2]<-'Not converged'
        } 
      }
        results[[variant]]$converge<-con.stat 
      }
  
    if ('best' %in% variants) {
      results$best$log_rr <- t(sapply(results$best$quantiles, getsdRR))
    }
  
    for (variant in variants) {
      results[[variant]]$rr_mean_intervals <-
        setNames(
          data.frame(
            makeInterval(results[[variant]]$rr_mean[, 2], results[[variant]]$rr_mean[, 3], results[[variant]]$rr_mean[, 1]),
            check.names = FALSE,
            row.names = analysis$groups
          ),
          c(
            paste(analysis$.private$variants[[variant]]$name, 'Estimate (95% CI)')
          )
        )
    }

    if ('time' %in% variants) {
      colnames(results$time$rr_mean) <-
        paste('Time_trend', colnames(results$time$rr_mean))
    }
  
    #Run a classic ITS analysis
    rr.its1 <-
      lapply(
        analysis$.private$data$time,
        its_func,
        post_period = analysis$post_period,
        eval_period = analysis$eval_period,
        time_points = analysis$time_points
      )
    rr.t <- sapply(rr.its1, `[[`, "rr.q.t", simplify = 'array')
    results$its = list()
    results$its$rr_end <-
      t(sapply(rr.its1, `[[`, "rr.q.post", simplify = 'array'))
    results$its$rr_mean_intervals <-
      data.frame(
        'Classic ITS (95% CI)' = makeInterval(
          results$its$rr_end[, 2],
          results$its$rr_end[, 3],
          results$its$rr_end[, 1]
        ),
        check.names = FALSE,
        row.names = analysis$groups
      )

    #Combine RRs into 1 for ease of plotting
    results$rr_mean_combo <- 
      do.call(rbind, lapply(seq_along(names(analysis$.private$variants)), function(idx) {
        variant = names(analysis$.private$variants)[[idx]]
        setNames(
          cbind(data.frame(
            Model=rep(idx, nrow(results[[variant]]$rr_mean)),
            Model_tag=rep(variant, nrow(results[[variant]]$rr_mean)),
            groups=analysis$groups,
            group.index=seq(
              from = 1,
              by = 1,
              length.out = nrow(results[[variant]]$rr_mean)
            )
          ), results[[variant]]$rr_mean), 
          c('Model', 'Model_tag', 'groups', 'group.index', 'lcl', 'mean.rr', 'ucl')
        )
      }))

    results$point.weights <-
      as.data.frame(matrix(rep(1, nrow(
        results$rr_mean_combo
      )), ncol = 1))
    names(results$point.weights) <- 'value'
  
    results$rr_mean_combo$group.index <-
      as.numeric(as.character(results$rr_mean_combo$group.index))
    results$rr_mean_combo$mean.rr <-
      as.numeric(as.character(results$rr_mean_combo$mean.rr))
    results$rr_mean_combo$lcl <-
      as.numeric(as.character(results$rr_mean_combo$lcl))
    results$rr_mean_combo$ucl <-
      as.numeric(as.character(results$rr_mean_combo$ucl))
    results$rr_mean_combo$group.index <- results$rr_mean_combo$group.index + (results$rr_mean_combo$Model - 1) * 0.15 # <-- TODO I suspect this is a kludge for plotting and should probably not be here at all
    results$rr_mean_combo$Model <- as.factor(unlist(lapply(analysis$.private$variants[results$rr_mean_combo$Model_tag], function(variant) variant$name)))
    results$rr_mean_combo$est.index <-
      as.factor(1:nrow(results$rr_mean_combo))
    #Fix order for axis
    results$rr_mean_combo$Model <-
      factor(
        results$rr_mean_combo$Model,
        lapply(
          analysis$.private$variants[c('time', 'time_no_offset', 'pca', 'full')], 
          function(variant) variant$name
        )
      )
    #print(levels(rr_mean_combo$Model))
  
  
    analysis$results$impact <- results
    analysis = evaluatr.prune(analysis, what="impact")
    return(analysis)
  })
  
  return(analysis$results$impact)
}

#' Perform cross-validation
#'
#' @param analysis Analysis object, initialized by TODO.init. You must call TODO.impact before calling TODO.sensitivity
#' @return Cross-validation results, `results`, as described below
#'
#' `results$full` TODO
#' `results$time` TODO
#' `results$time_no_offset` TODO
#' `results$pca` TODO
#' `results$ann_pred_quantiles_stack` TODO
#' `results$cumsum_prevented_stack` TODO
#' `results$log_rr_quantiles_stack` TODO
#' `results$log_rr_samples.prec.post_stack` TODO
#' `results$point.weights` TODO
#' `results$pred_quantiles_stack` TODO
#' `results$quantiles_stack` TODO
#' `results$rr_mean_combo` TODO
#' `results$rr_mean_stack` TODO
#' `results$rr_mean_stack_intervals` TODO
#' `results$rr_roll_stack` TODO
#' `results$stacking_weights` TODO
#' `results$stacking_weights.all` TODO
#' `results$stacking_weights.all.m` TODO
#'
#' @export

evaluatr.crossval = function(analysis) {
  addProgress(analysis, sprintf("Cross-validation (%s)", lapply(analysis$.private$variants, function(variant) variant$name)))
  results = list()
  
  #Creates List of lists: 1 entry for each stratum; within this, there are CV datasets for each year left out, and within this, there are 2 lists, one with full dataset, and one with the CV dataset
  for (variant in names(analysis$.private$variants)) {
    analysis$.private$data.cv[[variant]] <-
      lapply(
        analysis$.private$data[[variant]],
        makeCV,
        analysis$time_points,
        analysis$intervention_date
      )
  }
  
  #Run the models on each of these datasets
  clusterEvalQ(cluster(analysis), {
    library(pogit, quietly = TRUE)
    
    library(lubridate, quietly = TRUE)
  })
  clusterExport(cluster(analysis), c('doCausalImpact'), environment())
  for (variant in names(analysis$.private$variants)) {
    progressStartPart(analysis)
    results[[variant]]$groups <- setNames(
      parLapply(
        cl = cluster(analysis),
        analysis$.private$data.cv[[variant]],
        fun = function(x)
          lapply(
            x,
            doCausalImpact,
            analysis$intervention_date,
            analysis$n_seasons,
            time_points = analysis$time_points,
            crossval.stage = TRUE,
            burnN=analysis$set.burnN,
            sampleN=analysis$set.sampleN,
            var.select.on = analysis$.private$variants[[variant]]$var.select.on
          )
      ),
      analysis$groups
    )
    progressEndPart(analysis)
  }
  stopCluster(analysis)
  
  ll.cv = list()
  
  #Calculate pointwise log likelihood for cross-val prediction sample vs observed
  #These are N_iter*N_obs*N_cross_val array
  for (variant in names(analysis$.private$variants)) {
    ll.cv[[variant]] <-
      lapply(results[[variant]]$groups, function(x)
        lapply(x, crossval.log.lik))
    ll.cv[[variant]] <- lapply(ll.cv[[variant]], reshape.arr)
  }
  
  #Create list that has model result for each stratum
  ll.compare <- vector("list", length(ll.cv$pca))
  results$stacking_weights.all <-
    matrix(NA, nrow = length(ll.cv$pca), ncol = 4)
  
  for (i in 1:length(ll.compare)) {
    #will get NAs if one of covariates is constant in fitting period (ie pandemic flu dummy)...should fix this above
    ll.compare[[i]] <-
      cbind(ll.cv$full[[i]],
            ll.cv$time_no_offset[[i]],
            ll.cv$time[[i]],
            ll.cv$pca[[i]])
    keep <- complete.cases(ll.compare[[i]])
    ll.compare[[i]] <- ll.compare[[i]][keep, ]
    #occasionally if there is a very poor fit, likelihood is very very small, which leads to underflow issue and log(0)...delete these rows to avoid this as a dirty solution. Better would be to fix underflow
    row.min <- apply(exp(ll.compare[[i]]), 1, min)
    ll.compare[[i]] <- ll.compare[[i]][!(row.min == 0), ]
    #if(min(exp(ll.compare[[i]]))>0){
    results$stacking_weights.all[i, ] <-
      stacking_weights(ll.compare[[i]])
    #}
  }
  results$stacking_weights.all <-
    as.data.frame(round(results$stacking_weights.all, 3))
  names(results$stacking_weights.all) <-
    lapply(analysis$.private$variants, function(v) {
      v$name
    })
  results$stacking_weights.all <-
    cbind.data.frame(data.frame(groups = analysis$groups),
                     results$stacking_weights.all)
  results$stacking_weights.all.m <-
    melt(results$stacking_weights.all, id.vars = 'groups')
  # results$stacking_weights.all.m<-results$stacking_weights.all.m[order(results$stacking_weights.all.m$groups),]
  
  stacked.ests <- mapply(
    FUN = stack.mean,
    group = analysis$groups,
    impact_full = analysis$results$impact$full$groups,
    impact_time = analysis$results$impact$time$groups,
    impact_time_no_offset = analysis$results$impact$time_no_offset$groups,
    impact_pca = analysis$results$impact$pca$groups,
    MoreArgs = list(
      stacking_weights.all = results$stacking_weights.all,
      outcome = analysis$outcome
    ),
    SIMPLIFY = FALSE
  )
  # plot.stacked.ests<-lapply(stacked.ests,plot.stack.est)
  results$quantiles_stack <-
    setNames(lapply(
      analysis$groups,
      FUN = function(group) {
        rrPredQuantiles(
          impact = stacked.ests[[group]],
          denom_data = analysis$.private$ds[[group]][, analysis$denom_name],
          analysis$eval_period,
          analysis$post_period,
          analysis$n_seasons,
          analysis$year_def,
          analysis$time_points
        )
      }
    ), analysis$groups)
  results$pred_quantiles_stack <-
    sapply(results$quantiles_stack, getPred, simplify = 'array')
  results$rr_roll_stack <-
    sapply(
      results$quantiles_stack,
      FUN = function(quantiles_stack) {
        quantiles_stack$roll_rr
      },
      simplify = 'array'
    )
  results$rr_mean_stack <-
    round(t(sapply(results$quantiles_stack, getRR)), 2)
  results$rr_mean_stack_intervals <-
    data.frame(
      'Stacking Estimate (95% CI)' = makeInterval(
        results$rr_mean_stack[, 2],
        results$rr_mean_stack[, 3],
        results$rr_mean_stack[, 1]
      ),
      check.names = FALSE,
      row.names = analysis$groups
    )
  results$cumsum_prevented_stack <-
    sapply(
      analysis$groups,
      FUN = cumsum_func,
      quantiles = results$quantiles_stack,
      outcome = analysis$outcome,
      analysis$time_points,
      analysis$post_period,
      simplify = 'array'
    )
  results$ann_pred_quantiles_stack <-
    sapply(results$quantiles_stack, getAnnPred, simplify = FALSE)
  #Preds: Compare observed and expected
  results$full$pred <-
    lapply(results$impact$full, function(x)
      sapply(x, pred.cv, simplify = 'array'))
  results$pca$pred <-
    lapply(results$impact$pca, function(x)
      sapply(x, pred.cv, simplify = 'array'))
  
  results$log_rr_quantiles_stack <-
    sapply(
      results$quantiles_stack,
      FUN = function(quantiles) {
        quantiles$log_rr_full_t_quantiles
      },
      simplify = 'array'
    )
  dimnames(results$log_rr_quantiles_stack)[[1]] <-
    analysis$time_points
  
  results$log_rr_samples.prec.post_stack <-
    sapply(
      results$quantiles_stack,
      FUN = function(quantiles) {
        quantiles$log_rr_full_t_samples.prec.post
      },
      simplify = 'array'
    )
  
  results$rr_mean_combo = analysis$results$impact$rr_mean_combo
  results$point.weights <- results$stacking_weights.all.m
  
  analysis$results$crossval <- results
  return(results)
}

#' Perform sensitivity analysis
#'
#' @param analysis Analysis object, initialized by TODO.init. You must call TODO.impact before calling TODO.sensitivity
#' @return Sensitivity analysis results, `results`, as described below
#'
#' `results$rr_table` TODO
#' `results$rr_table_intervals` TODO
#' `results$sensitivity_pred_quantiles` Fitted values and credible intervals from models where the top covariates were dropped
#' `results$sensitivity_table` Unformatted matrix with rate ratio estimates from the original model and those with top covariates dropped
#' `results$sensitivity_table_intervals` Data frame containing rate ratio estimates from the original synthetic controls model
#'        as well as models where the top-weighted 1,2, or 3 covariates were dropped.
#'
#' @export

evaluatr.sensitivity = function(analysis) {
  results = list()
  bad_sensitivity_groups <-
    sapply(analysis$covars$full, function (covar) {
      ncol(covar) <= analysis$n_seasons - 1 + 3
    })
  sensitivity_covars_full <-
    analysis$covars$full[!bad_sensitivity_groups]
  sensitivity_ds <- analysis$.private$ds[!bad_sensitivity_groups]
  sensitivity_impact_full <-
    analysis$results$impact$full$groups[!bad_sensitivity_groups]
  sensitivity_groups <- analysis$groups[!bad_sensitivity_groups]
  
  if (length(sensitivity_groups) != 0) {
    #Weight Sensitivity Analysis - top weighted variables are excluded and analysis is re-run.
    clusterEvalQ(cluster(analysis), {
      library(pogit, quietly = TRUE)
      library(lubridate, quietly = TRUE)
      library(RcppRoll, quietly = TRUE)
    })
    clusterExport(
      cluster(analysis),
      c(
        'sensitivity_ds',
        'weightSensitivityAnalysis',
        'sensitivity_groups'
      ),
      environment()
    )
    sensitivity_analysis_full <-
      setNames(
        parLapply(
          cl = cluster(analysis),
          sensitivity_groups,
          fun = weightSensitivityAnalysis,
          covars = sensitivity_covars_full,
          ds = sensitivity_ds,
          impact = sensitivity_impact_full,
          time_points = analysis$time_points,
          intervention_date = analysis$intervention_date,
          n_seasons = analysis$n_seasons,
          outcome = analysis$outcome,
          eval_period = analysis$eval_period,
          post_period = analysis$post_period,
          year_def = analysis$year_def
        ),
        sensitivity_groups
      )
    stopCluster(analysis)

    results$sensitivity_pred_quantiles <-
      lapply(
        sensitivity_analysis_full,
        FUN = function(sensitivity_analysis) {
          pred_list <-
            vector(mode = 'list',
                   length = length(sensitivity_analysis))
          for (sensitivity_index in 1:length(sensitivity_analysis)) {
            pred_list[[sensitivity_index]] <-
              getPred(sensitivity_analysis[[sensitivity_index]])
          }
          return(pred_list)
        }
      )
    
    #Table of rate ratios for each sensitivity analysis level
    results$sensitivity_table <-
      t(
        sapply(
          sensitivity_groups,
          sensitivityTable,
          sensitivity_analysis = sensitivity_analysis_full,
          original_rr = analysis$results$impact$full$rr_mean
        )
      )
    results$sensitivity_table_intervals <- data.frame(
      'Estimate (95% CI)' = makeInterval(
        results$sensitivity_table[, 2],
        results$sensitivity_table[, 3],
        results$sensitivity_table[, 1]
      ),
      'Top Control 1' = results$sensitivity_table[, 'Top Control 1'],
      'Inclusion Probability of Control 1' = results$sensitivity_table[, 'Inclusion Probability of Control 1'],
      'Control 1 Estimate (95% CI)' = makeInterval(
        results$sensitivity_table[, 7],
        results$sensitivity_table[, 8],
        results$sensitivity_table[, 6]
      ),
      'Top Control 2' = results$sensitivity_table[, 'Top Control 2'],
      'Inclusion Probability of Control 2' = results$sensitivity_table[, 'Inclusion Probability of Control 2'],
      'Control 2 Estimate (95% CI)' = makeInterval(
        results$sensitivity_table[, 12],
        results$sensitivity_table[, 13],
        results$sensitivity_table[, 11]
      ),
      'Top Control 3' = results$sensitivity_table[, 'Top Control 3'],
      'Inclusion Probability of Control 3' = results$sensitivity_table[, 'Inclusion Probability of Control 3'],
      'Control 3 Estimate (95% CI)' = makeInterval(
        results$sensitivity_table[, 17],
        results$sensitivity_table[, 18],
        results$sensitivity_table[, 16]
      ),
      check.names = FALSE
    )
    results$rr_table <-
      cbind.data.frame(round(analysis$results$impact$time$rr_mean[!bad_sensitivity_groups,], 2),
                       results$sensitivity_table)
    results$rr_table_intervals <-
      cbind(
        'Trend Estimate (95% CI)' = analysis$results$impact$time$rr_mean_intervals[!bad_sensitivity_groups,],
        results$sensitivity_table_intervals
      )
  } else {
    results$sensitivity_table_intervals <- NA
  }
  
  analysis$results$sensitivity <- results
  return(results)
}


#Formats the data
#' @importFrom plyr rbind.fill arrange
evaluatr.impact.pre = function(analysis, run.stl=TRUE) {
  clusterUpdateAnalysis(analysis, function(analysis) {
    # Setup data
    prelog_data <-
      analysis$input_data[!is.na(analysis$input_data[, analysis$outcome_name]), ]#If outcome is missing, delete
    analysis$groups <-
      as.character(unique(unlist(prelog_data[, analysis$group_name], use.names = FALSE)))
    analysis$groups <-
      analysis$groups[!(analysis$groups %in% analysis$.private$exclude_group)]
    
    # Identify variables and groups that have zero variance
    # Zero-variance groups within a variable are replaced with NA
    # Entire zero-variance variables are removed from analysis
    prelog_data.split <- split(prelog_data, prelog_data[[analysis$group_name]])
    variance.vars <- lapply(prelog_data.split, function(group) apply(group[1:nrow(group),], 2, function(row) { 
      suppressWarnings(var(row))
    }))
    for(i in 1:length(prelog_data.split)){
      var.floor <- 1e-6
      # Flag low-variance columns; don't touch data in group/date/denom/outcome columns
      bad.cols <- !is.na(variance.vars[[i]]) & variance.vars[[i]] < 1e-6 & !(names(prelog_data) %in% c(analysis$group_name, analysis$date_name, analysis$outcome_name, analysis$denom_name))
      for (name in names(prelog_data.split[[i]])[bad.cols]) {
        dataCheckWarning(paste0("Data for '", name, "' removed from group '", names(prelog_data.split)[i], "' due to zero variance.\n"))
      }
      prelog_data.split[[i]] <- prelog_data.split[[i]][,!bad.cols]
    }
    names.before = names(prelog_data)
    prelog_data <- rbind.fill(prelog_data.split)
    prelog_data <- arrange(prelog_data, prelog_data[[analysis$group_name]], prelog_data[[analysis$date_name]])
    names.after = names(prelog_data)
    for (name in setdiff(names.before, names.after)) {
      dataCheckWarning(paste0("Data for '", name, "' removed from all groups due to zero variance.\n"))
    }
  
    # Format covars
    #test<-split(prelog_data, factor(prelog_data[,analysis$group_name]))
    #outcome.na<-sapply(test, function(x) sum(is.na(x[,analysis$outcome_name])))
    prelog_data[, analysis$date_name] <-
      formatDate(prelog_data[, analysis$date_name])
    prelog_data <- setNames(
      lapply(
        analysis$groups,
        FUN = splitGroup,
        ungrouped_data = prelog_data,
        group_name = analysis$group_name,
        date_name = analysis$date_name,
        start_date = analysis$start_date,
        end_date = analysis$end_date,
        no_filter = c(
          analysis$group_name,
          analysis$date_name,
          analysis$outcome_name,
          analysis$denom_name
        ),
        sparse_threshold = analysis$sparse_threshold
      ),
      analysis$groups
    )
    #if (exists('exclude_group')) {prelog_data <- prelog_data[!(names(prelog_data) %in% exclude_group)]}
  
    #Log-transform all variables, adding 0.5 to counts of 0.
 if(analysis$log.covars){
    analysis$.private$ds <-
      setNames(lapply(
        prelog_data,
        FUN = logTransform,
        no_log = c(
          analysis$group_name,
          analysis$date_name,
          analysis$outcome_name
        )
      ),
      analysis$groups)
 }else{
   analysis$.private$ds <-
     setNames( prelog_data,
     analysis$groups)
   
 }
   
    analysis$time_points <-
      unique(analysis$.private$ds[[1]][, analysis$date_name])
  
    #Monthly dummies
    if (analysis$n_seasons == 4) {
      dt <- quarter(as.Date(analysis$time_points))
    }
    if (analysis$n_seasons == 12) {
      dt <- month(as.Date(analysis$time_points))
    }
    if (analysis$n_seasons == 3) {
      dt.m <- month(as.Date(analysis$time_points))
      dt <- dt.m
      dt[dt.m %in% c(1, 2, 3, 4)] <- 1
      dt[dt.m %in% c(5, 6, 7, 8)] <- 2
      dt[dt.m %in% c(9, 10, 11, 12)] <- 3
    }
    analysis$.private$season.dummies <- dummies::dummy(dt)
    analysis$.private$season.dummies <- as.data.frame(analysis$.private$season.dummies)
    names(analysis$.private$season.dummies) <- paste0('s', 1:analysis$n_seasons)
    analysis$.private$season.dummies <- analysis$.private$season.dummies[, -analysis$n_seasons]
  
    analysis$.private$ds <-
      lapply(analysis$.private$ds, function(ds) {
        if (!(analysis$denom_name %in% colnames(ds))) {
          ds[analysis$denom_name] <- 0
        }
        return(ds)
      })
  
    analysis$sparse_groups <-
      sapply(analysis$.private$ds, function(ds) {
        return(ncol(ds[!(
          colnames(ds) %in% c(
            analysis$date_name,
            analysis$group_name,
            analysis$denom_name,
            analysis$outcome_name,
            analysis$.private$exclude_covar
          )
        )]) == 0)
      })
    analysis$.private$ds <-
      analysis$.private$ds[!analysis$sparse_groups]
    if (length(analysis$.private$ds) == 0) {
      stop("Unable to proceed with analysis: all groups are sparse")
    }
    analysis$groups <- analysis$groups[!analysis$sparse_groups]
  
    #Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
    analysis$covars = list()
    analysis$covars$full <-
      setNames(lapply(analysis$.private$ds, function(group) {
        makeCovars(
          analysis$country,
          analysis$time_points,
          analysis$intervention_date,
          analysis$.private$season.dummies,
          group
        )
      }), analysis$groups)
    analysis$covars$full <-
      lapply(
        analysis$covars$full,
        FUN = function(covars) {
          covars[,!(colnames(covars) %in% analysis$.private$exclude_covar), drop = FALSE]
        }
      )
    analysis$covars$full <-
      lapply(
        analysis$covars$full, function(x){ 
          x[,analysis$outcome_name]<-NULL
          return(x)
        }
      )
    analysis$covars$time <-
      setNames(lapply(
        analysis$covars$full,
        FUN = function(covars) {
          as.data.frame(list(cbind(
            analysis$.private$season.dummies, time_index = 1:nrow(covars)
          )))
        }
      ),
      analysis$groups)
    analysis$covars$null <-
      setNames(lapply(
        analysis$covars$full,
        FUN = function(covars) {
          as.data.frame(list(cbind(analysis$.private$season.dummies)))
        }
      ),
      analysis$groups)
  
    #Standardize the outcome variable and save the original mean and SD for later analysis.
    analysis$outcome <-
      sapply(
        analysis$.private$ds,
        FUN = function(data) {
          data[, analysis$outcome_name]
        }
      )
    
    analysis
  })
    
    
  offset <-
    sapply(
      analysis$.private$ds,
      FUN = function(data)
        exp(data[, analysis$denom_name])
    ) #offset term on original scale; 1 column per age group
  
  ##SECTION 1: CREATING SMOOTHED VERSIONS OF CONTROL TIME SERIES AND APPENDING THEM ONTO ORIGINAL DATAFRAME OF CONTROLS
  #EXTRACT LONG TERM TREND WITH DIFFERENT LEVELS OF SMOOTHNESS USING STL
  # Set a list of parameters for STL
  if(run.stl==TRUE){
        stl.covars <-
          mapply(
            smooth_func,
            ds.list = analysis$.private$ds,
            covar.list = analysis$covars$full,
            SIMPLIFY = FALSE,
            MoreArgs = list(n_seasons = analysis$n_seasons)
          )
        post.start.index <-
          which(analysis$time_points == analysis$post_period[1])
        
        if (length(analysis$groups) > 1) {
          stl.data.setup <-
            mapply(
              stl_data_fun,
              covars = stl.covars,
              ds.sub = analysis$.private$ds ,
              SIMPLIFY = FALSE,
              MoreArgs = list(
                n_seasons = analysis$n_seasons,
                outcome_name = analysis$outcome_name,
                post.start.index = post.start.index
              )
            ) #list of lists that has covariates for each regression for each strata
        } else{
          stl.data.setup <-
            list(
              mapply(
                stl_data_fun,
                covars = stl.covars,
                ds.sub = analysis$.private$ds,
                MoreArgs = list(
                  n_seasons = analysis$n_seasons,
                  outcome_name = analysis$outcome_name,
                  post.start.index = post.start.index
                )
              )
            )
        }
        
        ##SECTION 2: run first stage models for STL
        addProgress(analysis, sprintf("STL first stage (group %s)", analysis$groups), after=0)
        glm.results <-
          vector("list", length = length(stl.data.setup)) #combine models into a list
        clusterEvalQ(cluster(analysis), {
          library(lme4, quietly = TRUE)
        })
        clusterExport(cluster(analysis),
                      c('stl.data.setup', 'glm.fun', 'post.start.index'),
                      environment())
        for (i in 1:length(stl.data.setup)) {
          progressStartPart(analysis)
          glm.results[[i]] <-
            parLapply(
              cl = cluster(analysis),
              stl.data.setup[[i]],
              fun = function(d) {
                glm.fun(d, post.start.index)
              }
            )
          progressEndPart(analysis)
        }
        stopCluster(analysis)
  }
  ######################
  
  clusterUpdateAnalysis(analysis, function(analysis) {
    # Combine data
    #Combine the outcome, covariates, and time point information.
    analysis$.private$data$full <-
      setNames(
        lapply(
          analysis$groups,
          makeTimeSeries,
          outcome = analysis$outcome,
          covars = analysis$covars$full
        ),
        analysis$groups
      )
    analysis$.private$data$time <-
      setNames(
        lapply(
          analysis$groups,
          makeTimeSeries,
          outcome = analysis$outcome,
          covars = analysis$covars$time,
          trend = TRUE,
          offset = offset
        ),
        analysis$groups
      )
    if(run.stl==TRUE){
    analysis$.private$data$pca <-
      mapply(
        FUN = pca_top_var,
        glm.results.in = glm.results,
        covars = stl.covars,
        ds.in = analysis$.private$ds,
        SIMPLIFY = FALSE,
        MoreArgs = list(
          outcome_name = analysis$outcome_name,
          season.dummies = analysis$.private$season.dummies
        )
      )
    names(analysis$.private$data$pca) <- analysis$groups
    }
    #Time trend model but without a denominator
    analysis$.private$data$time_no_offset <-
      setNames(
        lapply(
          analysis$groups,
          makeTimeSeries,
          outcome = analysis$outcome,
          covars = analysis$covars$time,
          trend = FALSE
        ),
        analysis$groups
      )
    analysis
  })
}

addProgress <- function(analysis, partNames, after=NULL) {
  if (is.null(after)) {
    after = length(analysis$.private$progress_parts)
  }
  analysis$.private$progress_parts = append(analysis$.private$progress_parts, partNames, after=after)
}

progressStartPart <- function(analysis) {
  analysis$.private$progress(analysis, analysis$.private$progress_done, analysis$.private$progress_parts)
}

progressEndPart <- function(analysis) {
  analysis$.private$progress_done = analysis$.private$progress_done + 1
  analysis$.private$progress(analysis, analysis$.private$progress_done, analysis$.private$progress_parts, before=FALSE)
}

showProgress = function(analysis, done, names, before=TRUE) {
  if (before) {
    write(sprintf("Starting analysis part %s of %s (%s)", done + 1, length(names), names[[done + 1]]), stdout())
  } else {
    write(sprintf("Finished analysis part %s of %s (%s)", done, length(names), names[[done]]), stdout())
  }
}


#' Perform analysis controling for 1 variable at a time
#' @param analysis Analysis object, initialized by TODO.init.
#' @return Univariate analysis results, `results`, as described below
#'
#' `results$rr` Rate ratio (RR) estimate (median) from single control model
#' 
#' `results$rr.ucl` RR Upper 95% CrI from single control model
#' 
#' `results$rr.lcl` RR Lower 95% CrI from single control model
#'  
#' `results$aic.wgt` Akaike weight of the univariate model
#' 
#' `results$covar` Name of the control variable used for adjustment
#'
#' @export
evaluatr.univariate <- function(analysis) {
  evaluatr.impact.pre(analysis,run.stl=FALSE) #formats the data
  #####
  # analysis$.private$data$full
  clusterUpdateAnalysis(analysis, function(analysis) {
    results<-lapply( analysis$.private$data$full, single.var.glmer, 
                     n_seasons=analysis$n_seasons,
                    intro.date=analysis$post_period[1],
                    time_points=analysis[['time_points']],
                    eval.period=analysis$eval_period)  
    univariate.aic <- lapply(results, '[[','aic.summary')
    covar.names <- lapply(results, '[[','covar.names')
    aic.weights<- lapply(univariate.aic, function(x) exp(-0.5*(x-min(x)))/sum(exp(-0.5*(x-min(x)))) )
    rr.post <- lapply(results, '[[','rr.post')
    summary.results<- vector("list", length(rr.post)) 
    for(i in 1:length(rr.post)){
      summary.results[[i]]<-cbind.data.frame(rr.post[[i]],round(aic.weights[[i]],3),covar.names[[i]])
     names(summary.results[[i]])<-c('rr.lcl','rr','rr.ucl','aic.wgt','covar')
     summary.results[[i]]<-summary.results[[i]][order(-summary.results[[i]]$aic.wgt),]
    }
    analysis$results$univariate<-summary.results
    return(analysis)
  })
  
  return(analysis$results$univariate)
}
  

dataCheckWarning = function(message) {
  warning(warningCondition(message, class="evaluatr.dataCheck"))
}

# Set up cluster operations
# If we are using an outside cluster (as we do in the Web UI) start and stop are functions used to start it or stop it
# If we are using our own cluster (as we do by default), then start and stop are NULL and ww will set them up here
cluster = function(analysis) {
  if (is.null(analysis$.private$cluster)) {
    analysis$.private$cluster = analysis$.private$startCluster()
  }
  analysis$.private$cluster
}

stopCluster = function(analysis) {
  analysis$.private$stopCluster(analysis$.private$cluster)
  analysis$.private$cluster = NULL
}

defaultStartCluster = function() {
  if (Sys.getenv("CI") != "") {
    # If running on GitLab, default to multi-session cluster using all cores
    n_cores <- availableCores(methods=c("system"))
  } else {
    # If running on someone's personal computer, default to multi-session cluster leaving one core free (if possible)
    n_cores <- max(availableCores(methods=c("system")) - 1, 1)
  }
  parallel::makeCluster(n_cores)
}

defaultStopCluster = function(cluster) {
  parallel::stopCluster(cluster)
}

# This evaluates ... on a single node in a cluster. This doesn't help with performance (in fact, it will slightly decrease it), but it allows WebUI to push computation off the CPU that's doing web stuff
clusterEval1 = function(analysis, func) {
  result = future::value(future::remote(func(analysis), workers=cluster(analysis)))
  stopCluster(analysis)
  result
}

# This evaluates ... on a single node in a cluster, then updates the analysis object with the result of ... evaluation. It allows WebUI to push updates to the analysis object off to another CPU
clusterUpdateAnalysis = function(analysis, func) {
  newanalysis = clusterEval1(analysis, func)
  # Sending .private$cluster through remote evaluation makes it sad
  newanalysis$.private$cluster = analysis$.private$cluster
  analysis[names(newanalysis)] = newanalysis
}

# This gets us compatibility with versions of R back to 3.0
.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}
