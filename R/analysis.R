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
  
  return(analysis)
}

# This is used by the web UI to set up parallel computation 
evaluatr.initParallel = function(analysis, startCluster, stopCluster, progress) {
  analysis$.private$startCluster = startCluster
  analysis$.private$stopCluster = stopCluster
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
#' `results$rr_mean_combo` Rate ratios during the evaluation period, as calculated in each of the models.
#'
#' @importFrom stats AIC as.formula cov dpois glm median poisson prcomp predict quantile rmultinom rnorm rpois sd setNames stl var vcov complete.cases
#' @importFrom lme4 glmer glmerControl fixef
#' @importFrom lubridate as_date %m+% year month quarter time_length interval yday
#' @importFrom splines bs
#' @importFrom reshape melt
#' @importFrom MASS mvrnorm
#' @importFrom HDInterval hdi
#' @importFrom RcppRoll roll_sum
#' @importFrom plyr rbind.fill
#' @importFrom pogit poissonBvs
#' @importFrom parallel makeCluster clusterEvalQ clusterExport stopCluster parLapply
#' @importFrom stats pnorm
#' @export

evaluatr.impact = function(analysis, variants=c('full','time')) {
  addProgress(analysis, sprintf("Impact analysis (%s)", lapply(analysis$.private$variants, function(variant) variant$name)))
  evaluatr.impact.pre(analysis, run.stl= ('pca' %in% variants))
  results1 = list()
  

  analysis$.private$variants = analysis$.private$variants[variants]
  for (variant in variants) {
    progressStartPart(analysis)
    results1[[variant]]$groups <- setNames(
      lapply(
        analysis$.private$data[['full']],
        FUN = inla_mods,
        model.variant=variant
      ),
      analysis$groups
    )
    progressEndPart(analysis)
  }
  

  
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

    for (variant in intersect(c('full'), variants)) {

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
          analysis$.private$variants[c('full','time')], 
          function(variant) variant$name
        )
      )
    #print(levels(rr_mean_combo$Model))
    
    
    analysis$results$impact <- results
    analysis = evaluatr.prune(analysis, what="impact")
    return( analysis$results$impact)
}  
  



#Formats the data
#' @importFrom plyr rbind.fill arrange
evaluatr.impact.pre = function(analysis, run.stl=TRUE) {
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
    

  
  
  offset <-
    sapply(
      analysis$.private$ds,
      FUN = function(data)
        exp(data[, analysis$denom_name])
    ) #offset term on original scale; 1 column per age group
  
  ######################
  
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
   return(analysis$results$univariate)
}


dataCheckWarning = function(message) {
  warning(warningCondition(message, class="evaluatr.dataCheck"))
}

# This gets us compatibility with versions of R back to 3.0
.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}
