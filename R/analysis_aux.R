#Rearrange date to YYYY-MM-DD format.
formatDate <- function(time_points) {
  time_points <- as_date(time_points)
  time_points <- as.Date(time_points, format = '%Y-%m-%d')
  return(time_points)
}

splitGroup <-
  function(ungrouped_data,
           group_name,
           group,
           date_name,
           start_date,
           end_date,
           no_filter = NULL,
           sparse_threshold = 5) {
    ds <- ungrouped_data[ungrouped_data[, group_name] == group,]
    ds <- ds[, colSums(is.na(ds)) == 0]
    ds <-
      ds[match(start_date, ds[, date_name]):match(end_date, ds[, date_name]),]
    ds <-
      cbind(ds[, colnames(ds) %in% no_filter], filterSparse(ds[,!(colnames(ds) %in% no_filter), drop =
                                                                 FALSE], threshold=sparse_threshold))
    return(ds)
  }

#Log-transform the covariate
logTransform <- function(prelog_data, no_log = NULL) {
  prelog_data[,!(colnames(prelog_data) %in% no_log)][prelog_data[,!(colnames(prelog_data) %in% no_log)] == 0] <-
    0.5
  prelog_data[,!(colnames(prelog_data) %in% no_log)] <-
    log(prelog_data[,!(colnames(prelog_data) %in% no_log)])
  return(prelog_data)
}

filterSparse <- function(dataset, threshold = 5) {
  return(dataset[, colMeans(dataset) > threshold, drop = FALSE])
}

#Used to adjust the Brazil data for a code shift in 2008.
getTrend <- function(covar_vector, data) {
  new_data <- data
  new_data[c('bs1', 'bs2', 'bs3', 'bs4')] <- 0
  new_data$month_i <- as.factor(1)
  trend <-
    predict(
      glm(covar_vector ~ month_i + ., family = 'gaussian', data = data),
      type = 'response',
      newdata = new_data
    ) #month_i is first to be the reference.
  names(trend) <- NULL
  return(trend)
}

makeCovars <-
  function(country,
           time_points,
           intervention_date,
           season.dummies,
           ds_group) {
    if (country == "BrazilADJ2008") {
      #Eliminates effects from 2008 coding change
      covars <- ds_group[, 4:ncol(ds_group)]
      month_i <- as.factor(as.numeric(format(time_points, '%m')))
      spline <-
        setNames(as.data.frame(bs(
          1:nrow(covars), knots = 5, degree = 3
        )), c('bs1', 'bs2', 'bs3', 'bs4'))
      year_2008 <- numeric(nrow(covars))
      year_2008[1:nrow(covars) >= match(as.Date('2008-01-01'), time_points)] <-
        1
      data <- cbind.data.frame(year_2008, spline, month_i)
      trend <- lapply(covars, getTrend, data = data)
      covars <- covars - trend
    } else {
      covars <- ds_group[, 4:ncol(ds_group), drop = FALSE]
    }
    # if (intervention_date > as.Date('2009-09-01')) {
    #   covars$pandemic <-
    #     ifelse(time_points == '2009-08-01',
    #            1,
    #            ifelse(time_points == '2009-09-01', 1, 0))
    # }
    covars <-
      as.data.frame(lapply(covars[, apply(covars, 2, var) != 0, drop = FALSE], scale), check.names = FALSE)
    covars <- cbind(season.dummies, covars)
    return(covars)
  }

#Combine the outcome and covariates.
makeTimeSeries <-
  function(group,
           outcome,
           covars,
           trend = FALSE,
           offset = NA) {
    if (trend == FALSE) {
      return(cbind(outcome = outcome[, group],  covars[[group]]))
    }
    if (trend == TRUE) {
      return(cbind(
        outcome = outcome[, group],
        log.offset = log(offset[, group] + 0.5),
        covars[[group]]
      ))
    }
  }


#Main analysis function.
doCausalImpact <-
  function(zoo_data,
           intervention_date,
           n_seasons,
           ri.select = TRUE,
           time_points,
           crossval.stage = FALSE,
           var.select.on = TRUE,
           n_iter = 10000,
           burnN,
           sampleN,
           trend = FALSE) {
    #Format outcome and covariates for regular and cross-validations
    if (crossval.stage) {
      #Data for cross-validation
      y.pre <- zoo_data$cv.data[, 1]
      y.full <- zoo_data$full.data[, 1]
      exclude.indices <- zoo_data$exclude.indices
      if (trend) {
        x <-
          as.matrix(zoo_data$full.data[,-c(1, 2)]) #Removes outcome column and offset from dataset
        offset.t <- as.vector(exp(as.matrix(zoo_data$full.data[, 2])))
        offset.t.pre <- as.vector(exp(as.matrix(zoo_data$cv.data[, 2])))
        x.pre <- as.matrix(zoo_data$cv.data[, -c(1, 2)])
      } else{
        x <-
          as.matrix(zoo_data$full.data[,-1]) #Removes outcome column from dataset
        x.pre <- as.matrix(zoo_data$cv.data[, -1])
        x.pre.var <-
          apply(x.pre, 2, var) #test if covariate changes at all in pre-period; if not delete (ie pandemic)
        x <- x[, x.pre.var > 0] #eliminate from full matrix
        x.pre <- x.pre[, x.pre.var > 0] #eliminate from cv matrix
      }
    } else{
      ##Data for non-cross-validation
      y.pre <- zoo_data[time_points < as.Date(intervention_date), 1]
      y.full <- zoo_data[, 1] #all y
      exclude.indices <- NA
      if (trend) {
        x <-
          as.matrix(zoo_data[,-c(1, 2)]) #Removes outcome column and offset from dataset
        offset.t <- as.vector(exp(as.matrix(zoo_data[, 2])))
        offset.t.pre <-
          offset.t[time_points < as.Date(intervention_date)]
      } else{
        x <-
          as.matrix(zoo_data[,-c(1)]) #Removes outcome column from dataset
      }
      x.pre <-
        as.matrix(x[time_points < as.Date(intervention_date),])
    }
    
    post_period_response <- y.full
    post_period_response <-
      as.vector(post_period_response[time_points >= as.Date(intervention_date)])
    
    covars <- x
    cID <- seq_along(y.pre) #used for observation-level random effect
    
    #Which variables are fixed in the analysis (not estimated)
    if (trend) {
      deltafix.mod <-
        c(rep(1, times = (ncol(x.pre) - 1)), 0) #monthly dummies, offset, fixed
      bsts_model.pois  <-
        poissonBvs(
          y = y.pre ,
          X = x.pre,
          offset = offset.t.pre,
          BVS = var.select.on,
          model = list(
            deltafix = deltafix.mod,
            ri = ri.select,
            clusterID = cID
          ),
          mcmc=list(
            burnin=burnN,
            M=sampleN
          )
        )
    } else{
      if (var.select.on) {
        deltafix.mod <- rep(0, times = (ncol(x.pre)))
        deltafix.mod[1:(n_seasons - 1)] <- 1 #fix  monthly dummies
        bsts_model.pois  <-
          poissonBvs(
            y = y.pre ,
            X = x.pre,
            BVS = TRUE,
            model = list(
              deltafix = deltafix.mod,
              ri = TRUE,
              clusterID = cID
            ),
            mcmc=list(
              burnin=burnN,
              M=sampleN
              
            )
          )
      } else{
        if (ri.select) {
          bsts_model.pois  <-
            poissonBvs(
              y = y.pre ,
              X = x.pre,
              BVS = FALSE,
              model = list(ri = TRUE, clusterID = cID),
              mcmc=list(
                burnin=burnN,
                M=sampleN
                
              )
            )
        } else{
          bsts_model.pois  <- poissonBvs(y = y.pre , X = x.pre, BVS = FALSE,
                                         mcmc=list(
                                           burnin=burnN,
                                           M=sampleN
                                           
                                         ))
        }
      }
    }
    
    beta.mat <- bsts_model.pois$samplesP$beta[-c(1:burnN), ]
    x.fit <- cbind(rep(1, nrow(x)), x)
    #Generate  predictions with prediction interval
    if (ri.select) {
      disp <-
        bsts_model.pois$samplesP$thetaBeta[-c(1:burnN), 1] ##note theta beta is signed--bimodal dist---take abs value
      disp.mat <-
        rnorm(
          n = length(disp) * length(y.full),
          mean = 0,
          sd = abs(disp)
        ) #Note: confirmed that median(abs(disp)) is same as sd(rand.eff)
      disp.mat <-
        t(matrix(disp.mat, nrow = length(disp), ncol = length(y.full)))
    } else{
      disp.mat = 0 #if no random effect in model, just set to 0.
    }
    if (trend) {
      reg.mean <-   exp((x.fit %*% t(beta.mat)) + disp.mat)  * offset.t
    } else{
      reg.mean <-   exp((x.fit %*% t(beta.mat)) + disp.mat)
    }
    predict.bsts <- rpois(length(reg.mean), lambda = reg.mean)
    predict.bsts <-
      matrix(predict.bsts,
             nrow = nrow(reg.mean),
             ncol = ncol(reg.mean))
    #predict.bsts.q<-t(apply(predict.bsts,1,quantile, probs=c(0.025,0.5,0.975)))
    #matplot(predict.bsts.q, type='l', col=c('gray','black','gray'), lty=c(2,1,2), bty='l', ylab="N hospitalizations")
    #points(y.full)
    
    #Inclusion probabilities Poisson model
    incl.probs.mat <- t(bsts_model.pois$samplesP$pdeltaBeta[-c(1:burnN), ])
    inclusion_probs <- apply(incl.probs.mat, 1, mean)
    summary.pois <- summary(bsts_model.pois)
    covar.names <- dimnames(x.pre)[[2]]
    if (ri.select) {
      rand.eff <- bsts_model.pois$samplesP$bi[-c(1:burnN), ]
    } else{
      rand.eff = 0
    }
    inclusion_probs <- cbind.data.frame(covar.names, inclusion_probs)
    if (trend) {
      impact <-
        list(
          reg.mean,
          exclude.indices,
          rand.eff,
          offset.t,
          covars,
          beta.mat,
          predict.bsts,
          inclusion_probs,
          post.period.response = post_period_response,
          observed.y = y.full
        )
      names(impact) <-
        c(
          'reg.mean',
          'exclude.indices',
          'rand.eff',
          'offset.t.pre',
          'covars',
          'beta.mat',
          'predict.bsts',
          'inclusion_probs',
          'post_period_response',
          'observed.y'
        )
    } else{
      impact <-
        list(
          reg.mean,
          exclude.indices,
          rand.eff,
          covars,
          beta.mat,
          predict.bsts,
          inclusion_probs,
          post.period.response = post_period_response,
          observed.y = y.full
        )
      names(impact) <-
        c(
          'reg.mean',
          'exclude.indices' ,
          'rand.eff',
          'covars',
          'beta.mat',
          'predict.bsts',
          'inclusion_probs',
          'post_period_response',
          'observed.y'
        )
    }
    return(impact)
  }

crossval.log.lik <- function(cv.impact) {
  exclude.id <- cv.impact$exclude.indices
  pred.exclude <-
    cv.impact$reg.mean[exclude.id, , drop = FALSE] #use predicted mean (which incorporates random effect)
  obs.exclude <- cv.impact$observed.y[exclude.id]
  point.ll1 <-
    matrix(NA,
           nrow = nrow(pred.exclude),
           ncol = ncol(pred.exclude))
  for (i in 1:ncol(pred.exclude)) {
    point.ll1[, i] <-
      dpois(obs.exclude, lambda = pred.exclude[, i], log = TRUE)
  }
  
  return(point.ll1)
}

# cross validated predictions vs observed
pred.cv <- function(cv.impact) {
  exclude.id <- cv.impact$exclude.indices
  pred.exclude <-
    cv.impact$reg.mean[exclude.id, , drop = FALSE] #use predicted mean (which incorporates random effect)
  cv.pred.q <-
    t(apply(pred.exclude, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  cv.pred.q <-
    cbind(cv.impact$observed.y[exclude.id], cv.pred.q) #combine observed and expected
  return(cv.pred.q)
}

stack.mean <-
  function(group,
           impact_full,
           impact_time,
           impact_time_no_offset,
           impact_pca,
           stacking_weights.all,
           outcome) {
    #Averaged--multiply each log(mean) by weight, then add, then exponentiate and draw from Poisson
    weights <-
      as.numeric(as.vector(stacking_weights.all[stacking_weights.all$groups ==
                                                  group, ]))
    
    #df[sample(nrow(df), 3), ]
    samp.probs <- weights[2:5]
    n.iter <- ncol(impact_full$reg.mean)
    n.samps <- rmultinom(n = 1, size = n.iter, prob = samp.probs)
    
    rm.full <- log(impact_full$reg.mean[, sample(n.iter, n.samps[1])])
    rm.time <- log(impact_time$reg.mean[, sample(n.iter, n.samps[2])])
    rm.time_no_offset <-
      log(impact_time_no_offset$reg.mean[, sample(n.iter, n.samps[3])])
    rm.pca <- log(impact_pca$reg.mean[, sample(n.iter, n.samps[4])])
    
    pred.full <- apply(impact_full$reg.mean, 1, median)
    pred.time <- apply(impact_time$reg.mean, 1, median)
    pred.time.no_offset <-
      apply(impact_time_no_offset$reg.mean, 1, median)
    pred.pca <- apply(impact_pca$reg.mean, 1, median)
    all.preds <-
      cbind(pred.full, pred.time, pred.time.no_offset, pred.pca)
    
    stack <- cbind(rm.full, rm.time, rm.time_no_offset, rm.pca)
    pred.stack.count <- rpois(n = length(stack), lambda = exp(stack))
    pred.stack.count <-
      matrix(pred.stack.count,
             nrow = nrow(stack),
             ncol = ncol(stack))
    pred.stack.q <-
      t(apply(pred.stack.count, 1, quantile, probs = c(0.025, 0.5, 0.975)))
    # log.rr.stack.q<-log((outcome[,group]+0.5)/pred.stack.q)
    # log.rr.iter<- log((outcome[,group]+0.5)/pred.stack.count)
    # log_rr_stack.cov<-cov(t(log.rr.iter))
    # log_rr_stack.prec<-solve(log_rr_stack.cov) #NOT INVERTIBLE?
    # #log_rr_stack.prec=log_rr_stack.cov
    stacked.est <-
      list(pred.stack.count, pred.stack.q, outcome[, group])
    names(stacked.est) <-
      list('predict.bsts', 'pred.stack.q', 'observed.y')
    return(stacked.est)
  }

plot.stack.est <- function(stacked.ests) {
  matplot(
    stacked.ests$all.preds,
    lty = 1,
    type = 'l',
    col = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c'),
    bty = 'l',
    lwd = 0.5,
    ylim = c(0, range(stacked.ests$pred.stack.q)[2])
  )
  matplot(
    stacked.ests$pred.stack.q,
    type = 'l',
    col = c('gray', 'black', 'gray'),
    lwd = c(1, 2, 1),
    lty = c(2, 1, 2),
    bty = 'l',
    add = TRUE
  )
  points(stacked.ests$y)
  legend(
    'bottomleft',
    inset = 0.02,
    legend = c(
      "Synthetic controls",
      "Time trend",
      'Time trend, no offset',
      'STL+PCA'
    ),
    col = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c'),
    lty = 1,
    cex = 0.8,
    box.lty = 0
  )
}

reshape.arr <- function(point.ll3) {
  mat1 <- do.call(rbind, point.ll3)
  # arr1=sapply(point.ll3, function(x) x, simplify='array')
  # arr2<-aperm(arr1,perm=c(1,3,2))
  # mat1=matrix(arr2,nrow=dim(arr2)[1]*dim(arr2)[2],ncol=dim(arr2)[3] ) #f want to consider them ll together
  mean.mat <-
    apply(mat1, 1, logmeanexp)  #want to take mean of the likelihood, not log-likelihood..so use logmean exp to avoid underflow
  return(mean.mat)
}

#Save inclusion probabilities.
inclusionProb <- function(impact) {
  return(impact$inclusion_probs)
}

##K-fold cross-validation
#Create K-fold datasets
makeCV <- function(ds, time_points, intervention_date) {
  impact <- ds
  impact.pre <- ds[time_points < intervention_date, ]
  year.pre <- year(time_points[time_points < intervention_date])
  year.pre.vec <- unique(year.pre)
  N.year.pre <- length(year.pre.vec)
  impact.list.cv <- vector("list",  length = N.year.pre)
  cv.data <- vector("list",  length = N.year.pre)
  for (i in 1:N.year.pre) {
    impact.list.cv <- impact.pre[!(year.pre == year.pre.vec[i]), ]
    exclude.indices <- which((year.pre == year.pre.vec[i]))
    cv.data[[i]] <-
      list(impact.list.cv, impact, exclude.indices) #combine full data (pre and post) and CV data
    names(cv.data[[i]]) <- c('cv.data', 'full.data', 'exclude.indices')
  }
  return(cv.data)
}

#Estimate the rate ratios during the evaluation period and return to the original scale of the data.
rrPredQuantiles <-
  function(impact,
           denom_data = NULL,
           eval_period,
           post_period,
           n_seasons,
           year_def,
           time_points) {
    pred_samples <- impact$predict.bsts
    
    pred <-
      t(apply(
        pred_samples,
        1,
        quantile,
        probs = c(0.025, 0.5, 0.975),
        na.rm = TRUE
      ))
    pred.hdi <- cbind( pred[,'50%'],t(hdi(t(pred_samples), credMass = 0.95)) )
      
    eval_indices <-
      match(which(time_points == eval_period[1]), (1:length(impact$observed.y))):match(which(time_points ==
                                                                                               eval_period[2]), (1:length(impact$observed.y)))
    
    pred_eval_sum <- colSums(pred_samples[eval_indices,])
    
    eval_obs <- sum(impact$observed.y[eval_indices])
    
    # Time values
    # * If year_def is cal_year, then times are (numeric) years.
    #   For example, a data point in year 2000 will have numeric time value of 2000
    # * If year_def is epi_year, then times are (numeric) years corresponding to the first year in the jul-jun epi year
    #   For example, a data point in June 2000 will have numeric time value of 1999.
    #   On the other hand, a data point in July 2000 will have numeric time value of 2000.
    
    if (year_def == 'epi_year') {
      yr <- year(time_points)
      month <- month(time_points)
      yr[month <= 6] <- yr[month <= 6] - 1
      year <- yr
    } else{
      year <- year(time_points)
    }

        #aggregate observed and predicted
    obs.y.year <- tapply(impact$observed.y, year, sum)
    pred.yr.spl <- split(as.data.frame(pred_samples), year)
    pred.yr.spl.sum <-
      lapply(pred.yr.spl, function(x)
        apply(x, 2, sum))
    pred.yr.spl.sum.hdi<-lapply(pred.yr.spl.sum, hdi, credMass = 0.95) 
    pred.yr.spl.sum.q <-
      lapply(pred.yr.spl.sum,
             quantile,
             probs = c(0.025, 0.50, 0.975),
             na.rm = TRUE)
    pred.yr.sum.q <-
      as.data.frame(matrix(
        unlist(pred.yr.spl.sum.q),
        ncol = 3,
        byrow = TRUE
      ))
    names(pred.yr.sum.q) <- c('2.5%', '50%', '97.5%')
    pred.yr.sum.q$obs <- obs.y.year
    pred.yr.sum.q$year <- unique(year)
    
    pred.yr.sum.hdi <-
      cbind.data.frame(pred.yr.sum.q[,'50%'] ,matrix(
        unlist(pred.yr.spl.sum.hdi),
        ncol = 2,
        byrow = TRUE
      ))
    names(pred.yr.sum.hdi) <- c('median','lcl', 'ucl')
    pred.yr.sum.hdi$obs <- obs.y.year
    pred.yr.sum.hdi$year <- unique(year)
    # matplot(unique(year),pred.yr.sum.q[,1:3], type='l', col='lightgray', lty=c(3,1,3), bty='l', xlab="", ylab='Cases (N)', ylim=c(0,max(pred.yr.sum.q) ))
    # points(unique(year),pred.yr.sum.q$obs, pch=16)
    # abline(v=year(intervention_date )+0.5, col='gray', lty=2)
    
    eval_rr_sum <- eval_obs / pred_eval_sum
    rr.iter<-eval_rr_sum
    rr <-
      quantile(eval_rr_sum,
               probs = c(0.025, 0.5, 0.975),
               na.rm = TRUE)
    names(rr) <- c('Lower CI', 'Point Estimate', 'Upper CI')
    rr.hdi <- c(rr['Point Estimate'], hdi(eval_rr_sum, credMass = 0.95)) 
    mean_rr <- mean(eval_rr_sum)
    sd_log_rr <- sd(log(eval_rr_sum))
    
    #Calculate RR for the N months prior to vaccine introduction as a bias corrrection factor
    pre_indices <-
      which(time_points == (eval_period[1] %m+% months(12))):which(time_points ==
                                                                     (eval_period[1] %m+% months(1)))
    pred_pre_sum <- colSums(pred_samples[pre_indices,])
    pre_obs <- sum(impact$observed.y[pre_indices])
    rr_sum_pre <- pre_obs / pred_pre_sum  #Should be 0!
    
    unbias_rr <- eval_rr_sum / rr_sum_pre # same as log_rr - log_rr_pre=log(A/B)
    unbias_rr_q <- quantile(unbias_rr, probs = c(0.025, 0.5, 0.975))
    
    plot_rr_start <- which(time_points == post_period[1]) - n_seasons
    roll_rr_indices <-
      match(plot_rr_start, (1:length(impact$observed.y))):match(which(time_points ==
                                                                        eval_period[2]), (1:length(impact$observed.y)))
    
    obs_full <- impact$observed.y
    
    roll_sum_pred <-
      roll_sum(pred_samples[roll_rr_indices,], n_seasons)
    roll_sum_obs <- roll_sum(obs_full[roll_rr_indices], n_seasons)
    roll_rr_est <- roll_sum_obs / roll_sum_pred
    roll_rr <-
      t(apply(
        roll_rr_est,
        1,
        quantile,
        probs = c(0.025, 0.5, 0.975),
        na.rm = TRUE
      ))
    
    pred_samples_post <- pred_samples[eval_indices,]
    
    # obs_full[obs_full==0]<-0.5 #continuity correction for small sampls
    #  pred_quantiles<-t(apply(pred_samples, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
    #  matplot(impact$observed.y/pred_quantiles, type='l')
    log_rr_full_t_samples <-
      t(log((obs_full + 0.5) / (pred_samples + 0.5))) #Continuity correction for 0s
    log_rr_full_t_quantiles <-
      t(apply(
        log_rr_full_t_samples,
        2,
        quantile,
        probs = c(0.025, 0.5, 0.975),
        na.rm = TRUE
      ))
    log_rr_full_t_hdi <- cbind(log_rr_full_t_quantiles[,'50%'],t(hdi(log_rr_full_t_samples, credMass=0.95) ))
    log_rr_full_t_sd <-
      t(apply(log_rr_full_t_samples, 2, sd, na.rm = TRUE))
    #
    #   #Covariance matrix for pooled analysis
    #log_rr_full_t_samples.covar <- cov(log_rr_full_t_samples)
    post.indices <-
      which(time_points == post_period[1]):which(time_points == post_period[2])
   # log_rr_full_t_samples.prec.post <-
     # solve(log_rr_full_t_samples.covar) #NOT INVERTIBLE?
    #
    # quantiles <- list(pred_samples_post_full = pred_samples_post,roll_rr=roll_rr, log_rr_full_t_samples.prec=log_rr_full_t_samples.prec, log_rr_full_t_samples=log_rr_full_t_samples,log_rr_full_t_quantiles=log_rr_full_t_quantiles,log_rr_full_t_sd=log_rr_full_t_sd, plot_pred = plot_pred,log_plot_pred=log_plot_pred, log_plot_pred_SD=log_plot_pred_SD, rr = rr, mean_rate_ratio = mean_rate_ratio,rr.iter=rr.iter)
    # quantiles <- list(pred_samples = pred_samples, pred = pred, rr = rr, roll_rr = roll_rr, mean_rr = mean_rr)
    
    quantiles <-
      list(
        log_rr_full_t_hdi=log_rr_full_t_hdi,
        rr.hdi=rr.hdi,
        pred.yr.sum.hdi=pred.yr.sum.hdi,
        pred.hdi=pred.hdi,
        unbias_rr_q = unbias_rr_q,
        pred.yr.sum.q = pred.yr.sum.q,
       # log_rr_full_t_samples.prec.post = log_rr_full_t_samples.prec.post,
        pred_samples = pred_samples,
        pred = pred,
        rr = rr,
        roll_rr = roll_rr,
        mean_rr = mean_rr,
        pred_samples_post_full = pred_samples_post,
        roll_rr = roll_rr,
        log_rr_full_t_quantiles = log_rr_full_t_quantiles,
        log_rr_full_t_sd = log_rr_full_t_sd,
        rr = rr,
        rr.iter=rr.iter
      )
    return(quantiles)
  }

getPred <- function(quantiles) {
  return(quantiles$pred)
}

getPredHDI <- function(quantiles) {
  return(quantiles$pred.hdi)
}

getAnnPred <- function(quantiles) {
  return(quantiles$pred.yr.sum.q)
}

getAnnPredHDI <- function(quantiles) {
  return(quantiles$pred.yr.sum.hdi)
}

getRR <- function(quantiles) {
  return(quantiles$rr)
}
getRRiter <- function(quantiles) {
  return(quantiles$rr.iter)
}
getRRHDI <- function(quantiles) {
  return(quantiles$rr.hdi)
}
getmeanRR <- function(quantiles) {
  return(quantiles$mean_rr)
}

getsdRR <- function(quantiles) {
  return(quantiles$sd_log_rr)
}


makeInterval <-
  function(point_estimate,
           upper_interval,
           lower_interval,
           digits = 2) {
    return(paste(
      round(as.numeric(point_estimate), digits),
      ' (',
      round(as.numeric(lower_interval), digits),
      ', ',
      round(as.numeric(upper_interval), digits),
      ')',
      sep = ''
    ))
  }

#Plot predictions.
plotPred <-
  function(pred_quantiles,
           time_points,
           post_period,
           ylim,
           outcome_plot,
           title = NULL,
           sensitivity_pred_quantiles = NULL,
           sensitivity_title = 'Sensitivity Plots',
           plot_sensitivity = FALSE) {
    post_period_start <- which(time_points == post_period[1])
    post_period_end <- which(time_points == post_period[2])
    post_dates <-
      c(time_points[post_period_start:post_period_end], rev(time_points[post_period_start:post_period_end]))
    
    if (!plot_sensitivity) {
      pred_plot <- ggplot() +
        #geom_polygon(data = data.frame(time = c(post_dates, rev(post_dates)), pred_bound = c(pred_quantiles[which(time_points %in% post_dates), 3], rev(pred_quantiles[which(time_points %in% post_dates), 1]))), aes_string(x = 'time', y = 'pred_bound', color='variable'), alpha = 0.3) +
        geom_ribbon(aes(
          x = time_points[post_period_start:post_period_end],
          ymin = pred_quantiles[which(time_points %in% post_dates), 1],
          ymax = pred_quantiles[which(time_points %in% post_dates), 3]
        ),
        alpha = 0.2) +
        
        geom_line(
          data = data.frame(time = time_points, outcome = outcome_plot),
          aes_string(x = 'time', y = 'outcome')
        ) +
        geom_line(
          data = data.frame(time = time_points[1:(post_period_start - 1)], pred_outcome = pred_quantiles[1:(post_period_start -
                                                                                                              1), 2]),
          aes_string(x = 'time', y = 'pred_outcome'),
          linetype = 'dashed',
          color = 'red'
        ) +
        geom_line(
          data = data.frame(time = time_points[post_period_start:post_period_end], pred_outcome = pred_quantiles[post_period_start:post_period_end, 2]),
          aes_string(x = 'time', y = 'pred_outcome'),
          linetype = 'dashed',
          color = 'white'
        ) +
        
        labs(x = 'Time', y = 'Number of Cases') +
        #scale_colour_manual(values = c('black', 'white')) +
        scale_fill_hue(guide = 'none') +
        ggtitle(title) +
        theme_bw() +
        theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
        ) +
        theme(
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      return(pred_plot)
    } else if (!is.null(sensitivity_pred_quantiles)) {
      sensitivity_df <-
        data.frame(
          'Outcome' = outcome_plot,
          'Estimate' = pred_quantiles[, 2],
          'Sensitivity 1' = sensitivity_pred_quantiles[[1]][, 2],
          'Sensitivity 2' = sensitivity_pred_quantiles[[2]][, 2],
          'Sensitivity 3' = sensitivity_pred_quantiles[[3]][, 2],
          check.names = FALSE
        )
      sensitivity_bound <-
        data.frame(
          'Sensitivity 1' = c(
            sensitivity_pred_quantiles[[1]][which(time_points %in% post_dates), 3],
            rev(sensitivity_pred_quantiles[[1]][which(time_points %in% post_dates), 1])
          ),
          'Sensitivity 2' = c(
            sensitivity_pred_quantiles[[2]][which(time_points %in% post_dates), 3],
            rev(sensitivity_pred_quantiles[[2]][which(time_points %in% post_dates), 1])
          ),
          'Sensitivity 3' = c(
            sensitivity_pred_quantiles[[3]][which(time_points %in% post_dates), 3],
            rev(sensitivity_pred_quantiles[[3]][which(time_points %in% post_dates), 1])
          ),
          check.names = FALSE
        )
      
      pred_plot <- ggplot() +
        geom_polygon(
          data = melt(sensitivity_bound, id.vars = NULL),
          aes_string(
            x = rep(post_dates, ncol(sensitivity_bound)),
            y = 'value',
            fill = 'variable'
          ),
          alpha = 0.3
        ) +
        geom_line(
          data = melt(sensitivity_df, id.vars = NULL),
          aes_string(
            x = rep(time_points, ncol(sensitivity_df)),
            y = 'value',
            color = 'variable'
          )
        ) +
        scale_colour_manual(values = c('black', '#e41a1c', '#4daf4a', '#4daf4a', '#984ea3')) +
        scale_fill_hue(guide = 'none') +
        labs(x = 'Time', y = 'Number of Cases') +
        ggtitle(sensitivity_title) +
        theme_bw() +
        theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
        ) +
        theme(
          legend.title = element_blank(),
          legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.background = element_rect(colour = NA, fill = 'transparent'),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      return(pred_plot)
    }
  }

#Plot aggregated predictions.
#' @importFrom ggplot2 scale_x_date waiver
plotPredAgg <-
  function(ann_pred_quantiles,
           time_points,
           year_def,
           intervention_date,
           post_period,
           ylim,
           outcome_plot,
           title = NULL,
           sensitivity_pred_quantiles = NULL,
           sensitivity_title = 'Sensitivity Plots',
           plot_sensitivity = FALSE) {
    
    # See rrPredQuantiles about the meaning of time values. Short version, numeric X for cal_years and factor "X/X+1" for epi_years.

    # Intervention marker is a vertical line placed midway between 
    # a. the last complete year preceding the post period
    # b. the first year overlapping the post period
    
    # Convert intervention date (a date) into the corresponding X-axis value
    if (year_def == 'epi_year') {
      # epi_year. Years span from July of year N to Jun of year N+1, and values are N from rrPredQuantiles.
      # Marker goes 0.5 years before the beginning of the first year overlapping the post period (identified from intervention_date+1), 
      # which is (approximately) the same as being at the beginning of year N
      month.post <- month(intervention_date + 1)
      year.post <- year(intervention_date + 1)
      year.post = ifelse(month.post <= 6, year.post - 1, year.post)
      # Time axis labels are "year/year" for epi years
      date.labels = function(date) { 
        paste0(year(date), '/', year(date) + 1)
      }
    } else {
      # cal_year. Time axis is numeric. Years span from N to N+1. Marker goes 0.5 years before beginning of year containing start of post period, which is at intervention_date+1
      year.post = year(intervention_date + 1)
      # Default time axis labels are fine for calendar years
      date.labels = waiver()
    }
    # Put the intervention marker at 0.5 years before the first post-intervention year
    year.intervention = as.Date(paste0(year.post - 1, "-07-01"))

    #how mny time points in each aggregation period?
    if (year_def == 'epi_year') {
      yrvec <- year(time_points)
      monthvec <- month(time_points)
      epiyrvec = yrvec
      epiyrvec[monthvec <= 6] = yrvec[monthvec <= 6] - 1
      n.months.year <- as.vector(table(epiyrvec))
    } else{
      n.months.year <- as.vector(table(year(time_points)))
    }
    
    year.date <- function(y) {
      as.Date(paste0(y, "-01-01"), "%Y-%m-%d")
    }
    
    ann_pred_quantiles[, c('2.5%', '50%', '97.5%', 'obs')] <-
      ann_pred_quantiles[, c('2.5%', '50%', '97.5%', 'obs')] * 12 / n.months.year # for partial years, inflate the counts proportional to N months
    pred_plot <- ggplot() +
      #geom_polygon(data = data.frame(time = c(post_dates, rev(post_dates)), pred_bound = c(pred_quantiles[which(time_points %in% post_dates), 3], rev(pred_quantiles[which(time_points %in% post_dates), 1]))), aes_string(x = 'time', y = 'pred_bound', color='variable'), alpha = 0.3) +
      geom_ribbon(
        data = ann_pred_quantiles,
        aes(
          x = year.date(year),
          ymin = ann_pred_quantiles$'2.5%',
          ymax = ann_pred_quantiles$'97.5%'
        ),
        alpha = 0.5,
        fill = 'lightgray'
      ) +
      geom_point(data = ann_pred_quantiles, aes_string(x = year.date(ann_pred_quantiles$year), y = ann_pred_quantiles$obs)) +
      geom_line(
        data = ann_pred_quantiles,
        aes_string(x = year.date(ann_pred_quantiles$year), y = ann_pred_quantiles$'50%'),
        linetype = 'solid',
        color = 'black'
      ) +
      labs(x = 'Year', y = 'Number of Cases') +
      geom_vline(xintercept = year.intervention,
                 linetype = 'dashed',
                 color = 'gray') +
      ylim(0, max(ann_pred_quantiles$'97.5%')) +
      ggtitle(title) +
      theme_bw() +
      scale_x_date(labels = date.labels) +
      #scale_x_continuous(breaks = 1:nrow(ann_pred_quantiles), labels = levels(ann_pred_quantiles$year)) +
      theme(
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = -45, hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    return(pred_plot)
  }

#Sensitivity analysis by dropping the top weighted covariates.
weightSensitivityAnalysis <-
  function(group,
           covars,
           ds,
           impact,
           time_points,
           intervention_date,
           n_seasons,
           outcome,
           n_iter = 10000,
           eval_period = NULL,
           post_period = NULL,
           year_def,
           burnN=2000,
           sampleN=8000) {
    par(mar = c(5, 4, 1, 2) + 0.1)
    covar_df <- as.matrix(covars[[group]])
    #colnames(covar_df)<-substring(colnames(covar_df), 2)
    
    incl_prob <- impact[[group]]$inclusion_probs[-c(1:(n_seasons - 1)), ]
    incl_prob <- incl_prob[order(incl_prob$inclusion_probs), ]
    max_var <- as.character(incl_prob$covar.names[nrow(incl_prob)])
    max_prob <- round(incl_prob$inclusion_probs[nrow(incl_prob)], 2)
    sensitivity_analysis <- vector('list', 3)
    
    for (i in 1:(min(nrow(incl_prob), 3))) {
      covar_df <- covar_df[, colnames(covar_df) != max_var, drop = FALSE]
      covar_df.pre <-
        covar_df[time_points < as.Date(intervention_date), ]
      
      #Combine covars, outcome, date
      y <- outcome[, group]
      y.pre <- outcome[time_points < as.Date(intervention_date), group]
      post_period_response <- outcome[, group]
      post_period_response <-
        as.vector(post_period_response[time_points >= as.Date(intervention_date)])
      cID <-
        seq_along(y.pre) #used for observation-level random effect
      deltafix.mod <- rep(0, times = (ncol(covar_df.pre)))
      deltafix.mod[1:(n_seasons - 1)] <- 1 #fix  monthly dummies
      bsts_model  <-
        poissonBvs(
          y = y.pre ,
          X = covar_df.pre,
          model = list(
            deltafix = deltafix.mod,
            ri = TRUE,
            clusterID = cID
          ),
          mcmc=list(
            burnin=burnN,
            M=sampleN
        ))
      
      beta.mat <- bsts_model$samplesP$beta[-c(1:burnN), ]
      #Generate  predictions with prediction interval
      disp <- bsts_model$samplesP$thetaBeta[-c(1:burnN), ]
      disp.mat <- rnorm(n = length(disp) * length(y),
                        mean = 0,
                        sd = abs(disp))
      disp.mat <- t(matrix(disp.mat, nrow = length(disp), ncol = length(y)))
      x.fit <- cbind(rep(1, nrow(covar_df)), covar_df)
      reg.mean <-   exp((x.fit %*% t(beta.mat)) + disp.mat)
      predict.bsts <- rpois(length(reg.mean), lambda = reg.mean)
      predict.bsts <-
        matrix(predict.bsts,
               nrow = nrow(reg.mean),
               ncol = ncol(reg.mean))
      
      incl.probs.mat <- t(bsts_model$samplesP$pdeltaBeta)
      inclusion_probs <- apply(incl.probs.mat, 1, mean)
      # summary.pois<- summary(bsts_model)
      covar.names <- dimnames(covar_df.pre)[[2]]
      inclusion_probs <- cbind.data.frame(covar.names, inclusion_probs)
      
      impact_sens <-
        list(
          predict.bsts,
          inclusion_probs,
          post.period.response = post_period_response,
          observed.y = outcome[, group]
        )
      names(impact_sens) <-
        c('predict.bsts',
          'inclusion_probs',
          'post_period_response',
          'observed.y')
      sensitivity_analysis[[i]] <-
        list(removed_var = max_var, removed_prob = max_prob)
      quantiles <-
        rrPredQuantiles(
          impact = impact_sens,
          eval_period = eval_period,
          post_period = post_period,
          n_seasons = n_seasons,
          year_def = year_def,
          time_points = time_points
        )
      sensitivity_analysis[[i]]$rr <- round(quantiles$rr, 2)
      sensitivity_analysis[[i]]$pred <- quantiles$pred
      
      #Set up for next cycle  tt
      incl_prob <-
        incl_prob[incl_prob$covar.names != max_var,] #EXCLUDE TOP VAR HERE from cycle i here
      incl_prob <- incl_prob[order(incl_prob$inclusion_probs), ]
      max_var <- as.character(incl_prob$covar.names[nrow(incl_prob)])
      max_prob <- round(incl_prob$inclusion_probs[nrow(incl_prob)], 2)
    }
    return(sensitivity_analysis)
  }

# predSensitivityAnalysis <- function(group, ds, zoo_data, denom_name, outcome_mean, outcome_sd, intervention_date, eval_period, post_period, time_points, n_seasons , n_pred, year_def) {
#   impact <- doCausalImpact(zoo_data[[group]], intervention_date, time_points, n_seasons, n_pred = n_pred)
#   quantiles <- lapply(group, FUN = function(group) {rrPredQuantiles(impact = impact,  eval_period = eval_period, post_period = post_period, year_def, time_points)})
#   rr_mean <- t(sapply(quantiles, getRR))
#   return(rr_mean)
# }
#
sensitivityTable <-
  function(group,
           sensitivity_analysis,
           original_rr = NULL) {
    top_controls <-
      lapply(
        1:length(sensitivity_analysis[[group]]),
        FUN = function(i) {
          top_control <-
            c(
              sensitivity_analysis[[group]][[i]]$removed_var,
              sensitivity_analysis[[group]][[i]]$removed_prob,
              sensitivity_analysis[[group]][[i]]$rr
            )
          names(top_control) <-
            c(
              paste('Top Control', i),
              paste('Inclusion Probability of Control', i),
              paste(names(sensitivity_analysis[[group]][[i]]$rr), i)
            )
          return(top_control)
        }
      )
    sensitivity_table <-
      c(original_rr[group,], c(top_controls, recursive = TRUE))
    return(sensitivity_table)
  }

###########################
#Functions for STL+PCA piece
#from STL package:
nextodd <- function(x) {
  x <- round(x)
  if (x %% 2 == 0)
    x <- x + 1
  as.integer(x)
}
DoSTL_trend <- function(new, t.windows, s.windows, n_seasons) {
  trend <- as.data.frame(matrix(NA, nrow = nrow(new), ncol = ncol(new)))
  for (j in 1:ncol(new)) {
    ts <- ts(new[, j], frequency = n_seasons)
    trend[, j] <-
      as.vector(stl(ts, s.window = s.windows, t.window = t.windows)[[1]][, 2])
  }
  colnames(trend) <-
    c(paste(colnames(new), ".trend.", t.windows, sep = ""))
  return(trend)
}

smooth_func <- function(ds.list, covar.list, n_seasons) {
  t.windows <-
    c(nextodd(0.04 * nrow(ds.list)),
      nextodd(0.2 * nrow(ds.list)),
      nextodd(0.5 * nrow(ds.list)))
  covars.raw.compile <- vector("list", length(t.windows))
  s.windows <- "periodic"
  # STL
  for (value_t in 1:length(t.windows)) {
    for (value_s in 1:length(s.windows)) {
      t <- t.windows[value_t]
      s <- s.windows[value_s]
      covar.noseason <- covar.list[, -c(1:(n_seasons - 1))]
      stl.covars <- DoSTL_trend(covar.noseason, t, s, n_seasons)
      covars.raw.compile[[value_t]] <- stl.covars
    }
  }
  covars.raw2 <- do.call(cbind, covars.raw.compile)
  covars.raw2 <- cbind(covar.list, covars.raw2)
  covars.stl <-
    covars.raw2 #COMBINE ALL VARIABLES WITH DIFFERENT SMOOTHING LEVELS, FROM RAW TO VERY SMOOTH
}

stl_data_fun <-
  function(covars,
           ds.sub,
           n_seasons,
           outcome_name,
           post.start.index) {
    aic.test <- vector(mode = "numeric", length = ncol(covars))
    V <-
      vector("list",  length = ncol(covars)) #combine models into a list
    pred.mean <-
      vector("list",  length = ncol(covars)) #combine models into a list
    pred.coefs <-
      vector("list",  length = ncol(covars)) #combine models into a list
    covar.test <-
      vector("list",  length = ncol(covars)) #combine models into a list
    #coef1<- vector("list",  length=ncol(covars)) #combine models into a list
    preds.stage2 <-
      vector("list",  length = ncol(covars)) #combine models into a list
    mod1 <-
      vector("list",  length = ncol(covars)) #combine models into a list
    pred.mean <-
      matrix(NA, nrow = nrow(covars), ncol = ncol(covars)) #combine models into a list
    aic.test[] <- NA
    covars.no.season <- covars[, -c(1:(n_seasons - 1))]
    covars.season <- covars[, c(1:(n_seasons - 1))]
    combos3 <- dimnames(covars.no.season)[[2]]
    ds.fit <-
      vector("list",  length = length(combos3)) #combine models into a list
    
    data.fit <- cbind.data.frame(ds.sub[, outcome_name], covars)
    names(data.fit)[1] <- 'y'
    data.fit <- data.fit[1:(post.start.index - 1), ]
    for (p in 1:length(combos3)) {
      incl.names <- c('y', names(covars.season), combos3[p])
      keep.cols <- which(names(data.fit) %in% incl.names)
      ds.fit[[p]] <- data.fit[, keep.cols]
      comment(ds.fit[[p]]) <- combos3[p]
    }
    return(ds.fit)
  }

modelsize_func <- function(ds, n_seasons) {
  mod1 <- ds$beta.mat
  nonzero <-
    apply(mod1, 1, function(x)
      sum(x != 0) - n_seasons) #How many non-zero covariates, aside from seasonal dummies and intercept?
  incl_prob <-
    apply(mod1, 2, function(x)
      mean(x != 0)) #How many non-zero covariates, aside from seasonal dummies and intercept?
  modelsize <-
    round(mean(nonzero), 2) #takes the mean model size among the 8000 MCMC iterations
  return(modelsize)
}

glm.fun <- function(ds.fit, post.start.index) {
  names(ds.fit) <- paste0('c', names(ds.fit))
  covars.fit <- ds.fit[, -1]
  pre.index <- 1:(post.start.index - 1)
  fixed.effects <- paste(names(covars.fit), collapse = "+")
  ds.fit$obs <- as.factor(1:nrow(ds.fit))
  form1 <- as.formula(paste0('cy~', fixed.effects, "+ (1|obs)"))
  mod1 <-
    glmer(
      form1,
      data = ds.fit[pre.index, ],
      family = 'poisson',
      control = glmerControl(optimizer = "bobyqa",
                             optCtrl =
                               list(maxfun = 2e5))
    )
  pred.mean <- predict(mod1, newdata = ds.fit, re.form = NA)
  aic.test <- AIC(mod1)
  test.var <-
    attributes(ds.fit)$comment  #THIS IS DIFFERENT FOR BIVARIATE
  glm.out <-
    list(pred.mean, ds.fit, mod1, aic.test, test.var) #save output in a named list
  names(glm.out) <-
    c('pred.mean', 'ds.fit.fun', 'mod1', 'aic.test', 'test.var')
  return(glm.out)
}

pca_top_var <-
  function(glm.results.in,
           covars,
           ds.in,
           outcome_name,
           season.dummies) {
    #Extract AICs from list into dataframe
    aics <-
      unlist(lapply(glm.results.in, '[[', 'aic.test'))  # This returns a vector with AIC score
    vars <-
      unlist(lapply(glm.results.in, '[[', 'test.var'))  # This returns a vector with the variable names
    pred.mean <-
      lapply(glm.results.in, '[[', 'pred.mean') # This returns a vector with the variable names
    pred.mean <- do.call(cbind, pred.mean)
    pred.mean <- exp(pred.mean)
    aic.df <- cbind.data.frame(vars, aics)
    names(aic.df) <- c('covars', 'aic')
    aic.df$model.index <- 1:nrow(aic.df)
    aic.df$grp <-
      as.numeric(as.factor(substr(aic.df$covars, 1, 3))) #for each smoothed or unsmoothed version of variable, assign it to a grouping
    aic.df$delta.aic <- aic.df$aic - min(aic.df$aic)
    aic.df$w_aic <-
      exp(-0.5 * aic.df$delta.aic) / sum(exp(-0.5 * aic.df$delta.aic))
    aic.df <- aic.df[order(-aic.df$w_aic), ]
    aic.df$cumsum <- cumsum(aic.df$w_aic)
    aic.df$keep.high.weight <-
      aic.df$cumsum <= 0.99  #only keep variables that contribute to 99% f weight
    aic.df$model.rank <-
      1:nrow(aic.df)  #only keep variables that contribute to 99% f weight
    
    #top.covar in each set
    aic.df2 <- aic.df
    aic.df2 <- aic.df2[order(-aic.df$w_aic), ]
    aic.df2$deseason = 0
    aic.df2$deseason[grep('trend', aic.df2$covars, fixed = TRUE)] <-
      1
    aic.df2 <-
      aic.df2[aic.df2$deseason == 1, ] #Only keep STL version of variable, not raw
    aic.df2$w_aic <- aic.df2$w_aic / sum(aic.df2$w_aic) #rescale weights
    top.covar.grp <- aic.df2[!duplicated(aic.df2$grp), ]
    remove <- c('t', 'nocovars')
    top.covars <-
      as.character(top.covar.grp$covars[!top.covar.grp$covars %in% remove])
    ###SETUP AND RUN MODELS WITH FIRST PC
    covars.keep.pca <- covars[, top.covars]
    #Run PCA
    pca <-
      prcomp(covars.keep.pca, scale = TRUE) # scale=TRUE should be added!!
    predictors2 <- as.data.frame(pca$x[, 1]) # First "1" PC
    names(predictors2) <- 'pca1'
    y = ds.in[, outcome_name]
    covar.matrix.pca <-
      cbind.data.frame(y, season.dummies, predictors2)
    #covar.matrix.pca$obs<-as.factor(1:nrow(covar.matrix.pca))
    return(covar.matrix.pca)
  }

cumsum_func <-
  function(group,
           quantiles,
           outcome,
           time_points,
           post_period, 
           hdi=FALSE) {
    is_post_period <- which(time_points >= post_period[1])
    is_pre_period <- which(time_points < post_period[1])
    
    #Cumulative sum of prevented cases
    cases_prevented <-
      quantiles[[group]]$pred_samples - outcome[, group]
    cumsum_cases_prevented_post <-
      apply(cases_prevented[is_post_period,], 2, cumsum)
    cumsum_cases_prevented_pre <-
      matrix(0,
             nrow = nrow(cases_prevented[is_pre_period,]),
             ncol = ncol(cases_prevented[is_pre_period,]))
    cumsum_cases_prevented <-
      rbind(cumsum_cases_prevented_pre, cumsum_cases_prevented_post)
    cumsum_prevented <-
      t(apply(
        cumsum_cases_prevented,
        1,
        quantile,
        probs = c(0.025, 0.5, 0.975),
        na.rm = TRUE
      ))
    if(hdi==TRUE){
    cumsum_prevented <- cbind(cumsum_prevented[,'50%'],  t(hdi( t(cumsum_cases_prevented),credMass = 0.95 )))
    }
    return(cumsum_prevented)
  }

#Classic its setup

its_func <- function(ds1,
                     time_points,
                     post_period,
                     eval_period) {
  ds1$post1 <- 0
  ds1$post1[time_points >= post_period[1] &
              time_points < eval_period[1]] <- 1
  ds1$post2 <- 0
  ds1$post2[time_points >= eval_period[1]] <- 1
  ds1$time_index <- ds1$time_index / max(ds1$time_index)
  ds1$obs <- 1:nrow(ds1)
  if (max(ds1$log.offset) > min(ds1$log.offset)) {
    ds1$log.offset <- scale(ds1$log.offset)
  }
  eval_indices <-
    match(which(time_points == eval_period[1]), (1:length(ds1$obs))):match(which(time_points ==
                                                                                   eval_period[2]), (1:length(ds1$obs)))
  #Fit classic ITS model
  mod1 <-
    glmer(
      outcome ~ s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8 + s9 + s10 + s11 + log.offset +
        time_index + post1 + post2 +
        time_index * post1 + time_index * post2 + (1 |
                                                     obs),
      data = ds1,
      family = poisson(link = log),
      control = glmerControl(optimizer = "bobyqa",
                             optCtrl =
                               list(maxfun = 2e5))
    )
  #GENERATE PREDICTIONS
  covars3 <-
    cbind(ds1[, c(
      's1',
      's2',
      's3',
      's4',
      's5',
      's6',
      's7',
      's8',
      's9',
      's10',
      's11',
      'log.offset',
      'time_index',
      'post1',
      'post2'
    )], ds1$time_index * ds1$post1,
    ds1$time_index * ds1$post2)
  # If log.offset is degenerate, then it was dropped by glmer and excluded from mod1, and
  # therefore it should also be excluded from covars3 (otherwise matrix multiplication will fail a few lines below)
  if (max(ds1$log.offset) == min(ds1$log.offset)) {
    covars3 = subset(covars3, select=-log.offset)
  }
  covars3 = as.matrix(covars3)
  covars3 <- cbind.data.frame(rep(1, times = nrow(covars3)), covars3)
  names(covars3)[1] <- "Intercept"
  pred.coefs.reg.mean <-
    mvrnorm(n = 10000,
            mu = fixef(mod1),
            Sigma = vcov(mod1))
  preds.stage1.regmean <-
    exp(as.matrix(covars3) %*% t(pred.coefs.reg.mean))
  #re.sd<-as.numeric(sqrt(VarCorr(mod1)[[1]]))
  #preds.stage1<-rnorm(n<-length(preds.stage1.regmean), mean=preds.stage1.regmean, sd=re.sd)
  #preds.stage1<-exp(matrix(preds.stage1, nrow=nrow(preds.stage1.regmean), ncol=ncol(preds.stage1.regmean)))
  
  #Then for counterfactual, set post-vax effects to 0.
  covars3.cf <-
    cbind(ds1[, c(
      's1',
      's2',
      's3',
      's4',
      's5',
      's6',
      's7',
      's8',
      's9',
      's10',
      's11',
      'log.offset',
      'time_index'
    )], matrix(
      0, nrow = nrow(ds1), ncol = 4
    ))
  if (max(ds1$log.offset) == min(ds1$log.offset)) {
    covars3.cf = subset(covars3.cf, select=-log.offset)
  }
  covars3.cf = as.matrix(covars3.cf)
  covars3.cf <-
    cbind.data.frame(rep(1, times = nrow(covars3.cf)), covars3.cf)
  preds.stage1.regmean.cf <-
    exp(as.matrix(covars3.cf) %*% t(pred.coefs.reg.mean))
  #preds.stage1.cf<-rnorm(n<-length(preds.stage1.regmean.cf), mean=preds.stage1.regmean.cf, sd=re.sd)
  #preds.stage1.cf<-exp(matrix(preds.stage1.cf, nrow=nrow(preds.stage1.regmean.cf), ncol=ncol(preds.stage1.regmean.cf)))
  
  rr.t <- preds.stage1.regmean / preds.stage1.regmean.cf
  rr.q.t <- t(apply(rr.t, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  
  preds.stage1.regmean.SUM <-
    apply(preds.stage1.regmean[eval_indices, ], 2, sum)
  preds.stage1.regmean.cf.SUM <-
    apply(preds.stage1.regmean.cf[eval_indices, ], 2, sum)
  rr.post <- preds.stage1.regmean.SUM / preds.stage1.regmean.cf.SUM
  rr.q.post <- quantile(rr.post, probs = c(0.025, 0.5, 0.975))
  #matplot(rr.q, type='l', bty='l', lty=c(2,1,2), col='gray')
  #abline(h=1)
  
  # rr.q.last<- rr.q.t[nrow(rr.q.t),]
  rr.out <- list(rr.q.post = rr.q.post, rr.q.t = rr.q.t)
  return(rr.out)
}

#' @importFrom lme4 VarCorr

single.var.glmer<-function(ds1, ds.labels, intro.date, time_points,n_seasons, eval.period){
  #GLMER
  outcome.pre<-ds1[,'outcome']
  outcome.pre[as.Date(time_points)>=intro.date] <-NA
  covars<-names(ds1)[-c(1:n_seasons)] #names of variables to test
  rr.post.q.glmer.manual <- vector("list", length(covars)) 
  aic.summary <- vector("list", length(covars)) 
  for(i in 1:length(covars)){
  covar1<-ds1[,covars[i]]
  eval.start.index<-which(time_points==eval.period[1])
  months<-month(time_points)
  season.dummies<-   dummies::dummy(months)[,1:(n_seasons-1)]
  dimnames(season.dummies)[[2]]<-paste0('s',1:(n_seasons-1))
  ds2<-cbind.data.frame(outcome.pre, season.dummies, scale(covar1))
  names(ds2)[ncol(ds2)]<-covars[i]
  ds2$obs<-as.factor(1:nrow(ds2))
  form1<-as.formula(paste0('outcome.pre~s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+' ,covars[i], '+(1|obs)')) 
                    mod1<-glmer(form1 , family=poisson(link=log) , data=ds2)
                    #Manually calculate CIs
                    aic.summary[[i]]<-AIC(mod1) #only for reference--tough to use for mixed model
                     covars3<-as.matrix(ds2[c('s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11',covars[i])])
                    covars3<-cbind.data.frame(rep(1, times=nrow(covars3)), covars3)
                    names(covars3)[1]<-"Intercept"
                    pred.coefs.reg.mean<- mvrnorm(n = 100, mu=fixef(mod1), Sigma=vcov( mod1))
                    re.sd<-as.numeric(sqrt(VarCorr(mod1)[[1]]))
                    preds.stage1.regmean<- as.matrix(covars3) %*% t(pred.coefs.reg.mean) 
                    re.int<-rnorm(n<-length(preds.stage1.regmean), mean=0, sd=re.sd) 
                    preds.stage1.regmean<-preds.stage1.regmean+re.int
             
                    preds.stage2<-rpois(n=length(preds.stage1.regmean)*100, exp(preds.stage1.regmean))
                    preds.stage2<-matrix(preds.stage2, nrow=nrow(preds.stage1.regmean), ncol=ncol(preds.stage1.regmean)*100)
                    
                    post.preds1.manual<-preds.stage2[eval.start.index:nrow(preds.stage1.regmean),]
                    post.preds.sums1.manual<-apply(post.preds1.manual,2,sum)
                    post.obs.sum<-  sum(ds1[ eval.start.index:nrow(preds.stage1.regmean),'outcome'])
                    post.rr1.manual<-post.obs.sum/post.preds.sums1.manual
                    rr.post.q.glmer.manual[[i]]<-quantile(post.rr1.manual,probs=c(0.025,0.5,0.975))
  }
  names(rr.post.q.glmer.manual)<-covars
  aic.summary<-unlist(aic.summary)
  rr.post.q.glmer.manual<-round(matrix(unlist(rr.post.q.glmer.manual), ncol=3, byrow=T),2)
  results<-list('aic.summary'=aic.summary,'rr.post'=rr.post.q.glmer.manual,'covar.names'=covars)
  return(results)
}

#' Generate plot of univariate results
#' @importFrom graphics plot axis abline arrows
#' @param ds Object imported from evaluatr.univariate
#' @return Univariate analysis plot, `results`, as described below
#' @export

evaluatr.univariate.plot<-function(ds, plot.labs='Univariate'){
  plot(y=1:nrow(ds), x=ds$rr, bty='l',yaxt='n', pch=16 , xlim=c(0.2,2), 
        ylab='' , xlab='Univariate Rate Ratio',ylim=c(nrow(ds),1), main=plot.labs)
  arrows(y0=1:nrow(ds), x0=ds$rr.lcl,x1=ds$rr.ucl,  length =0)
  axis(side=2, at=1:length(ds$covar), labels=ds$covar,las=1, cex.axis=0.6)
  abline(v=1, lty=2, col='gray')
}
