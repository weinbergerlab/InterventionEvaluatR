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
    #Filter columsn with 0 variations in the covariate in the pre-vax period
    x.var<-apply(covars,2, function(xx) var(xx))  
    x.var[is.na(x.var)]<-0
    covars<-covars[,x.var>0] 
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


#Main function
inla_mods<-function(zoo_data,
                    intervention_date,
                    n_seasons,
                    time_points,
                    error_dist,
                    denom_name,
                    eval_period,
                    post_period,
                    model.variant){
  
  y <- zoo_data[, 1] #all y
  y.pre <-y
  y.pre[time_points >= as.Date(intervention_date)]<-NA
  post_period_response<- y[time_points >= as.Date(intervention_date)]
  x.all <-as.matrix(zoo_data[,-c(1)]) #Removes outcome column from dataset
  
  y.pre[time_points>=intervention_date]<-NA
  
  #Filter columsn with 0 variations in the covariate in the pre-vax period
  #x.var<-apply(x,2, function(xx) var(xx))  
  #x.var[is.na(x.var)]<-0
  #x<-x[,x.var>0] 
  
  ##NEEDTO SEPRATE OUT SEASONAL DUMMIES FROM COVARS
  x<-x.all[,-(1:(n_seasons-1))]
  x.months<-x.all[,(1:(n_seasons-1))]
  dimnames(x.months)[[2]]<-paste0('season', 1:(n_seasons-1))
  
  #Should be already scaled, but doesn't hurt...
  x.scale<-apply(x,2, function(z) scale(z)) 
  
  y.aware.scale<- apply(x, 2, function(x1){
    log.y.pre.scale<- scale(log(y.pre+0.5))
    log.y.pre.scale<-log.y.pre.scale[!is.na(log.y.pre.scale)]
    reg<-lm(log.y.pre.scale~x1[time_points<intervention_date])
    slope<- reg$coefficients[2]
    x.scale<-x1*slope - mean(x1*slope)
    return(x.scale)
  })
  pca1<- prcomp(y.aware.scale, center = FALSE,scale. = FALSE)
  #n.pcs.keep<-sum(pca1$sdev>1)
  n.pcs.keep<-2
  pcs<-pca1$x
  pcs<-apply(pcs,2, scale) #SCALE THE PCS prior to regression!
  pcs.combo<-cbind.data.frame( 'month'=as.factor(month(time_points)),pcs[,1:min(n.pcs.keep, ncol(pcs)), drop=F] )
  form0<-as.formula(paste0('~', paste(names(pcs.combo), collapse='+')))
  pcs.df<-as.data.frame(model.matrix(form0, data=pcs.combo)[,-1])
  
  df1<-cbind.data.frame(y.pre, 'month'=as.factor(month(time_points)), y, x.scale )
  df1$t<-1:nrow(df1)
  covar.df.full<-cbind.data.frame( x.scale, 'month'=df1$month) 
  
  #Matrix for restriction random effects
  x.in.full<-as.data.frame(model.matrix(~. , data=covar.df.full))[,-1] #dummies for month, remove intercept
  mod.df.full<- cbind.data.frame(y,y.pre, x.in.full)
  
  if(model.variant=='full'){  
    
    A.full<- rbind(t(pcs.df[,names(pcs.df), drop=F]))
    A.full[,is.na(mod.df.full$y.pre)]<-0  #extraploation period shouldn't factor into constraint
    e.full=rep(0, nrow(A.full))
    
    n <- nrow(mod.df.full)
    X <- matrix(1,nrow = n, ncol= 1)
    Z <- as.matrix(x.in.full) 
    pX = ncol(X) 
    pZ = ncol(Z)
    idx.X = c(1:pX, rep(NA, pZ))
    idx.Z = c(rep(NA,pX), 1:pZ)
    hyper.fixed = list(prec = list(initial = log(0.00001), fixed=TRUE))
    param.beta = list(prec = list(param = c(1.0e-3, 1.0e-3))) 
    param.Z =  param.beta   
    param.data = list(prec = list(param = c(1.0e-3, 1.0e-3)))
    t=1:nrow(x.in.full)
    #http://www.r-inla.org/faq
    add.beta.cols<-as.data.frame(matrix(rep(1:ncol(x.in.full), each=nrow(x.in.full)), nrow=nrow(x.in.full)))
    names(add.beta.cols)<-paste0('beta', 1:ncol(x.in.full))
    mod.df.full<-cbind.data.frame(mod.df.full,add.beta.cols)
    #Setup for Ridge regression
    #https://bitbucket.org/hrue/r-inla/src/0d9389b414979bcf3507c60f69a50488d18e4c6c/r-inla.org/examples/stacks/stacks.R?at=default
    x.in.full.no.season<-x.in.full[,-grep('month', names(x.in.full))]
    x.in.full.season<-x.in.full[,grep('month', names(x.in.full))]
    first.beta<-paste0("f(beta1,", names(x.in.full.no.season)[1], ",model='iid',hyper = param.beta,", "values=c(", paste(1:ncol(x.in.full.no.season), collapse=',') ,"))" )
    next.beta.f<-function(x, covar.index){ 
      paste0("f(beta",covar.index,",", names(x.in.full.no.season)[covar.index], ",copy='beta1', fixed=T)")
    }
    betas.list<- mapply(next.beta.f, x=x.in.full.no.season[-1], covar.index=2:ncol(x.in.full.no.season), SIMPLIFY=T)
    next.betas<- paste(betas.list, collapse='+')  
    all.betas<-paste(first.beta,next.betas , sep='+')    
    if(error_dist=='ar1'){
      form1<- as.formula(paste0("y.pre ~", paste(names(x.in.full.season), collapse="+"),"+", all.betas,    "+ f(t, model = 'ar1', constr=T,extraconstr=list(A=A.full, e=e.full))") )
      #form1<- as.formula(paste0("y.pre ~", paste(names(x.in.full.season), collapse="+"),"+", all.betas,    "+ f(t, model = 'ar1')") )
    }else{
      form1<- as.formula(paste0("y.pre ~", paste(names(x.in.full.season), collapse="+"),"+", all.betas,    "+ f(t, model = 'iid')") )
      
    }
    mod1.full = inla(form1, data = mod.df.full, 
                     control.predictor = list(compute=TRUE, link = 1), 
                     family='poisson',
                     control.compute=list(config = TRUE,waic=TRUE))
    waic<-mod1.full$waic$waic
    pd<- mod1.full$waic$p.eff
    ds<-list('fitted.model'=mod1.full, 'mod.df'=mod.df.full, 'x.in'=x.in.full,'offset'=NA,'offset.used'=F,'ridge'=T, 'waic'=waic, 'pd'=pd)
  }
  
  if(model.variant %in% c('time','time_no_offset')){
    if(denom_name  %in% names(x)){
      log.offset<- x[,c(denom_name)]
    }else{
      log.offset<- rep(0, nrow(mod.df.full))
    }
    offset1<-exp(log.offset)
    x.in.time<-as.data.frame(model.matrix(~. , data=covar.df.full[,c('month'), drop=F]))[,-1] #dummies for month, remove intercept
    x.in.time$trend<-1:nrow(x.in.time)
    mod.df.time<- cbind.data.frame(y,y.pre, x.in.time)
    A.time<- rbind(t(x.in.time[,names(x.in.time), drop=F]))
    A.time[,is.na(mod.df.time$y.pre)]<-0  #extraploation period shouldn't factor into constraint
    e.time=rep(0, nrow(A.time))
    mod.df.time$t<-1:nrow(mod.df.time)
    #mod.df.time.offset<-cbind.data.frame(mod.df.time,'log.offset'=log.offset)
    if(error_dist=='ar1'){
      form2<- as.formula(paste0('y.pre ~',  paste(names(x.in.time), collapse='+'), "+  f(t, model = 'ar1', constr=T,extraconstr=list(A=A.time, e=e.time))") )
    }else{
      form2<- as.formula(paste0('y.pre ~',  paste(names(x.in.time), collapse='+'), "+  f(t, model = 'iid')") )
    }
    if(model.variant %in% c('time_no_offset')){
      mod2.time.no.offset = inla(form2, data = mod.df.time, control.predictor =   list(compute=TRUE, link = 1), family='poisson',control.compute=list(config = TRUE,waic=TRUE))
      waic<-mod2.time.no.offset$waic$waic
      pd<- mod2.time.no.offset$waic$p.eff
      ds<-list('fitted.model'=mod2.time.no.offset, 'mod.df'=mod.df.time, 'x.in'=x.in.time,'offset'=NA, 'offset.used'=F, 'ridge'=F, 'waic'=waic, 'pd'=pd)
    }else if(model.variant=='time'){
      #Time regression with offset
      mod3.time = inla(form2, data = mod.df.time, control.predictor =   list(compute=TRUE, link = 1), family='poisson',E=offset1, control.compute=list(config = TRUE,waic=TRUE))
      waic<-mod3.time$waic$waic
      pd<- mod3.time$waic$p.eff
      ds<-list('fitted.model'=mod3.time, 'mod.df'=mod.df.time, 'x.in'=x.in.time, 'offset'=offset1,'offset.used'=T,'ridge'=F, 'waic'=waic, 'pd'=pd)
    }
  }
  
  if(model.variant=='pca'){  
    x.pca <- pcs.df
    mod.df.pca<- cbind.data.frame(y,y.pre, x.pca)
    A.pca<- rbind(t(x.pca[,names(x.pca), drop=F]))
    A.pca[,is.na(mod.df.pca$y.pre)]<-0  #extraploation period shouldn't factor into constraint
    e.pca=rep(0, nrow(A.pca))
    t=1:nrow(x.pca)
    
    if(error_dist=='ar1'){
      form2<- as.formula(paste0('y.pre ~',  paste(names(x.pca), collapse='+'), "+  f(t, model = 'iid')") )
    }else{
      form2<- as.formula(paste0('y.pre ~',  paste(names(x.pca), collapse='+'), "+  f(t, model = 'ar1', constr=T,extraconstr=list(A=A.pca, e=e.pca))") )
      
    }
    mod1.pca = inla(form2, data = mod.df.pca, 
                    control.predictor = list(compute=TRUE, link = 1), 
                    family='poisson',
                    control.compute=list(config = TRUE,waic=TRUE))
    waic<-mod1.pca$waic$waic
    pd<- mod1.pca$waic$p.eff
    ds<-list('fitted.model'=mod1.pca, 'mod.df'=mod.df.pca, 'x.in'=x.pca,'offset'=NA,'offset.used'=F,'ridge'=F, 'waic'=waic, 'pd'=pd)
  }
  
  if(ds$ridge==T){  
    coefs<-ds$fitted.model$summary.random[['beta1']]
    coefs$covar<-names(ds$x.in)[-grep('month',names(ds$x.in)) ]
  }else{
    coefs<-ds$fitted.model$summary.fixed
  }
  coefs<-coefs[order(coefs$`0.5quant`),]
  
  posterior.list<-inla.posterior.sample(n=500, ds$fitted.model, seed=123)
  post.labels<-dimnames(posterior.list[[1]]$latent)[[1]]
  #post.labels<-gsub(":", "_", post.labels)
  posterior.samples<- sapply(posterior.list, '[[', 'latent')
  preds.select<-grep('Predictor',post.labels )
  rand.eff.select.t1<-which(substr(post.labels,1,2 )=='t:')
  
  if(ds$ridge==T){  
    covar.select.betas<-grep('beta1:', post.labels)
    covar.select.months<-grep('month', post.labels)
    covar.select<-c(covar.select.betas,covar.select.months)
  }else{
    covar.select<- which(names(ds$x.in)  %in% sub("\\:.*", "", post.labels))
  }
  
  beta.posterior<-(posterior.samples[covar.select,]) 
  beta.posterior.median<- apply(beta.posterior,1,median)
  beta.posterior.hdi<-cbind(beta.posterior.median,t(hdi(t(beta.posterior), credMass = 0.95)))
  row.names(beta.posterior.hdi)<-names(ds$x.in)
  
  beta.fix.posterior.hdi<-beta.posterior.hdi[-grep('month',names(ds$x.in)),]
  
  fixed.effect<- as.matrix(ds$x.in) %*% beta.posterior #fixed piece of regression, excluding intercept
  fixed.effect.hdi<-t(hdi(t(fixed.effect), credMass = 0.95))
  fixed.effect.median<-apply(fixed.effect,1, median)
  fixed.effect.hdi<-cbind.data.frame('median'=fixed.effect.median, fixed.effect.hdi)
  
  posterior.preds<-exp(posterior.samples[preds.select,]) #lambda
  #now take Poisson samples withmean of lambda
  posterior.preds.counts<- matrix(rpois(n=length(posterior.preds), lambda=posterior.preds), nrow=nrow(posterior.preds), ncol=ncol(posterior.preds))
  
  rand.eff.t1<-posterior.samples[rand.eff.select.t1,]
  rand.eff1.q<-t(apply(rand.eff.t1, 1, quantile, probs=c(0.025,0.5,0.975)))
  
  posterior.preds.q<-t(apply(posterior.preds.counts,1,quantile, probs=c(0.025, 0.5, 0.975)))
  posterior.median<-as.integer(round(t(apply(posterior.preds.counts,1,median))))
  ci<- t(hdi(t(posterior.preds.counts), credMass = 0.95))
  posterior.pred.hdi<- cbind.data.frame('median'=posterior.median, ci)
  
  rho1<-ds$fitted.model$summary.hyperpar[3,c('0.5quant','0.025quant', '0.975quant')]
  
  rand.eff.combined<-rand.eff.t1 #+rand.eff.t2
  rand.eff.combined.q<-t(apply(rand.eff.combined, 1, quantile, probs=c(0.025,0.5,0.975)))
  
  log.rr.pointwise<- apply(posterior.preds.counts, 2, function(x)log( (mod.df.full$y+1)/ (x+1)) )
  log.rr.pointwise.ci<- t(hdi(t(log.rr.pointwise), credMass = 0.95))
  log.rr.pointwise.median<- apply(log.rr.pointwise.ci,1,median)
  log_rr_full_t_hdi<-cbind(log.rr.pointwise.median,log.rr.pointwise.ci)
  # matplot(log.rr, type='l', col='gray', lty=c(2,1,2))
  # abline(h=0, col='red')
  log_rr_full_t_quantiles<- t(apply(log.rr.pointwise,1,quantile, probs=c(0.025,0.5,0.975)))
  log_rr_full_t_sd<- apply(log.rr.pointwise,1,sd)
  #log_rr_full_t_samples.prec.post<-1/log_rr_full_t_sd^2
  
  post.period<- which(time_points>= eval_period[1]  &time_points<= eval_period[2] )
  
  post.samples<-posterior.preds.counts[post.period,]
  
  post.samples.sum<-apply(post.samples,2,sum)
  obs.post.sum<- sum(mod.df.full$y[post.period])
  rr.agg<-obs.post.sum/post.samples.sum
  rr.q<-quantile(rr.agg, probs=c(0.025, 0.5, 0.975)) 
  rr.hdi<-c(rr.q['50%'],hdi(rr.agg, credMass = 0.95))
  mean_rr<-mean(rr.agg)
  
  post.samples.sum.q<-quantile(post.samples.sum, probs=c(0.025, 0.5, 0.975))
  pred.hdi <- cbind( post.samples.sum['50%'],hdi(post.samples.sum, credMass = 0.95))
  
  year<-year(time_points)  
  pred.sum.year.post<-apply(posterior.preds.counts,2, function(x) aggregate(x, by=list('year'=year), FUN=sum))
  pred.sum.year.post<-sapply(pred.sum.year.post,function(x) x[,2])
  pred.sum.year.ci<- t(hdi(t(pred.sum.year.post), credMass = 0.95))
  
  pred.sum.year<-apply(posterior.preds.counts,2, function(x) aggregate(x, by=list('year'=year), FUN=sum))
  y.agg<-aggregate(y, by=list('year'=year), FUN=sum)
  pred.sum.year<-sapply(pred.sum.year,function(x) x[,2])
  pred.yr.sum.q<- as.data.frame(t(apply(pred.sum.year,1,quantile, probs=c(0.025,0.5,0.975))))
  names(pred.yr.sum.q)<-c('2.5%','50%','97.5%')
  pred.yr.sum.q$obs <-y.agg[,2]
  pred.yr.sum.q$year <- unique(year)
  
  pred.yr.sum.hdi<-cbind.data.frame(pred.yr.sum.q[,'50%'],pred.sum.year.ci)
  pred_samples_post_full<-post.samples
  pred<-t(apply(posterior.preds.counts,1, quantile, probs=c(0.025, 0.5, 0.975)))
  
  #Cases averted
  is_post_period <- which(time_points >= post_period[1])
  is_pre_period <- which(time_points < post_period[1])
  cases_prevented<- apply(posterior.preds.counts, 2, function(x1)  x1 -y )
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
      na.rm = TRUE))
  cumsum_prevented_hdi <- cbind(cumsum_prevented[,'50%'],  t(hdi( t(cumsum_cases_prevented),credMass = 0.95 )))
  
  
  #out.list<- list('posterior.pred.hdi'=posterior.pred.hdi, 'rr.hdi'=rr.hdi,'inla.mod'=ds$fitted.model,'log.rr.pointwise.hdi'=log.rr.pointwise.hdi,'obs.y'=mod.df.full$y,'rand.eff1.q'=rand.eff1.q,'betas'=beta.posterior.hdi, 'post.period.start.index'=post.period[1], 'rand.eff.t1'=rand.eff.t1,'rand.eff.combined.q'=rand.eff.combined.q, 'fixed.effect'=fixed.effect.hdi,'beta.posterior.hdi'=beta.posterior.hdi, 'rho1'=rho1)
  #return(out.list)
  impact <-
    list(
      #   'reg.mean'= reg.mean,
      #   'exclude.indices'=exclude.indices,
      #   'rand.eff'=rand.eff,
      #   'covars'=covars,
      #   'beta.mat'= beta.mat,
      #   'predict.bsts'=predict.bsts,
      #   'inclusion_probs'=inclusion_probs,
      'post_period_response' = post_period_response,
      'observed.y' = y,
      'waic'= ds$waic,
      'p_D'=ds$pd,
      'rho'=rho1,
      'fixed.effect.hdi'=fixed.effect.hdi,
      'beta'=beta.fix.posterior.hdi,
      'rand.eff.combined.q'=rand.eff.combined.q
      
    )
  quantiles <-
    list(
      log_rr_full_t_hdi=log_rr_full_t_hdi,
      rr.hdi=rr.hdi,
      pred.yr.sum.hdi=pred.yr.sum.hdi,
      pred.hdi=pred.hdi,
      pred.yr.sum.q = pred.yr.sum.q,
      # log_rr_full_t_samples.prec.post = log_rr_full_t_samples.prec.post,
      #    pred_samples = pred_samples,
      pred = pred,
      #     rr = rr,
      #     roll_rr = roll_rr,
      mean_rr = mean_rr,
      pred_samples_post_full = pred_samples_post_full,
      #     roll_rr = roll_rr,
      log_rr_full_t_quantiles = log_rr_full_t_quantiles,
      log_rr_full_t_sd = log_rr_full_t_sd,
      rr = rr.q,
      rr.iter=rr.agg
    )
  results<-list('impact'=impact, 'quantiles'=quantiles,'cumsum_prevented_hdi'=cumsum_prevented_hdi,'cumsum_prevented'=cumsum_prevented)  
  return(results)
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
      match(which(time_points == eval_period[1]), (1:length( impact$observed.y))):match(which(time_points ==
                                                                                                eval_period[2]), (1:length( impact$observed.y)))
    
    pred_eval_sum <- colSums(pred_samples[eval_indices,])
    
    eval_obs <- sum( impact$observed.y[eval_indices])
    
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
    obs.y.year <- tapply( impact$observed.y, year, sum)
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
    pre_obs <- sum( impact$observed.y[pre_indices])
    rr_sum_pre <- pre_obs / pred_pre_sum  #Should be 0!
    
    #unbias_rr <- eval_rr_sum / rr_sum_pre # same as log_rr - log_rr_pre=log(A/B)
    #unbias_rr_q <- quantile(unbias_rr, probs = c(0.025, 0.5, 0.975))
    
    plot_rr_start <- which(time_points == post_period[1]) - n_seasons
    roll_rr_indices <-
      match(plot_rr_start, (1:length( impact$observed.y))):match(which(time_points ==
                                                                         eval_period[2]), (1:length( impact$observed.y)))
    
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
        #  unbias_rr_q = unbias_rr_q,
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
    pre_period <- which(time_points < post_period[1])
    post_period_start <- which(time_points == post_period[1])
    post_period_end <- which(time_points == post_period[2])
    post_dates <-
      c(time_points[post_period_start:post_period_end], rev(time_points[post_period_start:post_period_end]))
    
    if (!plot_sensitivity) {
      pred_plot <- ggplot() +
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
          data = data.frame(time = time_points[pre_period], pred_outcome = pred_quantiles[pre_period, 2]),
          aes_string(x = 'time', y = 'pred_outcome'),
          linetype = 'dotted',
          color = 'black'
        ) +
        geom_line(
          data = data.frame(time = time_points[post_period_start:post_period_end], pred_outcome = pred_quantiles[post_period_start:post_period_end, 2]),
          aes_string(x = 'time', y = 'pred_outcome'),
          linetype = 'dashed',
          color = 'black'
        ) +
        labs(x = 'Time', y = 'Number of Cases') +
        coord_cartesian(
          ylim=c(0, max(pred_quantiles[,3], outcome_plot))
        ) +
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
#' @importFrom ggplot2 scale_x_date waiver geom_pointrange
#' @importFrom lubridate date_decimal
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
    year.intervention = as.Date(date_decimal(year.post - 0.5))
    
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
      as.Date(date_decimal(y))
    }
    
    # ann_pred_quantiles has observed and predicted data for all years
    ann_pred_quantiles[, c('2.5%', '50%', '97.5%', 'obs')] <-
      ann_pred_quantiles[, c('2.5%', '50%', '97.5%', 'obs')] * 12 / n.months.year # for partial years, inflate the counts proportional to N months
    
    # We care about observed data from all years
    ann_obs <- ann_pred_quantiles
    
    # But we don't care about predicted data before intervention date
    ann_pred_quantiles_pre <- subset(ann_pred_quantiles, year < year.post)
    ann_pred_quantiles_post <- subset(ann_pred_quantiles, year >= year.post)
    
    # However, for display purposes, we force the predicted time series to converge on the observed time series at the point of intervention
    # last_pre_intervention = as.data.frame(t(colMeans(ann_obs[ann_obs$year %in% c(year.post - 1, year.post),])))
    last_pre_intervention = ann_obs[ann_obs$year == year.post - 1,]
    last_pre_intervention$'2.5%' = last_pre_intervention$'50%' = last_pre_intervention$'97.5%' = last_pre_intervention$obs
    ann_pred_quantiles_post <- rbind(last_pre_intervention, ann_pred_quantiles_post)
    
    pred_plot <- ggplot() +
      geom_ribbon(
        aes_(
          x = year.date(ann_pred_quantiles_post$year),
          ymin = ann_pred_quantiles_post$'2.5%',
          ymax = ann_pred_quantiles_post$'97.5%'
        ),
        alpha = 0.5,
        fill = 'lightgray'
      ) +
      geom_line(
        aes_(x = year.date(ann_obs$year), y = ann_obs$obs), 
        linetype = "solid", 
        color = "black"
      ) +
      geom_line(
        aes_(x = year.date(ann_pred_quantiles_pre$year), y = ann_pred_quantiles_pre$'50%'),
        linetype = 'dotted',
        color = 'black'
      ) +
      geom_line(
        aes_(x = year.date(ann_pred_quantiles_post$year), y = ann_pred_quantiles_post$'50%'),
        linetype = 'dashed',
        color = 'black'
      ) +
      labs(x = 'Year', y = 'Number of Cases') +
      coord_cartesian(
        ylim=c(0, max(ann_pred_quantiles_post$'97.5%', ann_obs$obs))
      ) +
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
                               list(maxfun = 2e6))
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

cumsum_func <-
  function(quantiles,
           outcome,
           time_points,
           post_period, 
           hdi=FALSE) {
    is_post_period <- which(time_points >= post_period[1])
    is_pre_period <- which(time_points < post_period[1])
    
    #Cumulative sum of prevented cases
    cases_prevented <-
      quantiles$pred_samples - outcome
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
                               list(maxfun = 2e6))
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
    months<- month(time_points)
    month.df <- cbind.data.frame('months'=factor(month(time_points)))
    #season.dummies <-   dummies::dummy(months)[,1:(n_seasons-1)]
    season.dummies <- model.matrix( ~ -1 + months, data=month.df)
    season.dummies <- season.dummies[,-1]
    dimnames(season.dummies)[[2]]<-paste0('s',1:(n_seasons-1))
    
    ds2<-cbind.data.frame(outcome.pre, season.dummies, scale(covar1))
    names(ds2)[ncol(ds2)]<-covars[i]
    ds2$obs<-as.factor(1:nrow(ds2))
    form1<-as.formula(paste0('outcome.pre~s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+' ,covars[i], '+(1|obs)')) 
    mod1<-glmer(form1 , family=poisson(link=log) , data=ds2,
                control = glmerControl(optimizer = "bobyqa",
                                       optCtrl =
                                         list(maxfun = 2e6)))
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
#' @param plot.labs Plot title
#' @return Univariate analysis plot, `results`, as described below
#' @export
#' @importFrom ggplot2 ggplot geom_segment geom_point coord_cartesian theme_minimal scale_y_discrete
#' @importFrom utils head tail

evaluatr.univariate.plot<-function(ds, plot.labs='Univariate'){
  ylevels = rev(ds$covar)
  ggplot(ds) +
    geom_point(aes(x=rr, y=covar)) + 
    geom_segment(aes(x=rr.lcl, xend=rr.ucl, y=covar, yend=covar)) +
    coord_cartesian(xlim=c(0.2, 2)) +
    geom_segment(aes(x=1, y=head(ylevels, 1), xend=1, yend=tail(ylevels, 1)), linetype="dashed", color="gray") +
    scale_y_discrete(limits = ylevels) +
    labs(title=plot.labs, x='Univariate Rate Ratio', y=NULL) +
    theme_minimal()
}