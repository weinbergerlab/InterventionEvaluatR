
##Process results

#Wrapper to call the model for each variant and loop through age groups
run_inla_mod<-function(model.variants.pass){
  mods<- lapply(unique(analysis$input_data[,analysis$group_name]),inla_mods, model.variant=model.variants.pass)
  return(mods)
}

 #run_inla_mod<-function(ages.pass){
#   mods<- lapply(model.variants,inla_mods, age.select=ages.pass)
#   return(mods)
# }

#Main function
inla_mods<-function(age.select, model.variant){
  
  x.all<-sa2[sa2$age==age.select,4:40]
  y<-  sa2[sa2$age==age.select,"Pneum"]
  dates<-unique(sa2$date)
  y.pre<-y
  y.pre[dates>=analysis$intervention_date]<-NA
  #Filter columsn with 0 variations in the covariate in the pre-vax period
  #x.var<-apply(x,2, function(xx) var(xx))  
  #x.var[is.na(x.var)]<-0
  #x<-x[,x.var>0] 
  
  ##NEEDTO SEPRATE OUT SEASONAL DUMMIES FROM COVARS
  x<-x.all[,-(1:(analysis$n_seasons-1))]
  x.months<-x[,(1:(analysis$n_seasons-1))]
  names(x.months)<-paste0('season', 1:(analysis$n_seasons-1))
  
  x.scale<-apply(x,2, function(z) scale(log(z+0.5))) 
  # x.scale$trend<-1:nrow(x.scale)
  
  y.aware.scale<- apply(x, 2, function(x1){
    log.y.pre.scale<- scale(log(y[dates<analysis$intervention_date]+0.5))
    reg<-lm(log.y.pre.scale~x1[dates<analysis$intervention_date])
    slope<- reg$coefficients[2]
    x.scale<-x1*slope - mean(x1*slope)
    return(x.scale)
  })
  pca1<- prcomp(y.aware.scale, center = FALSE,scale. = FALSE)
  n.pcs.keep<-sum(pca1$sdev>1)
  pcs<-pca1$x
  pcs<-apply(pcs,2, scale) #SCALE THE PCS prior to regression!
  pcs.combo<-cbind.data.frame( 'month'=as.factor(month(dates)),pcs[,1:n.pcs.keep, drop=F] )
  form0<-as.formula(paste0('~', paste(names(pcs.combo), collapse='+')))
  pcs.df<-as.data.frame(model.matrix(form0, data=pcs.combo)[,-1])
  
  df1<-cbind.data.frame(y.pre, x.months, y, x.scale )
  df1$t<-1:nrow(df1)
  covar.df.full<-cbind.data.frame( x.scale, x.months) 
  
 #Matrix for restriction random effects
  #x.in.full<-as.data.frame(model.matrix(~. , data=covar.df.full))[,-1] #dummies for month, remove intercept
  x.in.full<-covar.df.full
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
  x.in.full.no.season<-x.in.full[,-grep('season', names(x.in.full))]
  x.in.full.season<-x.in.full[,grep('season', names(x.in.full))]
  first.beta<-paste0("f(beta1,", names(x.in.full.no.season)[1], ",model='iid',hyper = param.beta,", "values=c(", paste(1:ncol(x.in.full.no.season), collapse=',') ,"))" )
  next.beta.f<-function(x, covar.index){ 
    paste0("f(beta",covar.index,",", names(x.in.full.no.season)[covar.index], ",copy='beta1', fixed=T)")
  }
  betas.list<- mapply(next.beta.f, x=x.in.full.no.season[-1], covar.index=2:ncol(x.in.full.no.season), SIMPLIFY=T)
  next.betas<- paste(betas.list, collapse='+')  
  all.betas<-paste(first.beta,next.betas , sep='+')    
  form1<- as.formula(paste0("y.pre ~", paste(names(x.in.full.season), collapse="+"),"+", all.betas,    "+ f(t, model = 'ar1', constr=T,extraconstr=list(A=A.full, e=e.full))") )
  
  mod1.full = inla(form1, data = mod.df.full, 
                   control.predictor = list(compute=TRUE, link = 1), 
                   family='poisson',
                   control.compute=list(config = TRUE))
  ds<-list('fitted.model'=mod1.full, 'mod.df'=mod.df.full, 'x.in'=x.in.full,'offset'=NA,'offset.used'=F,'ridge'=T)
  }
  
  if(model.variant %in% c('time','time_no_offset')){
    if(analysis$denom_name  %in% names(x)){
      log.offset<- x[,c(analysis$denom_name)]
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
   
    form2<- as.formula(paste0('y.pre ~',  paste(names(x.in.time), collapse='+'), "+  f(t, model = 'ar1', constr=T,extraconstr=list(A=A.time, e=e.time))") )
    
    if(model.variant=='time_no_offset'){
    mod2.time.no.offset = inla(form2, data = mod.df.time, control.predictor =   list(compute=TRUE, link = 1), family='poisson',control.compute=list(config = TRUE))
    ds<-list('fitted.model'=mod2.time.no.offset, 'mod.df'=mod.df.time, 'x.in'=x.in.time,'offset'=NA, 'offset.used'=F, 'ridge'=F)
    }else if(model.variant=='time'){
    #Time regression with offset
    mod3.time = inla(form2, data = mod.df.time, control.predictor =   list(compute=TRUE, link = 1), family='poisson',E=offset1, control.compute=list(config = TRUE))
    ds<-list('fitted.model'=mod3.time, 'mod.df'=mod.df.time, 'x.in'=x.in.time, 'offset'=offset1,'offset.used'=T,'ridge'=F)
    }
  }
  
  if(ds$ridge==T){  
    coefs<-ds$fitted.model$summary.random[['beta1']]
    coefs$covar<-names(ds$x.in)[-grep('season',names(ds$x.in)) ]
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
    covar.select.months<-grep('season', post.labels)
    covar.select<-c(covar.select.betas,covar.select.months)
  }else{
    covar.select<- which(names(ds$x.in)  %in% sub("\\:.*", "", post.labels))
  }
  
  beta.posterior<-(posterior.samples[covar.select,]) 
  beta.posterior.median<- apply(beta.posterior,1,median)
  beta.posterior.hdi<-cbind(beta.posterior.median,t(hdi(t(beta.posterior), credMass = 0.95)))
  row.names(beta.posterior.hdi)<-names(ds$x.in)
  
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
  log.rr.pointwise.hdi<-cbind(log.rr.pointwise.median,log.rr.pointwise.ci)
  # matplot(log.rr, type='l', col='gray', lty=c(2,1,2))
  # abline(h=0, col='red')
  
  post.period<- which(dates>= analysis$eval_period[1]  &dates<= analysis$eval_period[2] )
  
  post.samples<-posterior.preds.counts[post.period,]
  
  post.samples.sum<-apply(post.samples,2,sum)
  obs.post.sum<- sum(mod.df.full$y[post.period])
  rr.agg<-obs.post.sum/post.samples.sum
  rr.q<-quantile(rr.agg, probs=c(0.025, 0.5, 0.975)) 
  rr.hdi<-c(rr.q['50%'],hdi(rr.agg, credMass = 0.95))
  
  
  out.list<- list('posterior.pred.hdi'=posterior.pred.hdi, 'rr.hdi'=rr.hdi,'inla.mod'=ds$fitted.model,'log.rr.pointwise.hdi'=log.rr.pointwise.hdi,'obs.y'=mod.df.full$y,'rand.eff1.q'=rand.eff1.q,'betas'=beta.posterior.hdi, 'post.period.start.index'=post.period[1], 'rand.eff.t1'=rand.eff.t1,'rand.eff.combined.q'=rand.eff.combined.q, 'fixed.effect'=fixed.effect.hdi,'beta.posterior.hdi'=beta.posterior.hdi, 'rho1'=rho1)
  return(out.list)
}


