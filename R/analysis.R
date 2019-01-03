syncon.impact = function(data_full, data_time, data_time_no_offset, data_pca, intervention_date, time_points, n_seasons, crossval) {
  #Start Cluster for CausalImpact (the main analysis function).
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE)})
  clusterExport(cl, c('doCausalImpact',  'intervention_date', 'time_points', 'n_seasons','crossval'), environment())
  impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, var.select.on=TRUE, time_points = time_points), groups)
  impact_time <- setNames(parLapply(cl, data_time, doCausalImpact, intervention_date = intervention_date,  var.select.on=FALSE,time_points = time_points, trend = TRUE), groups)
  impact_time_no_offset <- setNames(parLapply(cl, data_time_no_offset, doCausalImpact, intervention_date = intervention_date,  var.select.on=FALSE,time_points = time_points,  trend = FALSE), groups)
  impact_pca <- setNames(parLapply(cl, data_pca, doCausalImpact, intervention_date = intervention_date, var.select.on=FALSE, time_points = time_points), groups)
  #No covariates, but with random intercept
  #impact_seas_only <- setNames(parLapply(cl, data_null, doCausalImpact, intervention_date = intervention_date, var.select.on=FALSE, time_points = time_points, n_seasons = n_seasons), groups)
  #No covariates except seasonal, no random intercept
  #impact_seas_only_no_re <- setNames(parLapply(cl, data_null, doCausalImpact, intervention_date = intervention_date, var.select.on=FALSE,ri.select=FALSE, time_points = time_points, n_seasons = n_seasons), groups)
  stopCluster(cl)
  
  list(impact_full=impact_full, impact_time=impact_time, impact_time_no_offset=impact_time_no_offset, impact_pca=impact_pca)
}

syncon.crossval = function(data_full, data_time, data_time_no_offset, data_pca, intervention_date, time_points, n_seasons, crossval) {
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
      cv_impact_full <-setNames(parLapply(cl, cv.data_full, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, intervention_date = intervention_date,  var.select.on=TRUE, time_points = time_points)), groups)
      cv_impact_time_no_offset <-setNames(parLapply(cl, cv.data_time_no_offset, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, trend=FALSE, intervention_date = intervention_date,  var.select.on=FALSE, time_points = time_points)), groups)
      cv_impact_time <-setNames(parLapply(cl, cv.data_time, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, trend=TRUE, intervention_date = intervention_date,  var.select.on=FALSE, time_points = time_points)), groups)
      cv_impact_pca <-setNames(parLapply(cl, cv.data_pca, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, intervention_date = intervention_date,  var.select.on=FALSE, time_points = time_points)), groups)
    stopCluster(cl)
  
    #Calculate pointwise log likelihood for cross-val prediction sample vs observed
    #These are N_iter*N_obs*N_cross_val array
    ll.cv.full<-lapply(cv_impact_full, function(x) lapply(x,crossval.log.lik))
    ll.cv.full2<-lapply(ll.cv.full, reshape.arr)
    #
    ll.cv.time_no_offset<-lapply(cv_impact_time_no_offset, function(x) lapply(x,crossval.log.lik))
    ll.cv.time_no_offset2<-lapply(ll.cv.time_no_offset, reshape.arr)
    #
    ll.cv.time<-lapply(cv_impact_time, function(x) lapply(x,crossval.log.lik))
    ll.cv.time2<-lapply(ll.cv.time, reshape.arr)
    #
    ll.cv.pca<-lapply(cv_impact_pca, function(x) lapply(x,crossval.log.lik))
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
    stacking_weights.all.m<-melt(stacking_weights.all, id.vars='groups')
   # stacking_weights.all.m<-stacking_weights.all.m[order(stacking_weights.all.m$groups),]
    
    stacked.ests<-mapply(  FUN=stack.mean,group=groups,impact_full=impact_full,impact_time=impact_time,impact_time_no_offset=impact_time_no_offset,impact_pca=impact_pca, SIMPLIFY=FALSE )
   # plot.stacked.ests<-lapply(stacked.ests,plot.stack.est)
    quantiles_stack <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = stacked.ests[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period)}), groups)
    pred_quantiles_stack <- sapply(quantiles_stack, getPred, simplify = 'array')
    rr_roll_stack <- sapply(quantiles_stack, FUN = function(quantiles_stack) {quantiles_stack$roll_rr}, simplify = 'array')
    rr_mean_stack <- round(t(sapply(quantiles_stack, getRR)),2)
    rr_mean_stack_intervals <- data.frame('Stacking Estimate (95% CI)'     = makeInterval(rr_mean_stack[, 2], rr_mean_stack[, 3], rr_mean_stack[, 1]), check.names = FALSE, row.names = groups)
    cumsum_prevented_stack <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_stack, simplify = 'array')
    ann_pred_quantiles_stack <- sapply(quantiles_stack, getAnnPred, simplify = FALSE)
    
    #Preds: Compare observed and expected
    pred.cv.full<-lapply(cv_impact_full, function(x) sapply(x,pred.cv,simplify='array'))
    pred.cv.pca<-lapply(cv_impact_pca, function(x) sapply(x,pred.cv,simplify='array'))
      
    # # par(mfrow=c(3,2))
    # plot.grp=9
    # for(i in 1:6){
    # matplot(pred.cv.full[[plot.grp]][,c(2:4),i], type='l', ylab='Count',col='#1b9e77', lty=c(2,1,2), bty='l', ylim=range(pred.cv.full[[plot.grp]][,c(1),i])*c(0.8,1.2))
    # points(pred.cv.full[[plot.grp]][,c(1),i], pch=16)
    # title("Synthetic controls: Cross validation")
    # matplot(pred.cv.pca[[plot.grp]][,c(2:4),i], type='l',ylab='Count', col='#d95f02', lty=c(2,1,2), bty='l', ylim=range(pred.cv.full[[plot.grp]][,c(1),i])*c(0.8,1.2))
    # points(pred.cv.pca[[plot.grp]][,c(1),i], pch=16)
    # title("STL+PCA: Cross validation")
    # }
    return list(ann_pred_quantiles_stack=ann_pred_quantiles_stack, pred_quantiles_stack=pred_quantiles_stack, rr_roll_stack=rr_roll_stack, rr_mean_stack=rr_mean_stack, rr_mean_stack_intervals=rr_mean_stack_intervals, cumsum_prevented_stack=cumsum_prevented_stack)
}