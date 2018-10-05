#This is the analysis file. The functions used in this file are cointained in synthetic_control_functions.R
#There are two model variants: 
# *_full - Full synthetic control model with all covariates (excluding user-specified covariates).
# *_time - Trend adjustment using the specified variable (e.g., non-respiratory hospitalization or population size) as the denominator.

#############################
#                           #
#    System Preparations    #
#                           #
#############################

source('synthetic_control_functions.R', local = TRUE)

#############################
#Automatically set working directory to desktop
#setwd('~/synthetic-control-master/main analysis components')  #directory where .Rmd file is saved
#Set working directory: default to desktop--different path for windows vs Mac
if(.Platform$OS.type == "windows") {
  desktop<-file.path(Sys.getenv("USERPROFILE"),"Desktop")
  desktop<-gsub(pattern='\\',replacement='/', desktop, fixed=TRUE)
} else {
  desktop<- "~/Desktop"
}
auto.wd<-file.path(paste0(desktop,'/synthetic-control-poisson-master/main analysis components/'))
#

packages <- c('parallel', 'splines', 'lubridate','loo', 'RcppRoll','pomp','lme4', 'BoomSpikeSlab', 'ggplot2', 'reshape','dummies')
packageHandler(packages, update_packages, install_packages)
sapply(packages, library, quietly = TRUE, character.only = TRUE)

#Detect if pogit package installed; if not download archive (no longer on cran)
if("BayesLogit" %in% rownames(installed.packages())==FALSE){
  if(.Platform$OS.type == "windows") {
  #url_BayesLogit<- "https://mran.microsoft.com/snapshot/2017-02-04/src/contrib/BayesLogit_0.6.tar.gz"
  install_github("jwindle/BayesLogit")
  }else{
    url_BayesLogit<- "https://github.com/weinbergerlab/synthetic-control-poisson/blob/master/packages/BayesLogit_0.6_mac.tgz?raw=true"
  }
  pkgFile_BayesLogit <- "BayesLogit.tar.gz"
  download.file(url = url_BayesLogit, destfile = pkgFile_BayesLogit)
  install.packages(url_BayesLogit, type="source", repos=NULL)
}
if("pogit" %in% rownames(installed.packages())==FALSE){
  url_pogit <- "https://cran.r-project.org/src/contrib/Archive/pogit/pogit_1.1.0.tar.gz"
  pkgFile_pogit <- "pogit_1.1.0.tar.gz"
  download.file(url = url_pogit, destfile = pkgFile_pogit)
  install.packages(pkgs=pkgFile_pogit, type="source", repos=NULL)
  install.packages('logistf')
}
library(pogit)

#Detects number of available cores on computers. Used for parallel processing to speed up analysis.
n_cores <- detectCores()
set.seed(1)

###################################################
#                                                 #
# Directory setup and initialization of constants #
#                                                 #
###################################################

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
groups <- as.character(unique(unlist(prelog_data[, group_name], use.names = FALSE)))
if (exists('exclude_group')) {groups <- groups[!(groups %in% exclude_group)]}

###############################################
#                                             #
# Data and covariate preparation for analysis #
#                                             #
###############################################

#Make sure we are in right format
prelog_data[,date_name]<-as.Date(as.character(prelog_data[,date_name]), tryFormats=c("%m/%d/%Y",'%Y-%m-%d' ))

#test<-split(prelog_data, factor(prelog_data[,group_name]))
#outcome.na<-sapply(test, function(x) sum(is.na(x[,outcome_name])))
prelog_data[, date_name] <- formatDate(prelog_data[, date_name])
prelog_data <- setNames(lapply(groups, FUN = splitGroup, ungrouped_data = prelog_data, group_name = group_name, date_name = date_name, start_date = start_date, end_date = end_date, no_filter = c(group_name, date_name, outcome_name, denom_name)), groups)
#if (exists('exclude_group')) {prelog_data <- prelog_data[!(names(prelog_data) %in% exclude_group)]}

#Log-transform all variables, adding 0.5 to counts of 0.
ds <- setNames(lapply(prelog_data, FUN = logTransform, no_log = c(group_name, date_name,outcome_name)), groups)
time_points <- unique(ds[[1]][, date_name])

#Monthly dummies
if(n_seasons==4){dt<-quarter(as.Date(time_points))}
if(n_seasons==12){dt<-month(as.Date(time_points))}
if(n_seasons==3){
  dt.m<-month(as.Date(time_points))
  dt<-dt.m
  dt[dt.m %in% c(1,2,3,4)]<-1
  dt[dt.m %in% c(5,6,7,8)]<-2
  dt[dt.m %in% c(9,10,11,12)]<-3
    }
season.dummies<-dummy(dt)
season.dummies<-as.data.frame(season.dummies)
names(season.dummies)<-paste0('s', 1:n_seasons)
season.dummies<-season.dummies[,-n_seasons]

ds <- lapply(ds, function(ds) {
	if (!(denom_name %in% colnames(ds))) {
		ds[denom_name] <- 0
	}
	return(ds)
})

sparse_groups <- sapply(ds, function(ds) {
	return(ncol(ds[!(colnames(ds) %in% c(date_name, group_name, denom_name, outcome_name, exclude_covar))]) == 0)
})
ds <- ds[!sparse_groups]
groups <- groups[!sparse_groups]

#Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
covars_full <- setNames(lapply(ds, makeCovars), groups)
covars_full <- lapply(covars_full, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_covar), drop = FALSE]})
covars_time <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies,time_index = 1:nrow(covars))))}), groups)
covars_null <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies)))}), groups)

#Standardize the outcome variable and save the original mean and SD for later analysis.
outcome      <- sapply(ds, FUN = function(data) {data[, outcome_name]})
outcome_plot=outcome
offset<- sapply(ds, FUN=function(data) exp(data[, denom_name]) )  #offset term on original scale; 1 column per age group
################################
#set up for STL+PCA
################################
##SECTION 1: CREATING SMOOTHED VERSIONS OF CONTROL TIME SERIES AND APPENDING THEM ONTO ORIGINAL DATAFRAME OF CONTROLS
#EXTRACT LONG TERM TREND WITH DIFFERENT LEVELS OF SMOOTHNESS USING STL
# Set a list of parameters for STL
stl.covars<-mapply(smooth_func,ds.list=ds,covar.list=covars_full, SIMPLIFY=FALSE) 
post.start.index<-which(time_points==post_period[1])

if (length(groups)>1){ 
  stl.data.setup<-mapply(stl_data_fun,covars=stl.covars, ds.sub=ds ,SIMPLIFY=FALSE)  #list of lists that has covariates for each regression for each strata
}else{
  stl.data.setup <- list(mapply(stl_data_fun,covars=stl.covars, ds.sub=ds ))
}

##SECTION 2: run first stage models
n_cores <- detectCores()-1
glm.results<- vector("list",  length=length(stl.data.setup)) #combine models into a list
cl1 <- makeCluster(n_cores)
clusterEvalQ(cl1, {library(lme4, quietly = TRUE)})
clusterExport(cl1, c('stl.data.setup',  'glm.fun', 'time_points', 'n_seasons','post.start.index'), environment())
for(i in 1:length(stl.data.setup)){
  print(i)
  glm.results[[i]]<-parLapply(cl=cl1 ,     stl.data.setup[[i]], fun=glm.fun )
}
stopCluster(cl1)
######################

#Combine the outcome, covariates, and time point information.
data_full <- setNames(lapply(groups, makeTimeSeries, outcome = outcome,       covars = covars_full), groups)
data_time <- setNames(lapply(groups, makeTimeSeries, outcome = outcome, covars = covars_time, trend=TRUE), groups)
data_pca<-mapply(FUN=pca_top_var,glm.results.in=glm.results, covars=stl.covars,ds.in=ds, SIMPLIFY=FALSE)
names(data_pca)<-groups
#Null model where we only include seasonal terms but no covariates
data_null <- setNames(lapply(groups, makeTimeSeries, outcome = outcome, covars = covars_null, trend=FALSE), groups)
#Time trend model but without a denominator
data_time_no_offset <- setNames(lapply(groups, makeTimeSeries, outcome = outcome, covars = covars_time, trend=FALSE), groups)

###############################
#                             #
#        Main analysis        #
#                             #
###############################

#Start Cluster for CausalImpact (the main analysis function).
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE)})
clusterExport(cl, c('doCausalImpact',  'intervention_date', 'time_points', 'n_seasons','crossval'), environment())
  impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, var.select.on=TRUE, time_points = time_points), groups)
  impact_time <- setNames(parLapply(cl, data_time, doCausalImpact, intervention_date = intervention_date,  var.select.on=FALSE,time_points = time_points, trend = TRUE), groups)
  impact_time_no_offset <- setNames(parLapply(cl, data_time_no_offset, doCausalImpact, intervention_date = intervention_date,  var.select.on=FALSE,time_points = time_points,  trend = FALSE), groups)
  impact_pca <- setNames(parLapply(cl, data_pca, doCausalImpact, intervention_date = intervention_date, var.select.on=FALSE, time_points = time_points), groups)
stopCluster(cl)

##Model size for SC model
model.size.sc<-sapply(impact_full,modelsize_func)

####################################################
####################################################
#CROSS VALIDATION
####################################################
    if(crossval){
      #Creates List of lists: 1 entry for each stratum; within this, there are CV datasets for each year left out, and within this, there are 2 lists, one with full dataset, and one with the CV dataset
        cv.data_full<-lapply(data_full, makeCV)
        cv.data_time<-lapply(data_time, makeCV)
        cv.data_time_no_offset<-lapply(data_time_no_offset, makeCV)
        cv.data_pca<-lapply(data_pca, makeCV)
      #zoo_data<-cv.data_time[[1]][[2]]
      #Run the models on each of these datasets
        # Start the clock!--takes ~45 minutes
          ptm <- proc.time()
        cl <- makeCluster(n_cores)
        clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE)})
        clusterExport(cl, c('doCausalImpact',  'intervention_date', 'time_points', 'n_seasons','crossval'), environment())
          cv_impact_full <-setNames(parLapply(cl, cv.data_full, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, intervention_date = intervention_date,  var.select.on=TRUE, time_points = time_points)), groups)
          cv_impact_time_no_offset <-setNames(parLapply(cl, cv.data_time_no_offset, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, trend=FALSE, intervention_date = intervention_date,  var.select.on=FALSE, time_points = time_points)), groups)
          cv_impact_time <-setNames(parLapply(cl, cv.data_time, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, trend=TRUE, intervention_date = intervention_date,  var.select.on=FALSE, time_points = time_points)), groups)
          cv_impact_pca <-setNames(parLapply(cl, cv.data_pca, function(x) lapply(x, doCausalImpact,crossval.stage=TRUE, intervention_date = intervention_date,  var.select.on=FALSE, time_points = time_points)), groups)
        stopCluster(cl)
        # Stop the clock
        proc.time() - ptm
      
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
          
         save.stack.est<-list(post_period,outcome_plot, time_points,ann_pred_quantiles_stack, pred_quantiles_stack,rr_roll_stack,rr_mean_stack,rr_mean_stack_intervals,cumsum_prevented_stack)
        names(save.stack.est)<-c('post_period','outcome_plot','time_points', 'ann_pred_quantiles_stack', 'pred_quantiles_stack','rr_roll_stack','rr_mean_stack','rr_mean_stack_intervals','cumsum_prevented_stack')
        saveRDS(save.stack.est, file=paste0(output_directory, country, "Stack estimates.rds"))
        
        #Pointwise RR and uncertainty for second stage meta analysis
        log_rr_quantiles_stack   <- sapply(quantiles_stack,   FUN = function(quantiles) {quantiles$log_rr_full_t_quantiles}, simplify = 'array')
        dimnames(log_rr_quantiles_stack)[[1]] <- time_points
        log_rr_full_t_samples.stack.prec<-sapply(quantiles_stack,   FUN = function(quantiles) {quantiles$log_rr_full_t_samples.prec.post}, simplify = 'array')
        #log_rr_sd.stack   <- sapply(quantiles_stack,   FUN = function(quantiles) {quantiles$log_rr_full_t_sd}, simplify = 'array')
        
        saveRDS(log_rr_quantiles_stack, file=paste0(output_directory, country, "_log_rr_quantiles_stack.rds"))
        saveRDS(log_rr_full_t_samples.stack.prec, file=paste0(output_directory, country, "_log_rr_full_t_samples.stack.prec.rds"))
      }
##########################################################################
##########################################################################

#Save the inclusion probabilities from each of the models
inclusion_prob_full <- setNames(lapply(impact_full, inclusionProb), groups)
inclusion_prob_time <- setNames(lapply(impact_time, inclusionProb), groups)

#All model results combined
quantiles_full <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_full[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period)}), groups)
quantiles_time <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_time[[group]], denom_data = ds[[group]][, denom_name],  eval_period = eval_period, post_period = post_period)}), groups)
quantiles_time_no_offset <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_time_no_offset[[group]], denom_data = ds[[group]][, denom_name],  eval_period = eval_period, post_period = post_period)}), groups)
quantiles_pca <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_pca[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period)}), groups)
#USE SC estimates if model size is at leatst 1, otherwise use STL+PCA
quantiles_best<-vector("list", length(quantiles_full)) 
quantiles_best[model.size.sc>=1]<-quantiles_full[model.size.sc>=1]
quantiles_best[model.size.sc<1]<-quantiles_pca[model.size.sc<1]
quantiles_best<-setNames(quantiles_best,groups)

#Model predicitons
pred_quantiles_full <- sapply(quantiles_full, getPred, simplify = 'array')
pred_quantiles_time <- sapply(quantiles_time, getPred, simplify = 'array')
pred_quantiles_time_no_offset <- sapply(quantiles_time_no_offset, getPred, simplify = 'array')
pred_quantiles_pca <- sapply(quantiles_pca, getPred, simplify = 'array')
pred_quantiles_best <- sapply(quantiles_best, getPred, simplify = 'array')

#Predictions, aggregated by year
ann_pred_quantiles_full <- sapply(quantiles_full, getAnnPred, simplify = FALSE)
ann_pred_quantiles_time <- sapply(quantiles_time, getAnnPred, simplify = FALSE)
ann_pred_quantiles_time_no_offset <- sapply(quantiles_time_no_offset, getAnnPred, simplify = FALSE)
ann_pred_quantiles_pca <- sapply(quantiles_pca, getAnnPred, simplify = FALSE)
ann_pred_quantiles_best<-sapply(quantiles_best, getAnnPred, simplify = FALSE)

#Pointwise RR and uncertainty for second stage meta analysis
log_rr_quantiles   <- sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_quantiles}, simplify = 'array')
dimnames(log_rr_quantiles)[[1]] <- time_points
log_rr_sd   <- sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_sd}, simplify = 'array')
log_rr_full_t_samples.prec<-sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_samples.prec}, simplify = 'array')
saveRDS(log_rr_quantiles, file=paste0(output_directory, country, "_log_rr_quantiles.rds"))
saveRDS(log_rr_sd, file=paste0(output_directory, country, "_log_rr_sd.rds"))
saveRDS(log_rr_full_t_samples.prec, file=paste0(output_directory, country, "_log_rr_full_t_samples.prec.rds"))

log_rr_quantiles_best   <- sapply(quantiles_best,   FUN = function(quantiles) {quantiles$log_rr_full_t_quantiles}, simplify = 'array')
dimnames(log_rr_quantiles_best)[[1]] <- time_points
log_rr_best_t_samples.prec<-sapply(quantiles_best,   FUN = function(quantiles) {quantiles$log_rr_full_t_samples.prec}, simplify = 'array')
saveRDS(log_rr_quantiles_best, file=paste0(output_directory, country, "_log_rr_quantiles_best.rds"))
saveRDS(log_rr_best_t_samples.prec, file=paste0(output_directory, country, "_log_rr_best_t_samples.prec.rds"))

#Rolling rate ratios
rr_roll_full <- sapply(quantiles_full, FUN = function(quantiles_full) {quantiles_full$roll_rr}, simplify = 'array')
rr_roll_time <- sapply(quantiles_time, FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array')
rr_roll_time_no_offset <- sapply(quantiles_time_no_offset, FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array')
rr_roll_pca <- sapply(quantiles_pca, FUN = function(quantiles_pca) {quantiles_pca$roll_rr}, simplify = 'array')
rr_roll_best <- sapply(quantiles_best, FUN = function(quantiles_best) {quantiles_best$roll_rr}, simplify = 'array')

#Rate ratios for evaluation period.
rr_mean_full <- t(sapply(quantiles_full, getRR))
rr_mean_time <- t(sapply(quantiles_time, getRR))
rr_mean_time_no_offset <- t(sapply(quantiles_time_no_offset, getRR))
rr_mean_pca <- t(sapply(quantiles_pca, getRR))
rr_mean_best <- t(sapply(quantiles_best, getRR))

rr_mean_full_intervals <- data.frame('SC Estimate (95% CI)'     = makeInterval(rr_mean_full[, 2], rr_mean_full[, 3], rr_mean_full[, 1]), check.names = FALSE, row.names = groups)
rr_mean_time_intervals <- data.frame('Time trend Estimate (95% CI)' = makeInterval(rr_mean_time[, 2], rr_mean_time[, 3], rr_mean_time[, 1]), check.names = FALSE, row.names = groups)
rr_mean_time_no_offset_intervals <- data.frame('Time trend (no offset) Estimate (95% CI)' = makeInterval(rr_mean_time_no_offset[, 2], rr_mean_time_no_offset[, 3], rr_mean_time_no_offset[, 1]), check.names = FALSE, row.names = groups)
rr_mean_pca_intervals <- data.frame('STL+PCA Estimate (95% CI)'     = makeInterval(rr_mean_pca[, 2], rr_mean_pca[, 3], rr_mean_pca[, 1]), check.names = FALSE, row.names = groups)
rr_mean_best_intervals <- data.frame('Best Estimate (95% CI)'     = makeInterval(rr_mean_best[, 2], rr_mean_best[, 3], rr_mean_best[, 1]), check.names = FALSE, row.names = groups)
colnames(rr_mean_time) <- paste('Time_trend', colnames(rr_mean_time))

#Combine RRs into 1 file for plotting
rr_mean_combo<- as.data.frame(rbind( cbind(rep(1, nrow(rr_mean_full)),groups,  seq(from=1, by=1, length.out=nrow(rr_mean_full)),rr_mean_full),
                       cbind(rep(2, nrow(rr_mean_time)),groups, seq(from=1, by=1, length.out=nrow(rr_mean_full)), rr_mean_time),
                       cbind(rep(3, nrow(rr_mean_time_no_offset)),groups, seq(from=1, by=1, length.out=nrow(rr_mean_full)), rr_mean_time_no_offset),
                    cbind(rep(4, nrow(rr_mean_pca)), groups, seq(from=1, by=1, length.out=nrow(rr_mean_full)),rr_mean_pca)))
        names(rr_mean_combo)<-c('Model', 'groups', 'group.index','lcl','mean.rr','ucl')
        if(crossval){
          point.weights2<-stacking_weights.all.m
        }else{
          point.weights2<-as.data.frame(matrix(rep(1,nrow(rr_mean_combo)), ncol=1))
          #point.weights2[]<-1.5
          names(point.weights2)<-'value'
        }
        rr_mean_combo$point.weights<-point.weights2$value
        rr_mean_combo$group.index<-as.numeric(as.character(rr_mean_combo$group.index))
        rr_mean_combo$mean.rr<-as.numeric(as.character(rr_mean_combo$mean.rr))
        rr_mean_combo$lcl<-as.numeric(as.character(rr_mean_combo$lcl))
        rr_mean_combo$ucl<-as.numeric(as.character(rr_mean_combo$ucl))
        rr_mean_combo$group.index[rr_mean_combo$Model==2]<-rr_mean_combo$group.index[rr_mean_combo$Model==2]+0.15
        rr_mean_combo$group.index[rr_mean_combo$Model==3]<-rr_mean_combo$group.index[rr_mean_combo$Model==3]+0.3
        rr_mean_combo$group.index[rr_mean_combo$Model==4]<-rr_mean_combo$group.index[rr_mean_combo$Model==4]+0.45
        rr_mean_combo$Model<-as.character(rr_mean_combo$Model)
        rr_mean_combo$Model[rr_mean_combo$Model=='1']<-"Synthetic Controls"
        rr_mean_combo$Model[rr_mean_combo$Model=='2']<-"Time trend"
        rr_mean_combo$Model[rr_mean_combo$Model=='3']<-"Time trend (No offset)"
        rr_mean_combo$Model[rr_mean_combo$Model=='4']<-"STL+PCA"
        cbPalette <- c("#1b9e77", "#d95f02", "#7570b3",'#e7298a')
        rr_mean_combo$est.index<-as.factor(1:nrow(rr_mean_combo))
        #Fix order for axis
        rr_mean_combo$Model<-as.factor(rr_mean_combo$Model)
        rr_mean_combo$Model = factor(rr_mean_combo$Model,levels(rr_mean_combo$Model)[c(2,3,4,1)])
        #print(levels(rr_mean_combo$Model))

cumsum_prevented <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_full, simplify = 'array')
cumsum_prevented_pca <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_pca, simplify = 'array')
cumsum_prevented_time <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_time, simplify = 'array')
cumsum_prevented_best <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_best, simplify = 'array')




################################
#                              #
#     Sensitivity Analyses     #
#                              #
################################

#Pred Sensitivity Analysis--tests effect of changing prior on Ncovars from 3 to 2 to 10
# cl <- makeCluster(n_cores)
# clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
# clusterExport(cl, c('doCausalImpact', 'predSensitivityAnalysis', 'inclusionProb', 'rrPredQuantiles', 'getPred', 'getRR', 'groups', 'ds', 'data_full', 'denom_name', 'outcome_mean', 'outcome_sd', 'intervention_date', 'eval_period', 'post_period', 'time_points', 'n_seasons'), environment())
# 
# sensitivity_analysis_pred_2  <- setNames(as.data.frame(t(parSapply(cl, groups, predSensitivityAnalysis, ds = ds, zoo_data = data_full, denom_name = denom_name, outcome_mean = outcome_mean, outcome_sd = outcome_sd, intervention_date = intervention_date, eval_period = eval_period, post_period = post_period, time_points = time_points, n_seasons = n_seasons, n_pred = 2 ))), c('Lower CI', 'Point Estimate', 'Upper CI'))
# sensitivity_analysis_pred_10 <- setNames(as.data.frame(t(parSapply(cl, groups, predSensitivityAnalysis, ds = ds, zoo_data = data_full, denom_name = denom_name, outcome_mean = outcome_mean, outcome_sd = outcome_sd, intervention_date = intervention_date, eval_period = eval_period, post_period = post_period, time_points = time_points, n_seasons = n_seasons, n_pred = 10))), c('Lower CI', 'Point Estimate', 'Upper CI'))
# 
# stopCluster(cl)
# 
# sensitivity_analysis_pred_2_intervals  <- data.frame('Estimate (95% CI)' = makeInterval(sensitivity_analysis_pred_2[, 2],  sensitivity_analysis_pred_2[, 3],  sensitivity_analysis_pred_2[, 1]),  row.names = groups, check.names = FALSE)
# sensitivity_analysis_pred_10_intervals <- data.frame('Estimate (95% CI)' = makeInterval(sensitivity_analysis_pred_10[, 2], sensitivity_analysis_pred_10[, 3], sensitivity_analysis_pred_10[, 1]), row.names = groups, check.names = FALSE)
# 


bad_sensitivity_groups <- sapply(covars_full, function (covar) {ncol(covar) <= n_seasons-1+3})
 sensitivity_covars_full <- covars_full[!bad_sensitivity_groups]
 sensitivity_ds <- ds[!bad_sensitivity_groups]
 sensitivity_impact_full <- impact_full[!bad_sensitivity_groups]
 sensitivity_groups <- groups[!bad_sensitivity_groups]

if (length(sensitivity_groups)!=0) {
 #Weight Sensitivity Analysis - top weighted variables are excluded and analysis is re-run.
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
clusterExport(cl, c('sensitivity_ds', 'doCausalImpact', 'year_def', 'weightSensitivityAnalysis', 'rrPredQuantiles', 'sensitivity_groups', 'intervention_date', 'outcome', 'time_points', 'n_seasons',  'eval_period', 'post_period','crossval'), environment())
  sensitivity_analysis_full <- setNames(parLapply(cl, sensitivity_groups, weightSensitivityAnalysis, covars = sensitivity_covars_full, ds = sensitivity_ds, impact = sensitivity_impact_full, time_points = time_points, intervention_date = intervention_date, n_seasons = n_seasons, outcome = outcome,  eval_period = eval_period, post_period = post_period), sensitivity_groups)
stopCluster(cl)

sensitivity_pred_quantiles  <- lapply(sensitivity_analysis_full, FUN = function(sensitivity_analysis) {
	pred_list <- vector(mode = 'list', length = length(sensitivity_analysis))
	for (sensitivity_index in 1:length(sensitivity_analysis)) {
		pred_list[[sensitivity_index]] <- getPred(sensitivity_analysis[[sensitivity_index]])
	}
	return(pred_list)
})

#Table of rate ratios for each sensitivity analysis level
sensitivity_table <- t(sapply(sensitivity_groups, sensitivityTable, sensitivity_analysis = sensitivity_analysis_full, original_rr = rr_mean_full))
sensitivity_table_intervals <- data.frame('Estimate (95% CI)' = makeInterval(sensitivity_table[, 2],  sensitivity_table[, 3],  sensitivity_table[, 1]),
																					'Top Control 1' = sensitivity_table[, 'Top Control 1'],
																					'Inclusion Probability of Control 1' = sensitivity_table[, 'Inclusion Probability of Control 1'],
																					'Control 1 Estimate (95% CI)' = makeInterval(sensitivity_table[, 7],  sensitivity_table[, 8],  sensitivity_table[, 6]),
																					'Top Control 2' = sensitivity_table[, 'Top Control 2'],
																					'Inclusion Probability of Control 2' = sensitivity_table[, 'Inclusion Probability of Control 2'],
																					'Control 2 Estimate (95% CI)' = makeInterval(sensitivity_table[, 12],  sensitivity_table[, 13],  sensitivity_table[, 11]),
																					'Top Control 3' = sensitivity_table[, 'Top Control 3'],
																					'Inclusion Probability of Control 3' = sensitivity_table[, 'Inclusion Probability of Control 3'],
																					'Control 3 Estimate (95% CI)' = makeInterval(sensitivity_table[, 17],  sensitivity_table[, 18],  sensitivity_table[, 16]), check.names = FALSE)
rr_table <- cbind.data.frame(round(rr_mean_time[!bad_sensitivity_groups, ],2), sensitivity_table)
rr_table_intervals <- cbind('ITS Estimate (95% CI)' = rr_mean_time_intervals[!bad_sensitivity_groups, ], sensitivity_table_intervals)
} else {
  sensitivity_table_intervals <- NA
}
