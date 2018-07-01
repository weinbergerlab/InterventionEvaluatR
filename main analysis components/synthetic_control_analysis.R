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

packages <- c('parallel', 'splines', 'lubridate','logistf','loo', 'RcppRoll','pomp','lme4', 'BoomSpikeSlab', 'ggplot2', 'reshape','dummies')
packageHandler(packages, update_packages, install_packages)
sapply(packages, library, quietly = TRUE, character.only = TRUE)

#Detect if pogit package installed; if not download archive (no longer on cran)
if("BayesLogit" %in% rownames(installed.packages())==FALSE){
  if(.Platform$OS.type == "windows") {
  url_BayesLogit<- "https://mran.microsoft.com/snapshot/2017-02-04/src/contrib/BayesLogit_0.6.tar.gz"
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
data_file <- paste(input_directory, file_name, sep = '')
prelog_data <- read.csv(data_file, check.names = FALSE)
groups <- as.character(unique(unlist(prelog_data[, group_name], use.names = FALSE)))
if (exists('exclude_group')) {groups <- groups[!(groups %in% exclude_group)]}

###############################################
#                                             #
# Data and covariate preparation for analysis #
#                                             #
###############################################


prelog_data[, date_name] <- formatDate(prelog_data[, date_name])
prelog_data <- setNames(lapply(groups, FUN = splitGroup, ungrouped_data = prelog_data, group_name = group_name, date_name = date_name, start_date = start_date, end_date = end_date, no_filter = c(group_name, date_name, outcome_name, denom_name)), groups)
#if (exists('exclude_group')) {prelog_data <- prelog_data[!(names(prelog_data) %in% exclude_group)]}

#Log-transform all variables, adding 0.5 to counts of 0.
ds <- setNames(lapply(prelog_data, FUN = logTransform, no_log = c(group_name, date_name,outcome_name)), groups)
time_points <- unique(ds[[1]][, date_name])


#Monthly dummies
if(n_seasons==4){x<-quarter(as.Date(time_points))}
if(n_seasons==12){x<-month(as.Date(time_points))}
if(n_seasons==3){
    x.m<-month(as.Date(time_points))
    x<-x.m
    x[x.m %in% c(1,2,3,4)]<-1
    x[x.m %in% c(5,6,7,8)]<-2
    x[x.m %in% c(9,10,11,12)]<-3
    }
season.dummies<-dummy(x)
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
covars_full <- setNames(lapply(ds, makeCovars, code_change = code_change,season.dummies=season.dummies,  intervention_date = intervention_date, time_points = time_points), groups)
covars_full <- sapply(covars_full, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_covar), drop = FALSE]})
covars_time <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies,time_index = 1:nrow(covars))))}), groups)

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
stl.covars<-mapply(smooth_func,ds.list=ds,covar.list=covars_full) 
post.start.index<-which(time_points==post_period[1])
stl.data.setup<-mapply(stl_data_fun,covars=stl.covars, ds.sub=ds )  #list of lists that has covariates for each regression for each strata

##SECTION 2: run first stage models
n_cores <- detectCores()-1
glm.results<- vector("list",  length=length(stl.data.setup)) #combine models into a list
cl1 <- makeCluster(n_cores)
clusterEvalQ(cl1, {library(lme4, quietly = TRUE)})
clusterExport(cl1, c('stl.data.setup',  'glm.fun', 'time_points', 'n_seasons','post.start.index'), environment())
for(i in 1:length(stl.data.setup)){
  glm.results[[i]]<-parLapply(cl=cl1 ,     stl.data.setup[[i]], fun=glm.fun )
}
stopCluster(cl1)
######################

#Combine the outcome, covariates, and time point information.
data_full <- setNames(lapply(groups, makeTimeSeries, outcome = outcome,       covars = covars_full), groups)
data_time <- setNames(lapply(groups, makeTimeSeries, outcome = outcome, covars = covars_time, trend=TRUE), groups)
data_pca<-mapply(FUN=pca_top_var,glm.results.in=glm.results, covars=stl.covars,ds.in=ds, SIMPLIFY=FALSE)
names(data_pca)<-groups

###############################
#                             #
#        Main analysis        #
#                             #
###############################

#Start Cluster for CausalImpact (the main analysis function).
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE)})
clusterExport(cl, c('doCausalImpact',  'intervention_date', 'time_points', 'n_seasons'), environment())

impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, var.select.on=TRUE, time_points = time_points, n_seasons = n_seasons), groups)
impact_time <- setNames(parLapply(cl, data_time, doCausalImpact, intervention_date = intervention_date,  var.select.on=FALSE,time_points = time_points, n_seasons = n_seasons, trend = TRUE), groups)
impact_pca <- setNames(parLapply(cl, data_pca, doCausalImpact, intervention_date = intervention_date, var.select.on=FALSE, time_points = time_points, n_seasons = n_seasons), groups)

stopCluster(cl)



#calculate WAIC
waic_full_all<-lapply(impact_full,waic_fun)
waic_time_all<-lapply(impact_time,waic_fun, trend=TRUE)
waic_pca_all<-lapply(impact_pca,waic_fun, trend=FALSE)
#Extract WAIC
waic_full<-sapply(waic_full_all, '[[', 'waic_2')
waic_time<-sapply(waic_time_all, '[[', 'waic_2')
waic_pca<-sapply(waic_pca_all, '[[', 'waic_2')
waic_combo<-cbind.data.frame(waic_full,waic_time, waic_pca)

#Extract PWAIC
p.waic_full<-sapply(waic_full_all, '[[', 'PWAIC_2')
p.waic_time<-sapply(waic_time_all, '[[', 'PWAIC_2')
p.waic_pca<-sapply(waic_pca_all, '[[', 'PWAIC_2')
pwaic_combo<-cbind.data.frame(p.waic_full,p.waic_time, p.waic_pca)

#Calculate WAIC weights for each stratum https://cran.r-project.org/web/packages/loo/vignettes/loo2-weights.html
waic_weights<- as.data.frame(round(t(apply(waic_combo,1, function(x) exp(-0.5*(x-min(x))) / sum(exp(-0.5*(x-min(x)))) ) ),2))
waic_weights$group<-groups
waic_weights.m<-melt(waic_weights, id='group')

#Not clear the log likelihoods are correctly calculated...?
log.lik.full<-lapply(waic_full_all, '[[', 'log.lik.mat')
log.lik.time<-lapply(waic_time_all, '[[', 'log.lik.mat')
log.lik.pca<-lapply(waic_pca_all, '[[', 'log.lik.mat')

#WEIGHTS FOR STACKING code--will probbaly need to do K-fold cross validation
# age.log.like<- vector("list",  length=length(log.lik.full)) #combine models into a list
# for(i in 1:length(log.lik.full)){
#   age.log.like[[i]]<-list(log.lik.full[[i]],log.lik.time[[i]],log.lik.pca[[i]])
# }
# r_eff_list.age<-lapply(age.log.like,r_eff_func1 )
# loo.weights<-round(t(mapply(FUN=loo_model_weights, x=age.log.like, r_eff_list=r_eff_list.age, method='stacking', cores=n_cores)),2)


#Save the inclusion probabilities from each of the models.
inclusion_prob_full <- setNames(lapply(impact_full, inclusionProb), groups)
inclusion_prob_time <- setNames(lapply(impact_time, inclusionProb), groups)

#All model results combined
quantiles_full <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_full[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period)}), groups)
quantiles_time <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_time[[group]], denom_data = ds[[group]][, denom_name],  eval_period = eval_period, post_period = post_period)}), groups)
quantiles_pca <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_pca[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period)}), groups)


#Model predicitons
pred_quantiles_full <- sapply(quantiles_full, getPred, simplify = 'array')
pred_quantiles_time <- sapply(quantiles_time, getPred, simplify = 'array')
pred_quantiles_pca <- sapply(quantiles_pca, getPred, simplify = 'array')

#Pointwise RR and uncertainty for second stage meta analysis
log_rr_quantiles   <- sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_quantiles}, simplify = 'array')
dimnames(log_rr_quantiles)[[1]] <- time_points
log_rr_sd   <- sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_sd}, simplify = 'array')
log_rr_full_t_samples.prec<-sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_samples.prec}, simplify = 'array')
saveRDS(log_rr_quantiles, file=paste0(output_directory, country, "_log_rr_quantiles.rds"))
saveRDS(log_rr_sd, file=paste0(output_directory, country, "_log_rr_sd.rds"))
saveRDS(log_rr_full_t_samples.prec, file=paste0(output_directory, country, "_log_rr_full_t_samples.prec.rds"))

#Rolling rate ratios
rr_roll_full <- sapply(quantiles_full, FUN = function(quantiles_full) {quantiles_full$roll_rr}, simplify = 'array')
rr_roll_time <- sapply(quantiles_time, FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array')
rr_roll_pca <- sapply(quantiles_pca, FUN = function(quantiles_pca) {quantiles_pca$roll_rr}, simplify = 'array')

#Rate ratios for evaluation period.
rr_mean_full <- t(sapply(quantiles_full, getRR))
rr_mean_time <- t(sapply(quantiles_time, getRR))
rr_mean_pca <- t(sapply(quantiles_pca, getRR))

rr_mean_full_intervals <- data.frame('SC Estimate (95% CI)'     = makeInterval(rr_mean_full[, 2], rr_mean_full[, 3], rr_mean_full[, 1]), check.names = FALSE, row.names = groups)
rr_mean_time_intervals <- data.frame('ITS Estimate (95% CI)' = makeInterval(rr_mean_time[, 2], rr_mean_time[, 3], rr_mean_time[, 1]), check.names = FALSE, row.names = groups)
rr_mean_pca_intervals <- data.frame('STL+PCA Estimate (95% CI)'     = makeInterval(rr_mean_pca[, 2], rr_mean_pca[, 3], rr_mean_pca[, 1]), check.names = FALSE, row.names = groups)

colnames(rr_mean_time) <- paste('ITS', colnames(rr_mean_time))

#Combine RRs into 1 file for plotting
rr_mean_combo<- as.data.frame(rbind( cbind(rep(1, nrow(rr_mean_full)),groups,  seq(from=1, by=1, length.out=nrow(rr_mean_full)),rr_mean_full),
                       cbind(rep(2, nrow(rr_mean_time)),groups, seq(from=1, by=1, length.out=nrow(rr_mean_full)), rr_mean_time),
                       cbind(rep(3, nrow(rr_mean_pca)), groups, seq(from=1, by=1, length.out=nrow(rr_mean_full)),rr_mean_pca)))
        names(rr_mean_combo)<-c('Model', 'groups', 'group.index','lcl','mean.rr','ucl')
        rr_mean_combo$waic_weights<-waic_weights.m$value
        rr_mean_combo$group.index<-as.numeric(as.character(rr_mean_combo$group.index))
        rr_mean_combo$mean.rr<-as.numeric(as.character(rr_mean_combo$mean.rr))
        rr_mean_combo$lcl<-as.numeric(as.character(rr_mean_combo$lcl))
        rr_mean_combo$ucl<-as.numeric(as.character(rr_mean_combo$ucl))
        rr_mean_combo$group.index[rr_mean_combo$Model==2]<-rr_mean_combo$group.index[rr_mean_combo$Model==2]+0.15
        rr_mean_combo$group.index[rr_mean_combo$Model==3]<-rr_mean_combo$group.index[rr_mean_combo$Model==3]+0.3
        rr_mean_combo$Model<-as.character(rr_mean_combo$Model)
        rr_mean_combo$Model[rr_mean_combo$Model=='1']<-"Synthetic Controls"
        rr_mean_combo$Model[rr_mean_combo$Model=='2']<-"Time trend"
        rr_mean_combo$Model[rr_mean_combo$Model=='3']<-"STL+PCA"
        cbPalette <- c("#1b9e77", "#d95f02", "#7570b3")
        rr_mean_combo$index<-as.factor(1:nrow(rr_mean_combo))
        #Fix order for axis
        rr_mean_combo$Model<-as.factor(rr_mean_combo$Model)
        rr_mean_combo$Model = factor(rr_mean_combo$Model,levels(rr_mean_combo$Model)[c(2,3,1)])
        #print(levels(rr_mean_combo$Model))

cumsum_prevented <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_full, simplify = 'array')
cumsum_prevented_pca <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_pca, simplify = 'array')
cumsum_prevented_time <- sapply(groups, FUN = cumsum_func, quantiles = quantiles_time, simplify = 'array')




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
 bad_sensitivity_groups <- sapply(covars_full, function (covar) {ncol(covar) <= 3})
 sensitivity_covars_full <- covars_full[!bad_sensitivity_groups]
 sensitivity_ds <- ds[!bad_sensitivity_groups]
 sensitivity_impact_full <- impact_full[!bad_sensitivity_groups]
 sensitivity_groups <- groups[!bad_sensitivity_groups]

 #Weight Sensitivity Analysis - top weighted variables are excluded and analysis is re-run.
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(pogit, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
clusterExport(cl, c('sensitivity_ds', 'doCausalImpact',  'weightSensitivityAnalysis', 'rrPredQuantiles', 'sensitivity_groups', 'intervention_date', 'outcome', 'time_points', 'n_seasons',  'eval_period', 'post_period'), environment())
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
