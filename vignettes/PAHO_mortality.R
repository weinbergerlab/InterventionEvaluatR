params <-
list(sensitivity = TRUE, crossval = FALSE)

## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 3,
  fig.width = 5,
  fig.align = "center", 
  dpi=300, 
	out.width="600px"
)

## ------------------------------------------------------------------------
#install.packages('devtools') 
#library(devtools) 
#devtools::install_github('https://github.com/weinbergerlab/InterventionEvaluatR') 

## ----setup_packages, include=FALSE, echo=TRUE----------------------------
library(xtable)
library(knitr)
library(InterventionEvaluatR)

## ---- echo=F-------------------------------------------------------------
#Save the Ecuador mortality data as part of the package
# library(RCurl)
# script <- getURL("https://raw.githubusercontent.com/weinbergerlab/paho-pneumonia-mortality/master/Data/PAHO%20all%20age%20cuts_SubChapters.csv", ssl.verifypeer = FALSE)
# 
# paho1 <-read.csv(text=getURL("https://raw.githubusercontent.com/weinbergerlab/paho-pneumonia-mortality/master/Data/PAHO%20all%20age%20cuts_SubChapters.csv"),  header=T)
# ec1<-paho1[paho1$age_group %in% c('ec 2-59m A','ec 2-23m A','ec 24-50m A'), ]
#ec1$age_group<-factor(ec1$age_group)
# saveRDS(ec1,'../data/paho_ecuador_mortality.rds')

## ----viewdata, include=TRUE----------------------------------------------
    data(pnas_brazil, package = "InterventionEvaluatR") #load the data
    head(pnas_brazil[,1:5]) #View first few rows and columns

## ------------------------------------------------------------------------

pnas_brazil$date<-as.Date(pnas_brazil$date, tryFormats=c('%Y-%m-%d',
                                                    '%m-%d-%Y',
                                                    '%m/%d/%Y',
                                                    '%Y/%m/%d',
                                                    '%d/%m/%Y'
                                                    ) )

## ------------------------------------------------------------------------
    pnas_brazil2<-pnas_brazil[pnas_brazil$age_group %in% c(9),] #Subset to age groups 8 and 9
    pnas_brazil2<-pnas_brazil2[order(pnas_brazil2$age_group, pnas_brazil2$date),] #Sort data by age group and month
    pnas_brazil2<-pnas_brazil2[as.Date(pnas_brazil2$date)>=as.Date('2004-01-01'),] #Ignore 2003

## ----setup_data, echo=TRUE-----------------------------------------------

analysis <- evaluatr.init(
  country = "Brazil", data = pnas_brazil2,
  post_period_start = "2010-01-01", #First 'post-intervention' month is Jan 2012
  eval_period_start = "2012-01-01", #We ignore first 2 years of data to allow for vaccine ramp up
  eval_period_end = "2013-12-01", #The evaluation period lasts 2 years
  n_seasons = 12, #This is monthly data, so select 12
  year_def = "cal_year", # we are in southern hemisphere, so aggregate results by calendar year (Jan-Dec)
  group_name = "age_group",  #Strata categry name
  date_name = "date", #Date variable name
  outcome_name = "J12_18", #Outcome variable name
  denom_name = "ach_noj" #Denominator variable name
)
set.seed(1)

## ----univariate, fig.width=3, fig.height=5-------------------------------
 glmer_results= evaluatr.univariate(analysis)
 lapply(glmer_results,evaluatr.univariate.plot)

## ----main analysis, include = FALSE--------------------------------------
impact_results = evaluatr.impact(analysis)

## ----sensitivity_analyses, include = FALSE-------------------------------
if (params$sensitivity) {
  sensitivity_results <- evaluatr.sensitivity(analysis)
}

## ----sparse, results="asis"----------------------------------------------
if (!is.null(names(analysis$sparse_groups[analysis$sparse_groups])) && length(names(analysis$sparse_groups[analysis$sparse_groups])) != 0) {
  print(xtable(data.frame("Sparse Groups" = names(analysis$sparse_groups[analysis$sparse_groups]), check.names = FALSE), align = "cc"), type="html")
}

## ----Comparison of estimates from different models, results="asis"-------
if (params$crossval) {
  print(xtable(cbind.data.frame(crossval_results$rr_mean_stack_intervals, impact_results$full$rr_mean_intervals, impact_results$time$rr_mean_intervals, impact_results$time_no_offset$rr_mean_intervals, impact_results$its$rr_mean_intervals, impact_results$pca$rr_mean_intervals), align = "ccccccc"), type="html")
} else {
  print(xtable(cbind.data.frame(impact_results$best$rr_mean_intervals, impact_results$full$rr_mean_intervals, impact_results$time$rr_mean_intervals, impact_results$time_no_offset$rr_mean_intervals, impact_results$its$rr_mean_intervals, impact_results$pca$rr_mean_intervals), align = "ccccccc"), type="html")
}

## ------------------------------------------------------------------------
last.point<-dim(impact_results$best$cumsum_prevented)[1]
cum.prevented<-impact_results$best$cumsum_prevented[last.point,,]
print(cum.prevented)

## ----mainplot1, echo=FALSE-----------------------------------------------
plots <- evaluatr.plots(analysis)
plots$summary

## ----modelsize, results="asis"-------------------------------------------
model_size = data.frame(t(analysis$model_size))
print(xtable(setNames(model_size, c("Model Size")), align=rep("c", ncol(model_size) + 1)), type="html")

## ----incl, include = FALSE-----------------------------------------------
incl_probs <- NULL
for (group in analysis$groups) {
  incl_prob <- impact_results$full$groups[[group]]$inclusion_probs[-c(1:(analysis$n_seasons - 1)), ]
  incl_prob <- incl_prob[order(-incl_prob$inclusion_probs), ]
  incl_prob <- incl_prob[c(1:3), ]
  incl_prob2 <- incl_prob[, 2]
  incl_prob_names <- incl_prob[, 1]
  incl_prob3 <- data.frame("Group" = group, "Greatest Inclusion Variable" = incl_prob_names[1], "Greatest Inclusion Probability" = incl_prob2[1], "Second Greatest Inclusion Variable" = incl_prob_names[2], "Second Greatest Inclusion Probability" = incl_prob2[2], "Third Greatest Inclusion Variable" = incl_prob_names[3], "Third Greatest Inclusion Probability" = incl_prob2[3], check.names = FALSE)
  incl_probs <- rbind(incl_probs, incl_prob3)
}
rownames(incl_probs) <- NULL

## ----incl_table, results="asis"------------------------------------------
print(xtable(incl_probs, align = rep("c", ncol(incl_probs) + 1)), type="html")

## ----sensitivity, results="asis"-----------------------------------------
if (exists("sensitivity_results")) {
  print(xtable(sensitivity_results$sensitivity_table_intervals, align = rep("c", ncol(sensitivity_results$sensitivity_table_intervals) + 1)), type="html")
}

## ----plots, results = 'asis', fig.width=5, fig.height=4------------------
for (group in names(plots$groups)) {
      par(mfrow=c(4,1))
      print(plots$groups[[group]]$pred_full )
      print(plots$groups[[group]]$pred_best )
      print(plots$groups[[group]]$pred_time )
      print(plots$groups[[group]]$pred_pca )
}

## ----plots2, results = 'asis', fig.width=5, fig.height=4-----------------
for (group in names(plots$groups)) {
      par(mfrow=c(4,1))
      print(plots$groups[[group]]$pred_full_agg )
      print(plots$groups[[group]]$pred_best_agg )
      print(plots$groups[[group]]$pred_time_agg )
      print(plots$groups[[group]]$pred_pca_agg )
}

## ----plots3, results = 'asis', fig.width=5, fig.height=4-----------------
for (group in names(plots$groups)) {
      par(mfrow=c(4,1))
      print(plots$groups[[group]]$cumsum_prevented )
}

## ------------------------------------------------------------------------
impact_results$best$rr_mean
impact_results$best$rr_mean_hdi


## ------------------------------------------------------------------------
matplot(impact_results$best$pred_quantiles_HDI[,,1], type='l', col='gray', lty=c(2,1,2), bty='l', ylim=c(0, max(impact_results$best$pred_quantiles_HDI[,,1])), ylab='Count')
matplot(impact_results$best$pred_quantiles[,,1], type='l', col='red', lty=c(2,1,2), bty='l', ylim=c(0, max(impact_results$best$pred_quantiles_HDI[,,1])), ylab='Count', add=T)
points(analysis$input_data[, analysis$outcome_name])

## ----save_results, echo=FALSE--------------------------------------------
output_file <- "Results" # Directory where results will be saved.
output_file <- paste0(output_file, "_", analysis$country, "_", format(Sys.time(), "%Y-%m-%d-%H%M%S"), ".Rds")
evaluatr.save(analysis, output_file)

