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

## ----setup_packages, include=FALSE, echo=TRUE----------------------------
library(xtable)
library(knitr)
library(InterventionEvaluatR)

## ----viewdata, include=TRUE----------------------------------------------
    data(pnas_brazil, package = "InterventionEvaluatR") #load the data
    head(pnas_brazil[,1:5]) #View first few rows and columns
    
    pnas_brazil2<-pnas_brazil[pnas_brazil$age_group %in% c(9,8),] #Subset to age groups 8 and 9
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

