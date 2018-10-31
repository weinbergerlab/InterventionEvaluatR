#This is the file used to set variables to be used in analysis, as well as to run the analysis.
#Make sure *_analysis.R, *_report.R, *_report.Rmd, *_functions.R, and *_plot.R are all in the same folder as this file.
#Model the setup shown in this file, then run this file from the console using source('This file's directory/This file's name'). 
#Clear the workspace
rm(list = ls(all = TRUE))
gc()
#source('synthetic_control_analysis.R', local = TRUE)
require(RCurl)
#Set the working directory
#This step only works as written if this file is run using the source() command. Otherwise, skip this step and set manually.

###WORKING DIRECTORY Should be set as the directory where .Rmd file is saved  ####
#setwd(auto.wd) ##automatically set working directory to '~desktop/synthetic-control-poisson-master/main analysis components/'

setwd('C:/Users/dmw63/Documents/synthetic-control-poisson/main analysis components')

#Used to check for relevant packages and update them if out of date or install them if not installed.
update_packages  <- TRUE #Whether to update outdated packages.
install_packages <- TRUE #Whether to install missing packages.
install_pandoc   <- TRUE #Whether to install pandoc, which requires an external installer, and rmarkdown, a package that depends on pandoc's successful installation.

#Assign variable values
country       <- 'br_monthly_SubChapter' #Country or region name.
n_seasons     <- 12       #Number of months (seasons) per year. 12 for monthly, 4 for quarterly, 3 for trimester data.
exclude_covar <- c('J09_J18_prim','J00_J99_prim' ,'J12_J18_any', 'J13_any',
                   'possible_pneumo_code',"population" ,"acm_noj_nodiarr_prim",
                   'A00_A09_prim', 'A30_49_prim', 'B95_B98_prim', 'G00_G09_prim', 'H65_H75_prim',
                   'F00_F99_prim', "O00_O99_prim")      #User-defined list of covariate columns to exclude from all analyses.
exclude_group <- c()      #User-defined list of groups to exclude from analyses.
if(country=="Brazil"){code_change   <- TRUE     #Used for Brazil data. Set to TRUE to adjust for year 2008 coding changes; otherwise, set to FALSE.
}else{
  code_change   <- FALSE
}

### ---------------- Import data ---------------- ###
input_directory  <- "C:/Users/dmw63/Your team Dropbox/PAHO mortality/Data/" #Directory or URL page containing input data file
file_name="PAHO all age cuts_SubChapters.csv"
output_directory <- "../Results/PAHO_Mortality"   #Directory where results will be saved.
output_directory <- paste(output_directory,'_', country,"_",format(Sys.time(),'%Y-%m-%d-%H%M%S'), '/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.
#Import direct from github
data_file <- paste0(input_directory, file_name)
prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM URL
prelog_data<-prelog_data[substr(prelog_data$age_group,1,2)=='br',]

table(prelog_data$age_group)[which(table(prelog_data$age_group)>0)]

### ---------------- For age groups with imcomplete data ---------------- ###

# Any age groups with incomplete data?
for (i in 1:length(unique(prelog_data$age_group))) {
  print(length(prelog_data$monthdate[which(prelog_data$age_group==unique(prelog_data$age_group)[i])]) / 12)
}

# For those age groups with incomplete data, add rows with zeros for missing months
CorrectNumMonths <- length(seq(range(as.Date(prelog_data$monthdate))[1], range(as.Date(prelog_data$monthdate))[2], by="month"))
CorrectMonths <- data.frame(monthdate = seq(range(as.Date(prelog_data$monthdate))[1], range(as.Date(prelog_data$monthdate))[2], by="month"))
CorrectMonths$monthdate <- as.Date(CorrectMonths$monthdate)
AgeGrp <- unique(prelog_data$age_group)
for (i in 1:length(AgeGrp)) {
  dt <- prelog_data[which(prelog_data$age_group==AgeGrp[i]),]
  dt$monthdate <- as.Date(dt$monthdate)
  dt<-dt[order(dt$monthdate),]
  if (length(dt$monthdate) != CorrectNumMonths) {
    FixMonths <- merge(CorrectMonths, dt, by = "monthdate", all.x=T)
    FixMonths$age_group <- AgeGrp[i]
    # Add zeros to missing months
    FixMonths[is.na(FixMonths)] <- 0
    # Replace with the corrected data
    prelog_data <- prelog_data[-c(which(prelog_data$age_group==AgeGrp[i])),] # Remove data for this specific age group
    prelog_data$monthdate <- as.Date(prelog_data$monthdate) 
    prelog_data <- rbind(prelog_data, FixMonths) # And merge it with corrected data
  } 
}

# Check
table(prelog_data$age_group)[which(table(prelog_data$age_group)>0)]


### ---------------- Filtering for sparse data ---------------- ###

# Sort data by age group and month before filtering 
prelog_data <- prelog_data[order(prelog_data$age_group, prelog_data$monthdate),] 

spl1<-split(prelog_data, as.character(prelog_data$age_group))
n.obs<-  sapply( spl1, function(x) nrow(x) )
n.deaths<-  sapply( spl1, function(x) sum(x$J12_J18_prim) )
n.obs2<-rep(n.obs, times=n.obs)
n.deaths2<-rep(n.deaths, times=n.obs)
n.deaths

# Filter prelog data to get rid of time series missing observations or with very few
nrow(prelog_data)
nrow(prelog_data[n.obs2==max(n.obs2) & n.deaths2>500,])
prelog_data<-prelog_data[n.obs2==max(n.obs2) & n.deaths2>500,]

# Check
table(prelog_data$age)[which(table(prelog_data$age)>0)]



prelog_data$age_group<-factor(prelog_data$age_group)
prelog_data$monthdate<-factor(prelog_data$monthdate)

group_name   <- 'age_group' #Name of column containing group labels.
date_name    <- 'monthdate'      #Name of column containing dates.
outcome_name <- 'J12_J18_prim'    #Name of column containing outcome.
denom_name   <- 'acm_noj_prim'   #Name of column containing denominator to be used in offset.

#MOST DATES MUST BE IN FORMAT "YYYY-MM-01", exception is end of pre period, which is 1 day before end of post period
start_date        <- as.Date('2005-01-01') #Indicates the date of the first data point.
intervention_date <- as.Date('2010-02-28') #Indicates the date of intervention in the data.
end_date          <- as.Date('2015-12-01') #Indicates the date of the last data point.
pre_period        <- as.Date(c('2005-01-01', '2010-02-28')) #Range over which the data is trained for the CausalImpact model.
post_period       <- as.Date(c('2010-03-01', '2015-12-01')) #Range from the intervention date to the end date.
eval_period       <- as.Date(c('2011-03-01', '2015-12-01')) #Range over which rate ratio calculation will be performed.
year_def   <-'cal_year'  #Can be cal_year to aggregate results by Jan-Dec; 'epi_year' to aggregate July-June

crossval=FALSE #run cross validation? Note this takes time...adds ~40 min with 10 age groups, 7 cores
sensitivity_switch=FALSE
#Run analysis, but don't generate HTML report
# source('synthetic_control_analysis.R', local = TRUE)
# source('synthetic_control_write_results.R', local = TRUE)
# source('synthetic_control_plot.R', local = TRUE)

###########################################################3
# To debug the analysis.R code
# prelog_data<-prelog_data[which(prelog_data$age_group=='br_1_1'|prelog_data$age_group=='br_1_2' | prelog_data$age_group=='br_2_2'),]

#Run analysis and generate HTML report
source('synthetic_control_report.R', local = TRUE)
source('synthetic_control_write_results.R', local = TRUE) #save .csv files with output tables

# #-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----#
# pca_ranking
# save(pca_ranking,file="ar_pca_ranking_monthly.RData")
# par(las=2)
# par(mar=c(10,4,4,4))
# groups
# par(mfrow=c(1,4))
# i=9
# barplot(pca_ranking[[i]][order(pca_ranking[[i]])], main=paste(groups[i]))
# #-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----#

library(beepr)
beep(3)
