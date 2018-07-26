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

setwd('C:/Users/dmw63/Documents/GitHub/synthetic-control-poisson/main analysis components')

#Used to check for relevant packages and update them if out of date or install them if not installed.
update_packages  <- TRUE #Whether to update outdated packages.
install_packages <- TRUE #Whether to install missing packages.
install_pandoc   <- TRUE #Whether to install pandoc, which requires an external installer, and rmarkdown, a package that depends on pandoc's successful installation.

#Assign variable values
country       <- 'Brazil_mortality_' #Country or region name.
n_seasons     <- 4       #Number of months (seasons) per year. 12 for monthly, 4 for quarterly, 3 for trimester data.
exclude_covar <- c()      #User-defined list of covariate columns to exclude from all analyses.
exclude_group <- c()      #User-defined list of groups to exclude from analyses.
if(country=="Brazil_state"){code_change   <- TRUE     #Used for Brazil data. Set to TRUE to adjust for year 2008 coding changes; otherwise, set to FALSE.
    }else{
    code_change   <- FALSE
  }

input_directory  <- "https://raw.githubusercontent.com/weinbergerlab/PCV-mortality-Brazil/master/" #Directory or URL page containing input data file
file_name="HDI.csv"
output_directory <- '../Results'   #Directory where results will be saved.
output_directory <- paste(output_directory,'_', country,"_",format(Sys.time(),'%Y-%m-%d-%H%M%S'), '/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.
#Import direct from github
data_file <- paste0(input_directory, file_name)
prelog_data <- read.csv(text=getURL(data_file), check.names = FALSE)# IF IMPORTING FROM URL

#prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM LOCAL

#Import from local file

group_name   <- 'age_group' #Name of column containing group labels.
date_name    <- 'date'      #Name of column containing dates.
outcome_name <- 'J12_J18'    #Name of column containing outcome.
denom_name   <- 'NoPCV'   #Name of column containing denominator to be used in offset.

#MOST DATES MUST BE IN FORMAT "YYYY-MM-01", exception is end of pre period, which is 1 day before end of post period
start_date        <- as.Date('2004-01-01') #Indicates the date of the first data point.
intervention_date <- as.Date('2009-12-31') #Indicates the date of intervention in the data.
end_date          <- as.Date('2013-10-01') #Indicates the date of the last data point.
pre_period        <- as.Date(c('2004-01-01', '2009-12-31')) #Range over which the data is trained for the CausalImpact model.
post_period       <- as.Date(c('2010-01-01', '2013-10-01')) #Range from the intervention date to the end date.
eval_period       <- as.Date(c('2012-01-01', '2013-10-01')) #Range over which rate ratio calculation will be performed.
year_def          <-'cal_year'  #Can be cal_year to aggregate results by Jan-Dec; 'epi_year' to aggregate July-June

crossval=TRUE #run cross validation? Note this takes time...adds ~40 min with 10 age groups, 7 cores
#Run analysis, but don't generate HTML report
# source('synthetic_control_analysis.R', local = TRUE)
# source('synthetic_control_write_results.R', local = TRUE)
# source('synthetic_control_plot.R', local = TRUE)

#Run analysis and generate HTML report
source('synthetic_control_report.R', local = TRUE)
source('synthetic_control_write_results.R', local = TRUE) #save .csv files with output tables
