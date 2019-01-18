#setwd("~/synthetic-control-poisson")
rmarkdown::render("./main analysis components/synthetic_control_report self contained.Rmd", output_file ='br.html', 
                  params = list(country = "Brazil", 
                                pre_period_start= '2004-01-01',
                                pre_period_end= '2009-12-31',
                                post_period_start= '2010-01-01', #when is vaccine introduced?
                                post_period_end= '2013-12-01',
                                eval_period_start= '2012-01-01',#when do we start measuring the RR? (ie 12 or 24 months post intro) 
                                eval_period_end= '2013-12-01',
                              year_def = 'cal_year',  #Can be cal_year to aggregate results by Jan-Dec; 'epi_year' to aggregate July-June
                              sensitivity=TRUE, #Run sensitivity analyses?
                              crossval=FALSE, #run cross validation? Note this takes time...adds ~40 min with 10 age groups, 7 cores
                              group_name = 'age_group', #Name of column containing group labels.
                              date_name = 'date',      #Name of column containing dates.
                              outcome_name = 'J12_18',    #Name of column containing outcome.
                              denom_name   = 'ach_noj',   #Name of column containing denominator to be used in offset.
                              n_seasons = 12, #12 for monthly, 4 for quarterly
                              input_directory  = 'https://raw.githubusercontent.com/weinbergerlab/synthetic-control/master/Datasets%20for%20PNAS/', #Directory (or URL) containing input data file. ,
                              file_name="Dataset%20S1%20Brazil.csv", #name of csv file with data
                              github.import=TRUE, #Import from github? or local
                              update_packages= FALSE ,
                              install_packages= TRUE
                      ))