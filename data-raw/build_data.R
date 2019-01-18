# Source this script in package root directory to (re)generate R data files from raw CSV files
datasets = list(
  pnas_us_ipd="PNAS/S6-US-IPD.csv", 
  pnas_us_pneumonia="PNAS/S5-US-pneumonia.csv", 
  pnas_brazil="PNAS/S1-Brazil.csv", 
  pnas_ecuador="PNAS/S3-Ecuador.csv", 
  pnas_mexico="PNAS/S4-Mexico.csv", 
  pnas_chile="PNAS/S2-Chile.csv"
)

for (dataset in names(datasets)) {
  assign(dataset, read.csv(paste0("data-raw/", datasets[[dataset]])))
  eval(substitute(use_data(d, overwrite=TRUE), list(d=as.name(dataset))))
}
