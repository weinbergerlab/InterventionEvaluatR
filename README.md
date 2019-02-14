# InterventionEvaluatR
This package is designed to run several different types of analyses to evaluate the effect of interventions using time series data of counts. The package implements synthetic controls analyses (Broderson; Bruhn), STL+PCA analysis (Shioda), a simple time trend analysis (with or without an offset term), and a classic interrupted time series analysis. In all approaches, a Poisson model is used with an observation-level random intercept to capture overdispersion. The synthetic controls analysis uses Bayesian variable selection (spike and slab priors, as implemented in the pogit package) to fit the model.

## Installing the package:
Install and load devtools in R. Then run:
 devtools::install_github('https://github.com/weinbergerlab/InterventionEvaluatR')
library(InterventionEvaluatR)

## Getting started
After loading the package, run: vignette('brazil')

## Sample Data

Sample data, which are described in Bruhn et al., are included with the package. There are five files, each containing time series of hospitalization for various conditions from five countries in the Americas. The data are subsetted by age group. 

##Citations
Brodersen, Kay H., et al. "Inferring causal impact using Bayesian structural time-series models." The Annals of Applied Statistics 9.1 (2015): 247-274.

Bruhn, Christian AW, et al. "Estimating the population-level impact of vaccines using synthetic controls." Proceedings of the National Academy of Sciences 114.7 (2017): 1524-1529.

Shioda, Kayoko, et al. "Challenges in estimating the impact of vaccination with sparse data." Epidemiology 30.1 (2019): 61-68.
