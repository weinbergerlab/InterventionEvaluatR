% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{evaluatr.init}
\alias{evaluatr.init}
\title{Initialize analysis}
\usage{
evaluatr.init(country, data, pre_period_start = "start",
  post_period_start, post_period_end = eval_period_end,
  eval_period_start, eval_period_end, n_seasons = 12,
  year_def = "cal_year", group_name, date_name, outcome_name,
  set.burnN = 5000, set.sampleN = 10000, denom_name,
  log.covars = TRUE, ridge = F, error_dist = "iid",
  sparse_threshold = 5)
}
\arguments{
\item{country}{A one-word label for output (eg country or state name).}

\item{data}{Input dataframe. There should be a row for each time point (e.g., month) and category (e.g., age group). There should be variables for the date, for the category (or a column of 1s if only 1 category), for the outcome variable (a count), for the denominator (or a column of 1s if no denominator), and columns for each control variable.}

\item{pre_period_start}{Date when analysis starts, YYYY-MM-01. defaults to first date in dataset.}

\item{post_period_start}{Month when intervention introduced. YYY-MM-01}

\item{post_period_end}{Date when analysis ends, YYYY-MM-01. defaults to first date in dataset. Defaults to end of evaluation period.}

\item{eval_period_start}{First month of the period when the effect of intervention is evaluated. YYYY-MM-01. typically 12-24 months after post_period_start.}

\item{eval_period_end}{Last month of the period when the effect of intervention is evaluated. YYYY-MM-01.}

\item{n_seasons}{How many observations per year? Defaults to 12 (monthly data) Change to 4 for quarterly}

\item{year_def}{Should results be aggregated by calendar year ('cal_year': the default) or epidemiological year ('epi_year'; July-June)}

\item{group_name}{Name of the stratification variable (e.g., age group). If only one age group present, add a column of 1s to the dataset}

\item{date_name}{Name of the variable with the date for the time series}

\item{outcome_name}{Name of the outcome (y) variable in the 'data' dataframe. Should be a count}

\item{set.burnN}{Number of MCMC iterations for burn in (default 5000),}

\item{set.sampleN}{Number of MCMC iterations post-burn-in to use for inference (default 10000),}

\item{denom_name}{Name of the denominator variable in the 'data' dataframe. if there is no denominator, include a column of 1s.}

\item{log.covars}{Should the covariate be log transformed? (default: TRUE)}

\item{ridge}{Run ridge regression with AR(1) random intercepts (faster) or spike and slab (with iid random intercept) for variable selection. Logical, Default TRUE.}

\item{error_dist}{For the INLA models: use an 'iid' or 'ar1' error on the random intercept. Defaults to iid. Use caution with AR(1) as it can introduce bias in some situations}

\item{sparse_threshold}{Threshold for filtering out control variables based on sparsity (mean number of cases per time period). Defaults to 5.}
}
\value{
Initialized analysis object, `analysis` as described below

`analysis$country` as passed to `country`

`analysis$input_data` as passed to `data`

`analysis$n_seasons` as passed to `n_seasons`

`analysis$year_def` as passed to `year_def`

`analysis$pre_period` Range of dates in the pre-intervention period

`analysis$post_period` Range of dates in the post-intervention period

`analysis$eval_period` Range of dates in the evaluation period

`analysis$start_date` First date of the pre-intervention period

`analysis$intervention_date` Last time point before the start of the post-period

`analysis$end_date` Last date in the evaluation period

`analysis$group_name` as passed to `group_name`

`analysis$date_name` as passed in in `date_name`

`analysis$outcome_name` as passed in in `outcome_name`

`analysis$denom_name` as passed in in `denom_name`

`analysis$time_points` Vector of time points in the dataset

`analysis$set.burnN` as passed in in `set.burnN`

`analysis$set.sampleN` as passed in in `set.sampleN`

`analysis$log.covars` as passed in in `log.covars`
  
`analysis$groups` Vector of groups analyzed

`analysis$sparse_groups` Vector indicating which groups were too sparse to analyze

`analysis$model_size` Average number of covariates included in the synthetic control model

`analysis$covars` Matrix of covariates used for analysis

`analysis$outcome` as passeed to `outcome_name`

`analysis$ridge` as passeed to `ridge`
`analysis$error_dist` as passeed to `error_dist`
`analysis$sparse_threshold` as passed to `sparse_threshold`
}
\description{
Initialize analysis
}
