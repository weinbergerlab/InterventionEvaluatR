params <-
list(sensitivity = TRUE, crossval = TRUE)

## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup_packages, include=FALSE---------------------------------------
library(knitr)
library(InterventionEvaluatR)

## ----setup_data, include=FALSE-------------------------------------------
data(pnas_brazil, package = "InterventionEvaluatR")
analysis <- evaluatr.init(
  country = "Brazil", data = pnas_brazil,
  pre_period_start = "2004-01-01", pre_period_end = "2009-12-31",
  post_period_start = "2010-01-01", post_period_end = "2013-12-01",
  eval_period_start = "2012-01-01", eval_period_end = "2013-12-01",
  n_seasons = 12, year_def = "cal_year",
  group_name = "age_group", date_name = "date", outcome_name = "J12_18", denom_name = "ach_noj"
)
set.seed(1)

## ----main analysis, include = FALSE--------------------------------------
impact_results = evaluatr.impact(analysis)

## ----crossval, include = FALSE-------------------------------------------
if (params$crossval) {
  crossval_results <- evaluatr.crossval(analysis)
}

## ----sensitivity_analyses, include = FALSE-------------------------------
if (params$sensitivity) {
  sensitivity_results <- evaluatr.sensitivity(analysis)
}

## ----sparse--------------------------------------------------------------
if (!is.null(names(analysis$sparse_groups[analysis$sparse_groups])) && length(names(analysis$sparse_groups[analysis$sparse_groups])) != 0) {
  kable(data.frame("Sparse Groups" = names(analysis$sparse_groups[analysis$sparse_groups]), check.names = FALSE), align = "c")
}

## ----Comparison of estimates from different models-----------------------
if (params$crossval) {
  kable(cbind.data.frame(crossval_results$rr_mean_stack_intervals, impact_results$full$rr_mean_intervals, impact_results$time$rr_mean_intervals, impact_results$time_no_offset$rr_mean_intervals, impact_results$its$rr_mean_intervals, impact_results$pca$rr_mean_intervals), align = "c")
} else {
  kable(cbind.data.frame(impact_results$best$rr_mean_intervals, impact_results$full$rr_mean_intervals, impact_results$time$rr_mean_intervals, impact_results$time_no_offset$rr_mean_intervals, impact_results$its$rr_mean_intervals, impact_results$pca$rr_mean_intervals), align = "c")
}

## ----fig.width=5, fig.height=3, fig.align = "center", dpi=300, echo=FALSE----
plots <- evaluatr.plots(analysis)
plots$summary

## ----Comparison of Cross validation Weights from different models--------
if (params$crossval) {
  kable(crossval_results$stacking_weights.all, align = "c")
} else {
  print("Cross-validation not performed")
}

## ----modelsize-----------------------------------------------------------
kable(analysis$model_size, col.names = c("Model Size"))

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

## ----incl_table----------------------------------------------------------
kable(incl_probs, align = "c")

## ----sensitivity---------------------------------------------------------
if (exists("sensitivity_results")) {
  kable(sensitivity_results$sensitivity_table_intervals, align = "c")
}

## ----plots,fig.height =3 , fig.width = 5, fig.align = "center", dpi=300,results = 'asis'----
for (group in names(plots$groups)) {
  for (group_plot in plots$groups[[group]]) {
    print(group_plot)
  }
}

## ----save_results, echo=FALSE--------------------------------------------
output_file <- "Results" # Directory where results will be saved.
output_file <- paste0(output_file, "_", analysis$country, "_", format(Sys.time(), "%Y-%m-%d-%H%M%S"), ".Rds")
evaluatr.save(analysis, output_file)

