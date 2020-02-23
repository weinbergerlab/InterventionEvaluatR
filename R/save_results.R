#' Save data for TODO
#'
#' @param analysis Analysis object, initialized by TODO.init.
#' @param output_file An Rds file in which analysis results will be saved.
#' @param prune If TRUE (default) then diagnostic data is not saved. If FALSE, diagnostic data is saved (making the output file much larger).
#'
#' @export

evaluatr.save <- function (analysis, output_file, prune=TRUE) {
  dir.create(dirname(output_file),
             recursive = TRUE,
             showWarnings = FALSE)
  
  if (prune) {
    analysis = evaluatr.prune(analysis)
  }
  saved = list(
    version = 2,
    input = list(
      country = analysis$country,
      data = analysis$input_data,
      pre_priod = analysis$pre_period,
      post_period = analysis$post_period,
      eval_period = analysis$eval_period,
      n_seasons = analysis$n_seasons,
      year_def = analysis$year_def,
      group_name = analysis$group_name,
      date_name = analysis$date_name,
      outcome_name = analysis$outcome_name,
      denom_name = analysis$denom_name
    ),
    results = analysis$results
  )
  
  saveRDS(saved, output_file)
}

#' Prune unneeded data from the analysis object. Used internally before saving to reduce size of saved data.

#' @export
#' @param what What to prune. 'impact' means prune intermediate results of impact analysis, but leave around data needed for sensitivity analysis and cross-validation. 'all' means prune everything.
#' @keywords internal

evaluatr.prune <- function(analysis, what='all') {
  if (what == 'all') {
    analysis$.private$ds = NULL
    analysis$.private$data = NULL
    analysis$.private$data.cv = NULL
  }
  
  for (variant in c("best", names(analysis$.private$variants))) {
    if (what == "all") {
      analysis$results$impact[[variant]]$log_rr_full_t_samples.prec = NULL
      analysis$results$impact[[variant]]$log_rr_quantiles = NULL
      analysis$results$impact[[variant]]$quantiles = NULL
      analysis$results$impact[[variant]]$groups = NULL
    } else if (what == "impact") {
      for (group in names(analysis$results$impact[[variant]]$groups)) {
        analysis$results$impact[[variant]]$groups[[group]]$reg.mean = NULL
        analysis$results$impact[[variant]]$groups[[group]]$rand.eff = NULL
        analysis$results$impact[[variant]]$groups[[group]]$predict.bsts = NULL
        analysis$results$impact[[variant]]$groups[[group]]$beta.mat = NULL
      }

      for (group in names(analysis$results$impact[[variant]]$quantiles)) {
        analysis$results$impact[[variant]]$quantiles[[group]]$pred_samples_post_full = NULL
        analysis$results$impact[[variant]]$quantiles[[group]]$pred_samples = NULL
      }

      analysis$results$impact[[variant]]$rr_iter = NULL
    }
  }
  
  analysis
}
