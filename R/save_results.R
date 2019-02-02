syncon.save <- function (analysis, output_file) {
  dir.create(dirname(output_file),
             recursive = TRUE,
             showWarnings = FALSE)
  
  saved = list(
    version = 1,
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
