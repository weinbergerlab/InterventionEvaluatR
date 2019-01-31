syncon.save <- function (syncon, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  saved = list(
    version=1,
    input=list(
      country=syncon$country, 
      data=syncon$input_data,
      pre_priod=syncon$pre_period,
      post_period=syncon$post_period,
      eval_period=syncon$eval_period,
      n_seasons=syncon$n_seasons, year_def=syncon$year_def, 
      group_name=syncon$group_name, date_name=syncon$date_name, 
      outcome_name=syncon$outcome_name, denom_name=syncon$denom_name
    ),
    results=syncon$results
  )
  
  saveRDS(saved, output_file)
}
