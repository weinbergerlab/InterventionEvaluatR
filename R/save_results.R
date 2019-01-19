sc_save_results <- function (syncon, output_file) {
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
  
  prefix = paste0(output_directory, prefix)

  with(syncon$results, {
    # Output the rate ratio estimates to a new file.
    write.csv(impact$full$rr_mean, paste0(prefix, '_rr_full.csv'))
    write.csv(impact$full$rr_mean_intervals, paste0(prefix, '_rr_full_intervals.csv'))
    write.csv(impact$full$rr_roll, paste0(prefix, '_rr_roll_full.csv'))
      
    write.csv(impact$time$rr_mean, paste0(prefix, '_rr_time_trend.csv'))
    write.csv(impact$time$rr_mean_intervals, paste0(prefix, '_rr_time_trend_intervals.csv'))
    write.csv(impact$time$rr_roll, paste0(prefix, '_rr_roll_time.csv'))
      
    write.csv(impact$pca$rr_mean, paste0(prefix, '_rr_pca.csv'))
    write.csv(impact$pca$rr_mean_intervals, paste0(prefix, '_rr_pca_intervals.csv'))
    write.csv(impact$pca$rr_roll, paste0(prefix, '_rr_roll_pca.csv'))
      
    write.csv(impact$best$rr_mean, paste0(prefix, '_rr_best.csv'))
    write.csv(impact$best$rr_mean_intervals, paste0(prefix, '_rr_best_intervals.csv'))
    write.csv(impact$best$rr_roll, paste0(prefix, '_rr_roll_best.csv'))

    # Stacking weights
    if(!is.na(cross_validation)) {
      write.csv(cross_validation$stacking_weights.all, paste0(prefix, '_stacking_weights.csv'))
    }
  
    # Tables for rate ratios.
    if (!is.na(sensitivity_analysis)) {
      write.csv(sensitivity_analysis$sensitivity, paste0(prefix, '_sensitivity_table.csv'))
      write.csv(sensitivity_analysis$sensitivity_intervals, paste0(prefix, '_sensitivity_table_intervals.csv'))
      write.csv(sensitivity_analysis$rr, paste0(prefix, '_rr_table.csv'))
      write.csv(sensitivity_analysis$rr_intervals, paste0(prefix, '_rr_table_intervals.csv'))
    }
  })
}

