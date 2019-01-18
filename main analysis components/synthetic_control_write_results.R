syncon.write_results <- function (syncon, output_directory, prefix) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

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

