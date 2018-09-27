#Use this script after running the analysis script to write important results to file.

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

#Output the rate ratio estimates to a new file.
write.csv(rr_mean_full, paste(output_directory, country, '_rr_full.csv', sep = ''))
write.csv(rr_mean_time, paste(output_directory, country, '_rr_time_trend.csv', sep = ''))
write.csv(rr_mean_pca, paste(output_directory, country, '_rr_pca.csv', sep = ''))

write.csv(rr_mean_full_intervals, paste(output_directory, country, '_rr_full_intervals.csv', sep = ''))
write.csv(rr_mean_time_intervals, paste(output_directory, country, '_rr_time_trend_intervals.csv', sep = ''))
write.csv(rr_mean_pca_intervals, paste(output_directory, country, '_rr_pca_intervals.csv', sep = ''))

write.csv(rr_roll_full, paste(output_directory, country, '_rr_roll_full.csv', sep = ''))
write.csv(rr_roll_time, paste(output_directory, country, '_rr_roll_time.csv', sep = ''))
write.csv(rr_roll_pca, paste(output_directory, country, '_rr_roll_pca.csv', sep = ''))

#STacking weights
if(crossval){
write.csv(stacking_weights.all, paste(output_directory, country, '_stacking_weights.csv', sep = ''))
}
  
#Tables for rate ratios.
if (exists('sensitivity_table_intervals')) {
write.csv(sensitivity_table, paste(output_directory, country, '_sensitivity_table.csv', sep = ''))
write.csv(sensitivity_table_intervals, paste(output_directory, country, '_sensitivity_table_intervals.csv', sep = ''))
write.csv(rr_table, paste(output_directory, country, '_rr_table.csv', sep = ''))
write.csv(rr_table_intervals, paste(output_directory, country, '_rr_table_intervals.csv', sep = ''))
}
