#' Generate plots for TODO
#'
#' @param analysis Analysis object, initialized by TODO.init.
#' @return Plots, `plots`, as described below
#'
#' `plots$summary` TODO
#' `plots$groups` TODO
#'
#' @importFrom graphics legend matplot par points
#' @importFrom ggplot2 aes aes_string aes_ geom_line scale_fill_manual geom_point geom_polygon geom_ribbon geom_vline ggplot ggtitle scale_colour_manual scale_fill_hue theme theme_bw element_blank element_line element_rect element_text labs geom_errorbar scale_shape_manual scale_size_manual guides geom_hline ylim coord_cartesian scale_y_continuous sec_axis
#'
#' @export

evaluatr.plots <- function(analysis) {
  plots = list(groups = list())
  cumsum_prevented_plot = list(groups = list())
  impact_results = analysis$results$impact
  cbPalette <- c("#1b9e77", "#d95f02", "#7570b3", 
                 "#e7298a")
  
  for (group in analysis$groups) {
    covars.sub <- analysis$covars$full[[group]][, -c(1:(analysis$n_seasons - 
                                                          1))]
    alpha1 = rep(impact_results$full$inclusion_probs[[group]][-c(1:(analysis$n_seasons - 
                                                                      1)), "inclusion_probs"], each = nrow(analysis$covars$full[[group]]))
    covar_plot <- ggplot(melt(covars.sub, id.vars = NULL), 
                         mapping = aes_string(x = rep(analysis$time_points, 
                                                      ncol(covars.sub)), y = "value", group = "variable", 
                                              alpha = alpha1)) + geom_line() + labs(x = "Time", 
                                                                                    y = "Scaled Covariates") + ggtitle(paste(group, 
                                                                                                                             "Scaled Covariates Weighted by Inclusion Probability")) + 
      theme_bw() + theme(axis.line = element_line(colour = "black"), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                         panel.border = element_blank(), panel.background = element_blank()) + 
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    min_max <- c(min(c(impact_results$full$pred_quantiles[, 
                                                          , group], analysis$outcome[, group])), max(c(impact_results$full$pred_quantiles[, 
                                                                                                                                          , group], analysis$outcome[, group])))
    pred_full_plot <- plotPred(impact_results$full$pred_quantiles[, 
                                                                  , group], analysis$time_points, analysis$post_period, 
                               min_max, analysis$outcome[, group], title = paste(group, 
                                                                                 "Synthetic controls estimate"))
    pred_time_plot <- plotPred(impact_results$time$pred_quantiles[, 
                                                                  , group], analysis$time_points, analysis$post_period, 
                               min_max, analysis$outcome[, group], title = paste(group, 
                                                                                 "Time trend adjustment estimate (with offset)"))
    if (length(impact_results$pca$pred_quantiles[, , group]) > 
        0) {
      pred_pca_plot <- plotPred(impact_results$pca$pred_quantiles[, 
                                                                  , group], analysis$time_points, analysis$post_period, 
                                min_max, analysis$outcome[, group], title = paste(group, 
                                                                                  "STL+PCA estimate"))
    }else {
      pred_pca_plot <- NA
    }
 
    
    pred_full_plot_agg <- plotPredAgg(impact_results$full$ann_pred_quantiles[[group]], 
                                      analysis$time_points, analysis$year_def, analysis$intervention_date, 
                                      analysis$post_period, min_max, analysis$outcome[, 
                                                                                      group], title = paste(group, "Synthetic controls estimate"))
    pred_time_plot_agg <- plotPredAgg(impact_results$time$ann_pred_quantiles[[group]], 
                                      analysis$time_points, analysis$year_def, analysis$intervention_date, 
                                      analysis$post_period, min_max, analysis$outcome[, 
                                                                                      group], title = paste(group, "Time trend adjustment estimate (with offset)"))
    if (length(impact_results$pca$ann_pred_quantiles[[group]]) > 
        0) {
      pred_pca_plot_agg <- plotPredAgg(impact_results$pca$ann_pred_quantiles[[group]], 
                                       analysis$time_points, analysis$year_def, analysis$intervention_date, 
                                       analysis$post_period, min_max, analysis$outcome[, 
                                                                                       group], title = paste(group, "STL+PCA estimate"))
    }else {
      pred_pca_plot_agg <- NA
    }
    cumsum.ds <- as.data.frame(impact_results$full$cumsum_prevented[, 
                                                                    , group])
    last.cumsum <- round(cumsum.ds[nrow(cumsum.ds), ])
    prevented.print <- paste0(last.cumsum["50%"], " (", 
                              last.cumsum["2.5%"], ", ", last.cumsum["97.5%"], 
                              ")")
    cumsum_prevented_plot[[group]] <- ggplot(cumsum.ds) + 
      geom_line(aes_(y = cumsum.ds[, "50%"], x = analysis$time_points, 
                     colour = "Cumulative Cases prevented")) + 
      geom_ribbon(aes_(ymin = cumsum.ds[, "2.5%"], 
                       ymax = cumsum.ds[, "97.5%"], x = analysis$time_points, 
                       fill = "band"), alpha = 0.3) + scale_colour_manual("", 
                                                                          values = "black") + scale_fill_manual("", 
                                                                                                                values = "grey12") + scale_y_continuous(position = "right") + 
      theme(legend.position = "none") + labs(x = "Time", 
                                             y = "Cumulative Cases Prevented") + ggtitle(paste(group, 
                                                                                               prevented.print, "Cases Prevented")) + theme_bw() + 
      theme(axis.line = element_line(colour = "black"), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_blank(), panel.background = element_blank()) + 
      theme(legend.title = element_blank(), legend.position = "none", 
            plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
    if (!is.list(impact_results$full$rr_roll)) {
      min_max <- c(min(impact_results$full$rr_roll[, , 
                                                   group], impact_results$time$rr_roll[, , group]), 
                   max(impact_results$full$rr_roll[, , group], impact_results$time$rr_roll[, 
                                                                                           , group]))

      if (length(impact_results$full$rr_roll[, , group]) > 
          0) {
        rr_roll_full_plot <- ggplot(melt(as.data.frame(impact_results$full$rr_roll[, 
                                                                                   , group]), id.vars = NULL), mapping = aes_string(x = rep(analysis$time_points[(length(analysis$time_points) - 
                                                                                                                                                                    nrow(impact_results$full$rr_roll[, , group]) + 
                                                                                                                                                                    1):length(analysis$time_points)], ncol(impact_results$full$rr_roll[, 
                                                                                                                                                                                                                                       , group])), y = "value", linetype = "variable")) + 
          geom_line() + geom_hline(yintercept = 1, linetype = 4) + 
          labs(x = "Time", y = "Rolling Rate Ratio") + 
          ggtitle(paste(group, "Synthetic Control Rolling Rate Ratio")) + 
          coord_cartesian(ylim = min_max) + theme_bw() + 
          theme(axis.line = element_line(colour = "black"), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.border = element_blank(), panel.background = element_blank()) + 
          theme(legend.title = element_blank(), legend.position = c(0, 
                                                                    1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, 
                                                                                                                                         fill = "transparent"), plot.title = element_text(hjust = 0.5), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
      else {
        rr_roll_full_plot <- NA
      }
      if (length(impact_results$time$rr_roll[, , group]) > 
          0) {
        rr_roll_time_plot <- ggplot(melt(as.data.frame(impact_results$time$rr_roll[, 
                                                                                   , group]), id.vars = NULL), mapping = aes_string(x = rep(analysis$time_points[(length(analysis$time_points) - 
                                                                                                                                                                    nrow(impact_results$time$rr_roll[, , group]) + 
                                                                                                                                                                    1):length(analysis$time_points)], ncol(impact_results$time$rr_roll[, 
                                                                                                                                                                                                                                       , group])), y = "value", linetype = "variable")) + 
          geom_line() + geom_hline(yintercept = 1, linetype = 4) + 
          labs(x = "Time", y = "Rolling Rate Ratio") + 
          ggtitle(paste(group, "TT Rolling Rate Ratio")) + 
          coord_cartesian(ylim = min_max) + theme_bw() + 
          theme(axis.line = element_line(colour = "black"), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.border = element_blank(), panel.background = element_blank()) + 
          theme(legend.title = element_blank(), legend.position = c(0, 
                                                                    1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, 
                                                                                                                                         fill = "transparent"), plot.title = element_text(hjust = 0.5), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
      else {
        rr_roll_time_plot <- NA
      }
      if (length(impact_results$pca$rr_roll[, , group]) > 
          0) {
        rr_roll_pca_plot <- ggplot(melt(as.data.frame(impact_results$pca$rr_roll[, 
                                                                                 , group]), id.vars = NULL), mapping = aes_string(x = rep(analysis$time_points[(length(analysis$time_points) - 
                                                                                                                                                                  nrow(impact_results$pca$rr_roll[, , group]) + 
                                                                                                                                                                  1):length(analysis$time_points)], ncol(impact_results$pca$rr_roll[, 
                                                                                                                                                                                                                                    , group])), y = "value", linetype = "variable")) + 
          geom_line() + geom_hline(yintercept = 1, linetype = 4) + 
          labs(x = "Time", y = "Rolling Rate Ratio") + 
          ggtitle(paste(group, "STL+PCA Rolling Rate Ratio")) + 
          coord_cartesian(ylim = min_max) + theme_bw() + 
          theme(axis.line = element_line(colour = "black"), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.border = element_blank(), panel.background = element_blank()) + 
          theme(legend.title = element_blank(), legend.position = c(0, 
                                                                    1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, 
                                                                                                                                         fill = "transparent"), plot.title = element_text(hjust = 0.5), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
      else {
        rr_roll_pca_plot <- NA
      }
    }  else {
      rr_roll_full_plot <- NA
      rr_roll_time_plot <- NA
      rr_roll_pca_plot <- NA
    }
    
    
    plots$groups[[group]] <- list(covar = covar_plot, 
                                  pred_full = pred_full_plot, 
                                  pred_time = pred_time_plot, pred_pca = pred_pca_plot, 
                                   pred_full_agg = pred_full_plot_agg, 
                                  pred_time_agg = pred_time_plot_agg, pred_pca_agg = pred_pca_plot_agg, 
                                  rr_roll_full = rr_roll_full_plot, 
                                  rr_roll_time = rr_roll_time_plot, rr_roll_pca = rr_roll_pca_plot, 
                                  cumsum_prevented = cumsum_prevented_plot[[group]])
    
  }
  
  
  return(plots)
}