syncon.plots <- function(analysis) {
  plots = list(groups = list())
  
  cbPalette  <-  c("#1b9e77", "#d95f02", "#7570b3", '#e7298a')
  #Compare rate ratios, with size of marker scaled to cross val weights
  if ("crossval" %in% names(analysis$results)) {
    point.weights = analysis$results$crossval$point.weights
  } else {
    point.weights = analysis$results$impact$point.weights
  }
  plots$summary <-
    ggplot(
      analysis$results$impact$rr_mean_combo,
      aes_(
        x =  ~ group.index,
        y =  ~ mean.rr,
        color =  ~ Model,
        group =  ~ Model
      )
    ) +
    geom_errorbar(aes_(ymin =  ~ lcl, ymax =  ~ ucl),
                  colour = "gray",
                  width = .0) +
    geom_point(aes_(
      shape =  ~ Model,
      size =  ~ analysis$results$impact$rr_mean_combo$est.index
    )) +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    scale_size_manual(values = c(point.weights$value * 2)) + #Scales area, which is optimal for bubbl plot
    #geom_errorbar(rr_mean_combo,aes_(ymin=~lcl, ymax=~ucl), colour="black", width=.1) +
    theme_bw() +
    guides(size = FALSE) + #turn off size axis
    scale_colour_manual(values = cbPalette) +
    labs(x = "Group", y = "Rate ratio") +
    geom_hline(yintercept = 1,
               colour = 'gray',
               linetype = 2) +
    theme(
      axis.line = element_line(colour = "black"),
      legend.position = c(0.2, 0.9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  
  for (group in analysis$groups) {
    #View scaled covariates
    covars.sub <-
      analysis$covars$full[[group]][, -c(1:(analysis$n_seasons - 1))]
    alpha1 = rep(
      analysis$results$impact$full$inclusion_probs[[group]][-c(1:(analysis$n_seasons -
                                                                    1)), 'inclusion_probs'],
      each = nrow(analysis$covars$full[[group]])
    )
    covar_plot <-
      ggplot(
        melt(covars.sub, id.vars = NULL),
        mapping = aes_string(
          x = rep(analysis$time_points, ncol(covars.sub)),
          y = 'value',
          group = 'variable',
          alpha = alpha1
        )
      ) +
      geom_line() +
      labs(x = 'Time', y = 'Scaled Covariates') +
      ggtitle(paste(group, 'Scaled Covariates Weighted by Inclusion Probability')) +
      theme_bw() +
      theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      ) +
      theme(
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    #Plot predictions
    min_max <-
      c(min(
        c(
          analysis$results$impact$full$pred_quantiles[, , group],
          analysis$outcome[, group]
        )
      ), max(
        c(
          analysis$results$impact$full$pred_quantiles[, , group],
          analysis$outcome[, group]
        )
      ))
    pred_best_plot <-
      plotPred(
        analysis$results$impact$best$pred_quantiles[, , group],
        analysis$time_points,
        analysis$post_period,
        min_max,
        analysis$outcome[, group],
        title = paste(group, 'Best estimate')
      )
    pred_full_plot <-
      plotPred(
        analysis$results$impact$full$pred_quantiles[, , group],
        analysis$time_points,
        analysis$post_period,
        min_max,
        analysis$outcome[, group],
        title = paste(group, 'Synthetic controls estimate')
      )
    pred_time_plot <-
      plotPred(
        analysis$results$impact$time$pred_quantiles[, , group],
        analysis$time_points,
        analysis$post_period,
        min_max,
        analysis$outcome[, group],
        title = paste(group, 'Interupted time series estimate')
      )
    pred_pca_plot <-
      plotPred(
        analysis$results$impact$pca$pred_quantiles[, , group],
        analysis$time_points,
        analysis$post_period,
        min_max,
        analysis$outcome[, group],
        title = paste(group, 'STL+PCA estimate')
      )
    if ("crossval" %in% names(analysis$results)) {
      pred_stack_plot <-
        plotPred(
          analysis$results$crossval$pred_quantiles_stack[, , group],
          analysis$time_points,
          analysis$post_period,
          min_max,
          analysis$outcome[, group],
          title = paste(group, 'Stacked estimate')
        )
      pred_stack_plot_agg <-
        plotPredAgg(
          analysis$results$crossval$ann_pred_quantiles_stack[[group]],
          analysis$time_points,
          analysis$year_def,
          analysis$intervention_date,
          analysis$post_period,
          min_max,
          analysis$outcome[, group],
          title = paste(group, 'Stacked estimate')
        )
      
    }
    pred_best_plot_agg <-
      plotPredAgg(
        analysis$results$impact$best$ann_pred_quantiles[[group]],
        analysis$time_points,
        analysis$year_def,
        analysis$intervention_date,
        analysis$post_period,
        min_max,
        analysis$outcome[, group],
        title = paste(group, 'Best estimate')
      )
    pred_full_plot_agg <-
      plotPredAgg(
        analysis$results$impact$full$ann_pred_quantiles[[group]],
        analysis$time_points,
        analysis$year_def,
        analysis$intervention_date,
        analysis$post_period,
        min_max,
        analysis$outcome[, group],
        title = paste(group, 'Synthetic controls estimate')
      )
    pred_time_plot_agg <-
      plotPredAgg(
        analysis$results$impact$time$ann_pred_quantiles[[group]],
        analysis$time_points,
        analysis$year_def,
        analysis$intervention_date,
        analysis$post_period,
        min_max,
        analysis$outcome[, group],
        title = paste(group, 'Interupted time series estimate')
      )
    pred_pca_plot_agg <-
      plotPredAgg(
        analysis$results$impact$pca$ann_pred_quantiles[[group]],
        analysis$time_points,
        analysis$year_def,
        analysis$intervention_date,
        analysis$post_period,
        min_max,
        analysis$outcome[, group],
        title = paste(group, 'STL+PCA estimate')
      )
    
    if ("sensitivity" %in% names(analysis$results)) {
      pred_sensitivity_plot <-
        plotPred(
          analysis$results$impact$full$pred_quantiles[, , group],
          analysis$time_points,
          analysis$post_period,
          min_max,
          analysis$outcome[, group],
          sensitivity_pred_quantiles = analysis$sensitivity_pred_quantiles[[group]],
          sensitivity_title = paste(group, 'Sensitivity Plots'),
          plot_sensitivity = TRUE
        )
    } else {
      pred_sensitivity_plot <- NA
    }
    #matplot(pred_quantiles_full[, , 10], ylim=c(0,22000), type='l')   ##Check
    #points(prelog_data[[10]]$J12_18)
    
    #Plot rolling rate ratio
    min_max <-
      c(
        min(
          analysis$results$impact$full$rr_roll[, , group],
          analysis$results$impact$time$rr_roll[, , group]
        ),
        max(
          analysis$results$impact$full$rr_roll[, , group],
          analysis$results$impact$time$rr_roll[, , group]
        )
      )
    rr_roll_best_plot <-
      ggplot(
        melt(
          as.data.frame(analysis$results$impact$best$rr_roll[, , group]),
          id.vars = NULL
        ),
        mapping = aes_string(
          x = rep(
            analysis$time_points[(length(analysis$time_points) - nrow(analysis$results$impact$best$rr_roll[, , group]) + 1):length(analysis$time_points)],
            ncol(analysis$results$impact$best$rr_roll[, , group])
          ),
          y = 'value',
          linetype = 'variable'
        )
      ) +
      geom_line() + geom_hline(yintercept = 1, linetype = 4) +
      labs(x = 'Time', y = 'Rolling Rate Ratio') +
      ggtitle(paste(group, 'Synthetic Control Rolling Rate Ratio')) +
      coord_cartesian(ylim = min_max) +
      theme_bw() +
      theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      ) +
      theme(
        legend.title = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(colour = NA, fill = 'transparent'),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    rr_roll_full_plot <-
      ggplot(
        melt(
          as.data.frame(analysis$results$impact$full$rr_roll[, , group]),
          id.vars = NULL
        ),
        mapping = aes_string(
          x = rep(
            analysis$time_points[(length(analysis$time_points) - nrow(analysis$results$impact$full$rr_roll[, , group]) + 1):length(analysis$time_points)],
            ncol(analysis$results$impact$full$rr_roll[, , group])
          ),
          y = 'value',
          linetype = 'variable'
        )
      ) +
      geom_line() + geom_hline(yintercept = 1, linetype = 4) +
      labs(x = 'Time', y = 'Rolling Rate Ratio') +
      ggtitle(paste(group, 'Synthetic Control Rolling Rate Ratio')) +
      coord_cartesian(ylim = min_max) +
      theme_bw() +
      theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      ) +
      theme(
        legend.title = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(colour = NA, fill = 'transparent'),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    rr_roll_time_plot <-
      ggplot(
        melt(
          as.data.frame(analysis$results$impact$time$rr_roll[, , group]),
          id.vars = NULL
        ),
        mapping = aes_string(
          x = rep(
            analysis$time_points[(length(analysis$time_points) - nrow(analysis$results$impact$time$rr_roll[, , group]) + 1):length(analysis$time_points)],
            ncol(analysis$results$impact$time$rr_roll[, , group])
          ),
          y = 'value',
          linetype = 'variable'
        )
      ) +
      geom_line() + geom_hline(yintercept = 1, linetype = 4) +
      labs(x = 'Time', y = 'Rolling Rate Ratio') +
      ggtitle(paste(group, 'TT Rolling Rate Ratio')) +
      coord_cartesian(ylim = min_max) +
      theme_bw() +
      theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      ) +
      theme(
        legend.title = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(colour = NA, fill = 'transparent'),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    rr_roll_pca_plot <-
      ggplot(
        melt(
          as.data.frame(analysis$results$impact$pca$rr_roll[, , group]),
          id.vars = NULL
        ),
        mapping = aes_string(
          x = rep(
            analysis$time_points[(length(analysis$time_points) - nrow(analysis$results$impact$pca$rr_roll[, , group]) + 1):length(analysis$time_points)],
            ncol(analysis$results$impact$pca$rr_roll[, , group])
          ),
          y = 'value',
          linetype = 'variable'
        )
      ) +
      geom_line() + geom_hline(yintercept = 1, linetype = 4) +
      labs(x = 'Time', y = 'Rolling Rate Ratio') +
      ggtitle(paste(group, 'STL+PCA Rolling Rate Ratio')) +
      coord_cartesian(ylim = min_max) +
      theme_bw() +
      theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      ) +
      theme(
        legend.title = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(colour = NA, fill = 'transparent'),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    if ("crossval" %in% names(analysis$results)) {
      rr_roll_stack_plot <-
        ggplot(
          melt(
            as.data.frame(analysis$results$crossval$rr_roll_stack[, , group]),
            id.vars = NULL
          ),
          mapping = aes_string(
            x = rep(
              analysis$time_points[(
                length(analysis$time_points) - nrow(analysis$results$crossval$rr_roll_stack[, , group]) + 1
              ):length(analysis$time_points)],
              ncol(analysis$results$crossval$rr_roll_stack[, , group])
            ),
            y = 'value',
            linetype = 'variable'
          )
        ) +
        geom_line() + geom_hline(yintercept = 1, linetype = 4) +
        labs(x = 'Time', y = 'Rolling Rate Ratio') +
        ggtitle(paste(group, 'Stacked Rolling Rate Ratio')) +
        coord_cartesian(ylim = min_max) +
        theme_bw() +
        theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
        ) +
        theme(
          legend.title = element_blank(),
          legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.background = element_rect(colour = NA, fill = 'transparent'),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
    }
    #Plot cumulative sums
    cumsum_prevented_plot <-
      ggplot(
        melt(
          as.data.frame(analysis$results$impact$best$cumsum_prevented[, , group]),
          id.vars = NULL
        ),
        mapping = aes_string(
          x = rep(
            analysis$time_points,
            ncol(analysis$results$impact$best$cumsum_prevented[, , group])
          ),
          y = 'value',
          linetype = 'variable'
        )
      ) +
      geom_line() + geom_hline(yintercept = 1, linetype = 4) +
      labs(x = 'Time', y = 'Cumulative Sum Prevented') +
      ggtitle(paste(group, 'Cumulative Number of Cases Prevented')) +
      theme_bw() +
      theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      ) +
      theme(
        legend.title = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(colour = NA, fill = 'transparent'),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    if ("crossval" %in% names(analysis$results)) {
      cumsum_prevented_stack_plot <-
        ggplot(
          melt(
            as.data.frame(analysis$results$crossval$cumsum_prevented_stack[, , group]),
            id.vars = NULL
          ),
          mapping = aes_string(
            x = rep(
              analysis$time_points,
              ncol(analysis$results$crossval$cumsum_prevented_stack[, , group])
            ),
            y = 'value',
            linetype = 'variable'
          )
        ) +
        geom_line() + geom_hline(yintercept = 1, linetype = 4) +
        labs(x = 'Time', y = 'Cumulative Sum Prevented') +
        ggtitle(paste(
          group,
          'Cumulative Number of Cases Prevented (Stacked model)'
        )) +
        theme_bw() +
        theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
        ) +
        theme(
          legend.title = element_blank(),
          legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.background = element_rect(colour = NA, fill = 'transparent'),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
    }
    
    
    if ("crossval" %in% names(analysis$results)) {
      plots$groups[[group]] <- list(
        covar = covar_plot,
        pred_stack = pred_stack_plot,
        pred_full = pred_full_plot,
        pred_time = pred_time_plot,
        pred_pca = pred_pca_plot,
        pred_stack_agg = pred_stack_plot_agg,
        pred_full_agg = pred_full_plot_agg,
        pred_time_agg = pred_time_plot_agg,
        pred_pca_agg = pred_pca_plot_agg,
        #pred_sensitivity = pred_sensitivity_plot,
        rr_roll_stack = rr_roll_stack_plot,
        rr_roll_full = rr_roll_full_plot,
        rr_roll_time = rr_roll_time_plot,
        rr_roll_pca = rr_roll_pca_plot,
        cumsum_prevented_stack = cumsum_prevented_stack_plot,
        cumsum_prevented = cumsum_prevented_plot
      )
    } else {
      plots$groups[[group]] <- list(
        covar = covar_plot,
        pred_best = pred_best_plot,
        pred_full = pred_full_plot,
        pred_time = pred_time_plot,
        pred_pca = pred_pca_plot,
        pred_best_agg = pred_best_plot_agg,
        pred_full_agg = pred_full_plot_agg,
        pred_time_agg = pred_time_plot_agg,
        pred_pca_agg = pred_pca_plot_agg,
        #  pred_sensitivity = pred_sensitivity_plot,
        rr_roll_best = rr_roll_best_plot,
        rr_roll_full = rr_roll_full_plot,
        rr_roll_time = rr_roll_time_plot,
        rr_roll_pca = rr_roll_pca_plot,
        cumsum_prevented = cumsum_prevented_plot
      )
    }
  }
  
  return(plots)
}
