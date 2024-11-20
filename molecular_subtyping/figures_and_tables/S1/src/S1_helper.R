### Helpers for S1 ###

# t-test at each split of cluster dendrogram
split_t_test <- function(branch_left_id, branch_right_id, n_factor, visit_max) {
  p_value <- matrix(NA, nrow = n_factor, ncol = visit_max)
  rownames(p_value) <- paste0("Factor", 1:n_factor)
  colnames(p_value) <- paste0("Visit", 1:visit_max)
  t_score <- matrix(NA, nrow = n_factor, ncol = visit_max)
  rownames(t_score) <- paste0("Factor", 1:n_factor)
  colnames(t_score) <- paste0("Visit", 1:visit_max)
  effect_size <- matrix(NA, nrow = n_factor, ncol = visit_max)
  rownames(effect_size) <- paste0("Factor", 1:n_factor)
  colnames(effect_size) <- paste0("Visit", 1:visit_max)
  n_sample <- matrix(NA, nrow = n_factor, ncol = visit_max)
  branch_left_idx = which(rownames(X) %in% branch_left_id)
  branch_right_idx = which(rownames(X) %in% branch_right_id)
  min_sample_size = 5
  for (i in 1:n_factor) {
    col_select <- seq(i, by = n_factor, length.out = visit_max)
    left_factor <- X[branch_left_idx, col_select]
    right_factor <- X[branch_right_idx, col_select]
    for (j in 1:visit_max) {
      if (sum(!is.na(left_factor[, j])) >= min_sample_size && 
          sum(!is.na(right_factor[, j])) >= min_sample_size) {
        t_test <- t.test(left_factor[, j], right_factor[, j])
        # remove sample size effect
        n1 <- sum(!is.na(left_factor[, j]))
        n2 <- sum(!is.na(right_factor[, j]))
        n_effect <- n1*n2/(n1+n2)
        p_value[i, j] <- t_test$p.value
        t_score[i, j] <- t_test$statistic
        effect_size[i, j] <- t_test$statistic/sqrt(n_effect)
        n_sample[i, j] <- sum(!is.na(left_factor[, j])) + 
          sum(!is.na(right_factor[, j]))
      }
    }
  }
  p_adjust = apply(p_value, 2, function(p) p.adjust(p, method = "BH"))
  out = list(p_value = p_value, p_adjust = p_adjust, t_score = t_score, 
             effect_size = effect_size, n_sample = n_sample)
  return(out)
}

# Merge subtype label and clinical data
cytof_mat_preprocess <- function(cytof_mat, subtype_hcl, clin_data) {
  cytof_mat <- cytof_mat[order(cytof_mat$event_id), ]
  cytof_mat <- left_join(cytof_mat, 
                         clin_data[,c("event_id", "event_type", 
                                      "enrollment_site",
                                      "sex", "discretized_admit_age_quantile",
                                      "participant_id", "event_date", 
                                      "symptom_date", "PASC_group")],  
                         by = "event_id")
  cytof_mat <- cytof_mat[!is.na(cytof_mat$event_type), ]
  cytof_mat$event_type <- gsub("Visit ", "V", cytof_mat$event_type)
  subtype_cytof <- subtype_hcl
  subtype_cytof$participant_id <- rownames(subtype_cytof)
  cytof_mat <- left_join(cytof_mat, subtype_cytof, by = "participant_id")
  cytof_mat <- cytof_mat[!is.na(cytof_mat$label), ]
  return(cytof_mat)
}

# Wrapper for plotting smoothing spline
helper_trajectory_data_smooth_spline <- 
  function(assay_mat, clin_colnames_select, endpoint = "label", 
           use_endpoint_levels, CI_interval = 0.95, max_date = 30, 
           subgroup_cr = c("A" = "#2ca02c", "B" = "#1f77b4", "C" = "#ff7f0e",
                           "D" = "#8c564b", "E" = "#9C3418", "F" = "#9467bd"), 
           legend_name = "Subtype", feature_select = FALSE, 
           limits = c(-Inf, Inf), feature_with_ylab = NULL,
           feature_with_xlab = NULL) {
    inputDF <- assay_mat
    inputDF <- inputDF %>%
      gather(key = "name", value = "value", -clin_colnames_select)
    inputDF$PASC_group <- factor(inputDF$PASC_group, 
                                 levels = c("minimal", "physical", 
                                            "cognitive", "multiple"))
    inputDF$name <- as.factor(inputDF$name)
    inputDF$sex <- factor(inputDF$sex, levels = c("Female", "Male"))
    inputDF$discretized_admit_age_quantile <- 
      as.factor(inputDF$discretized_admit_age_quantile)
    inputDF$participant_id <- as.factor(inputDF$participant_id)
    inputDF$enrollment_site <- as.factor(inputDF$enrollment_site)
    inputDF <- inputDF[!is.na(inputDF$value) & !is.na(inputDF[[endpoint]]),]
    inputDF <- inputDF[inputDF[[endpoint]] %in% use_endpoint_levels,]
    inputDF[[endpoint]] <- factor(inputDF[[endpoint]])
    inputDF$value <- as.numeric(inputDF$value)
    
    smooth_model_loop <- model_loop(inputDF, age_sex = T,
                                    modelType = "smoothSpline",
                                    endpoint = endpoint,
                                    old_p_corrections = TRUE)
    slope_select <- apply(combn(use_endpoint_levels, 2), 2,
                          function(x) paste0(x[1], "v", x[2]))
    slope_select <- paste0("p.slope_", slope_select)
    smooth_model_loop1 <- smooth_model_loop[, c("p.slope", "p.intercept",
                                                slope_select)]
    smooth_model_loop1_adj <- apply(smooth_model_loop1, 2,
                                    function(p) p.adjust(p, method = "BH"))
    pthr <- 0.001; particular_interest_features <- c()
    if (feature_select) {
      smooth_model_loop1 <- smooth_model_loop1[unique(c(
        particular_interest_features,
        which(apply(smooth_model_loop1_adj[,],1,min)<pthr))),]
      fnames <- rownames(smooth_model_loop1)
    } else {
      fnames <- unique(inputDF$name)
    }
    particular_interest_features <- fnames
    test_pval_list <- list(pval = c(), adj_pval = c())
    test_fig_list <- list()
    test_pval_list$pval <- smooth_model_loop1[fnames,]
    test_pval_list$adj_pval <- smooth_model_loop1_adj[fnames,]
    test_fig_list <- trajectory_plot(inputDF = inputDF, 
                                     model_DF = smooth_model_loop1_adj,
                                     endpoint = endpoint, 
                                     feature_names = fnames,
                                     legend_name = legend_name,
                                     legend_label = use_endpoint_levels,
                                     endpoint_order = NULL, 
                                     max_date = max_date,
                                     colors = subgroup_cr,
                                     limits = limits,
                                     feature_with_ylab = feature_with_ylab,
                                     feature_with_xlab = feature_with_xlab)
    return(list(test_pval_list = test_pval_list, test_fig_list = test_fig_list,
                particular_interest_features = particular_interest_features))
  }


# Longitudinal trajectory plot
trajectory_plot <- function(inputDF, model_DF, endpoint, feature_names,
                            legend_name, legend_label, endpoint_order = NULL,
                            max_date = 28, xlabel = "Days from admission",
                            CI_interval = 0.95, group_trendline_dropout = TRUE,
                            colors = c("#639A21", "#39828C",
                                       "#6371AD", "#BD7D31", "#9C3418"),
                            limits = c(-Inf, Inf), feature_with_ylab = NULL,
                            feature_with_xlab = NULL) {
  plotDF <- inputDF
  plotDF <- plotDF[!is.na(plotDF[[endpoint]]),]
  if(!is.null( endpoint_order)){
    plotDF[[endpoint]] <- factor(plotDF[[endpoint]], levels = endpoint_order)
  }
  plotDF[[endpoint]] <- ordered(plotDF[[endpoint]])
  plot_list <- vector("list", length(feature_names))
  names(plot_list) <- feature_names
  for(feature in names(plot_list)){
    plotExample <- plotDF[plotDF$name == feature,]
    
    # smoothSpline from plot_model in data_analysis_template_codebase.R
    formula_use <- formula(paste0("value~ s(event_date, bs = 'cr')+
              s(event_date, bs = 'cr', by =", endpoint, ")+", endpoint,
                                  "+ sex + discretized_admit_age_quantile"))
    fit <- gamm4::gamm4(formula_use, data = plotExample, 
                        random = ~(1|enrollment_site/participant_id))
    plotExample$yhat <- predict(fit$mer)
    pred <- ggpredict(fit, c("event_date", endpoint), type = "fixed",
                      ci.lvl = CI_interval)
    colnames(pred) <- gsub("x", "event_date", 
                           gsub("group", endpoint, colnames(pred)))
    # group_trendline dropout
    if(group_trendline_dropout) {
      pred <- pred %>%
        left_join(plotExample %>%
                    mutate({{endpoint}} := as.character(.data[[endpoint]])) %>%
                    group_by(.data[[endpoint]]) %>% 
                    dplyr::summarize(max = max(event_date)), 
                  by = c({{endpoint}})) %>%
        filter(event_date <= max)
    }
    
    plotExample <- plotExample[plotExample$value >= limits[1] &
                                 plotExample$value <= limits[2],]
    plot_list[[feature]] <- ggplot(data = plotExample, 
                                   mapping = aes(x = event_date, y = value)) +
      geom_jitter(mapping = aes_string(color = endpoint),
                  width = 0, height = 0, size = 1.0, alpha = 0.2) +
      geom_line(inherit.aes = F, data = pred, linewidth = 1.1,
                mapping = aes_string(group = endpoint, y = "predicted", 
                                     x = "event_date", color = endpoint)) +
      geom_ribbon(inherit.aes = F, data = pred, alpha = 0.2,
                  aes_string(x = "event_date", ymin = "conf.low", 
                             ymax = "conf.high", fill = endpoint)) +
      scale_fill_manual(values = colors, name = legend_name,
                        labels = legend_label) +
      scale_color_manual(values = colors, name = legend_name,
                         labels = legend_label) + 
      xlim(c(0, max_date)) + # max_date+1
      ylim(c(min(plotExample$value), max(plotExample$value))) +
      theme_minimal() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.x = element_line(linewidth = 1),
            axis.line.y = element_line(linewidth = 1),
            axis.ticks = element_line(linewidth = 1),
            axis.text.x = element_text(color = "black", size = 15),  
            axis.text.y = element_text(color = "black", size = 15),  
            axis.title = element_text(color = "black", size = 15),
            axis.ticks.length = unit(0.15, "cm"),
            plot.subtitle = element_text(color = "black", hjust = 0.5, 
                                         size = 15),
            legend.text = element_text(color = "black", size = 15),
            legend.title = element_text(color = "black", size = 15),
            legend.position = "bottom")
    
    if (is.null(feature_with_ylab) || feature %in% feature_with_ylab) {
      plot_list[[feature]] <- plot_list[[feature]] +
        ylab(paste0("normalized value")) +
        theme(plot.margin = margin(r=-5))
    } else {
      plot_list[[feature]] <- plot_list[[feature]] +
        ylab("") +
        theme(plot.margin = margin(l=-5))
    }
    
    if (is.null(feature_with_xlab) || feature %in% feature_with_xlab) {
      plot_list[[feature]] <- plot_list[[feature]] +
        xlab(xlabel)
    } else {
      plot_list[[feature]] <- plot_list[[feature]] +
        xlab("")
    }
    
    p_shape <- model_DF[feature,"p.slope"]
    p_average <- model_DF[feature,"p.intercept"]
    
    plot_list[[feature]] <- plot_list[[feature]] +
      labs(subtitle = feature) +
      annotate("text", x = 10,  # x = max_date/3
               y = max(plotExample$value, na.rm = TRUE),
               vjust = 0.75, hjust = 0, size = 5,
               label = paste0("shape p.val = ",
                              ifelse(p_shape < 0.001,
                                     format(p_shape, scientific = TRUE, 
                                            digits = 2),
                                     format(p_shape, digits = 2)), 
                              "\naverage p.val = ",
                              ifelse(p_average < 0.001,
                                     format(p_average, scientific = TRUE, 
                                            digits = 2),
                                     format(p_average, digits = 2))))
  }
  return(plot_list)
}

# Pie chart of top factors
pie_plot <- function(split_t, idx_factor, capped_effect_size = 1) {
  combined_data <- data.frame()
  for (idx in idx_factor) {
    adj_p = split_t$p_adjust[idx, ]
    t_score = split_t$t_score[idx, ]
    effect_size = split_t$effect_size[idx, ]
    n_sample = split_t$n_sample[idx, ]
    effect_size[effect_size > capped_effect_size] <- capped_effect_size
    effect_size[effect_size < -capped_effect_size] <- -capped_effect_size
    data_heatmap <- data.frame(
      effect_size = effect_size,
      adj_p = adj_p,
      n_sample = n_sample,
      visit = paste0("V", 1:visit_max)
    )
    combined_data <- rbind(combined_data, data_heatmap)
  }
  combined_data$factor <- factor(rep(idx_factor, each = visit_max))
  levels(combined_data$factor) = paste("Factor", idx_factor)
  
  ggplot(combined_data, aes(x = "", y = n_sample, fill = effect_size)) +
    geom_bar(stat = "identity") +
    coord_polar(theta = "y") +
    facet_wrap(~factor, ncol = 3) +
    geom_text(aes(x = 1.25, label = visit),
              size = 2.5, position = position_stack(vjust = 0.5)) +
    geom_text(aes(x = 0.9, 
                  label = ifelse(adj_p <= 0.001, "***",
                                 ifelse(adj_p <= 0.01, "**",
                                        ifelse(adj_p <= 0.05, "*", "")))),
              size = 3, position = position_stack(vjust = 0.42)) +
    scale_fill_gradient2(low = "#284AA0", mid = "white", high = "#BC4A4A", 
                         midpoint = 0, 
                         limits = c(-capped_effect_size, capped_effect_size)) +
    labs(fill = expression("Effect Size")) +
    theme_void() +
    theme(strip.text = element_text(size = 15))
}