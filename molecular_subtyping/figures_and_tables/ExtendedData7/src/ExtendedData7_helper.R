### Helpers for Extended Data 7 ###

# Aggregate enrichment p-values - gsea
enrichment_aggregation_subtype_gsea = 
  function(factor_enrichment, adj.p.thr = 0.05,
           negation = 1, thr = .2, type = "gsea"){
    tmp0 =list()
    if (type == "gsea") {
      for(database_id in 1:length(factor_enrichment)){
        tmp0[[database_id]] <- try(enrichment_aggregation1(
          pos_list = factor_enrichment[[database_id]],
          thr = thr, type = type))
      }
    }
    
    pos_value=negation
    database_names = names(factor_enrichment)
    tmp = list()
    for(database_id in 1:length(factor_enrichment)){
      if(class(tmp0[[database_id]])!="try-error"){
        aa =  tmp0[[database_id]]$direction_mat[!is.na(tmp0[[database_id]]$direction_mat)]
        idx = which(!is.na(tmp0[[database_id]]$direction_mat))
        ww = as.matrix(tmp0[[database_id]]$direction_mat)
        if (type == "gsea") {
          ww[idx] = ifelse(aa > 0, pos_value, -pos_value)
        }
        tmp0[[database_id]]$direction_mat[,] = ww 
        tmp0[[database_id]]$direction_mat = apply(tmp0[[database_id]]$direction_mat,
                                                  c(1,2), as.numeric)
        tmp0[[database_id]]$pval_mat = -log(tmp0[[database_id]]$pval_mat, base = 10)
        tmp0[[database_id]]$pval_mat[,!(colnames(tmp0[[database_id]]$pval_mat) %in% 
                                          c("joint", "joint_adj"))] = 
          tmp0[[database_id]]$pval_mat[,!(colnames(tmp0[[database_id]]$pval_mat) %in%
                                            c("joint", "joint_adj"))]*
          sign(tmp0[[database_id]]$direction_mat)
        
        tmp[[database_id]] = data.frame(matrix(NA, 
                                               nrow = nrow(tmp0[[database_id]]$pval_mat), 
                                               ncol = length(omic_short_names)+2))
        colnames(tmp[[database_id]]) = c(omic_short_names, c("joint", "joint_adj"))
        tmp[[database_id]][,colnames(tmp0[[database_id]]$pval_mat)] = 
          tmp0[[database_id]]$pval_mat
        tmp[[database_id]]$pathway = rownames(tmp0[[database_id]]$pval_mat)
        rownames( tmp[[database_id]]) = rownames(tmp0[[database_id]]$pval_mat)
        tmp[[database_id]]$database = database_names[database_id]
      }
    }
    aggregated_mat = NULL
    for(database_id in 1:length(tmp)){
      aggregated_mat=rbind( aggregated_mat, tmp[[database_id]])
    }
    aggregated_mat =  aggregated_mat[aggregated_mat$joint_adj>abs(log(adj.p.thr,base=10)),]
    aggregated_mat=aggregated_mat[order(-aggregated_mat$joint),]
    return(aggregated_mat)
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
           feature_with_xlab = NULL, multiple_correction = TRUE) {
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
    if (multiple_correction) {
      smooth_model_use <- smooth_model_loop1_adj
    } else {
      smooth_model_use <- smooth_model_loop1
    }
    test_fig_list <- trajectory_plot(inputDF = inputDF, 
                                     model_DF = smooth_model_use,
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
      annotate("text", x = max_date/3,  # x = max_date/3
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

# Two-group comparison controlling basic demographics, admission date,
# participant ID, and enrollment site
two_pasc_mixed_test <- function(feature_name, data_mat, visit_numbers) {
  est = array(NA, length(feature_name))
  names(est) = feature_name
  pval = array(NA, length(feature_name))
  names(pval) = feature_name
  
  for(j in 1:length(feature_name)) {
    fname = feature_name[j]
    tdata = data_mat[data_mat$event_type %in% visit_numbers, ]
    tdata = tdata[!is.na(tdata[[fname]]),]
    pasc = factor(tdata$PASC)
    
    df_gamm = data.frame(fname = tdata[[fname]],
                         event_date = tdata$event_date,
                         participant_id = tdata$participant_id,
                         enrollment_site = tdata$enrollment_site,
                         pasc = tdata$PASC,
                         age = tdata$discretized_admit_age_quantile,
                         sex = tdata$sex)
    form = as.formula(paste0("fname ~ s(event_date) + age + sex + pasc"))
    fit = gamm4::gamm4(formula = form, data = df_gamm,
                       random = ~(1|enrollment_site/participant_id))
    coef_tab = summary(fit$gam)$p.table["pascOther Deficits",]
    est[j] = coef_tab[1]
    pval[j] = coef_tab[4]
  }
  return(list(est = est, pval = pval))
}

# Heatmap of hallmark pathway correlation with Factor 30
heatmap_hallmark_cor <- function(hallmark_cor, hallmark_pval, hallmark_feature) {
  df_cor <- as.data.frame(hallmark_cor)
  df_cor <- pivot_longer(df_cor, cols = everything())
  pval_cor <- as.data.frame(hallmark_pval)
  pval_cor <- pivot_longer(pval_cor, cols = everything())
  df_cor$pval <- pval_cor$value
  df_cor$feature <- rep(hallmark_feature, each = 2)
  df_cor$feature <- factor(df_cor$feature)
  df_cor$feature <- factor(df_cor$feature, levels = rev(hallmark_feature))
  
  df_cor %>% ggplot(mapping = aes(x = name, y = feature, fill = value)) +
    geom_tile(linewidth = 0.2) +
    scale_fill_gradient2(low = "#284AA0", mid = "white", 
                         high = "#BC4A4A", midpoint = 0,
                         name = expression("Correlation with Factor 30"),
                         limits = c(-1, 1), breaks = c(-1, 0, 1)) +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    ylab("") +
    xlab("") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15, hjust = 0),
          axis.text.x = element_text(size = 15, , vjust = 0.5, angle = 90),
          legend.position = "left",
          legend.title = element_text(size = 15, angle = 90),
          legend.text = element_text(size = 15, angle = 90, hjust = 0.5),
          plot.margin = margin(b = 50, t = 10)) +
    geom_text(aes(label = ifelse(pval <= 0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 4, vjust = 0.76, color = "black")
}