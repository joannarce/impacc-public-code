### Helpers for Extended Data 4 ###

# Child population preprocessing (no granulo)
preprocess_bld_cytof_child_wo_granulo <- function(counts) {
  bld_cytof_child_nongran <- counts %>%
    rownames_to_column(var = "event_id") %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Tier1_Undefined`,
                  -`Tier2_undefined_Undefined`,
                  -`Neutrophil (CD16low)`,
                  -`Neutrophil (CD16hi)`,
                  -`Eosinophil`,
                  -`Basophil`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_child_nongran <- bld_cytof_child_nongran/rowSums(bld_cytof_child_nongran)
  bld_cytof_child_nongran <- log1p(bld_cytof_child_nongran)
  
  bld_cytof_child_nongran <- bld_cytof_child_nongran %>%
    apply(2, scale) %>%
    data.frame(check.names = F) %>%
    `rownames<-`(row.names(bld_cytof_child_nongran))
  
  return(value = bld_cytof_child_nongran)
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

# Child population preprocessing (granulo)
preprocess_bld_cytof_child_seperate_granulo <- function(counts) {
  bld_cytof_child_counts_processed <- preprocess_bld_cytof_child_wo_granulo(counts)
  
  bld_cytof_child_gran <- counts %>%
    rownames_to_column(var = "event_id") %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Tier1_Undefined`,
                  -`Tier2_undefined_Undefined`#,
                  #-`Neutrophil (CD16low)`,
                  #-`Neutrophil (CD16hi)`,
                  #-`Eosinophil`,
                  #-`Basophil`
    ) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_child_gran <- bld_cytof_child_gran/rowSums(bld_cytof_child_gran)
  bld_cytof_child_gran <- log1p(bld_cytof_child_gran)
  
  bld_cytof_child_gran <- bld_cytof_child_gran %>%
    apply(2, scale) %>%
    data.frame(check.names = F) %>%
    `rownames<-`(row.names(bld_cytof_child_gran))
  
  final <- cbind.data.frame(bld_cytof_child_counts_processed,
                            bld_cytof_child_gran[,!colnames(bld_cytof_child_gran) %in% colnames(bld_cytof_child_counts_processed)])
  return(final)
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

# Quantile normalization
quantileNormalizeNoTies <- function(data) {
  # Calculate ranks for each column
  ranks <- apply(data, 2, rank, ties.method = "average")
  ranks_norm = ranks/(apply(ranks,2,max)+1/2)
  normalizedData = qnorm(ranks_norm)
  return(normalizedData)
}

# Muti-group comparison controlling basic demographics, admission date,
# participant ID, and enrollment site
multi_group_mixed_test <- function(feature_name, data_mat, visit_numbers,
                                   subtype_remove = "D") {
  if (!"label" %in% names(data_mat)) {
    stop("Error: 'label'(subtype) column not found in 'data_mat'")
  }
  if (!is.null(subtype_remove)) {
    data_mat = data_mat[data_mat$label != subtype_remove, ]
  }
  
  est = matrix(NA, nrow = 5, ncol = length(feature_name))
  rownames(est) = c("EF vs ABC", "F vs E", "A vs BC", "B vs AC", "C vs AB")
  colnames(est) = feature_name
  pval = matrix(NA, nrow = 5, ncol = length(feature_name))
  rownames(pval) = c("EF vs ABC", "F vs E", "A vs BC", "B vs AC", "C vs AB")
  colnames(pval) = feature_name
  
  for(j in 1:length(feature_name)) {
    fname = feature_name[j]
    tdata = data_mat[data_mat$event_type %in% visit_numbers, ]
    tdata = tdata[!is.na(tdata[[fname]]),]
    subtype = factor(tdata$label)
    
    x1 = model.matrix(~subtype)
    colnames(x1) = c("Intercept", "EFvsABC", "FvsE", "AvsBC", "BvsAC")
    x1[, "EFvsABC"] = -1
    x1[subtype == "E" | subtype == "F", "EFvsABC"] = 1
    x1[, "FvsE"] = 0
    x1[subtype == "E", "FvsE"] = -1
    x1[subtype == "F", "FvsE"] = 1
    x1[, "AvsBC"] = 0
    x1[subtype == "A", "AvsBC"] = 1
    x1[subtype == "C", "AvsBC"] = -1
    x1[, "BvsAC"] = 0
    x1[subtype == "B", "BvsAC"] = 1
    x1[subtype == "C", "BvsAC"] = -1
    
    x2 = model.matrix(~subtype)
    colnames(x2) = c("Intercept", "EFvsABC", "FvsE", "BvsAC", "CvsAB")
    x2[, "EFvsABC"] = -1
    x2[subtype == "E" | subtype == "F", "EFvsABC"] = 1
    x2[, "FvsE"] = 0
    x2[subtype == "E", "FvsE"] = -1
    x2[subtype == "F", "FvsE"] = 1
    x2[, "BvsAC"] = 0
    x2[subtype == "B", "BvsAC"] = 1
    x2[subtype == "A", "BvsAC"] = -1
    x2[, "CvsAB"] = 0
    x2[subtype == "C", "CvsAB"] = 1
    x2[subtype == "A", "CvsAB"] = -1
    
    df_gamm1 = cbind(data.frame(fname = tdata[[fname]],
                                event_date = tdata$event_date,
                                participant_id = tdata$participant_id,
                                enrollment_site = tdata$enrollment_site,
                                age = tdata$discretized_admit_age_quantile,
                                sex = tdata$sex), x1)
    form1 = as.formula(
      paste0("fname ~ s(event_date) + age + sex + EFvsABC + FvsE + AvsBC + BvsAC"))
    fit1 = gamm4::gamm4(formula = form1, data = df_gamm1,
                        random =~(1|enrollment_site/participant_id))
    coef_tab = summary(fit1$gam)$p.table[c("EFvsABC", "FvsE", "AvsBC", "BvsAC"),]
    est[1:4,j] = coef_tab[,1]
    pval[1:4,j] = coef_tab[,4]
    
    df_gamm2 = cbind(data.frame(fname = tdata[[fname]],
                                event_date = tdata$event_date,
                                participant_id = tdata$participant_id,
                                enrollment_site = tdata$enrollment_site,
                                age = tdata$discretized_admit_age_quantile,
                                sex = tdata$sex), x2)
    form2 = as.formula(
      paste0("fname ~ s(event_date) + age + sex + EFvsABC + FvsE + BvsAC + CvsAB"))
    fit2 = gamm4::gamm4(formula = form2, data = df_gamm2,
                        random = ~(1|enrollment_site/participant_id))
    coef_tab = summary(fit2$gam)$p.table[c("EFvsABC", "FvsE", "BvsAC", "CvsAB"),]
    est[5,j] = coef_tab["CvsAB",1]
    pval[5,j] = coef_tab["CvsAB",4]
  }
  return(list(est = est, pval = pval))
}

### x1
# S <- matrix(c(1, -1, 0, 1, 0,
#               1, -1, 0, 0, 1,
#               1, -1, 0, -1, -1,
#               1, 1, -1, 0, 0,
#               1, 1, 1, 0, 0), byrow = T, nrow = 5, ncol = 5)
# S
# L <- ginv(S)
# L %>% fractions()
# solve(t(S)%*%S)%*%t(S) %>% fractions()

### x2
# S <- matrix(c(1, -1, 0, -1, -1,
#               1, -1, 0, 1, 0,
#               1, -1, 0, 0, 1,
#               1, 1, -1, 0, 0,
#               1, 1, 1, 0, 0), byrow = T, nrow = 5, ncol = 5)
# S
# L <- ginv(S)
# L %>% fractions()
# solve(t(S)%*%S)%*%t(S) %>% fractions()

# two group mixed test
two_group_mixed_test = function(fname_vec, data_mat, visit_numbers,
                                grp_col, grp0, grp1){
  res_tab = data.frame(feature = fname_vec)
  res_tab$coef = NA; res_tab$pval = NA
  for(j in 1:length(fname_vec)){
    fname = fname_vec[j]
    tdata = data.frame(fname = data_mat[[fname]],
                       event_date = data_mat$event_date,
                       event_id = data_mat$event_id,
                       label = data_mat[[grp_col]],
                       event_type = data_mat$event_type, 
                       participant_id = data_mat$participant_id,
                       age = data_mat$discretized_admit_age_quantile,
                       sex = data_mat$sex)
    tdata = tdata[tdata$event_type %in% visit_numbers & 
                    ((tdata$label %in% grp0)|(tdata$label %in% grp1)),]
    tdata = tdata[!is.na(tdata$fname),]
    tdata$grp = ifelse(tdata$label %in% grp1, 1, 0)
    formula1 = as.formula("fname ~ s(event_date) + grp + age + sex")
    fitted1 = gamm4::gamm4(formula = formula1, data = tdata, 
                           random = ~(1|participant_id))
    tab1 = summary(fitted1$gam)
    # tab1$p.coeff
    res_tab$coef[j] = tab1$p.coeff[2]
    res_tab$pval[j] = tab1$p.pv[2]
  }
  return(res_tab)
}

