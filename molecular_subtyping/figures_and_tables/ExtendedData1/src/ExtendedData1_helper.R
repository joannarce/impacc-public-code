### Helpers for Extended Data 1 ###

# Ordinal regression adjusted for demographics (two-group comparison)
test_ordinal_demo_two_group <- function(subtype, response, demo_data, 
                                        base_group, test_group) {
  id <- (!is.na(subtype)) & (!is.na(response)) & 
    (rowSums(is.na(demo_data)) == 0) &
    ((subtype == base_group) | (subtype == test_group))
  subtype <- factor(subtype[id])
  response <- factor(response[id])
  demo_data <- demo_data[id,]
  nclass <- length(unique(subtype))
  nlevel_resp <- nlevels(response)
  
  subtype <- relevel(subtype, ref = base_group)
  df_ordinal <- cbind(data.frame(response = response), subtype, demo_data)
  model1 <- ordinal::clm(response ~ ., data = df_ordinal)
  coef_ordinal <- head(summary(model1)$coefficients, 
                       ((nlevel_resp-1)+(nclass-1)))
  coef_ordinal <- tail(coef_ordinal, 1)
  rownames(coef_ordinal) <- paste0(test_group, " vs ", base_group)
  est_ordinal <- coef_ordinal[,1]
  p_ordinal <- coef_ordinal[,4]
  return(coef_ordinal)
}

# Heatmap for ordinal regression before merging
heatmap_ordinal_merge <- function(test_ord) {
  comparison_name <- rownames(test_ord)
  df_ord <- as.data.frame(test_ord[,1])
  colnames(df_ord) <- c("TG")
  df_ord <- pivot_longer(df_ord, cols = everything())
  
  pval_ord <- as.data.frame(test_ord[,4])
  colnames(pval_ord) <- c("TG")
  pval_ord <- pivot_longer(pval_ord, cols = everything())
  df_ord$pval <- pval_ord$value
  
  df_ord$subtype <- comparison_name
  df_ord$subtype <- factor(df_ord$subtype)
  df_ord$subtype <- factor(df_ord$subtype,
                           levels = rev(levels(df_ord$subtype)))
  
  df_ord %>% ggplot(mapping = aes(x = name, y = subtype, color = value)) +
    geom_point(size=10) +
    geom_point(aes(x = name, y = subtype), shape = 1, 
               size = 10, color = "grey30") +
    scale_color_gradient2(low = "#284AA0", mid="white", 
                          high = "#BC4A4A", midpoint = 0,
                          name = expression("Coefficient"),
                          limits = c(-max(abs(df_ord$value)),
                                     max(abs(df_ord$value)))) +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    ylab("") +
    xlab("") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15, hjust = 0),
          axis.text.x = element_text(size = 15, angle = 90),
          legend.position = "left",
          legend.title = element_text(size = 15, angle = 90),
          legend.text = element_text(size = 15, angle = 90, hjust = 0.5),
          legend.key.size = unit(2, "lines"),
          plot.margin = margin(b = 46, t = 10),
          legend.spacing.y = unit(1.5, 'cm')) +
    geom_text(aes(label = ifelse(pval <= 0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 7, vjust = 0.8, color = "black")
}