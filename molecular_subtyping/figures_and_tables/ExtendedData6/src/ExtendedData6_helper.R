### Helpers for Figure 6 ###

# Fisher's exact test
test_clinical_fisher_exact = function(subtype, response, 
                                      response_binary, seed = 123){
  nclass = sort(unique(subtype))
  nresponse = unique(response)
  fisher_table_base = matrix(NA, nrow = length(nclass), 
                             ncol = length(nresponse))
  odds_ratio_base = matrix(NA, nrow = length(nclass), 
                           ncol = length(nresponse))
  fisher_table_rest = matrix(NA, nrow = length(nclass), 
                             ncol = length(nresponse))
  odds_ratio_rest = matrix(NA, nrow = length(nclass), 
                           ncol = length(nresponse))
  fisher_table_row = matrix(NA, nrow = length(nclass), ncol = 1)
  fisher_table_min = matrix(NA, nrow = length(nclass), ncol = 1)
  odds_ratio_min = matrix(NA, nrow = length(nclass), ncol = 1)
  tt1 = table(factor(subtype), response)
  set.seed(seed)
  global_pvalue_multiresponse = fisher.test(tt1, simulate.p.value = T, 
                                            B = 2E5)$p.value
  tt1 = table(factor(subtype), factor(response_binary))
  set.seed(seed)
  global_pvalue_binary = fisher.test(tt1, simulate.p.value = T, 
                                     B = 2E5)$p.value
  for(j in 1:length(nclass)){
    tt1 <-  table(factor(ifelse(subtype == nclass[j],1,0)), response)
    fisher_table_row[j] = fisher.test(tt1,simulate.p.value=TRUE, 
                                      B = 2E5)$p.value
    tt1 <- table(factor(ifelse(subtype == nclass[j],1,0)), response_binary)
    ft = fisher.test(tt1,simulate.p.value=TRUE, B = 2E5)
    fisher_table_min[j] = ft$p.value
    odds_ratio_min[j] = ft$estimate
    for(l in 1:length(nresponse)){
      tt = data.frame(subtype = factor(ifelse(subtype == nclass[j],1,0)), 
                      response = factor(ifelse(response==nresponse[l], 1, 0)))
      tt1 <- table(tt$subtype , tt$ response)
      ft = fisher.test(tt1)
      fisher_table_rest[j,l] <- ft$p.value
      odds_ratio_rest[j,l] <- ft$estimate
      tt=tt[response %in%c(nresponse[c(1,l)]),]
      tt1 <-  table(tt$subtype , tt$response)
      ft = fisher.test(tt1)
      fisher_table_base[j,l] <- ft$p.value
      odds_ratio_base[j,l] <- ft$estimate
    }
  }
  rownames(fisher_table_row) = rownames(fisher_table_min) = 
    rownames(fisher_table_rest) = rownames(fisher_table_base) = 
    paste0("Subtype", nclass)
  colnames(fisher_table_rest) = colnames(fisher_table_base) = nresponse
  global = c(global_pvalue_multiresponse, global_pvalue_binary)
  names(global) = c("multi", "binary")
  return(list(global = global, fisher_table_rest = fisher_table_rest, 
              fisher_table_base = fisher_table_base,
              fisher_table_row = fisher_table_row, 
              fisher_table_min = fisher_table_min,
              odds_ratio_rest = odds_ratio_rest,
              odds_ratio_base = odds_ratio_base,
              odds_ratio_min = odds_ratio_min))
}

# Heatmap for Fisher's exact test - post acute status
heatmap_fisher_post_acute <- function(test_fisher, feature_name,
                                      plot_margin_b = 0) {
  df_fisher <- as.data.frame(log(test_fisher$odds_ratio_min))
  colnames(df_fisher) = feature_name
  df_fisher <- pivot_longer(df_fisher, cols = everything())
  df_fisher$subtype <- c("A", "B", "C", "E", "F")
  df_fisher$subtype <- factor(df_fisher$subtype)
  df_fisher$subtype <- factor(df_fisher$subtype,
                              levels = rev(levels(df_fisher$subtype)))
  
  pval_fisher <- as.data.frame(test_fisher$fisher_table_min)
  colnames(pval_fisher) <- feature_name
  pval_fisher <- pivot_longer(pval_fisher, cols = everything())
  df_fisher$pval <- pval_fisher$value
  
  df_fisher %>% ggplot(mapping = aes(x = name, y = subtype, fill = value)) +
    geom_tile(linewidth = 0.2) +
    scale_fill_gradient2(low = "#284AA0", mid = "white", 
                         high = "#BC4A4A", midpoint = 0,
                         name = expression("log"*"("*"odds ratio"*")"),
                         limits = c(-1, 1)) +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    ylab("") +
    xlab("") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
          legend.position = "left",
          legend.title = element_text(size = 15, angle = 90, hjust = 0.5),
          legend.text = element_text(size = 15, angle = 90, hjust = 0.5),
          legend.spacing.y = unit(1, 'cm'),
          legend.key.size = unit(1.8, 'lines'),
          plot.margin = margin(b = plot_margin_b, t = 65)) +
    geom_text(aes(label = ifelse(pval <= 0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 9, vjust = 0.72, color = "black")
}