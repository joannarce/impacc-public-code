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

# Heatmap for Fisher's exact test
heatmap_fisher <- function(test_fisher) {
  df_fisher <- cbind(test_fisher$odds_ratio_min, 
                     test_fisher$odds_ratio_base[, -1])
  df_fisher <- log(df_fisher)
  colnames(df_fisher) <- c("Any Deficit", "Physical", "Cognitive", "Multiple")
  df_fisher <- as.data.frame(df_fisher)
  df_fisher <- pivot_longer(df_fisher, cols = everything())
  
  pval_fisher <- cbind(test_fisher$fisher_table_min, 
                       test_fisher$fisher_table_base[, -1])
  colnames(pval_fisher) <- c("Any Deficit", "Physical", "Cognitive", "Multiple")
  pval_fisher <- as.data.frame(pval_fisher)
  pval_fisher <- pivot_longer(pval_fisher, cols = everything())
  
  df_fisher$pval <- pval_fisher$value
  
  df_fisher$subtype <- rep(c("Subtype A", "Subtype B", "Subtype C", 
                             "Subtype E", "Subtype F"),
                           each = ncol(test_fisher$fisher_table_base))
  df_fisher$subtype <- factor(df_fisher$subtype)
  df_fisher$subtype <- factor(df_fisher$subtype,
                              levels = rev(levels(df_fisher$subtype)))
  df_fisher$name <- factor(df_fisher$name, 
                           levels = c("Any Deficit", "Physical",
                                      "Cognitive", "Multiple"))
  
  df_fisher %>% ggplot(mapping = aes(x = name, y = subtype, fill = value)) +
    geom_tile(linewidth = 0.2) +
    scale_fill_gradient2(low = "#284AA0", mid = "white", 
                         high = "#BC4A4A", midpoint = 0,
                         name = expression("log"*"("*"odds ratio"*")"),
                         limits = c(-max(abs(df_fisher$value)),
                                    max(abs(df_fisher$value)))) +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    ylab("") +
    xlab("") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 15, angle = 90),
          legend.position = "left",
          legend.title = element_text(size = 15, angle = 90),
          legend.text = element_text(size = 15, angle = 90, hjust = 0.5),
          legend.spacing.y = unit(1, 'cm'),
          legend.key.size = unit(2, 'lines'))+
    geom_text(aes(label = ifelse(pval <= 0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 9, vjust = 0.72, color = "black")
}

# Heatmap for Fisher's exact test - global test
heatmap_fisher_global <- function(test_fisher) {
  df_global <- data.frame("Global" = -log10(test_fisher$fisher_table_row))
  df_global <- pivot_longer(df_global, cols = everything())
  df_global$subtype <- c("A", "B", "C", "E", "F")
  df_global$subtype <- factor(df_global$subtype)
  df_global$subtype <- factor(df_global$subtype,
                              levels = rev(levels(df_global$subtype)))
  
  df_global %>% ggplot(mapping = aes(x = name, y = subtype, size = value)) +
    geom_point(color = "grey30") +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    ylab("") +
    xlab("") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15, hjust = 0),
          axis.text.x = element_text(size = 15, angle = 90),
          legend.position = "left",
          legend.title = element_text(size = 15, angle = 90),
          legend.text = element_text(size = 15, angle = 90),
          plot.margin = margin(b = 20)) +
    labs(size = expression(-log[10]*"("*italic("P")*"-value"*")")) +
    geom_text(aes(label = ifelse(abs(value) >= -log10(0.001), "***",
                                 ifelse(abs(value) >= -log10(0.01), "**",
                                        ifelse(abs(value) >= -log10(0.05), "*", "")))),
              size = 7, vjust = 0.78, color = "white") +
    scale_size(range = c(3, 12))
}

# Cox proportional hazards regression adjusted for demographics
# (one vs other comparison)
coxph_coef <- function(surv_object, subtype, demo_data) {
  subtype = factor(subtype)
  nclass = nlevels(subtype)
  x_coxph = model.matrix(~subtype-1)
  x_coxph = data.frame(x_coxph[,1:(nclass-1)]-x_coxph[,nclass])
  colnames(x_coxph) = paste0("Subtype", levels(subtype)[1:(nclass-1)])
  df_coxph = data.frame(x_coxph)
  df_coxph = cbind(data.frame(x_coxph), demo_data)
  model1 <- coxph(surv_object ~., data = df_coxph)
  coef_coxph <- summary(model1)$coefficients[1:(nclass-1), ]
  
  x_coxph = model.matrix(~subtype-1)
  x_coxph = data.frame(x_coxph[,2:(nclass)]-x_coxph[,1])
  colnames(x_coxph) = paste0("Subtype", levels(subtype)[2:nclass])
  df_coxph = data.frame(x_coxph)
  df_coxph = cbind(data.frame(x_coxph), demo_data)
  model2 <- coxph(surv_object ~., data = df_coxph)
  coef_coxph <- rbind(coef_coxph, summary(model2)$coefficients[(nclass-1),])
  rownames(coef_coxph) <- paste0("Subtype", levels(subtype)[1:nclass])
  return(coef_coxph)
}