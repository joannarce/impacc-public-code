### Helpers for Extended Data 2 ###

# Complication adjusted by demographics and comorbidities (group-wise comparison)
# Predicted probability isolating the effect of subtype
comp_comparison_adjust_pred <- function(df_comp, df_adjust, subtype_label, 
                                        subtype_remove = "D") {
  if (!is.null(subtype_remove)) {
    id = (subtype_label != subtype_remove) & (!is.na(subtype_label)) & 
      (!rowSums(is.na(df_comp))) & (!rowSums(is.na(df_adjust)))
    df_comp = df_comp[id, ]
    df_adjust = df_adjust[id, ]
    subtype_label = subtype_label[id]
  }
  subtype_label = factor(subtype_label)
  
  est_logit = matrix(NA, nrow = 5, ncol = ncol(df_comp))
  rownames(est_logit) = c("EF vs ABC", "F vs E",
                          "A vs BC", "B vs AC", "C vs AB")
  colnames(est_logit) = colnames(df_comp)
  p_logit = matrix(NA, nrow = 5, ncol = ncol(df_comp))
  rownames(p_logit) = c("EF vs ABC", "F vs E",
                        "A vs BC", "B vs AC", "C vs AB")
  colnames(p_logit) = colnames(df_comp)
  pred_mean_logit = matrix(NA, nrow = 5, ncol = ncol(df_comp))
  rownames(pred_mean_logit) = c("A", "B", "C", "E", "F")
  colnames(pred_mean_logit) = colnames(df_comp)
  
  x1_comp = model.matrix(~subtype_label)
  colnames(x1_comp) = c("Intercept", "EFvsABC", "FvsE", "AvsBC", "BvsAC")
  x1_comp[, "EFvsABC"] = -1
  x1_comp[subtype_label == "E" | subtype_label == "F", "EFvsABC"] = 1
  x1_comp[, "FvsE"] = 0
  x1_comp[subtype_label == "E", "FvsE"] = -1
  x1_comp[subtype_label == "F", "FvsE"] = 1
  x1_comp[, "AvsBC"] = 0
  x1_comp[subtype_label == "A", "AvsBC"] = 1
  x1_comp[subtype_label == "C", "AvsBC"] = -1
  x1_comp[, "BvsAC"] = 0
  x1_comp[subtype_label == "B", "BvsAC"] = 1
  x1_comp[subtype_label == "C", "BvsAC"] = -1
  
  x2_comp = model.matrix(~subtype_label)
  colnames(x2_comp) = c("Intercept", "EFvsABC", "FvsE", "BvsAC", "CvsAB")
  x2_comp[, "EFvsABC"] = -1
  x2_comp[subtype_label == "E" | subtype_label == "F", "EFvsABC"] = 1
  x2_comp[, "FvsE"] = 0
  x2_comp[subtype_label == "E", "FvsE"] = -1
  x2_comp[subtype_label == "F", "FvsE"] = 1
  x2_comp[, "BvsAC"] = 0
  x2_comp[subtype_label == "B", "BvsAC"] = 1
  x2_comp[subtype_label == "A", "BvsAC"] = -1
  x2_comp[, "CvsAB"] = 0
  x2_comp[subtype_label == "C", "CvsAB"] = 1
  x2_comp[subtype_label == "A", "CvsAB"] = -1
  
  for (i in 1:ncol(df_comp)) {
    response = df_comp[[i]]
    df_logit = cbind(data.frame(response = response), df_adjust, x1_comp)
    model1 = glm(response ~ ., family = binomial, data = df_logit)
    coef_logit = summary(model1)$coef[c("EFvsABC", "FvsE", "AvsBC", "BvsAC"),]
    est_logit[1:4,i] <- coef_logit[,1]
    p_logit[1:4,i] <- coef_logit[,4]
    
    df_new_A <- df_logit
    df_new_A$Intercept <- 1
    df_new_A$EFvsABC <- -1
    df_new_A$FvsE <- 0
    df_new_A$AvsBC <- 1
    df_new_A$BvsAC <- 0
    pred_A <- predict(model1, newdata = df_new_A, type = "response")
    pred_mean_logit["A", i] <- mean(pred_A)
    
    df_new_B <- df_logit
    df_new_B$Intercept <- 1
    df_new_B$EFvsABC <- -1
    df_new_B$FvsE <- 0
    df_new_B$AvsBC <- 0
    df_new_B$BvsAC <- 1
    pred_B <- predict(model1, newdata = df_new_B, type = "response")
    pred_mean_logit["B", i] <- mean(pred_B)
    
    df_new_C <- df_logit
    df_new_C$Intercept <- 1
    df_new_C$EFvsABC <- -1
    df_new_C$FvsE <- 0
    df_new_C$AvsBC <- -1
    df_new_C$BvsAC <- -1
    pred_C <- predict(model1, newdata = df_new_C, type = "response")
    pred_mean_logit["C", i] <- mean(pred_C)
    
    df_new_E <- df_logit
    df_new_E$Intercept <- 1
    df_new_E$EFvsABC <- 1
    df_new_E$FvsE <- -1
    df_new_E$AvsBC <- 0
    df_new_E$BvsAC <- 0
    pred_E <- predict(model1, newdata = df_new_E, type = "response")
    pred_mean_logit["E", i] <- mean(pred_E)
    
    df_new_F <- df_logit
    df_new_F$Intercept <- 1
    df_new_F$EFvsABC <- 1
    df_new_F$FvsE <- 1
    df_new_F$AvsBC <- 0
    df_new_F$BvsAC <- 0
    pred_F <- predict(model1, newdata = df_new_F, type = "response")
    pred_mean_logit["F", i] <- mean(pred_F)
    
    df_logit = cbind(data.frame(response = response), df_adjust, x2_comp)
    model2 = glm(response ~ ., family = binomial, data = df_logit)
    est_logit[5,i] <- summary(model2)$coef["CvsAB",1]
    p_logit[5,i] <- summary(model2)$coef["CvsAB",4]
  }
  return(list(est_logit = est_logit, p_logit = p_logit, 
              pred_mean_logit = pred_mean_logit))
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

### predicted probability
# solve(solve(t(S)%*%S)%*%t(S) %>% fractions())
# logit(p) = b0 + b1*I_EFvsABC + b2*I_FvsE + b3*I_AvsBC + b4*I_BvsAC
# A: b0-b1+b3
# B: b0-b1+b4
# C: b0-b1-b3-b4
# E: b0+b1-b2
# F: b0+b1+b2

# linear regression adjusted for demographics
test_clinical_linear_demo_group_comparison <- function(
    subtype_label, clin_data, demo_data) {
  subtype_label = factor(subtype_label)
  nclass = length(unique(subtype_label))
  
  est_linear = matrix(NA, nrow = 5, ncol = ncol(clin_data))
  rownames(est_linear) = c("EF vs ABC", "F vs E",
                           "A vs BC", "B vs AC", "C vs AB")
  colnames(est_linear) = colnames(clin_data)
  p_linear = matrix(NA, nrow = 5, ncol = ncol(clin_data))
  rownames(p_linear) = c("EF vs ABC", "F vs E",
                         "A vs BC", "B vs AC", "C vs AB")
  colnames(p_linear) = colnames(clin_data)
  p_lr = array(NA, ncol(clin_data))
  names(p_lr) = colnames(clin_data)
  est_comorb = matrix(NA, nrow = sum(grepl("^comorb", names(demo_data))), 
                      ncol = ncol(clin_data))
  rownames(est_comorb) = names(demo_data)[grepl("^comorb", names(demo_data))]
  colnames(est_comorb) = colnames(clin_data)
  p_comorb = matrix(NA, nrow = sum(grepl("^comorb", names(demo_data))), 
                    ncol = ncol(clin_data))
  rownames(p_comorb) = names(demo_data)[grepl("^comorb", names(demo_data))]
  colnames(p_comorb) = colnames(clin_data)
  
  x1_clin = model.matrix(~subtype_label)
  colnames(x1_clin) = c("Intercept", "EFvsABC", "FvsE", "AvsBC", "BvsAC")
  x1_clin[, "EFvsABC"] = -1
  x1_clin[subtype_label == "E" | subtype_label == "F", "EFvsABC"] = 1
  x1_clin[, "FvsE"] = 0
  x1_clin[subtype_label == "E", "FvsE"] = -1
  x1_clin[subtype_label == "F", "FvsE"] = 1
  x1_clin[, "AvsBC"] = 0
  x1_clin[subtype_label == "A", "AvsBC"] = 1
  x1_clin[subtype_label == "C", "AvsBC"] = -1
  x1_clin[, "BvsAC"] = 0
  x1_clin[subtype_label == "B", "BvsAC"] = 1
  x1_clin[subtype_label == "C", "BvsAC"] = -1
  
  x2_clin = model.matrix(~subtype_label)
  colnames(x2_clin) = c("Intercept", "EFvsABC", "FvsE", "BvsAC", "CvsAB")
  x2_clin[, "EFvsABC"] = -1
  x2_clin[subtype_label == "E" | subtype_label == "F", "EFvsABC"] = 1
  x2_clin[, "FvsE"] = 0
  x2_clin[subtype_label == "E", "FvsE"] = -1
  x2_clin[subtype_label == "F", "FvsE"] = 1
  x2_clin[, "BvsAC"] = 0
  x2_clin[subtype_label == "B", "BvsAC"] = 1
  x2_clin[subtype_label == "A", "BvsAC"] = -1
  x2_clin[, "CvsAB"] = 0
  x2_clin[subtype_label == "C", "CvsAB"] = 1
  x2_clin[subtype_label == "A", "CvsAB"] = -1
  
  for (i in 1:ncol(clin_data)) {
    response = clin_data[,i]
    id_single = (!is.na(subtype_label)) & (!is.na(response))
    subtype = factor(subtype_label[id_single])
    response = response[id_single]
    demo = demo_data[id_single,]
    names(response) = "response"
    df_lr = cbind(data.frame(response = response), demo, subtype)
    model1 = lm(response ~ ., data = df_lr)
    p_lr[i] <- anova(model1)[ncol(demo)+1,5]
    
    df_logit = cbind(data.frame(response = response), x1_clin, demo)
    model2 = lm(response ~ ., data = df_logit)
    coef_logit = summary(model2)$coef[c("EFvsABC", "FvsE", "AvsBC", "BvsAC"),]
    est_linear[1:4,i] <- coef_logit[,1]
    p_linear[1:4,i] <- coef_logit[,4]
    coef_comorb <- summary(model2)$coef[names(demo_data)[grepl("^comorb", 
                                                               names(demo_data))],]
    est_comorb[,i] <- coef_comorb[,1]
    p_comorb[,i] <- coef_comorb[,4]
    
    df_logit = cbind(data.frame(response = response), x2_clin, demo)
    model2 = lm(response ~ ., data = df_logit)
    est_linear[5,i] <- summary(model2)$coef["CvsAB",1]
    p_linear[5,i] <- summary(model2)$coef["CvsAB",4]
  }
  return(list(est_linear = est_linear, p_linear = p_linear, p_lr = p_lr,
              est_comorb = est_comorb, p_comorb = p_comorb))
}

# Heatmap for complication category - subtype
heatmap_comp_category_subtype <- function(comp_test_category, 
                                          comp_test_category_comorb,
                                          comp_names) {
  df1 <- as.data.frame(comp_test_category$est_linear)
  n_feature <- ncol(df1)
  n_subtype <- nrow(df1)
  rowname_heatmap <- rownames(df1)
  df1 <- pivot_longer(df1, cols = everything())
  colnames(df1) <- c("name", "est")
  df1$group <- "Vanilla"
  df1$subtype <- rep(rowname_heatmap, each = n_feature)
  df2 <- as.data.frame(comp_test_category_comorb$est_linear)
  df2 <- pivot_longer(df2, cols = everything())
  colnames(df2) <- c("name", "est")
  df2$group <- "Comorbidity Adjusted"
  df2$subtype <- rep(rowname_heatmap, each = n_feature)
  df_heatmap <- rbind(df1, df2)
  
  pval1 <- as.data.frame(comp_test_category$p_adj)
  pval1 <- pivot_longer(pval1, cols = everything())
  pval2 <- as.data.frame(comp_test_category_comorb$p_adj)
  pval2 <- pivot_longer(pval2, cols = everything())
  df_pval <- rbind(pval1, pval2)
  df_heatmap$pval <- df_pval$value
  df_heatmap$name <- comp_names[df_heatmap$name]
  df_heatmap$group <- factor(df_heatmap$group)
  df_heatmap$group <- relevel(df_heatmap$group, ref = "Vanilla")
  
  df_heatmap$subtype <- factor(df_heatmap$subtype)
  df_heatmap$subtype <- factor(df_heatmap$subtype,
                               levels = rev(levels(df_heatmap$subtype)))
  
  df_heatmap %>% ggplot(aes(x = name, y = subtype, fill = est)) +
    geom_tile(color = "grey50", linewidth = 0.1) +
    scale_fill_gradient2(low = "#284AA0", mid="white", high = "#BC4A4A", midpoint=0,
                         name = expression("Coefficient"),
                         limits = c(-max(abs(df_heatmap$est)),
                                    max(abs(df_heatmap$est)))) +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("") +
    xlab("") +
    facet_grid(.~group, scales = "free", space = "free") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15, hjust = 1),
          axis.text.x = element_text(size = 15),
          legend.position = "top",
          legend.title = element_text(size = 15, vjust = 0.8),
          legend.text = element_text(size = 15),
          legend.key.size = unit(2, 'lines'),
          panel.spacing = unit(0, "lines"),
          legend.spacing.y = unit(5, 'mm'),
          legend.spacing.x = unit(5, 'mm'),
          strip.text.x = element_text(size = 15, face = "bold"),
          strip.background = element_rect(fill = "white")) +
    geom_text(aes(label = ifelse(pval <= 0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 5, vjust = 0.72)
}

# Heatmap for complication category - comorbidity
heatmap_comp_category_comorb <- function(test_logit, comp_names, 
                                         clin_names, adjust_p = TRUE) {
  df_heatmap = as.data.frame(test_logit$est_comorb)
  n_feature = ncol(df_heatmap)
  n_subtype = nrow(df_heatmap)
  rowname_heatmap = rownames(df_heatmap)
  if (adjust_p) {
    df_pval <- as.data.frame(test_logit$p_comorb_adj)
  } else {
    df_pval <- as.data.frame(test_logit$p_comorb)
  }
  
  df_heatmap <- pivot_longer(df_heatmap, cols = everything())
  colnames(df_heatmap) <- c("name", "est")
  df_pval <- pivot_longer(df_pval, cols = everything())
  df_heatmap$pval <- df_pval$value
  df_heatmap$subtype <- rep(rowname_heatmap,
                            each = n_feature)
  df_heatmap$name <- comp_names[df_heatmap$name]
  df_heatmap$subtype <- clin_names[df_heatmap$subtype]
  df_heatmap$subtype <- factor(df_heatmap$subtype)
  df_heatmap$subtype <- factor(df_heatmap$subtype,
                               levels = rev(levels(df_heatmap$subtype)))
  
  df_heatmap %>% ggplot(aes(x=name, y = subtype, fill = est)) +
    geom_tile(color = "grey50", size=0.1) +
    scale_fill_gradient2(low = "#284AA0", mid = "white", high = "#BC4A4A", 
                         midpoint = 0, name = expression("Coefficient"),
                         limits = c(-max(abs(df_heatmap$est)),
                                    max(abs(df_heatmap$est)))) +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("") +
    xlab("") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15, hjust = 1),
          axis.text.x = element_text(size = 15),
          legend.position = "top",
          legend.title = element_text(size = 15, vjust = 0.8),
          legend.text = element_text(size = 15, hjust = 0.5),
          legend.key.size = unit(2, 'lines'),
          panel.spacing = unit(0, "lines"),
          legend.spacing.y = unit(5, 'mm'),
          legend.spacing.x = unit(5, 'mm'),
          strip.text.x = element_text(size = 15, face = "bold"),
          strip.background = element_rect(fill = "white")) +
    geom_text(aes(label = ifelse(pval <= 0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 5, vjust = 0.72)
}