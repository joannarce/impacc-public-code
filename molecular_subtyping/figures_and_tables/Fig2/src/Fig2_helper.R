### Helpers for Figure 2 ###

# Logistic regression adjusted for demographics (one vs other comparison)
test_clinical_logit_demo <- function(subtype_label, clin_data, demo_data) {
  subtype_label = factor(subtype_label)
  nclass = length(unique(subtype_label))
  
  est_logit = matrix(NA, nrow = nclass, ncol = ncol(clin_data))
  rownames(est_logit) = paste0("Subtype", levels(subtype_label))
  colnames(est_logit) = colnames(clin_data)
  p_logit = matrix(NA, nrow = nclass, ncol = ncol(clin_data))
  rownames(p_logit) = paste0("Subtype", levels(subtype_label))
  colnames(p_logit) = colnames(clin_data)
  p_lr = array(NA, ncol(clin_data))
  names(p_lr) = colnames(clin_data)
  
  for (i in 1:ncol(clin_data)) {
    response = clin_data[,i]
    id_single = (!is.na(subtype_label)) & (!is.na(response))
    subtype = factor(subtype_label[id_single])
    response = response[id_single]
    demo_data = demo_data[id_single,]
    names(response) = "response"
    df_lr = cbind(data.frame(response = response), demo_data, subtype)
    model1 = glm(response ~ ., family = binomial, data = df_lr)
    p_lr[i] <- anova(model1, test = "Chisq")[ncol(demo_data)+2,5]
    
    x_clin = model.matrix(~subtype-1)
    x_clin = data.frame(x_clin[,1:(nclass-1)]-x_clin[,nclass])
    colnames(x_clin) = paste0("Subtype", levels(subtype)[1:(nclass-1)])
    df_logit = cbind(data.frame(response = response), x_clin, demo_data)
    model2 = glm(response ~ ., family = binomial, data = df_logit)
    est_logit[1:(nclass-1),i] <- summary(model2)$coef[2:nclass, 1]
    p_logit[1:(nclass-1),i] <- summary(model2)$coef[2:nclass, 4]
    
    x_clin = model.matrix(~subtype-1)
    x_clin = data.frame(x_clin[,2:(nclass)]-x_clin[,1])
    colnames(x_clin) = paste0("Subtype", levels(subtype)[2:nclass])
    df_logit = cbind(data.frame(response = response), x_clin, demo_data)
    model2 = glm(response ~ ., family = binomial, data = df_logit)
    est_logit[nclass,i] <- summary(model2)$coef[nclass, 1]
    p_logit[nclass,i] <- summary(model2)$coef[nclass, 4]
  }
  return(list(est_logit = est_logit, p_logit = p_logit, p_lr = p_lr))
}

# Logistic regression (one vs other comparison)
test_clinical_logit <- function(subtype_label, clin_data) {
  subtype_label = factor(subtype_label)
  nclass = length(unique(subtype_label))
  
  est_logit = matrix(NA, nrow = nclass, ncol = ncol(clin_data))
  rownames(est_logit) = paste0("Subtype", LETTERS[1:nclass])
  colnames(est_logit) = colnames(clin_data)
  p_logit = matrix(NA, nrow = nclass, ncol = ncol(clin_data))
  rownames(p_logit) = paste0("Subtype", LETTERS[1:nclass])
  colnames(p_logit) = colnames(clin_data)
  p_lr = array(NA, ncol(clin_data))
  names(p_lr) = colnames(clin_data)
  
  for (i in 1:ncol(clin_data)) {
    response = clin_data[,i]
    id_single = (!is.na(subtype_label)) & (!is.na(response))
    subtype = factor(subtype_label[id_single])
    response = response[id_single]
    names(response) = "response"
    df_lr = cbind(data.frame(response = response), subtype)
    model1 = glm(response ~ subtype, family = binomial, data = df_lr)
    p_lr[i] <- anova(model1, test = "Chisq")[2,5]
    
    x_clin = model.matrix(~subtype-1)
    x_clin = data.frame(x_clin[,1:(nclass-1)]-x_clin[,nclass])
    colnames(x_clin) = paste0("Subtype", levels(subtype)[1:(nclass-1)])
    df_logit = cbind(data.frame(response = response), x_clin)
    model2 = glm(response ~ ., family = binomial, data = df_logit)
    est_logit[1:(nclass-1),i] <- summary(model2)$coef[-1, 1]
    p_logit[1:(nclass-1),i] <- summary(model2)$coef[-1, 4]
    
    x_clin = model.matrix(~subtype-1)
    x_clin = data.frame(x_clin[,2:(nclass)]-x_clin[,1])
    colnames(x_clin) = paste0("Subtype", levels(subtype)[2:nclass])
    df_logit = cbind(data.frame(response = response), x_clin)
    model2 = glm(response ~ ., family = binomial, data = df_logit)
    est_logit[nclass,i] <- summary(model2)$coef[nclass, 1]
    p_logit[nclass,i] <- summary(model2)$coef[nclass, 4]
  }
  return(list(est_logit = est_logit, p_logit = p_logit, p_lr = p_lr))
}

# Ordinal regression (one vs other comparison)
test_ordinal <- function(subtype, response) {
  id = (!is.na(subtype)) & (!is.na(response))
  subtype = factor(subtype[id])
  response = factor(response[id])
  nclass = length(unique(subtype))
  x_ordinal = model.matrix(~subtype-1)
  x_ordinal = data.frame(x_ordinal[,1:(nclass-1)]-x_ordinal[,nclass])
  colnames(x_ordinal) = paste0("Subtype", levels(subtype)[1:(nclass-1)])
  df_ordinal = cbind(data.frame(response = response),x_ordinal)
  model1 <- clm(response ~ ., data = df_ordinal)
  coef_ordinal <- tail(summary(model1)$coefficients, nclass-1)
  est_ordinal <- coef_ordinal[,1]
  p_ordinal <- coef_ordinal[,4]
  
  x_ordinal = model.matrix(~subtype-1)
  x_ordinal = data.frame(x_ordinal[,2:(nclass)]-x_ordinal[,1])
  colnames(x_ordinal) = paste0("Subtype", levels(subtype)[2:nclass])
  df_ordinal = cbind(data.frame(response = response),x_ordinal)
  model1 <- clm(response ~ ., data = df_ordinal)
  coef_ordinal1 <- tail(summary(model1)$coefficients, 1)
  est_ordinal1 <- coef_ordinal1[,1]
  p_ordinal1 <- coef_ordinal1[,4]
  est_ordinal = append(est_ordinal, est_ordinal1)
  names(est_ordinal) = paste0("Subtype", levels(subtype)[1:nclass])
  p_ordinal = append(p_ordinal, p_ordinal1)
  names(p_ordinal) = paste0("Subtype", levels(subtype)[1:nclass])
  return(list(est_ordinal = est_ordinal, p_ordinal = p_ordinal))
}

# Ordinal regression adjusted for demographics (one vs other comparison)
test_ordinal_demo <- function(subtype, response, demo_data) {
  id = (!is.na(subtype)) & (!is.na(response)) & (!rowSums(is.na(demo_data)))
  subtype = factor(subtype[id])
  response = factor(response[id])
  demo_data = demo_data[id,]
  nclass = length(unique(subtype))
  nlevel_resp = nlevels(response)
  
  x_ordinal = model.matrix(~subtype-1)
  x_ordinal = data.frame(x_ordinal[,1:(nclass-1)]-x_ordinal[,nclass])
  colnames(x_ordinal) = paste0("Subtype", levels(subtype)[1:(nclass-1)])
  df_ordinal = cbind(data.frame(response = response), x_ordinal, demo_data)
  model1 <- clm(response ~ ., data = df_ordinal)
  model_AIC <- AIC(model1)
  model_loglik <- logLik(model1)
  coef_ordinal <- head(summary(model1)$coefficients, 
                       ((nlevel_resp-1)+(nclass-1)))
  coef_ordinal <- tail(coef_ordinal, nclass-1)
  est_ordinal <- coef_ordinal[,1]
  p_ordinal <- coef_ordinal[,4]
  
  x_ordinal = model.matrix(~subtype-1)
  x_ordinal = data.frame(x_ordinal[,2:(nclass)]-x_ordinal[,1])
  colnames(x_ordinal) = paste0("Subtype", levels(subtype)[2:nclass])
  df_ordinal = cbind(data.frame(response = response), x_ordinal, demo_data)
  model1 <- clm(response ~ ., data = df_ordinal)
  coef_ordinal1 <- head(summary(model1)$coefficients, 
                        ((nlevel_resp-1)+(nclass-1)))
  coef_ordinal1 <- coef_ordinal1[(nlevel_resp-1)+(nclass-1),]
  est_ordinal1 <- coef_ordinal1[1]
  p_ordinal1 <- coef_ordinal1[4]
  est_ordinal = append(est_ordinal, est_ordinal1)
  names(est_ordinal) = paste0("Subtype", levels(subtype)[1:nclass])
  p_ordinal = append(p_ordinal, p_ordinal1)
  names(p_ordinal) = paste0("Subtype", levels(subtype)[1:nclass])
  return(list(est_ordinal = est_ordinal, p_ordinal = p_ordinal,
              model_AIC = model_AIC, model_loglik = model_loglik))
}

# Heatmap for clinical characteristics
heatmap_clinical <- function(test_logit, clin_names, select_columns,
                             adjust_p = TRUE, capped_est = 2) {
  if (adjust_p) {
    pval_heatmap <- as.data.frame(test_logit$p_adj)
  } else {
    pval_heatmap <- as.data.frame(test_logit$p_logit)
  }
  pval_heatmap <- pivot_longer(pval_heatmap, cols = everything())
  
  df_heatmap <- data.frame(test_logit$est_logit)
  df_heatmap[df_heatmap < -capped_est] = -capped_est
  df_heatmap[df_heatmap > capped_est] = capped_est
  n_feature = ncol(df_heatmap)
  n_subtype = nrow(df_heatmap)
  
  df_heatmap <- pivot_longer(df_heatmap, cols = everything())
  df_heatmap$pval <- pval_heatmap$value
  df_heatmap$subtype <- rep(sapply(c("A","B", "C", "E", "F"), 
                                   function(i) paste0("Subtype ", i)), 
                            each = n_feature)
  df_heatmap$subtype <- factor(df_heatmap$subtype)
  df_heatmap$subtype <- factor(df_heatmap$subtype,
                               levels = rev(levels(df_heatmap$subtype)))
  
  df_heatmap <- df_heatmap %>%
    mutate(category = case_when(name == "ever_icu" ~ " ",
                                name == "ever_ever_esc" ~ " ",
                                str_starts(name, "baseline") ~ "Baseline",
                                str_starts(name, "comorb") ~ "Comorbidity",
                                str_starts(name, "comp") ~ "Complication"))
  logit_categories <- c("Baseline", "Comorbidity", "Complication", " ")
  df_heatmap$category[!(df_heatmap$category %in% logit_categories)] <- "Demographic"
  df_heatmap$category <- factor(df_heatmap$category, 
                                levels = c("Demographic", "Comorbidity", "Baseline", 
                                           "Complication", " "))
  df_heatmap$name <- clin_names[df_heatmap$name]
  
  df_heatmap %>% ggplot(aes(x=name, y = subtype, fill = value)) +
    geom_tile(color = "grey50", linewidth=0.1) +
    scale_fill_gradient2(low = "#284AA0", mid="white", high = "#BC4A4A", midpoint=0, 
                         name = expression("Coefficient"),
                         limits = c(-capped_est, capped_est)) +
    cowplot::theme_cowplot() + 
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylab("") +
    xlab("") +
    facet_grid(.~category, scales = "free", space = "free") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 18, hjust = 0),
          axis.text.x = element_text(size = 15),
          legend.position = "top",
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          panel.spacing = unit(0, "lines"),
          strip.text.x = element_text(size = 18, face = "bold"),
          strip.background = element_rect(fill = "white"),
          legend.key.size = unit(7.5, "mm")) +
    geom_text(aes(label = ifelse(pval <= 0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 5, vjust = 0.78)
}

# Donut plot for clinical characteristics
clin_donut_plot <- function(clin_var, subtype, palette_name, label_position = 3,
                            label_vjust = 0.5) {
  df_donut = data.frame(data.frame(Subtype = factor(subtype_hcl$label),
                                   Var = clin_var))
  
  df_donut = df_donut %>%
    group_by(Subtype, Var) %>%
    summarise(Count = n(), .groups = "drop_last") %>%
    mutate(Total = sum(Count)) %>%
    mutate(Percentage = (Count/Total)*100) %>%
    mutate(Percentage = round(Percentage, 0))
  
  ggplot(df_donut[df_donut$Subtype == subtype, ], 
         aes(x = 2, y = Percentage, fill = Var)) +
    geom_col(color = "black") +
    geom_text(aes(label = paste0(Percentage, "%")), size = 4,
              position = position_stack(vjust = 0.5)) +
    geom_text(aes(x = label_position, label = Var),
              position = position_stack(vjust = label_vjust), size = 5, 
              color = "black") +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = palette_name) +
    xlim(0.5, 3.4) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.margin = margin(b = -500, t = -500)) +
    annotate("text", x = 0.55, y = max(df_donut$Percentage)/2, 
             label = paste(subtype), size = 5)
}

# Mosaic plot for clinical characteristics
mosaic_plot_clinical <- function(subtype, feature, name, cex_title = 1.2) {
  subtype <- paste0(subtype)
  feature <- ifelse(feature == 1, "Yes", "No")
  df_mosaic <- data.frame(subtype=subtype, feature=factor(feature))
  original_par <- par(no.readonly = TRUE)
  par(mgp = c(0.2, 0, 0), mar = c(1.5,1.5,1,1))
  mosaicplot(table(df_mosaic$subtype, df_mosaic$feature), 
             main = "", xlab = "", ylab = "",
             color = c("No" = "#00BFC480", "Yes" = "#F8766D80"), 
             cex.axis = 1.2)
  title(xlab = "Subtype", cex.lab = 1.2)
  title(ylab = name, cex.lab = cex_title)
  par(original_par)
}

# Complication adjusted by demographics (group-wise comparison)
comp_comparison_adjust <- function(df_comp, df_adjust, subtype, 
                                   subtype_remove = "D") {
  if (!is.null(subtype_remove)) {
    idx = which(subtype != subtype_remove)
    df_comp = df_comp[idx, ]
    df_adjust = df_adjust[idx, ]
    subtype = subtype[idx]
  }
  subtype = factor(subtype)
  
  est_logit = matrix(NA, nrow = 5, ncol = ncol(df_comp))
  rownames(est_logit) = c("EF vs ABC", "F vs E",
                          "A vs BC", "B vs AC", "C vs AB")
  colnames(est_logit) = colnames(df_comp)
  p_logit = matrix(NA, nrow = 5, ncol = ncol(df_comp))
  rownames(p_logit) = c("EF vs ABC", "F vs E",
                        "A vs BC", "B vs AC", "C vs AB")
  colnames(p_logit) = colnames(df_comp)
  
  x1_comp = model.matrix(~subtype)
  colnames(x1_comp) = c("Intercept", "EFvsABC", "FvsE", "AvsBC", "BvsAC")
  x1_comp[, "EFvsABC"] = -1
  x1_comp[subtype == "E" | subtype == "F", "EFvsABC"] = 1
  x1_comp[, "FvsE"] = 0
  x1_comp[subtype == "E", "FvsE"] = -1
  x1_comp[subtype == "F", "FvsE"] = 1
  x1_comp[, "AvsBC"] = 0
  x1_comp[subtype == "A", "AvsBC"] = 1
  x1_comp[subtype == "C", "AvsBC"] = -1
  x1_comp[, "BvsAC"] = 0
  x1_comp[subtype == "B", "BvsAC"] = 1
  x1_comp[subtype == "C", "BvsAC"] = -1
  
  x2_comp = model.matrix(~subtype)
  colnames(x2_comp) = c("Intercept", "EFvsABC", "FvsE", "BvsAC", "CvsAB")
  x2_comp[, "EFvsABC"] = -1
  x2_comp[subtype == "E" | subtype == "F", "EFvsABC"] = 1
  x2_comp[, "FvsE"] = 0
  x2_comp[subtype == "E", "FvsE"] = -1
  x2_comp[subtype == "F", "FvsE"] = 1
  x2_comp[, "BvsAC"] = 0
  x2_comp[subtype == "B", "BvsAC"] = 1
  x2_comp[subtype == "A", "BvsAC"] = -1
  x2_comp[, "CvsAB"] = 0
  x2_comp[subtype == "C", "CvsAB"] = 1
  x2_comp[subtype == "A", "CvsAB"] = -1
  
  for (i in 1:ncol(df_comp)) {
    response = df_comp[[i]]
    df_logit = cbind(data.frame(response = response), df_adjust, x1_comp)
    model1 = glm(response ~ ., family = binomial, data = df_logit)
    coef_logit = summary(model1)$coef[c("EFvsABC", "FvsE", "AvsBC", "BvsAC"),]
    est_logit[1:4,i] <- coef_logit[,1]
    p_logit[1:4,i] <- coef_logit[,4]
    
    df_logit = cbind(data.frame(response = response), df_adjust, x2_comp)
    model2 = glm(response ~ ., family = binomial, data = df_logit)
    est_logit[5,i] <- summary(model2)$coef["CvsAB",1]
    p_logit[5,i] <- summary(model2)$coef["CvsAB",4]
  }
  return(list(est_logit = est_logit, p_logit = p_logit))
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