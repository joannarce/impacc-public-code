### Helpers for Extended Data 8 ###

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

# Heatmap for medications
heatmap_med <- function(test_logit, med_names, capped_est = 2) {
  pval_heatmap <- as.data.frame(test_logit$p_adj)
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
  df_heatmap$name <- med_names[match(df_heatmap$name, names(med_names))]
  
  df_heatmap %>% ggplot(aes(x = name, y = subtype, fill = value)) +
    geom_tile(color = "grey50", size=0.1) +
    scale_fill_gradient2(low = "#284AA0", mid="white", high = "#BC4A4A", midpoint=0, 
                         name = expression("Coefficient"),
                         limits = c(-capped_est, capped_est)) +
    cowplot::theme_cowplot() + 
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylab("") +
    xlab("") +
    # facet_grid(.~category, scales = "free", space = "free") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15, hjust = 0),
          axis.text.x = element_text(size = 15),
          legend.position = "top",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.spacing = unit(0, "lines"),
          # strip.text.x = element_text(size = 9, face = "bold"),
          # strip.background = element_rect(fill = "white"),
          legend.key.size = unit(7.5, "mm")) +
    geom_text(aes(label = ifelse(pval <= 0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 7, vjust = 0.78)
}

# Heatmap for medications - global test
heatmap_med_global <- function(test_logit, med_names) {
  df_global <- data.frame("Global" = -log10(test_logit$p_lr_adj))
  df_global <- pivot_longer(df_global, cols = everything())
  df_global$medication <- names(test_logit$p_lr_adj)
  df_global$medication <- factor(df_global$medication)
  df_global$medication <- med_names[match(df_global$medication, names(med_names))]
  
  df_global %>% ggplot(mapping = aes(x = medication, y = name, size = value)) +
    geom_point(color = "grey30") +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    ylab("") +
    xlab("") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15, hjust = 0),
          axis.text.x = element_blank(),
          legend.position = c(-0.1, 1), # legend position at top and slightly left
          legend.justification = c(0, 1), 
          legend.direction = "horizontal",
          legend.text = element_text(size = 15),
          plot.margin = margin(l = 32)) +
    labs(size = expression(-log[10]*"("*italic("P")*"-value"*")")) +
    geom_text(aes(label = ifelse(abs(value) >= -log10(0.001), "***",
                                 ifelse(abs(value) >= -log10(0.01), "**",
                                        ifelse(abs(value) >= -log10(0.05), "*", "")))),
              size = 5, vjust = 0.78, color = "white") +
    scale_size(range = c(3, 10))
}

# Binomial regression adjusted for demographics
test_clinical_binomial_demo = function(subtype, response_binary, demo_data){
  id = (!is.na(subtype)) & (!is.na(response_binary)) & (!rowSums(is.na(demo_data)))
  subtype = factor(as.character(subtype[id]))
  response_binary = response_binary[id]
  demo_data = demo_data[id,]
  nclass = length(unique(subtype))
  
  x_pasc = model.matrix(~subtype-1)
  x_pasc = data.frame(x_pasc[,1:(nclass-1)]-x_pasc[,nclass])
  colnames(x_pasc) = paste0("Subtype", levels(subtype)[1:(nclass-1)])
  df_binom = cbind(data.frame(PASC = response_binary), x_pasc, demo_data)
  model1 <- glm(PASC ~ ., family = binomial, data = df_binom)
  p_trees_binary <- coef(summary(model1))[2:nclass, ]
  
  x_pasc = model.matrix(~subtype-1)
  x_pasc = data.frame(x_pasc[,2:(nclass)]-x_pasc[,1])
  colnames(x_pasc) = paste0("Subtype", levels(subtype)[-1])
  df_binom = cbind(data.frame(PASC = response_binary), x_pasc, demo_data)
  model1 <- glm(PASC ~ ., family = binomial, data = df_binom)
  p_trees_binary1 <- coef(summary(model1))[2:nclass, ]
  p_trees_binary = rbind(p_trees_binary, p_trees_binary1[nclass-1,])
  rownames(p_trees_binary)[nclass] = paste0("Subtype", levels(subtype)[nclass])
  
  return(p_binary = p_trees_binary)
}

# Heatmap for binomial regression
heatmap_binomial <- function(test_binomial, val_capped = 1, col_name) {
  pval_heatmap <- as.data.frame(t(test_binomial[, "Pr(>|z|)"]))
  pval_heatmap <- pivot_longer(pval_heatmap, cols = everything())
  
  df_heatmap <- data.frame(t(test_binomial[, "Estimate"]))
  df_heatmap[df_heatmap > val_capped] = val_capped
  df_heatmap[df_heatmap < -val_capped] = -val_capped
  n_subtype = nrow(df_heatmap)
  
  df_heatmap <- pivot_longer(df_heatmap, cols = everything())
  df_heatmap$pval <- pval_heatmap$value
  df_heatmap$subtype <- c("A","B", "C", "E", "F")
  
  df_heatmap$subtype <- factor(df_heatmap$subtype)
  df_heatmap$subtype <- factor(df_heatmap$subtype,
                               levels = rev(levels(df_heatmap$subtype)))
  df_heatmap$name <- col_name
  
  df_heatmap %>% ggplot(mapping = aes(x = name, y = subtype, fill = value)) +
    geom_tile(size = 0.2) +
    scale_fill_gradient2(low = "#284AA0", mid = "white", 
                         high = "#BC4A4A", midpoint = 0,
                         name = expression("Coefficient"),
                         limits = c(-val_capped, val_capped),
                         breaks = c(-val_capped, 0, val_capped)) +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    ylab("") +
    xlab("") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5),
          legend.position = "left",
          legend.title = element_text(size = 15, angle = 90, vjust = 1),
          legend.text = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 0.5),
          legend.spacing.y = unit(1, 'cm'),
          plot.margin = margin(b = -5, t = 8))+
    geom_text(aes(label = ifelse(pval <= -0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 9, vjust = 0.75, color = "black")
}

# Binomial regression adjusted for demographics (multiple responses)
test_multi_pasc_demo <- function(subtype, response, demo_data, 
                                 subtype_remove = "D") {
  id = (!is.na(subtype)) & (!is.na(response)) & (!rowSums(is.na(demo_data))) &
    (subtype != subtype_remove)
  subtype = factor(as.character(subtype[id]))
  response = response[id]
  demo_data = demo_data[id,]
  nclass = length(unique(subtype))
  
  response = data.frame(model.matrix(~response-1))
  colnames(response) = gsub("response", "", colnames(response))
  pasc_names = setdiff(colnames(response), "minimal")
  
  pval_pasc = matrix(NA, nrow = nclass, ncol = length(pasc_names))
  colnames(pval_pasc) = pasc_names
  rownames(pval_pasc) = levels(subtype)
  coef_pasc = matrix(NA, nrow = nclass, ncol = length(pasc_names))
  colnames(coef_pasc) = pasc_names
  rownames(coef_pasc) = levels(subtype)
  
  for (fname in pasc_names) {
    idx_test = which(response[[fname]] == 1 |
                       response[["minimal"]] == 1)
    response_binary = response[idx_test, fname]
    
    x_pasc = model.matrix(~subtype-1)
    x_pasc = data.frame(x_pasc[,1:(nclass-1)]-x_pasc[,nclass])
    colnames(x_pasc) = paste0("Subtype", levels(subtype)[1:(nclass-1)])
    df_binom = cbind(data.frame(PASC = response_binary), 
                     x_pasc[idx_test, ], demo_data[idx_test, ])
    model1 <- glm(PASC ~ ., family = binomial, data = df_binom)
    coef_pasc[1:(nclass-1), fname] <- coef(summary(model1))[2:nclass, "Estimate"]
    pval_pasc[1:(nclass-1), fname] <- coef(summary(model1))[2:nclass, "Pr(>|z|)"]
    
    x_pasc = model.matrix(~subtype-1)
    x_pasc = data.frame(x_pasc[,2:(nclass)]-x_pasc[,1])
    colnames(x_pasc) = paste0("Subtype", levels(subtype)[-1])
    df_binom = cbind(data.frame(PASC = response_binary), 
                     x_pasc[idx_test, ], demo_data[idx_test, ])
    model2 <- glm(PASC ~ ., family = binomial, data = df_binom)
    coef_pasc[nclass, fname] <- coef(summary(model2))[nclass, "Estimate"]
    pval_pasc[nclass, fname] <- coef(summary(model2))[nclass, "Pr(>|z|)"]
  }
  return(list(pval_pasc = pval_pasc, coef_pasc = coef_pasc))
}

# Heatmap for binomial regression (multiple responses)
heatmap_pasc_multi <- function(test_multi, plot_title, val_capped = 1) {
  pval_multi <- as.data.frame(test_multi$pval_pasc)
  pval_multi <- pivot_longer(pval_multi, cols = everything())
  
  df_multi <- as.data.frame(test_multi$coef_pasc)
  colnames(df_multi) <- c("Physical", "Cognitive", "Multiple")
  df_multi[df_multi > val_capped] = val_capped
  df_multi[df_multi < -val_capped] = -val_capped
  df_multi <- pivot_longer(df_multi, cols = everything())
  df_multi$pval <- pval_multi$value
  df_multi$subtype <- rep(c("A","B", "C", "E", "F"),
                          each = ncol(test_multi$coef_pasc))
  df_multi$subtype <- factor(df_multi$subtype)
  df_multi$subtype <- factor(df_multi$subtype,
                             levels = rev(levels(df_multi$subtype)))
  df_multi$name <- factor(df_multi$name)
  df_multi$name <- factor(df_multi$name, 
                          levels = c("Physical", "Cognitive", "Multiple"))
  
  df_multi %>% ggplot(mapping = aes(x = name, y = subtype, fill = value)) +
    geom_tile(size=0.2) +
    scale_fill_gradient2(low = "#284AA0", mid = "white", 
                         high = "#BC4A4A", midpoint = 0,
                         name = expression("Coefficient"),
                         limits = c(-val_capped, val_capped),
                         breaks = c(-val_capped, 0, val_capped)) +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    ylab("") +
    xlab("") +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5),
          legend.position = "left",
          legend.title = element_text(size = 15, angle = 90, vjust = 1),
          legend.text = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 0.5),
          legend.spacing.y = unit(1, 'cm'),
          plot.title = element_text(size = 15, face = "plain"))+
    geom_text(aes(label = ifelse(pval <= -0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 9, vjust = 0.75, color = "black") +
    ggtitle(plot_title)
}

# Global test adjusted for demographics (multiple responses)
test_global_multi_pasc <- function(subtype, response, demo_data,
                                   subtype_remove = "D") {
  id = (!is.na(subtype)) & (!is.na(response)) & (!rowSums(is.na(demo_data))) &
    (subtype != subtype_remove)
  subtype = factor(as.character(subtype[id]))
  response = response[id]
  demo_data = demo_data[id,]
  nclass = length(unique(subtype))
  
  x_pasc = model.matrix(~subtype-1)
  global_pval = array(NA, nclass)
  names(global_pval) = c("A","B", "C", "E", "F")
  for (i in 1:nclass) {
    x_binary = x_pasc[ ,i]
    df_global = cbind(data.frame(subtype = x_binary),
                      response, demo_data)
    demo_names = colnames(demo_data)
    form1 = as.formula(paste0("subtype ~ response + ",
                              paste0(demo_names, collapse = "+")))
    full_model = glm(form1, family = binomial, data = df_global)
    form2 = as.formula(paste0("subtype ~ ",
                              paste0(demo_names, collapse = "+")))
    reduced_model = glm(form2, family = binomial, data = df_global)
    global_test = anova(reduced_model, full_model, test = "Chisq")
    global_pval[i] = global_test[["Pr(>Chi)"]][2]
  }
  return(global_pval)
}

# Heatmap for global test (multiple responses)
heatmap_multi_global <- function(test_global) {
  df_global <- data.frame("Global" = -log10(test_global))
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
          axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
          legend.position = "left",
          legend.title = element_text(size = 15, angle = 90),
          legend.text = element_text(size = 15, angle = 90),
          plot.margin = margin(b = 20)) +
    labs(size = expression(-log[10]*"("*italic("P")*"-value"*")")) +
    geom_text(aes(label = ifelse(abs(value) >= -log10(0.001), "***",
                                 ifelse(abs(value) >= -log10(0.01), "**",
                                        ifelse(abs(value) >= -log10(0.05), "*", "")))),
              size = 7, vjust = 0.78, color = "white") +
    scale_size(range = c(3, 12), limits = c(0, 1.6),
               breaks = c(0.4, 0.8, 1.2, 1.6), 
               labels = c("0.4", "0.8", "1.2", "1.6"))
}