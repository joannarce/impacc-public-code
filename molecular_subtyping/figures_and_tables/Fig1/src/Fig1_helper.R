### Helpers for Figure 1 ###

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

# Heatmap for ordinal regression
heatmap_ordinal <- function(test_ord) {
  df_ord <- as.data.frame(test_ord$est_ordinal)
  colnames(df_ord) <- c("TG")
  df_ord <- pivot_longer(df_ord, cols = everything())
  pval_ord <- as.data.frame(test_ord$p_ordinal)
  colnames(pval_ord) <- c("TG")
  pval_ord <- pivot_longer(pval_ord, cols = everything())
  df_ord$pval <- pval_ord$value
  df_ord$subtype <- rep(sapply(1:length(test_ord$p_ordinal),
                               function(i) paste0(LETTERS[i])),
                        each = 1)
  df_ord$subtype <- factor(df_ord$subtype)
  df_ord$subtype <- factor(df_ord$subtype,
                           levels = rev(levels(df_ord$subtype)))
  
  df_ord %>% ggplot(mapping = aes(x = name, y = subtype, fill = value)) +
    geom_tile(linewidth = 0.2) +
    scale_fill_gradient2(low = "#284AA0", mid = "white", 
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
          plot.margin = margin(b = 50, t = 10),
          legend.spacing.y = unit(0.5, 'cm')) +
    geom_text(aes(label = ifelse(pval <= 0.001, "***",
                                 ifelse(pval <= 0.01, "**",
                                        ifelse(pval <= 0.05, "*", "")))),
              size = 7, vjust = 0.72, color = "black")
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