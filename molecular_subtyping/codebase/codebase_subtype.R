### Codebase for Molecular Subtyping ###

# Arrange data to long format
make_data_long = function(factor_scores, sample_meta, visit_max, mean_removal = F){
  participant_ids_visit1 = unique(sample_meta$participant_id[sample_meta$event_type =="Visit 1"])
  #summarize
  sample_meta_participant_level =  unique(sample_meta[,c("participant_id", "trajectory_group", "sex", "discharge_dates", "death_date" ,"admit_age", "discretized_admit_age_quantile")])
  sample_meta_participant_level$discharge_dates[is.na(sample_meta_participant_level$discharge_dates)] = 999
  sample_meta_participant_level$discharge_dates[sample_meta_participant_level$discharge_dates==""] = 999
  sample_meta_participant_level$death_date[is.na(sample_meta_participant_level$death_date)] = 999
  sample_meta_participant_level$death_date[is.infinite(sample_meta_participant_level$death_date)] = 999
  sample_meta_participant_level$death_date[sample_meta_participant_level$death_date==""] = 999
  #get the max discharging date, min discharging date and discharging date
  sample_meta_participant_level$discharge_min_date = as.numeric(sapply( sample_meta_participant_level$discharge_dates,
                                                                        function(z){min(as.integer(unlist(strsplit(z, split=","))))}))
  sample_meta_participant_level$discharge_min_date[is.infinite(sample_meta_participant_level$discharge_min_date)] = 999
  sample_meta_participant_level$discharge_max_date = as.numeric(sapply( sample_meta_participant_level$discharge_dates,
                                                                        function(z){max(as.integer(unlist(strsplit(z, split=","))))}))
  sample_meta_participant_level$discharge_max_date[is.infinite(sample_meta_participant_level$discharge_max_date)] = 999
  sample_meta_participant_level$discharge_max_acute_date = sapply( sample_meta_participant_level$discharge_dates,
                                                                   function(z){
                                                                     ll = unlist(strsplit(z, split=","))
                                                                     ll = as.integer(ll)
                                                                     ll[is.infinite(ll)] = 999
                                                                     ll1 = ll[ll<=42]
                                                                     if(length(ll1) > 0){
                                                                       max(ll1)
                                                                     }else{
                                                                       return(999)
                                                                     }
                                                                   })
  sample_meta_participant_level$discharge_dates = NULL
  sample_meta_participant_level <-  sample_meta_participant_level %>% arrange(participant_id)%>%
    group_by(participant_id) %>%
    mutate(discharge_min_date = min(discharge_min_date)) %>%
    mutate(discharge_max_date = max(discharge_max_date)) %>%
    mutate(discharge_max_acute_date = min(discharge_max_acute_date))
  sample_meta_participant_level = data.frame(sample_meta_participant_level)
  sample_meta_participant_level = unique(sample_meta_participant_level)
  
  rownames(sample_meta_participant_level) = sample_meta_participant_level$participant_id
  sample_meta_participant_level = sample_meta_participant_level[participant_ids_visit1,]
  X = data.frame(matrix(NA, nrow = length(participant_ids_visit1), ncol =  visit_max*ncol(factor_scores)))
  rownames(X) =  participant_ids_visit1
  for(v in 1:visit_max){
    ll = which(sample_meta$event_type ==paste0("Visit ",v))
    sub_factors = factor_scores[ll,]
    sub_sample_meta = sample_meta[ll,]
    rownames(sub_factors) =sub_sample_meta$participant_id
    ll1 = intersect(rownames(sub_factors),rownames(X))
    ll2 = (ncol(factor_scores)*(v-1)+1):(ncol(factor_scores)*v)
    X[ll1,ll2] =  sub_factors[ll1,]
  }
  if(mean_removal){
    for(j in 1:ncol(X)){
      ll = !is.na(X[,j])
      X[ll,j] = X[ll,j] - mean(X[ll,j])
    }
  }
  return(list(X = X, sample_meta_participant_level=sample_meta_participant_level))
  ###back trimming until 
}

# Distance matrix calculation
calculate_distance_matrix <- function(X, n_factor,visit_max){
  n_participant <- nrow(X)
  dist_mat <- matrix(0, nrow = n_participant, ncol = n_participant)
  count_mat <- matrix(0, nrow = n_participant, ncol = n_participant)
  for (k in 1:visit_max) {
    start_col <- (k - 1) * n_factor + 1
    end_col <- k * n_factor
    subset_X <- X[, start_col:end_col]
    visit_dist <- as.matrix(dist(subset_X, method = "euclidean"))
    visit_count <- ifelse(is.na(visit_dist), 0, 1)
    visit_dist[is.na(visit_dist)] <- 0
    dist_mat <- dist_mat + visit_dist
    count_mat <- count_mat + visit_count
  }
  dist_mat <- dist_mat / count_mat
  dist_mat <- as.dist(dist_mat)
  return(dist_mat)
}

# Pairwise survival comparison
pairwise_survival_comparisons=function(surv_object, subtype){
  subtype0 = sort(unique(subtype))
  pvals = matrix(NA, length(subtype0), length(subtype0))
  for(i in 1:(length(subtype0)-1)){
    for(j in (i+1):length(subtype0)){
      data_tmp = data.frame(subtype = factor(subtype[subtype%in%subtype0[c(i,j)]]))
      surv_object1 = surv_object[subtype%in%subtype0[c(i,j)]]
      tmp1 <- survdiff( surv_object1 ~ subtype, data =    data_tmp)
      pvals[i,j] = pchisq(tmp1$chisq, 1, lower.tail = F)
      pvals[j,i] = pvals[i,j]
    }
  }
  colnames(pvals) <-rownames(pvals) <- paste0("subtype",1:length( subtype0))
  return(pvals)
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

# Heatmap for ordinal regression
heatmap_ordinal <- function(test_ord) {
  df_ord <- as.data.frame(test_ord$est_ordinal)
  colnames(df_ord) <- c("TG")
  df_ord <- pivot_longer(df_ord, cols = everything())
  pval_ord <- as.data.frame(test_ord$p_ordinal)
  colnames(pval_ord) <- c("TG")
  pval_ord <- pivot_longer(pval_ord, cols = everything())
  df_ord$pval <- pval_ord$value
  df_ord$subtype <- gsub("Subtype", "", names(test_ord$est_ordinal))
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
          legend.title = element_text(size = 15, angle = 90, vjust = 0.8),
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

# Linear regression (one vs other comparison)
test_clinical_linear <- function(subtype, response) {
  id = (!is.na(subtype)) & (!is.na(response))
  subtype = factor(subtype[id])
  response = response[id]
  nclass = length(unique(subtype))
  x_linear = model.matrix(~subtype-1)
  x_linear = data.frame(x_linear[,1:(nclass-1)]-x_linear[,nclass])
  colnames(x_linear) = paste0("Subtype", levels(subtype)[1:(nclass-1)])
  df_linear = cbind(data.frame(response = response),x_linear)
  model1 <- lm(response ~ ., data = df_linear)
  coef_linear <- tail(summary(model1)$coefficients, nclass-1)
  est_linear <- coef_linear[,1]
  p_linear <- coef_linear[,4]
  
  x_linear = model.matrix(~subtype-1)
  x_linear = data.frame(x_linear[,2:(nclass)]-x_linear[,1])
  colnames(x_linear) = paste0("Subtype", levels(subtype)[2:nclass])
  df_linear = cbind(data.frame(response = response),x_linear)
  model1 <- lm(response ~ ., data = df_linear)
  coef_linear1 <- tail(summary(model1)$coefficients, 1)
  est_linear1 <- coef_linear1[,1]
  p_linear1 <- coef_linear1[,4]
  est_linear = append(est_linear, est_linear1)
  names(est_linear) = paste0("Subtype", levels(subtype)[1:nclass])
  p_linear = append(p_linear, p_linear1)
  names(p_linear) = paste0("Subtype", levels(subtype)[1:nclass])
  
  return(list(est_linear = est_linear, p_linear = p_linear))
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
heatmap_clinical <- function(test_logit, clin_names,
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
  df_heatmap$subtype <- rep(
    sapply(gsub("Subtype", "", rownames(test_logit$est_logit)),
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
          legend.title = element_text(size = 18, vjust = 0.8),
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
           use_endpoint_levels, CI_interval = 0.95, 
           max_date = 30, min_date = 0,
           subgroup_cr = c("A" = "#2ca02c", "B" = "#1f77b4", "C" = "#ff7f0e",
                           "D" = "#8c564b", "E" = "#9C3418", "F" = "#9467bd"), 
           legend_name = "Subtype", feature_select = FALSE, 
           limits = c(-Inf, Inf), feature_with_ylab = NULL,
           feature_with_xlab = NULL, multiple_correction = TRUE,
           xlabel = "Days from admission", text_position = 3) {
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
                                     min_date = min_date,
                                     colors = subgroup_cr,
                                     limits = limits,
                                     feature_with_ylab = feature_with_ylab,
                                     feature_with_xlab = feature_with_xlab,
                                     xlabel = xlabel,
                                     text_position = text_position)
    return(list(test_pval_list = test_pval_list, test_fig_list = test_fig_list,
                particular_interest_features = particular_interest_features))
  }

# Longitudinal trajectory plot
trajectory_plot <- function(inputDF, model_DF, endpoint, feature_names,
                            legend_name, legend_label, endpoint_order = NULL,
                            max_date = 28, min_date = 0, 
                            xlabel = "Days from admission",
                            CI_interval = 0.95, group_trendline_dropout = TRUE,
                            colors = c("#639A21", "#39828C",
                                       "#6371AD", "#BD7D31", "#9C3418"),
                            limits = c(-Inf, Inf), feature_with_ylab = NULL,
                            feature_with_xlab = NULL, text_position) {
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
      xlim(c(min_date, max_date)) + # max_date+1
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
      annotate("text", x = min_date+(max_date-min_date)/text_position,
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

# Serum Olink Core Distance Preprocessing
serum_olink_core_dist_preprocess <- function(serum_olink, so_idx, visit_max) {
  event_type <- serum_olink$event_type
  participant_id <- serum_olink$participant_id
  unique_id <- unique(participant_id)
  n_participant <- length(unique(unique_id))
  serum_olink <- serum_olink[, so_idx]
  serum_olink_core <- matrix(NA, nrow = ncol(serum_olink),
                             ncol = visit_max*n_participant)
  rownames(serum_olink_core) <- colnames(serum_olink)
  for(v in 1:visit_max){
    ll <- which(event_type == paste0("V",v))
    sub_so <- serum_olink[ll,]
    rownames(sub_so) <- participant_id[ll]
    ll1 <- which(unique_id %in% rownames(sub_so))
    ll2 <- (n_participant*(v-1)+1):(n_participant*v)
    ll2 <- ll2[ll1]
    serum_olink_core[,ll2] <- t(sub_so)
  }
  return(serum_olink_core)
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
                        random = ~(1|enrollment_site/participant_id))
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

# Pathway shortlist selection
pathway_shortlist_selection = function(aggregated_mat, data_base_names, size_max, KEGG_metabo_max){
  idx_short = c()
  for(database_id in 1:length( data_base_names)){
    idx = which(aggregated_mat$database == data_base_names[database_id])
    if(length(idx) > 0){
      idx1 = idx[1:min(length(idx),size_max[database_id])]
    }else{
      idx1 = c()
    }
    if(data_base_names[database_id] == "kegg"){
      idx = which((aggregated_mat$database == data_base_names[database_id]) & aggregated_mat$pathway %in% kegg_metabopathnames)
      if(length(idx) > 0){
        idx1_B = idx[1:min(length(idx),KEGG_metabo_max)]
        idx1 = unique(c(idx1, idx1_B))
      }
    }
    idx_short = c(idx_short, idx1)
  }
  aggregated_mat_short = aggregated_mat[idx_short,]
  return(aggregated_mat_short)
}

# Enrichment aggregation
enrichment_aggregation_subtype = function(factor_enrichment, adj.p.thr = 0.05, negation = 1, thr = .2, type = "mhg"){
  tmp0 =list()
  for(database_id in 1:length(factor_enrichment)){
    tmp0[[database_id]] <- try(enrichment_aggregation1(pos_list=factor_enrichment[[database_id]][["pos"]], neg_list =factor_enrichment[[database_id]][["neg"]], thr =thr, type = type))
    if(class(tmp0[[database_id]]) == "try-error"){
      tmp0[[database_id]] <- try(enrichment_aggregation1(pos_list=factor_enrichment[[database_id]][["pos"]], thr = thr, type =type))
      if(class(tmp0[[database_id]])!="try-error"){
        tmp0[[database_id]]$direction_mat[,] = array("+", dim=dim(tmp0[[database_id]]$pval_mat))
      }
      if(class(tmp0[[database_id]]) == "try-error"){
        tmp0[[database_id]] <- try(enrichment_aggregation1(pos_list=factor_enrichment[[database_id]][["neg"]], thr = thr, type =type))
        if(class(tmp0[[database_id]])!="try-error"){
          tmp0[[database_id]]$direction_mat[,] = array("-", dim=dim(tmp0[[database_id]]$pval_mat))
        }
      }
    }
  }
  pos_value=negation
  database_names = names(factor_enrichment)
  tmp = list()
  for(database_id in 1:length(factor_enrichment)){
    if(class(tmp0[[database_id]])!="try-error"){
      aa =  tmp0[[database_id]]$direction_mat[!is.na(tmp0[[database_id]]$direction_mat)]
      idx = which(!is.na(tmp0[[database_id]]$direction_mat))
      ww=as.matrix(tmp0[[database_id]]$direction_mat)
      ww[idx] = ifelse(aa=="+", pos_value, -pos_value)
      tmp0[[database_id]]$direction_mat[,] =ww 
      tmp0[[database_id]]$direction_mat=apply(tmp0[[database_id]]$direction_mat,c(1,2), as.numeric)
      tmp0[[database_id]]$pval_mat = -log(tmp0[[database_id]]$pval_mat, base = 10)
      tmp0[[database_id]]$pval_mat[,!(colnames(tmp0[[database_id]]$pval_mat)%in%c("joint", "joint_adj"))] = tmp0[[database_id]]$pval_mat[,!(colnames(tmp0[[database_id]]$pval_mat)%in%c("joint", "joint_adj"))]*sign(tmp0[[database_id]]$direction_mat)
      
      tmp[[database_id]] = data.frame(matrix(NA, nrow = nrow(tmp0[[database_id]]$pval_mat), ncol = length(omic_short_names)+2))
      colnames(tmp[[database_id]]) = c(omic_short_names,c("joint", "joint_adj"))
      tmp[[database_id]][,colnames(tmp0[[database_id]]$pval_mat)] = tmp0[[database_id]]$pval_mat
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
  
  df_fisher$subtype <- rep(
    rownames(full_clust_test[["PASC"]][["fisher"]][["fisher_table_rest"]]),
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
          legend.title = element_text(size = 15, angle = 90, vjust = 0.8),
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
  df_global$subtype <- gsub(
    "Subtype", "",
    rownames(full_clust_test[["PASC"]][["fisher"]][["fisher_table_row"]]))
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
          legend.title = element_text(size = 15, angle = 90, vjust = 0.8),
          legend.text = element_text(size = 15, angle = 90),
          plot.margin = margin(b = 20)) +
    labs(size = expression(-log[10]*"("*italic("P")*"-value"*")")) +
    geom_text(aes(label = ifelse(abs(value) >= -log10(0.001), "***",
                                 ifelse(abs(value) >= -log10(0.01), "**",
                                        ifelse(abs(value) >= -log10(0.05), "*", "")))),
              size = 7, vjust = 0.78, color = "white") +
    scale_size(range = c(3, 12))
}

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

# t-test compared to mean
cytof_t_test <- function(cytof_mat, idx, min_sample_size = 5, 
                         unique_group, unique_visit) {
  n_group = length(unique_group)
  n_visit = length(unique_visit)
  p_value = matrix(NA, nrow = length(idx), ncol = n_group*n_visit)
  p_adjust = matrix(NA, nrow = length(idx), ncol = n_group*n_visit)
  t_score = matrix(NA, nrow = length(idx), ncol = n_group*n_visit)
  effect_size = matrix(NA, nrow = length(idx), ncol = n_group*n_visit)
  for (i in idx) {
    for (g in 1:n_group) {
      d1 = cytof_mat[cytof_mat$label == unique_group[g],]
      d2 = cytof_mat[cytof_mat$label != unique_group[g],]
      for (v in 1:n_visit) {
        s1 = d1[d1$event_type == unique_visit[v],]
        s2 = d2[d2$event_type == unique_visit[v],]
        if (sum(!is.na(s1[,i])) >= min_sample_size && 
            sum(!is.na(s2[,i])) >= min_sample_size) {
          t_test = t.test(s1[,i], s2[,i])
          # remove sample size effect
          n1 <- sum(!is.na(s1[,i]))
          n2 <- sum(!is.na(s2[,i]))
          n_effect <- n1*n2/(n1+n2)
          p_value[i, (n_visit*(g-1)+v)] = t_test$p.value
          t_score[i, (n_visit*(g-1)+v)] = t_test$statistic
          effect_size[i, (n_visit*(g-1)+v)] = t_test$statistic/sqrt(n_effect)
        }
      }
    }
  }
  comb = as.vector(t(outer(unique_group, unique_visit, FUN = paste, sep = "_")))
  rownames(p_value) = colnames(cytof_mat)[idx]
  colnames(p_value) = comb
  rownames(t_score) = colnames(cytof_mat)[idx]
  colnames(t_score) = comb
  rownames(effect_size) = colnames(cytof_mat)[idx]
  colnames(effect_size) = comb
  p_adjust = p_value
  p_adjust = apply(p_adjust, 2, function(p) p.adjust(p, method = "BH"))
  return(list(p_value = p_value, p_adjust = p_adjust, 
              t_score = t_score, effect_size = effect_size))
}

# Baseline comparison (group-wise comparison)
baseline_comparison <- function(anc_norm, demo_data, subtype_label, 
                                subtype_remove = "D") {
  if (!is.null(subtype_remove)) {
    idx = which(subtype_label != subtype_remove)
    anc_norm = anc_norm[idx, ]
    demo_data = demo_data[idx, ]
    subtype_label = subtype_label[idx]
  }
  subtype_label = factor(subtype_label)
  
  est_lm = matrix(NA, nrow = 5, ncol = ncol(anc_norm))
  rownames(est_lm) = c("EF vs ABC", "F vs E",
                       "A vs BC", "B vs AC", "C vs AB")
  colnames(est_lm) = colnames(anc_norm)
  p_lm = matrix(NA, nrow = 5, ncol = ncol(anc_norm))
  rownames(p_lm) = c("EF vs ABC", "F vs E",
                     "A vs BC", "B vs AC", "C vs AB")
  colnames(p_lm) = colnames(anc_norm)
  
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
  
  for (i in 1:ncol(anc_norm)) {
    response = anc_norm[[i]]
    df_lm = cbind(data.frame(response = response), demo_data, x1_clin)
    model1 = lm(response ~ ., data = df_lm)
    coef_lm = summary(model1)$coef[c("EFvsABC", "FvsE", "AvsBC", "BvsAC"),]
    est_lm[1:4,i] <- coef_lm[,1]
    p_lm[1:4,i] <- coef_lm[,4]
    
    df_lm = cbind(data.frame(response = response), demo_data, x2_clin)
    model2 = lm(response ~ ., data = df_lm)
    est_lm[5,i] <- summary(model2)$coef["CvsAB",1]
    p_lm[5,i] <- summary(model2)$coef["CvsAB",4]
  }
  return(list(est_lm = est_lm, p_lm = p_lm))
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

# Assess interaction between subtype and factor/pathway
factor_subtype_interaction <- function(datDF0, fname0, sig_comorb) {
  subtype_label = factor(datDF0$label)
  x1 = model.matrix(~subtype_label)
  colnames(x1) = c("Intercept", "EFvsABC", "FvsE", "AvsBC", "BvsAC")
  x1[, "EFvsABC"] = -1
  x1[subtype_label == "E" | subtype_label == "F", "EFvsABC"] = 1
  x1[, "FvsE"] = 0
  x1[subtype_label == "E", "FvsE"] = -1
  x1[subtype_label == "F", "FvsE"] = 1
  x1[, "AvsBC"] = 0
  x1[subtype_label == "A", "AvsBC"] = 1
  x1[subtype_label == "C", "AvsBC"] = -1
  x1[, "BvsAC"] = 0
  x1[subtype_label == "B", "BvsAC"] = 1
  x1[subtype_label == "C", "BvsAC"] = -1
  
  x2 = model.matrix(~subtype_label)
  colnames(x2) = c("Intercept", "EFvsABC", "FvsE", "BvsAC", "CvsAB")
  x2[, "EFvsABC"] = -1
  x2[subtype_label == "E" | subtype_label == "F", "EFvsABC"] = 1
  x2[, "FvsE"] = 0
  x2[subtype_label == "E", "FvsE"] = -1
  x2[subtype_label == "F", "FvsE"] = 1
  x2[, "BvsAC"] = 0
  x2[subtype_label == "B", "BvsAC"] = 1
  x2[subtype_label == "A", "BvsAC"] = -1
  x2[, "CvsAB"] = 0
  x2[subtype_label == "C", "CvsAB"] = 1
  x2[subtype_label == "A", "CvsAB"] = -1
  
  # vanilla (x1)
  dataDF = cbind(datDF0, x1)
  formula0 = paste0(fname0," ~ s(event_date) + sex + discretized_admit_age_quantile")
  formula0 = formula(paste0(formula0, "+ EFvsABC + FvsE + AvsBC + BvsAC"))
  fit10 = gamm4::gamm4(formula0, random = ~(1|participant_id/enrollment_site),
                       data = dataDF)
  coef_fac0 = summary(fit10$gam)$p.table[c("EFvsABC", "FvsE", "AvsBC", "BvsAC"),
                                         "Estimate"]
  pval_fac0 = summary(fit10$gam)$p.table[c("EFvsABC", "FvsE", "AvsBC", "BvsAC"),
                                         "Pr(>|t|)"]
  
  # comorbidity adjusted (x1)
  formula1 = paste0(fname0," ~s(event_date) + sex+discretized_admit_age_quantile +", 
                    paste(sig_comorb, collapse = "+"))
  formula1 = formula(paste0(formula1,"+ EFvsABC + FvsE + AvsBC + BvsAC"))
  fit11 = gamm4::gamm4(formula1, random = ~(1|participant_id/enrollment_site),
                       data = dataDF)
  coef_fac1 = summary(fit11$gam)$p.table[c("EFvsABC", "FvsE", "AvsBC", "BvsAC"), 
                                         "Estimate"]
  pval_fac1 = summary(fit11$gam)$p.table[c("EFvsABC", "FvsE", "AvsBC", "BvsAC"), 
                                         "Pr(>|t|)"]
  
  # vanilla (x2)
  dataDF = cbind(datDF0, x2)
  formula0 = paste0(fname0," ~ s(event_date) + sex + discretized_admit_age_quantile")
  formula0 = formula(paste0(formula0,"+ EFvsABC + FvsE + BvsAC + CvsAB"))
  fit10 = gamm4::gamm4(formula0, random = ~(1|participant_id/enrollment_site),
                       data = dataDF)
  coef_fac0 = c(coef_fac0, summary(fit10$gam)$p.table["CvsAB", "Estimate"])
  pval_fac0 = c(pval_fac0, summary(fit10$gam)$p.table["CvsAB", "Pr(>|t|)"])
  
  # comorbidity adjusted (x2)
  formula1 = paste0(fname0,"~s(event_date)+sex+discretized_admit_age_quantile+", 
                    paste(sig_comorb, collapse = "+"))
  formula1 = formula(paste0(formula1,"+ EFvsABC + FvsE + BvsAC + CvsAB"))
  fit11 = gamm4::gamm4(formula1, random = ~(1|participant_id/enrollment_site),
                       data = dataDF)
  coef_fac1 = c(coef_fac1, summary(fit11$gam)$p.table["CvsAB", "Estimate"])
  pval_fac1 = c(pval_fac1, summary(fit11$gam)$p.table["CvsAB", "Pr(>|t|)"])
  
  # vanilla (residual)
  dataDF = cbind(datDF0, subtype_label)
  formula0 = paste0(fname0," ~ s(event_date) + sex + discretized_admit_age_quantile")
  formula0 = formula(paste0(formula0,"+ subtype_label"))
  fit_full = gamm4::gamm4(formula0, random = ~(1|participant_id/enrollment_site),
                          data = dataDF)
  coef_tmp = coef(fit_full$gam)
  new_data = model.matrix(~sex+discretized_admit_age_quantile, data = dataDF)[,-1]
  select_col = c("sexMale", names(coef_tmp)[grepl("^discretized", names(coef_tmp))])
  val = new_data %*% coef_tmp[select_col]
  res0 = dataDF[[fname0]] - val
  
  # comorbidity adjusted (residual)
  formula1 = paste0(fname0,"~s(event_date) + sex + discretized_admit_age_quantile + ", 
                    paste(sig_comorb, collapse = "+"))
  formula1 = formula(paste0(formula1,"+ subtype_label"))
  fit_full = gamm4::gamm4(formula1, random = ~(1|participant_id/enrollment_site),
                          data = dataDF)
  coef_comorb = summary(fit_full$gam)$p.table[sig_comorb, "Estimate"]
  pval_comorb = summary(fit_full$gam)$p.table[sig_comorb, "Pr(>|t|)"]
  coef_tmp = coef(fit_full$gam)
  formula2 = paste0("~ sex + discretized_admit_age_quantile + ", 
                    paste(sig_comorb, collapse = "+"))
  new_data = model.matrix(as.formula(formula2), data = dataDF)[,-1]
  select_col = c("sexMale", names(coef_tmp)[grepl("^discretized", names(coef_tmp))],
                 names(coef_tmp)[grepl("^comorb_", names(coef_tmp))])
  val = new_data %*% coef_tmp[select_col]
  res1 = dataDF[[fname0]] - val
  
  return(list(coef_fac0 = coef_fac0, pval_fac0 = pval_fac0,
              coef_fac1 = coef_fac1, pval_fac1 = pval_fac1, 
              res0 = res0, res1 = res1, label_use = subtype_label,
              coef_comorb = coef_comorb, pval_comorb = pval_comorb))
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
  df_heatmap$subtype <- c("A", "B", "C", "E", "F")
  
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