### Helpers for Extended Data 3 ###

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