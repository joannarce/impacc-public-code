### Helpers for Extended Data 5 ###

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

# Quantile normalization
quantileNormalizeNoTies <- function(data) {
  # Calculate ranks for each column
  ranks <- apply(data, 2, rank, ties.method = "average")
  ranks_norm = ranks/(apply(ranks,2,max)+1/2)
  normalizedData = qnorm(ranks_norm)
  return(normalizedData)
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
