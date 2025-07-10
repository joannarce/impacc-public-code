organized_data = readRDS("/scratch/data-integration/molecular_subtyping/Data Preprocessing/phases123_combined_jeremy_072023.rds")
base_obj = readRDS(file ="/scratch/data-integration/cleanedup/mod1_interpretation_base.rds")
ct_train = data.frame(
  "event_id" = rownames(base_obj$others$nasal_viralload$train$nasal_viralload_counts),
  "N1_CT" = base_obj$others$nasal_viralload$train$nasal_viralload_counts$N1_CT,
  "RP_CT" = base_obj$others$nasal_viralload$train$nasal_viralload_counts$RP_CT)
ct_test = data.frame(
  "event_id" = rownames(base_obj$others$nasal_viralload$test$nasal_viralload_counts),
  "N1_CT" = base_obj$others$nasal_viralload$test$nasal_viralload_counts$N1_CT,
  "RP_CT" = base_obj$others$nasal_viralload$test$nasal_viralload_counts$RP_CT)
auc_train = data.frame(
  "event_id" = rownames(base_obj$others$serum_rbd_abtiters$train$serum_rbd_abtiters_counts),
  "log2_AUC" = base_obj$others$serum_rbd_abtiters$train$serum_rbd_abtiters_counts$`AUC RBD IgG`)
auc_test = data.frame(
  "event_id" = rownames(base_obj$others$serum_rbd_abtiters$test$serum_rbd_abtiters_counts),
  "log2_AUC" = base_obj$others$serum_rbd_abtiters$test$serum_rbd_abtiters_counts$`AUC RBD IgG`)

ncol(organized_data$datasets$plasma_proteomics_targeted) +
  ncol(organized_data$datasets$plasma_proteomics_global_dda) +
  ncol(organized_data$datasets$serum_olink) +
  ncol(organized_data$datasets$plasma_metabolomics_global) +
  length(unique(c(colnames(organized_data$datasets$nasal_transcriptomics),
                  colnames(organized_data$datasets$pbmc_transcriptomics))))
  
sum(complete.cases(organized_data$datasets$plasma_proteomics_targeted)) +
  sum(complete.cases(organized_data$datasets$plasma_proteomics_global_dda)) +
  sum(complete.cases(organized_data$datasets$serum_olink)) +
  sum(complete.cases(organized_data$datasets$plasma_metabolomics_global)) +
  sum(complete.cases(organized_data$datasets$nasal_transcriptomics)) +
  sum(complete.cases(organized_data$datasets$pbmc_transcriptomics)) +
  sum(complete.cases(ct_train)) + # nasal viral load
  sum(complete.cases(ct_test)) +
  sum(complete.cases(auc_train)) + # antibody
  sum(complete.cases(auc_test)) +
  sum(complete.cases(base_obj$others$bld_cytof$train$bld_cytof_counts)) + # cytof
  sum(complete.cases(base_obj$others$bld_cytof$test$bld_cytof_counts))



