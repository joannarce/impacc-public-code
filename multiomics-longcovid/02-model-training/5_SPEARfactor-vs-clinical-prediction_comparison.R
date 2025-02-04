
# Install SPEAR package
#remotes::install_bitbucket("kleinstein/SPEAR@main")

# Install necessary packages
packages <- c("yardstick","tidyverse", "SPEAR", "glmnet", "mltools", "caret")

new_pkg <- packages[!(packages %in% installed.packages())]
if (length(new_pkg)>0) {install.packages(new_pkg)}

# list (and install if needed) Bioconductor packages
packages_bioconductor <- c("MultiAssayExperiment")
new_pkg_bioconductor <- packages_bioconductor[!(packages_bioconductor %in% installed.packages())]
if (length(new_pkg_bioconductor)>0) {BiocManager::install(new_pkg_bioconductor)}

# load packages
# load packages
packages <- sort(append(packages, packages_bioconductor)) # so the list of versions will be in order
for (n in seq_along(packages)) {
 suppressPackageStartupMessages(library(packages[n], character.only = TRUE))
 cat(paste0(packages[n], ": ", packageVersion(packages[n]), "\n")) # print simplified package versions
}

select <- dplyr::select
predict <- stats::predict


# Function for CV model ---------------

# CV scores lasso model (glmnet)

# fold_ids: sample ids (should match rownames) with fold numbers assigned to all the rows in X and Y.
# X feature matrix (samples in rows, features in columns)
# Y true class matrix (samples in rows, class in column, encoded as numbers)
# Y_full class matrix (samples in rows, class in column, as strings)

cv_scores_lasso_model_min_other <- function( X,
                                             model_name = "model1",
                                             Y=physpro_groups_encoded, 
                                             Y_full=physpro_groups, 
                                             n_folds = 10, 
                                             bootstrap_sample = 1
                                             ){  
  if (dim(X)[1] != length(Y)){
    stop("X matrix first dimension does not match Y matrix first dimension")
  }
  if (!length(Y)==length(Y_full)){
    stop("Y matrix dimensions do not match Y_full matrix dimensions")
  }
  
  random_seeds = seq(1,bootstrap_sample)
  
  print(paste0("Training model ",model_name))
  
  model.preds.te <- data.frame()
  imp_all_fold <- data.frame()
  res.df <- data.frame()
  
  for (random_seed in random_seeds){
    print(paste0("Bootstrap ",random_seed," of ",bootstrap_sample))
    set.seed(random_seed)
    fold_ids <- createFolds(as.factor(Y_full), k = n_folds, list = F, returnTrain = FALSE)

    for(f in 1:n_folds){
      # Get ids
      train_ids <- f != fold_ids
      test_ids <- f == fold_ids
      f.ids.tmp <- fold_ids[train_ids]
      f.ids.tmp[f.ids.tmp > f] = f.ids.tmp[f.ids.tmp > f] - 1
      # Subset:
      if (dim(X)[2] == 1) {
        X.tr <- cbind(X[train_ids,], X[train_ids,])
        X.te <- cbind(X[test_ids,], X[test_ids,])
      } else {
        X.tr <- X[train_ids,]
        X.te <- X[test_ids,]
      }
      
      Y.tr <- Y[train_ids]
      Y.te <- Y[test_ids]
      Y.tr.full <- Y_full[train_ids]
      Y.te.full <- Y_full[test_ids]
      # Train
      lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = f.ids.tmp, family = "multinomial")
      
      #Get feature importance
      imp_fold <- caret::getModelInfo("glmnet")$glmnet$varImp(lasso_model, lambda = "lambda.min")
      imp_all_fold <- rbind(imp_all_fold, data.frame(
        Model = model_name,
        Seed = random_seed,
        Fold = f,
        features = rownames(imp_fold),
        coef = imp_fold[,1]
      ))
      
      # Evaluate
      lasso.preds.tr = predict(lasso_model, X.tr, s = "lambda.min")[,,1]
      lasso.preds.te = predict(lasso_model, X.te, s = "lambda.min")[,,1]
      
      # Convert to df:
      lasso.preds.tr <- as.data.frame(lasso.preds.tr)
      colnames(lasso.preds.tr) <- c("MIN", "OTHER")
      lasso.preds.tr$estimate = c("MIN", "OTHER")[apply(lasso.preds.tr, 1, which.max)]
      lasso.preds.tr$truth = Y.tr.full
      lasso.preds.tr <- lasso.preds.tr %>% select(truth, estimate, everything())
      lasso.preds.te <- as.data.frame(lasso.preds.te)
      colnames(lasso.preds.te) <- c("MIN", "OTHER")
      lasso.preds.te$estimate = c("MIN", "OTHER")[apply(lasso.preds.te, 1, which.max)]
      lasso.preds.te$truth = Y.te.full
      lasso.preds.te <- lasso.preds.te %>% select(truth, estimate, everything())
      lasso.preds.te$Fold <- f
      lasso.preds.te$Seed <- random_seed
      lasso.preds.te$Model <- model_name
      model.preds.te <- rbind(model.preds.te, lasso.preds.te)
      # Assess:
      in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
      cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate
      
      res.df <- rbind(res.df, data.frame(
        Model = model_name,
        Seed = random_seed,
        Fold = f,
        Train = in.sample,
        CV = cv.sample
      ))
    }

  }
    
  # Compute mean coef of features over folds
  imp_all_fold_mean <- imp_all_fold %>%
    group_by(Model, Seed, features) %>%
    summarize(mean_coef = mean(coef)) %>%
    arrange(desc(abs(mean_coef)))
    
  return(list(resdf = res.df, testscores = model.preds.te, var_importance_cv_mean = imp_all_fold_mean, var_importance_all = imp_all_fold))

}

cv_scores_lasso_model_regression <- function( X,
                                             model_name = "model1",
                                             Y=physpro_groups_encoded, 
                                             n_folds = 10, 
                                             fold_ids = fphys.ids){  
  if (dim(X)[1] != length(Y)){
    stop("X matrix first dimension does not match Y matrix first dimension")
  }
  if (dim(X)[1] != length(fold_ids)){
    stop("X matrix first dimension does not match fold_ids matrix first dimension")
  }
  
  model.preds.te <- data.frame()
  imp_all_fold <- data.frame()
  res.df <- data.frame()
  for(f in 1:n_folds){
    print(paste0("Fold", f, "..."))
    # Get ids
    train_ids <- f != fold_ids
    test_ids <- f == fold_ids
    f.ids.tmp <- fold_ids[train_ids]
    f.ids.tmp[f.ids.tmp > f] = f.ids.tmp[f.ids.tmp > f] - 1
    # Subset:
    if (dim(X)[2] == 1) {
      X.tr <- cbind(X[train_ids,], X[train_ids,])
      X.te <- cbind(X[test_ids,], X[test_ids,])
    } else {
      X.tr <- X[train_ids,]
      X.te <- X[test_ids,]
    }
    
    Y.tr <- Y[train_ids]
    Y.te <- Y[test_ids]

    # Train
    lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = f.ids.tmp)
    
    #Get feature importance
    imp_fold <- caret::getModelInfo("glmnet")$glmnet$varImp(lasso_model, lambda = "lambda.min")
    imp_all_fold <- rbind(imp_all_fold, data.frame(
      Model = model_name,
      Fold = f,
      features = rownames(imp_fold),
      coef = imp_fold[,1]
    ))
    
    # Evaluate
    lasso.preds.tr = predict(lasso_model, X.tr, s = "lambda.min")
    lasso.preds.te = predict(lasso_model, X.te, s = "lambda.min")
    
    # Convert to df:
    lasso.preds.tr <- as.data.frame(lasso.preds.tr)
    lasso.preds.tr$truth <- Y.tr
    colnames(lasso.preds.tr) <- c("estimate", "truth")
    lasso.preds.tr <- lasso.preds.tr %>% select(truth, estimate)
    lasso.preds.te <- as.data.frame(lasso.preds.te)
    lasso.preds.te$truth <- Y.te
    colnames(lasso.preds.te) <- c("estimate", "truth")
    lasso.preds.te <- lasso.preds.te %>% select(truth, estimate)
    model.preds.te <- rbind(model.preds.te, lasso.preds.te)
    
    # Assess:
    in.sample <- yardstick::rmse(data = lasso.preds.tr, truth, estimate)$.estimate
    cv.sample <- yardstick::rmse(data = lasso.preds.te, truth, estimate)$.estimate
    
    res.df <- rbind(res.df, data.frame(
      Model = model_name,
      Fold = f,
      Train = in.sample,
      CV = cv.sample
    ))
  }
  model.preds.te$Model <- model_name
  
  # Compute mean coef of features over folds
  imp_all_fold_mean <- imp_all_fold %>%
    group_by(Model, features) %>%
    summarize(mean_coef = mean(coef)) %>% #geometric mean
    arrange(desc(abs(mean_coef)))
  
  return(list(resdf = res.df, testscores = model.preds.te, var_importance_cv_mean = imp_all_fold_mean, var_importance_all = imp_all_fold))
  
}

# CV Data: AUROC for all Models

# Loading MOFA and SPEAR objects...

# # Load MOFA object
# MOFAobj <- readRDS("/scratch/data-fullcohort/MOFA_factors/MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.rds")
# MOFA.fs <- MOFA2::get_factors(MOFAobj)$group1
# 
# # Load SPEAR binomial object
# SPEARobj <- readRDS("/scratch/data-fullcohort/spear_model/results/020624/021624_binomial_SPEAR.rds")
# SPEARobj$set.weights(method = "min")
# SPEAR.fs <- SPEARobj$get.factor.scores()
# SPEAR.cv.fs <- SPEARobj$get.factor.scores(cv = TRUE)
# f.ids <- SPEARobj$params$fold.ids
# pro_groups <- SPEARobj$data$train$pro_group_labels
# pro_groups_encoded <- ifelse(SPEARobj$data$train$pro_group_labels == "MIN", 0, 1)
# 
# # Load SPEAR Physical object
# SPEARobj.phys <- readRDS("/scratch/data-fullcohort/spear_model/results/020624/021924_physical_gaussianNEW_SPEAR.rds")
# SPEARobj.phys$set.weights(method = "min")
# SPEARphys.fs <- SPEARobj.phys$get.factor.scores()
# SPEARphys.cv.fs <- SPEARobj.phys$get.factor.scores(cv = TRUE)
# fphys.ids <- SPEARobj.phys$params$fold.ids
# physpro_groups <- factor(ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", "MIN", "OTHER"), levels = c("MIN", "OTHER"))
# physpro_groups_encoded <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", 0, 1)
# tmp <- SPEARobj.phys$data$test@ExperimentList$PPG
# tmp[is.na(tmp)] <- 0
# SPEARobj.phys$data$test@ExperimentList$PPG <- tmp


# Prepare metadata for clinical model based on integration manuscript ------

# Load clinical data
data_env_train <- readRDS("/scratch/data-fullcohort/batch-effects/impacc_convalescent_data_env_train_Feb_11_2024.RDS")
clinical_data_train <- data_env_train$clinical_data
data_env_test <- readRDS("/scratch/data-fullcohort/batch-effects/impacc_convalescent_data_env_test_Feb_11_2024.RDS")
clinical_data_test <- data_env_test$clinical_data

# # Clinical model integration --------------
# clinical_vars <- c("sex","race", "ethnicity", "discretized_admit_age_quantile", "bmi")
# include_vars_comorb <- c("comorb_htn", "comorb_isaric_dm", "comorb_anyresp_noasthma", "comorb_isaric_asthma", "comorb_isaric_cardiac", 
#                          "comorb_isaric_ckd", "comorb_isaric_neoplasm", "comorb_neuro", "comorb_isaric_liver", "comorb_hxtrans", 
#                          "comorb_hiv2", "comorb_eversmkvape", "comorb_anysubstance", "comorb_count", "comorb_countcat")
# include_vars_baseline <- c("baseline_img_infil", "baseline_lymph_abn", "baseline_platelets_abn", "baseline_alt_abn", "baseline_cr_abn",    
#                            "baseline_crp_abn", "baseline_ddimer_abn", "baseline_trop_abn" )
# all_clinical_vars <- c("event_id", clinical_vars, include_vars_comorb, include_vars_baseline)
# 
# # Clinical data train
# clinical_vars_train <- clinical_data_train %>%
#   dplyr::select(all_of(all_clinical_vars)) %>%
#   distinct()
# 
# # Add mean BMI for missing values, set missing baselines to 0.
# clinical_vars_train$bmi[is.na(clinical_vars_train$bmi)]<-mean(clinical_vars_train$bmi,na.rm=TRUE)
# clinical_vars_train$baseline_img_infil[is.na(clinical_vars_train$baseline_img_infil)] = 0
# clinical_vars_train$baseline_img_infil[clinical_vars_train$baseline_img_infil==999]<-0
# 
# # One-hot encoding multi-variate data
# temp <- data.frame(clinical_vars_train$race)
# colnames(temp) <-c("race")
# temp$race<-as.factor(temp$race)
# newtemp <- one_hot(as.data.table(temp), dropUnusedLevels = F)
# clinical_vars_train <- subset(clinical_vars_train, select = -c(race) )
# clinical_vars_train <- cbind(clinical_vars_train,newtemp)
# 
# temp <- data.frame(clinical_vars_train$ethnicity)
# colnames(temp)<-c("ethnicity")
# temp$ethnicity<-as.factor(temp$ethnicity)
# newtemp <- one_hot(as.data.table(temp))
# clinical_vars_train <- subset(clinical_vars_train, select = -c(ethnicity) )
# clinical_vars_train <- cbind(clinical_vars_train, newtemp)
# 
# temp <- data.frame(clinical_vars_train$discretized_admit_age_quantile)
# colnames(temp)<-c("discretized_admit_age_quantile")
# temp$discretized_admit_age_quantile<-as.factor(temp$discretized_admit_age_quantile)
# newtemp <- one_hot(as.data.table(temp))
# clinical_vars_train <- subset(clinical_vars_train, select = -c(discretized_admit_age_quantile) )
# clinical_vars_train <-cbind(clinical_vars_train, newtemp)
# 
# clinical_vars_train$sex[clinical_vars_train$sex=="Male"]=0
# clinical_vars_train$sex[clinical_vars_train$sex=="Female"]=1
# 
# # Move event_id to rownames and standardize
# row_names <- clinical_vars_train$event_id
# clinical_vars_train$event_id <- NULL
# clinical_vars_train <- as.data.frame(sapply(clinical_vars_train, as.numeric))
# clinical_vars_train <- scale(clinical_vars_train)
# rownames(clinical_vars_train) <- row_names
# 
# # Order in same order as SPEAR object
# clinical_vars_ordered <- clinical_vars_train[rownames(SPEARphys.cv.fs),]
# 
# # Clinical data test
# clinical_vars_test <- clinical_data_test %>%
#   dplyr::select(all_of(all_clinical_vars)) %>%
#   distinct()
# 
# # Add mean BMI for missing values, set missing baselines to 0.
# clinical_vars_test$bmi[is.na(clinical_vars_test$bmi)]<-mean(clinical_vars_test$bmi,na.rm=TRUE)
# clinical_vars_test$baseline_img_infil[is.na(clinical_vars_test$baseline_img_infil)] = 0
# clinical_vars_test$baseline_img_infil[clinical_vars_test$baseline_img_infil==999]<-0
# 
# # One-hot encoding multi-variate data
# temp <- data.frame(clinical_vars_test$race)
# colnames(temp) <-c("race")
# temp$race<- factor(temp$race, levels = unique(clinical_data_train$race))
# newtemp <- one_hot(as.data.table(temp), dropUnusedLevels = F)
# clinical_vars_test <- subset(clinical_vars_test, select = -c(race) )
# clinical_vars_test <- cbind(clinical_vars_test,newtemp)
# 
# temp <- data.frame(clinical_vars_test$ethnicity)
# colnames(temp) <- c("ethnicity")
# temp$ethnicity <- factor(temp$ethnicity, levels = unique(clinical_data_train$ethnicity))
# newtemp <- one_hot(as.data.table(temp), dropUnusedLevels = F)
# clinical_vars_test <- subset(clinical_vars_test, select = -c(ethnicity) )
# clinical_vars_test <- cbind(clinical_vars_test, newtemp)
# 
# temp <- data.frame(clinical_vars_test$discretized_admit_age_quantile)
# colnames(temp)<-c("discretized_admit_age_quantile")
# temp$discretized_admit_age_quantile<-factor(temp$discretized_admit_age_quantile, levels = unique(clinical_data_train$discretized_admit_age_quantile))
# newtemp <- one_hot(as.data.table(temp), dropUnusedLevels = F)
# clinical_vars_test <- subset(clinical_vars_test, select = -c(discretized_admit_age_quantile) )
# clinical_vars_test <-cbind(clinical_vars_test, newtemp)
# 
# clinical_vars_test$sex[clinical_vars_test$sex=="Male"]=0
# clinical_vars_test$sex[clinical_vars_test$sex=="Female"]=1
# 
# # Move event_id to rownames and standardize
# row_names <- clinical_vars_test$event_id
# clinical_vars_test$event_id <- NULL
# clinical_vars_test <- as.data.frame(sapply(clinical_vars_test, as.numeric))
# clinical_vars_test <- scale(clinical_vars_test)
# rownames(clinical_vars_test) <- row_names
# 
# # Order columns in same order as train data
# clinical_vars_test <- clinical_vars_test[,colnames(clinical_vars_ordered)]
# stopifnot(all(colnames(clinical_vars_ordered) == colnames(clinical_vars_test)))


# Prepare metadata for clinical model based on Naresh paper and subtyping manuscript ---

# Load clinical data (acute phase)
data_env_acute_train <- readRDS("/scratch/data-fullcohort/batch-effects/impacc_acute_data_env_train_Mar_23_2024.RDS")
clinical_data_acute_train <- data_env_acute_train$clinical_data
data_env_acute_test <- readRDS("/scratch/data-fullcohort/batch-effects/impacc_acute_data_env_test_Mar_23_2024.RDS")
clinical_data_acute_test <- data_env_acute_test$clinical_data

# Read assays that are included in clinical model
antibody_titers_acute_env <- readRDS("/scratch/data-fullcohort/batch-effects/preprocessed_assay_antibody_titers_acute_20240701.rds")
nasal_viral_load_acute <- readRDS("/scratch/data-fullcohort/batch-effects/nasal_viral_load_acute_20240701.rds")
  

# # Clinical model subtyping ------------
# clinical_vars_subtyping <- c("sex","discretized_admit_age_quantile", "bmi", "los", "baseline_sofa_score_new")
# include_vars_comorb_subtyping <- c("comorb_htn", "comorb_isaric_dm", "comorb_anyresp_noasthma", "comorb_isaric_asthma", "comorb_isaric_cardiac", 
#                          "comorb_isaric_ckd", "comorb_isaric_neoplasm", "comorb_neuro", "comorb_isaric_liver", "comorb_hxtrans", 
#                          "comorb_hiv2", "comorb_eversmkvape", "comorb_anysubstance", "comorb_count")
# 
# all_clinical_vars_subtyping <- c("event_id", clinical_vars_subtyping, include_vars_comorb_subtyping)
# 
# # Clinical data train
# clinical_vars_subtyping_train <- clinical_data_acute_train %>%
#   dplyr::select(all_of(all_clinical_vars_subtyping)) %>%
#   distinct()
# 
# # Add mean BMI for missing values, set missing baselines to 0.
# clinical_vars_subtyping_train$bmi[is.na(clinical_vars_subtyping_train$bmi)]<-mean(clinical_vars_subtyping_train$bmi,na.rm=TRUE)
# 
# temp <- data.frame(clinical_vars_subtyping_train$discretized_admit_age_quantile)
# colnames(temp)<-c("discretized_admit_age_quantile")
# temp$discretized_admit_age_quantile<-as.factor(temp$discretized_admit_age_quantile)
# newtemp <- one_hot(as.data.table(temp))
# clinical_vars_subtyping_train <- subset(clinical_vars_subtyping_train, select = -c(discretized_admit_age_quantile) )
# clinical_vars_subtyping_train <-cbind(clinical_vars_subtyping_train, newtemp)
# 
# clinical_vars_subtyping_train$sex[clinical_vars_subtyping_train$sex=="Male"]=0
# clinical_vars_subtyping_train$sex[clinical_vars_subtyping_train$sex=="Female"]=1
# 
# #Include vars assays SPIKE IgG, N1-CT
# nasal_viral_load_acute_train <- nasal_viral_load_acute %>%
#   select(event_id, N1_CT) %>%
#   filter(event_id %in% clinical_data_acute_train$event_id )
# 
# antibody_titers_acute_train <- as.data.frame(antibody_titers_acute_env$preprocessed_data$serum_rbd_abtiters) %>%
#   rownames_to_column("sample_id") %>%
#   left_join(clinical_data_acute_train %>% 
#   select(sample_id,event_id)) %>%
#   select(-c(sample_id)) %>%
#   filter(!is.na(event_id)) %>%
#   select(event_id, `AUC Spike IgG`)
# 
# # Join all variables
# clinical_vars_subtyping_train <- clinical_vars_subtyping_train %>%
#   left_join(nasal_viral_load_acute_train) %>%
#   left_join(antibody_titers_acute_train) %>%
#   na.omit()
# 
# # Move event_id to rownames and standardize
# row_names <- clinical_vars_subtyping_train$event_id
# clinical_vars_subtyping_train$event_id <- NULL
# clinical_vars_subtyping_train <- as.data.frame(sapply(clinical_vars_subtyping_train, as.numeric))
# clinical_vars_subtyping_train <- scale(clinical_vars_subtyping_train)
# rownames(clinical_vars_subtyping_train) <- row_names

# # TODO update clinical data test
# # Clinical data test
# clinical_vars_test <- clinical_data_test %>%
#   dplyr::select(all_of(all_clinical_vars)) %>%
#   distinct()
# 
# # Add mean BMI for missing values, set missing baselines to 0.
# clinical_vars_test$bmi[is.na(clinical_vars_test$bmi)]<-mean(clinical_vars_test$bmi,na.rm=TRUE)
# clinical_vars_test$baseline_img_infil[is.na(clinical_vars_test$baseline_img_infil)] = 0
# clinical_vars_test$baseline_img_infil[clinical_vars_test$baseline_img_infil==999]<-0
# 
# # One-hot encoding multi-variate data
# temp <- data.frame(clinical_vars_test$race)
# colnames(temp) <-c("race")
# temp$race<- factor(temp$race, levels = unique(clinical_data_train$race))
# newtemp <- one_hot(as.data.table(temp), dropUnusedLevels = F)
# clinical_vars_test <- subset(clinical_vars_test, select = -c(race) )
# clinical_vars_test <- cbind(clinical_vars_test,newtemp)
# 
# temp <- data.frame(clinical_vars_test$ethnicity)
# colnames(temp) <- c("ethnicity")
# temp$ethnicity <- factor(temp$ethnicity, levels = unique(clinical_data_train$ethnicity))
# newtemp <- one_hot(as.data.table(temp), dropUnusedLevels = F)
# clinical_vars_test <- subset(clinical_vars_test, select = -c(ethnicity) )
# clinical_vars_test <- cbind(clinical_vars_test, newtemp)
# 
# temp <- data.frame(clinical_vars_test$discretized_admit_age_quantile)
# colnames(temp)<-c("discretized_admit_age_quantile")
# temp$discretized_admit_age_quantile<-factor(temp$discretized_admit_age_quantile, levels = unique(clinical_data_train$discretized_admit_age_quantile))
# newtemp <- one_hot(as.data.table(temp), dropUnusedLevels = F)
# clinical_vars_test <- subset(clinical_vars_test, select = -c(discretized_admit_age_quantile) )
# clinical_vars_test <-cbind(clinical_vars_test, newtemp)
# 
# clinical_vars_test$sex[clinical_vars_test$sex=="Male"]=0
# clinical_vars_test$sex[clinical_vars_test$sex=="Female"]=1
# 
# # Move event_id to rownames and standardize
# row_names <- clinical_vars_test$event_id
# clinical_vars_test$event_id <- NULL
# clinical_vars_test <- as.data.frame(sapply(clinical_vars_test, as.numeric))
# clinical_vars_test <- scale(clinical_vars_test)
# rownames(clinical_vars_test) <- row_names
# 
# # Order columns in same order as train data
# clinical_vars_test <- clinical_vars_test[,colnames(clinical_vars_ordered)]
# stopifnot(all(colnames(clinical_vars_ordered) == colnames(clinical_vars_test)))


# Clinical subtyping IgSpike viral load acute average ----------------

clinical_vars_subtyping <- c("sex","discretized_admit_age_quantile", "bmi", "los", "baseline_sofa_score_new")
include_vars_comorb_subtyping <- c("comorb_htn", "comorb_isaric_dm", "comorb_anyresp_noasthma", "comorb_isaric_asthma", "comorb_isaric_cardiac",
                         "comorb_isaric_ckd", "comorb_isaric_neoplasm", "comorb_neuro", "comorb_isaric_liver", "comorb_hxtrans",
                         "comorb_hiv2", "comorb_eversmkvape", "comorb_anysubstance", "comorb_count")

all_clinical_vars_subtyping <- c("event_id", clinical_vars_subtyping, include_vars_comorb_subtyping)


# Clinical data train
clinical_vars_subtyping_acute_train <- clinical_data_acute_train %>%
  dplyr::select("participant_id", all_of(clinical_vars_subtyping), all_of(include_vars_comorb_subtyping)) %>%
  distinct()

# Add mean BMI for missing values, set missing baselines to 0.
clinical_vars_subtyping_acute_train$bmi[is.na(clinical_vars_subtyping_acute_train$bmi)]<-mean(clinical_vars_subtyping_acute_train$bmi,na.rm=TRUE)

# One-hot encoding of discretized age
temp <- data.frame(clinical_vars_subtyping_acute_train$discretized_admit_age_quantile)
colnames(temp)<-c("discretized_admit_age_quantile")
temp$discretized_admit_age_quantile<-as.factor(temp$discretized_admit_age_quantile)
newtemp <- one_hot(as.data.table(temp))
clinical_vars_subtyping_acute_train <- subset(clinical_vars_subtyping_acute_train, select = -c(discretized_admit_age_quantile) )
clinical_vars_subtyping_acute_train <-cbind(clinical_vars_subtyping_acute_train, newtemp)

# Encoding of sex
clinical_vars_subtyping_acute_train$sex[clinical_vars_subtyping_acute_train$sex=="Male"]=0
clinical_vars_subtyping_acute_train$sex[clinical_vars_subtyping_acute_train$sex=="Female"]=1

#Include vars average of acute time point of SPIKE IgG titers, N1-CT viral load
nasal_viral_load_acute_participant_mean <- nasal_viral_load_acute %>%
  select(participant_id, event_id, N1_CT) %>%
  pivot_longer(cols = -c(event_id, participant_id), names_to = "feature", values_to = "value") %>%
  group_by(participant_id, feature) %>%
  summarize(mean_value = mean(value, na.rm=T)) %>%
  pivot_wider(id_cols = participant_id, names_from = feature, values_from = mean_value)

antibody_titers_acute_participant_mean <- as.data.frame(antibody_titers_acute_env$preprocessed_data$serum_rbd_abtiters) %>%
  select(`AUC Spike IgG`) %>%
  rownames_to_column("sample_id") %>%
  left_join(clinical_data_acute_train %>%
              select(sample_id,event_id,participant_id)) %>%
  select(-c(sample_id)) %>%
  filter(!is.na(event_id)) %>%
  pivot_longer(cols = -c(event_id, participant_id), names_to = "feature", values_to = "value") %>%
  group_by(participant_id, feature) %>%
  summarize(mean_value = mean(value, na.rm=T)) %>%
  pivot_wider(id_cols = participant_id, names_from = feature, values_from = mean_value)

clinical_vars_subtyping_acute_train <- clinical_vars_subtyping_acute_train %>%
  left_join(nasal_viral_load_acute_participant_mean, by = join_by(participant_id)) %>%
  left_join(antibody_titers_acute_participant_mean, by = join_by(participant_id)) %>%
  na.omit()
rownames(clinical_vars_subtyping_acute_train) <- clinical_vars_subtyping_acute_train$participant_id
clinical_vars_subtyping_acute_train$participant_id <- NULL


# Clinical subtyping model with additional variables ---------------------

clinical_additional <- read.csv("/scratch/data-fullcohort/clinical_data_supplements/longitudinal_multiomics_lab_clinical_supplement.csv")
clinical_additional <- clinical_additional %>%
  dplyr::rename(participant_id = studyid)

clinical_vars_additional <- clinical_data_acute_train %>%
  select(participant_id, baseline_crp_abn, comp_anemia, ever_steroids) %>%
  distinct() %>%
  left_join(clinical_additional %>% select(participant_id, baseline_lab_hgb, #baseline_lab_ferritin, 
                                           baseline_lab_hct, ever_nsaids, #baseline_lab_bili
  ))

summary(clinical_vars_additional)

clinical_acute_with_additional <- as.data.frame(clinical_vars_subtyping_acute_train) %>%
  rownames_to_column("participant_id") %>%
  inner_join(clinical_vars_additional, by = join_by(participant_id)) %>%
  na.omit() %>%
  remove_rownames() %>%
  tibble::column_to_rownames("participant_id")

clinical_acute_with_additional_labels <- clinical_data_acute_train %>%
  select(participant_id, pro_2_groups) %>%
  mutate(label_encoded = if_else(pro_2_groups == "MIN", 0 , 1)) %>%
  distinct() %>%
  filter(participant_id %in% rownames(clinical_acute_with_additional)) %>%
  column_to_rownames("participant_id")

clinical_acute_with_additional_labels_full <- clinical_acute_with_additional_labels[rownames(clinical_acute_with_additional), "pro_2_groups"]
clinical_acute_with_additional_labels_encoded <- clinical_acute_with_additional_labels[rownames(clinical_acute_with_additional), "label_encoded"]

set.seed(1)
clinical_acute_with_additional_train_fold <- createFolds(as.factor(clinical_acute_with_additional_labels_full), k = 10, list = F, returnTrain = FALSE)
clinical_acute_with_additional_train_fold_partid <- data.frame(fold = clinical_acute_with_additional_train_fold)
rownames(clinical_acute_with_additional_train_fold_partid) <- rownames(clinical_acute_with_additional_labels_full)

# SPEAR acute ----------------
## SPEAR acute average ----------------------------------------------

spear_object <- readRDS("/scratch/data-fullcohort/spear_model/results/020624/032124_spearphysical_finalscores_allsamples.rds")

spear_acute_train_participant_mean <- as.data.frame(spear_object$Acute) %>%
  rownames_to_column("event_id") %>%
  inner_join(clinical_data_acute_train %>% select(event_id, participant_id, event_type) %>% distinct(), by="event_id") %>%
  pivot_longer(cols = -c(event_id, participant_id, event_type), names_to = "feature", values_to = "value") %>%
  group_by(participant_id, feature) %>%
  summarize(mean_value = mean(value)) %>%
  filter(feature == "Factor1") %>%
  pivot_wider(id_cols = participant_id, names_from = feature, values_from = mean_value) %>%
  column_to_rownames("participant_id")

spear_acute_train_labels <- clinical_data_acute_train %>% 
  select(participant_id, pro_2_groups) %>%
  distinct() %>%
  filter(participant_id %in% rownames(spear_acute_train_participant_mean)) %>%
  mutate(label_encoded = if_else(pro_2_groups == "MIN", 0 , 1)) %>%
  column_to_rownames("participant_id")

spear_acute_train_encoded <- spear_acute_train_labels[rownames(spear_acute_train_participant_mean),c("label_encoded"), drop=F]
spear_acute_train_full <- spear_acute_train_labels[rownames(spear_acute_train_participant_mean),c("pro_2_groups"), drop=F]

set.seed(1)
spear_acute_train_fold <- createFolds(as.factor(spear_acute_train_full$pro_2_groups), k = 10, list = F, returnTrain = FALSE)
spear_acute_train_fold_partid <- data.frame(fold = spear_acute_train_fold)
rownames(spear_acute_train_fold_partid) <- rownames(spear_acute_train_full)


## Ensemble SPEAR acute average with clinical subtyping average -----------

# Order clinical subtyping data according to spear acute
spear_acute_train_participant_mean_clinical_subset <- spear_acute_train_participant_mean %>%
  rownames_to_column("participant_id") %>%
  filter(participant_id %in% rownames(clinical_vars_subtyping_acute_train))

clinical_vars_subtyping_acute_train_ordered <- clinical_vars_subtyping_acute_train[spear_acute_train_participant_mean_clinical_subset$participant_id,]
spear_acute_train_average_labels_encoded_clinical_subset <- spear_acute_train_encoded[spear_acute_train_participant_mean_clinical_subset$participant_id,]
spear_acute_train_average_labels_full_clinical_subset <- spear_acute_train_full[spear_acute_train_participant_mean_clinical_subset$participant_id,]

set.seed(1)
clinical_vars_subtyping_acute_train_fold <- createFolds(as.factor(spear_acute_train_average_labels_full_clinical_subset), k = 10, list = F, returnTrain = FALSE)
clinical_vars_subtyping_acute_train_fold_partid <- data.frame(fold = clinical_vars_subtyping_acute_train_fold)
rownames(clinical_vars_subtyping_acute_train_fold_partid) <- rownames(clinical_vars_subtyping_acute_train_ordered)

ensemble_acute_average_subtyping_train <- bind_cols(spear_acute_train_participant_mean_clinical_subset %>%
                                          column_to_rownames("participant_id"),
                                        clinical_vars_subtyping_acute_train_ordered)
ensemble_acute_average_subtyping_train <- sapply(ensemble_acute_average_subtyping_train, as.numeric)
ensemble_acute_average_subtyping_train <- scale(ensemble_acute_average_subtyping_train)


## Ensemble SPEAR acute average with clinical subtyping additional ---------

# Order clinical subtyping data according to spear acute
spear_acute_train_participant_mean_clinical_subset_additional <- spear_acute_train_participant_mean %>%
  rownames_to_column("participant_id") %>%
  filter(participant_id %in% rownames(clinical_acute_with_additional))

clinical_vars_subtyping_additional_train_ordered <- clinical_acute_with_additional[spear_acute_train_participant_mean_clinical_subset_additional$participant_id,]
spear_acute_train_average_labels_encoded_clinical_subset_additional <- spear_acute_train_encoded[spear_acute_train_participant_mean_clinical_subset_additional$participant_id,]
spear_acute_train_average_labels_full_clinical_subset_additional <- spear_acute_train_full[spear_acute_train_participant_mean_clinical_subset_additional$participant_id,]

set.seed(1)
clinical_vars_subtyping_acute_train_additional_fold <- createFolds(as.factor(spear_acute_train_average_labels_full_clinical_subset_additional), k = 10, list = F, returnTrain = FALSE)
clinical_vars_subtyping_acute_train_additional_fold_partid <- data.frame(fold = clinical_vars_subtyping_acute_train_additional_fold)
rownames(clinical_vars_subtyping_acute_train_additional_fold_partid) <- rownames(clinical_vars_subtyping_additional_train_ordered)

ensemble_acute_average_subtyping_additional_train <- bind_cols(spear_acute_train_participant_mean_clinical_subset_additional %>%
                                                      column_to_rownames("participant_id"),
                                                    clinical_vars_subtyping_additional_train_ordered)
ensemble_acute_average_subtyping_additional_train <- sapply(ensemble_acute_average_subtyping_additional_train, as.numeric)
ensemble_acute_average_subtyping_additional_train <- scale(ensemble_acute_average_subtyping_additional_train)


# ## Clinical subtyping visit 1 --------------
# 
# clinical_data_acute_train_visit1 <- clinical_data_acute_train %>%
#   filter(event_type == "Visit 1")
# 
# clinical_vars_subtyping_train_visit1 <- as.data.frame(clinical_vars_subtyping_train) %>%
#   rownames_to_column("event_id") %>%
#   left_join(clinical_data_acute_train_visit1 %>% select(event_id, participant_id, event_type) %>% distinct()) %>%
#   filter(event_type == "Visit 1") %>%
#   select(-c(event_id, event_type)) %>%
#   column_to_rownames("participant_id")
#   
# clinical_vars_subtyping_train_encoded <- spear_acute_train_labels[rownames(clinical_vars_subtyping_train_visit1),c("label_encoded"), drop=F]
# clinical_vars_subtyping_train_full <- spear_acute_train_labels[rownames(clinical_vars_subtyping_train_visit1),c("pro_2_groups"), drop=F]
# 
# set.seed(1)
# clinical_vars_subtyping_train_fold <- createFolds(as.factor(clinical_vars_subtyping_train_full$pro_2_groups), k = 10, list = F, returnTrain = FALSE)
# clinical_vars_subtyping_train_fold_partid <- data.frame(fold = clinical_vars_subtyping_train_fold)
# rownames(clinical_vars_subtyping_train_fold_partid) <- rownames(clinical_vars_subtyping_train_full)

# ensemble_visit1_subtyping_train <- cbind(spear_visit1_train_clinical_subset_mat, 
#                                          clinical_vars_subtyping_visit1_train_ordered)
# ensemble_visit1_subtyping_train <- scale(ensemble_visit1_subtyping_train)


## SPEAR acute per visit ----------------

acute_visits <- c("Visit 1", "Visit 2", "Visit 3", "Visit 4", "Visit 5", "Visit 6")

spear_acute_visits <- list()
for (visit in acute_visits) {
  spear_acute_train_visit <- as.data.frame(spear_object$Acute) %>%
    rownames_to_column("event_id") %>%
    inner_join(clinical_data_acute_train %>% 
                 select(event_id, participant_id, event_type) %>% 
                 distinct(), by="event_id") %>%
    filter(event_type == visit) %>% # Filter by visit number
    pivot_longer(cols = -c(event_id, participant_id, event_type), names_to = "feature", values_to = "value") %>%
    filter(feature == "Factor1") %>%
    pivot_wider(id_cols = participant_id, names_from = feature, values_from = value) %>%
    column_to_rownames("participant_id")
  
  spear_acute_train_visit_labels <- clinical_data_acute_train %>% 
    select(participant_id, pro_2_groups) %>%
    distinct() %>%
    filter(participant_id %in% rownames(spear_acute_train_visit)) %>%
    mutate(label_encoded = if_else(pro_2_groups == "MIN", 0 , 1)) %>%
    column_to_rownames("participant_id")
  
  spear_acute_train_visit_encoded <- spear_acute_train_visit_labels[rownames(spear_acute_train_visit),c("label_encoded"), drop=F]
  spear_acute_train_visit_full <- spear_acute_train_visit_labels[rownames(spear_acute_train_visit),c("pro_2_groups"), drop=F]
  
  set.seed(1)
  spear_acute_train_visit_fold <- createFolds(as.factor(spear_acute_train_visit_full$pro_2_groups), k = 10, list = F, returnTrain = FALSE)
  spear_acute_train_visit_fold_partid <- data.frame(fold = spear_acute_train_visit_fold)
  rownames(spear_acute_train_visit_fold_partid) <- rownames(spear_acute_train_visit_full)
  
  spear_acute_visits[[visit]][["matrix"]] <- spear_acute_train_visit
  spear_acute_visits[[visit]][["labels"]] <- spear_acute_train_visit_labels
  spear_acute_visits[[visit]][["labels_encoded"]] <- spear_acute_train_visit_encoded
  spear_acute_visits[[visit]][["labels_full"]] <- spear_acute_train_visit_full
  spear_acute_visits[[visit]][["folds"]] <- spear_acute_train_visit_fold_partid
  
  # Ensemble with clinical subtyping data
  spear_acute_train_visit_clinical_subset <- spear_acute_train_visit %>%
    rownames_to_column("participant_id") %>%
    inner_join(clinical_vars_subtyping_acute_train %>% rownames_to_column("participant_id"))
  
  spear_acute_visit_labels_encoded_clinical_subset <- spear_acute_train_visit_encoded[spear_acute_train_visit_clinical_subset$participant_id,]
  spear_acute_visit_labels_full_clinical_subset <- spear_acute_train_visit_full[spear_acute_train_visit_clinical_subset$participant_id,]
  
  set.seed(1)
  spear_acute_visit_labels_encoded_clinical_subset_fold <- createFolds(as.factor(spear_acute_visit_labels_full_clinical_subset), k = 10, list = F, returnTrain = FALSE)

  spear_acute_train_visit_clinical_subset_mat <- spear_acute_train_visit_clinical_subset %>% column_to_rownames("participant_id")
  ensemble_acute_visit_subtyping_train <- sapply(spear_acute_train_visit_clinical_subset_mat, as.numeric)
  ensemble_acute_visit_subtyping_train <- scale(ensemble_acute_visit_subtyping_train)
  
  spear_acute_visits[[visit]][["ensemble_clinical_matrix"]] <- ensemble_acute_visit_subtyping_train
  spear_acute_visits[[visit]][["ensemble_clinical_labels_encoded"]] <- spear_acute_visit_labels_encoded_clinical_subset
  spear_acute_visits[[visit]][["ensemble_clinical_labels_full"]] <- spear_acute_visit_labels_full_clinical_subset
  spear_acute_visits[[visit]][["ensemble_clinical_folds"]] <- spear_acute_visit_labels_encoded_clinical_subset_fold
  
  # Ensemble with clinical subtyping additional
  spear_acute_train_visit_clinical_subset <- spear_acute_train_visit %>%
    rownames_to_column("participant_id") %>%
    inner_join(clinical_acute_with_additional %>% rownames_to_column("participant_id"))
  
  spear_acute_visit_labels_encoded_clinical_subset <- spear_acute_train_visit_encoded[spear_acute_train_visit_clinical_subset$participant_id,]
  spear_acute_visit_labels_full_clinical_subset <- spear_acute_train_visit_full[spear_acute_train_visit_clinical_subset$participant_id,]
  
  set.seed(1)
  spear_acute_visit_labels_encoded_clinical_subset_fold <- createFolds(as.factor(spear_acute_visit_labels_full_clinical_subset), k = 10, list = F, returnTrain = FALSE)
  
  spear_acute_train_visit_clinical_subset_mat <- spear_acute_train_visit_clinical_subset %>% column_to_rownames("participant_id")
  ensemble_acute_visit_subtyping_train <- sapply(spear_acute_train_visit_clinical_subset_mat, as.numeric)
  ensemble_acute_visit_subtyping_train <- scale(ensemble_acute_visit_subtyping_train)
  
  spear_acute_visits[[visit]][["ensemble_clinical_additional_matrix"]] <- ensemble_acute_visit_subtyping_train
  spear_acute_visits[[visit]][["ensemble_clinical_additional_labels_encoded"]] <- spear_acute_visit_labels_encoded_clinical_subset
  spear_acute_visits[[visit]][["ensemble_clinical_additional_labels_full"]] <- spear_acute_visit_labels_full_clinical_subset
  spear_acute_visits[[visit]][["ensemble_clinical_additional_folds"]] <- spear_acute_visit_labels_encoded_clinical_subset_fold
  
}

# SPEAR convalescent --------------------------

## SPEAR convalescent average ------------------

spear_convalescent_train <- as.data.frame(rbind(spear_object$Train, spear_object$TrainNoPhysical_Overlap, spear_object$TrainNoPhysical_NoOverlap))

spear_convalescent_train_participant_mean <- spear_convalescent_train %>%
  #select(Factor1) %>%
  rownames_to_column("event_id") %>%
  left_join(clinical_data_train %>% select(event_id, participant_id, event_type) %>% distinct(), by="event_id") %>%
  pivot_longer(cols = -c(event_id, participant_id, event_type), names_to = "feature", values_to = "value") %>%
  group_by(participant_id, feature) %>%
  summarize(mean_value = mean(value)) %>%
  filter(feature == "Factor1") %>%
  pivot_wider(id_cols = participant_id, names_from = feature, values_from = mean_value) %>%
  column_to_rownames("participant_id")

spear_convalescent_train_labels <- clinical_data_train %>% 
  select(participant_id, pro_2_groups) %>%
  distinct() %>%
  filter(participant_id %in% rownames(spear_convalescent_train_participant_mean)) %>%
  mutate(label_encoded = if_else(pro_2_groups == "MIN", 0 , 1)) %>%
  column_to_rownames("participant_id")

spear_convalescent_train_encoded <- spear_convalescent_train_labels[rownames(spear_convalescent_train_participant_mean),c("label_encoded"), drop=F]
spear_convalescent_train_full <- spear_convalescent_train_labels[rownames(spear_convalescent_train_participant_mean),c("pro_2_groups"), drop=F]

set.seed(1)
spear_convalescent_train_fold <- createFolds(as.factor(spear_convalescent_train_full$pro_2_groups), k = 10, list = F, returnTrain = FALSE)
spear_convalescent_train_fold_partid <- data.frame(fold = spear_convalescent_train_fold)
rownames(spear_convalescent_train_fold_partid) <- rownames(spear_convalescent_train_full)


## Ensemble SPEAR convalescent clinical acute ----------------

# Order clinical subtyping data according to spear convalescent
spear_convalescent_train_participant_mean_clinical_subset <- spear_convalescent_train_participant_mean %>%
  rownames_to_column("participant_id") %>%
  filter(participant_id %in% rownames(clinical_vars_subtyping_acute_train))

clinical_vars_subtyping_acute_train_ordered_spearconv <- clinical_vars_subtyping_acute_train[spear_convalescent_train_participant_mean_clinical_subset$participant_id,]
spear_convalescent_train_average_labels_encoded_clinical_subset <- spear_convalescent_train_encoded[spear_convalescent_train_participant_mean_clinical_subset$participant_id,]
spear_convalescent_train_average_labels_full_clinical_subset <- spear_convalescent_train_full[spear_convalescent_train_participant_mean_clinical_subset$participant_id,]
spear_convalescent_train_average_labels_full_clinical_subset_fold <- clinical_vars_subtyping_acute_train_fold_partid[spear_convalescent_train_participant_mean_clinical_subset$participant_id,]

ensemble_conv_average_subtyping_train <- bind_cols(spear_convalescent_train_participant_mean_clinical_subset %>%
                                                      column_to_rownames("participant_id"),
                                                   clinical_vars_subtyping_acute_train_ordered_spearconv)
ensemble_conv_average_subtyping_train <- sapply(ensemble_conv_average_subtyping_train, as.numeric)
ensemble_conv_average_subtyping_train <- scale(ensemble_conv_average_subtyping_train)


## Ensemble SPEAR convalescent clinical acute additional --------

# Order clinical subtyping data according to spear convalescent
spear_convalescent_train_participant_mean_clinical_subset_additional <- spear_convalescent_train_participant_mean %>%
  rownames_to_column("participant_id") %>%
  filter(participant_id %in% rownames(clinical_acute_with_additional))

clinical_vars_subtyping_additional_train_ordered_spearconv <- clinical_acute_with_additional[spear_convalescent_train_participant_mean_clinical_subset_additional$participant_id,]
spear_convalescent_train_average_labels_encoded_clinical_subset_additional <- spear_convalescent_train_encoded[spear_convalescent_train_participant_mean_clinical_subset_additional$participant_id,]
spear_convalescent_train_average_labels_full_clinical_subset_additional <- spear_convalescent_train_full[spear_convalescent_train_participant_mean_clinical_subset_additional$participant_id,]
spear_convalescent_train_average_labels_full_clinical_subset_additional_fold <- clinical_vars_subtyping_acute_train_additional_fold_partid[spear_convalescent_train_participant_mean_clinical_subset_additional$participant_id,]

ensemble_conv_average_subtyping_additional_train <- bind_cols(spear_convalescent_train_participant_mean_clinical_subset_additional %>%
                                                     column_to_rownames("participant_id"),
                                                     clinical_vars_subtyping_additional_train_ordered_spearconv)
ensemble_conv_average_subtyping_additional_train <- sapply(ensemble_conv_average_subtyping_additional_train, as.numeric)
ensemble_conv_average_subtyping_additional_train <- scale(ensemble_conv_average_subtyping_additional_train)


## SPEAR convalescent per visit --------------

conv_visits <- c("Visit 7", "Visit 8", "Visit 9", "Visit 10")

spear_conv_visits <- list()
for (visit in conv_visits) {
  spear_conv_train_visit <- as.data.frame(spear_convalescent_train) %>%
    rownames_to_column("event_id") %>%
    inner_join(clinical_data_train %>% 
                 select(event_id, participant_id, event_type) %>% 
                 distinct(), by="event_id") %>%
    filter(event_type == visit) %>% # Filter by visit number
    pivot_longer(cols = -c(event_id, participant_id, event_type), names_to = "feature", values_to = "value") %>%
    filter(feature == "Factor1") %>%
    pivot_wider(id_cols = participant_id, names_from = feature, values_from = value) %>%
    column_to_rownames("participant_id")
  
  spear_conv_train_visit_labels <- clinical_data_train %>% 
    select(participant_id, pro_2_groups) %>%
    distinct() %>%
    filter(participant_id %in% rownames(spear_conv_train_visit)) %>%
    mutate(label_encoded = if_else(pro_2_groups == "MIN", 0 , 1)) %>%
    column_to_rownames("participant_id")
  
  spear_conv_train_visit_encoded <- spear_conv_train_visit_labels[rownames(spear_conv_train_visit),c("label_encoded"), drop=F]
  spear_conv_train_visit_full <- spear_conv_train_visit_labels[rownames(spear_conv_train_visit),c("pro_2_groups"), drop=F]
  
  set.seed(1)
  spear_conv_train_visit_fold <- createFolds(as.factor(spear_conv_train_visit_full$pro_2_groups), k = 10, list = F, returnTrain = FALSE)
  spear_conv_train_visit_fold_partid <- data.frame(fold = spear_conv_train_visit_fold)
  rownames(spear_conv_train_visit_fold_partid) <- rownames(spear_conv_train_visit_full)
  
  spear_conv_visits[[visit]][["matrix"]] <- spear_conv_train_visit
  spear_conv_visits[[visit]][["labels"]] <- spear_conv_train_visit_labels
  spear_conv_visits[[visit]][["labels_encoded"]] <- spear_conv_train_visit_encoded
  spear_conv_visits[[visit]][["labels_full"]] <- spear_conv_train_visit_full
  spear_conv_visits[[visit]][["folds"]] <- spear_conv_train_visit_fold_partid
  
  # Ensemble with clinical subtyping data
  spear_conv_train_visit_clinical_subset <- spear_conv_train_visit %>%
    rownames_to_column("participant_id") %>%
    inner_join(clinical_vars_subtyping_acute_train %>% rownames_to_column("participant_id"))
  
  spear_conv_visit_labels_encoded_clinical_subset <- spear_conv_train_visit_encoded[spear_conv_train_visit_clinical_subset$participant_id,]
  spear_conv_visit_labels_full_clinical_subset <- spear_conv_train_visit_full[spear_conv_train_visit_clinical_subset$participant_id,]
  
  set.seed(1)
  spear_conv_visit_labels_encoded_clinical_subset_fold <- createFolds(as.factor(spear_conv_visit_labels_full_clinical_subset), k = 10, list = F, returnTrain = FALSE)
  
  spear_conv_train_visit_clinical_subset_mat <- spear_conv_train_visit_clinical_subset %>% column_to_rownames("participant_id")
  ensemble_conv_visit_subtyping_train <- sapply(spear_conv_train_visit_clinical_subset_mat, as.numeric)
  ensemble_conv_visit_subtyping_train <- scale(ensemble_conv_visit_subtyping_train)
  
  spear_conv_visits[[visit]][["ensemble_clinical_matrix"]] <- ensemble_conv_visit_subtyping_train
  spear_conv_visits[[visit]][["ensemble_clinical_labels_encoded"]] <- spear_conv_visit_labels_encoded_clinical_subset
  spear_conv_visits[[visit]][["ensemble_clinical_labels_full"]] <- spear_conv_visit_labels_full_clinical_subset
  spear_conv_visits[[visit]][["ensemble_clinical_folds"]] <- spear_conv_visit_labels_encoded_clinical_subset_fold
  
  # Ensemble with clinical subtyping additional
  spear_conv_train_visit_clinical_subset <- spear_conv_train_visit %>%
    rownames_to_column("participant_id") %>%
    inner_join(clinical_acute_with_additional %>% rownames_to_column("participant_id"))
  
  spear_conv_visit_labels_encoded_clinical_subset <- spear_conv_train_visit_encoded[spear_conv_train_visit_clinical_subset$participant_id,]
  spear_conv_visit_labels_full_clinical_subset <- spear_conv_train_visit_full[spear_conv_train_visit_clinical_subset$participant_id,]
  
  set.seed(1)
  spear_conv_visit_labels_encoded_clinical_subset_fold <- createFolds(as.factor(spear_conv_visit_labels_full_clinical_subset), k = 10, list = F, returnTrain = FALSE)
  
  spear_conv_train_visit_clinical_subset_mat <- spear_conv_train_visit_clinical_subset %>% column_to_rownames("participant_id")
  ensemble_conv_visit_subtyping_train <- sapply(spear_conv_train_visit_clinical_subset_mat, as.numeric)
  ensemble_conv_visit_subtyping_train <- scale(ensemble_conv_visit_subtyping_train)
  
  spear_conv_visits[[visit]][["ensemble_clinical_additional_matrix"]] <- ensemble_conv_visit_subtyping_train
  spear_conv_visits[[visit]][["ensemble_clinical_additional_labels_encoded"]] <- spear_conv_visit_labels_encoded_clinical_subset
  spear_conv_visits[[visit]][["ensemble_clinical_additional_labels_full"]] <- spear_conv_visit_labels_full_clinical_subset
  spear_conv_visits[[visit]][["ensemble_clinical_additional_folds"]] <- spear_conv_visit_labels_encoded_clinical_subset_fold
  
}

# Analyte models ---------------------

mofa_model_train <- readRDS("/scratch/data-fullcohort/MOFA_factors/MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.rds")
mofa_model_test <- readRDS("/scratch/data-fullcohort/MOFA_factors/MOFA_test_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.rds")

train_data <- MOFA2::get_data(mofa_model_train)
test_data <- MOFA2::get_data(mofa_model_test)


feature_data_train <- list(
  SO = as.data.frame(t(train_data$SO$group1)),
  PMG = as.data.frame(t(train_data$PMG$group1)),
  PPT = as.data.frame(t(train_data$PPT$group1)),
  PPG = as.data.frame(t(train_data$PPG$group1)),
  PGX = as.data.frame(t(train_data$PGX$group1)),
  BCT = as.data.frame(t(train_data$BCT$group1))
)

path_imputed <- "/scratch/data-fullcohort/spear/imputed_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212/"

feature_data_train_imputed <- list(
  SO = read.csv(paste0(path_imputed, "spear_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_SO.csv")) %>% column_to_rownames("Features"),
  PMG = read.csv(paste0(path_imputed, "spear_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_PMG.csv")) %>% column_to_rownames("Features"),
  PPT = read.csv(paste0(path_imputed, "spear_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_PPT.csv")) %>% column_to_rownames("Features"),
  PPG = read.csv(paste0(path_imputed, "spear_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_PPG.csv")) %>% column_to_rownames("Features"),
  PGX = read.csv(paste0(path_imputed, "spear_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_PGX.csv")) %>% column_to_rownames("Features"),
  BCT = read.csv(paste0(path_imputed, "spear_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_BCT.csv")) %>% column_to_rownames("Features")
)

feature_data_test <- list(
  SO = as.data.frame(t(test_data$SO$group1)),
  PMG = as.data.frame(t(test_data$PMG$group1)),
  PPT = as.data.frame(t(test_data$PPT$group1)),
  PPG = as.data.frame(t(test_data$PPG$group1)),
  PGX = as.data.frame(t(test_data$PGX$group1)),
  BCT = as.data.frame(t(test_data$BCT$group1))
)

metabolite_row_DF = data_env_train$plasma_metabolomics_global_rowfeature
metabolite_row_DF$CHEMICAL_NAME <- str_replace_all(metabolite_row_DF$CHEMICAL_NAME," ","_")
metabolite_row_DF$CHEMICAL_NAME <- str_replace_all(metabolite_row_DF$CHEMICAL_NAME,"[[:punct:]]","")

metabolomics_translate <- data.frame(metabolite_id = row.names(metabolite_row_DF),
                                     CHEMICAL_NAME = str_replace_all(metabolite_row_DF$CHEMICAL_NAME," ","_"))


gene_ensembl_PGX <- data_env_train$pbmc_transcriptomics_rowfeature %>%
  rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, gene_name)

sample_event_train <- clinical_data_train %>% dplyr::select(sample_id, event_id) %>% distinct()


## Hallmark convalescent average ---------------

GSEA_result_hallmark_Factor1_PGX <- read.table("../enrichment_spear_factors_promis_physical_Factor1/GSEA_hallmark_Factor1_PGX_results.tsv", sep = "\t", header=T)
top_analytes_core_list <- GSEA_result_hallmark_Factor1_PGX %>% 
  filter(ID == "HALLMARK_HEME_METABOLISM") %>%
  pull(core_enrichment)

top_analytes_core <- str_split(top_analytes_core_list, "/")[[1]]

feat_data_assay <- feature_data_train_imputed[["PGX"]]
feat_data_assay_translate <- as.data.frame(t(feat_data_assay)) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_ensembl_PGX, by="ensembl_gene_id") %>%
  dplyr::select(!ensembl_gene_id) %>%
  column_to_rownames("gene_name")
feat_data_assay_translated <- as.data.frame(t(feat_data_assay_translate))

top_feat_data_df  <- feat_data_assay_translated %>% 
  dplyr::select(all_of(top_analytes_core))

heme_convalescent_train_participant_mean <- top_feat_data_df %>%
  rownames_to_column("event_id") %>%
  left_join(clinical_data_train %>% select(event_id, participant_id, event_type) %>% distinct(), by="event_id") %>%
  pivot_longer(cols = -c(event_id, participant_id, event_type), names_to = "feature", values_to = "value") %>%
  group_by(participant_id, feature) %>%
  summarize(mean_value = mean(value, na.rm=T)) %>%
  pivot_wider(id_cols = participant_id, names_from = feature, values_from = mean_value) %>%
  na.omit() %>%
  column_to_rownames("participant_id")

## Androgenic steroids convalescent average -----------------

GSEA_result_androgenic_steroids_Factor1_PMG <- read.table("../enrichment_spear_factors_promis_physical_Factor1/GSEA_subpath_Factor1_PMG_results.tsv", sep = "\t", header = T)
top_analytes_core_list <- GSEA_result_androgenic_steroids_Factor1_PMG %>%
  filter(ID == "Androgenic Steroids") %>%
  pull(core_enrichment)
top_analytes_core <- str_split(top_analytes_core_list, "/")[[1]]
top_analytes_core_simplified <- gsub("[[:punct:]]", "", top_analytes_core)
top_analytes_core_simplified <- gsub(" ","",top_analytes_core_simplified)

if (!all(top_analytes_core_simplified %in% metabolomics_translate$CHEMICAL_NAME)) {
  stop("ERROR: There are some metabolites in the list that are not found in the translation table!")
}

feat_data_assay <- feature_data_train_imputed[["PMG"]]
feat_data_assay_translate <- as.data.frame(t(feat_data_assay)) %>%
  rownames_to_column("metabolite_id") %>%
  left_join(metabolomics_translate, by = "metabolite_id") %>%
  dplyr::select(!metabolite_id) %>%
  column_to_rownames("CHEMICAL_NAME")
feat_data_assay_translated <- as.data.frame(t(feat_data_assay_translate))

top_feat_data_df <- feat_data_assay_translated %>%
  dplyr::select(all_of(top_analytes_core_simplified))
if (ncol(top_feat_data_df) != length(top_analytes_core_simplified)){
  stop("ERROR: some metabolites in the core enrichment were not found in the features table.")
}


androgenic_convalescent_train_participant_mean <- top_feat_data_df %>%
  rownames_to_column("event_id") %>%
  left_join(clinical_data_train %>% select(event_id, participant_id, event_type) %>% distinct(), by="event_id") %>%
  pivot_longer(cols = -c(event_id, participant_id, event_type), names_to = "feature", values_to = "value") %>%
  group_by(participant_id, feature) %>%
  summarize(mean_value = mean(value, na.rm=T)) %>%
  pivot_wider(id_cols = participant_id, names_from = feature, values_from = mean_value) %>%
  na.omit() %>%
  column_to_rownames("participant_id")

## SPEAR significant analytes convalescent average ----------

top_analyte_scores <- spear_object$analyte.scores %>% 
  dplyr::filter(Factor == "Factor1") %>%
  dplyr::mutate(joint.probability.rounded = round(joint.probability, 2)) %>%
  dplyr::filter(joint.probability.rounded >= 0.95)

top_feat_data <- list()
for (assay in unique(top_analyte_scores$Assay)) {
  print(assay)
  
  selected_features_assay <- top_analyte_scores %>% filter(Assay == assay)
  
  top_feature_data <- list()
  
  if (nrow(selected_features_assay > 0)) {
    top_feat_data_assay <- feature_data_train_imputed[[assay]] %>%
      dplyr::select(all_of(selected_features_assay$Analyte))
    
    if (assay == "PGX"){
    
      top_feat_data_assay_translate <- as.data.frame(t(top_feat_data_assay)) %>%
        rownames_to_column("ensembl_gene_id") %>%
        left_join(gene_ensembl_PGX, 
                  by="ensembl_gene_id") %>%
        dplyr::select(!ensembl_gene_id) %>%
        column_to_rownames("gene_name")
      top_feat_data_assay_translated <- as.data.frame(t(top_feat_data_assay_translate))
      
    } else if (assay == "PMG") {
      
      top_feat_data_assay_translate <- as.data.frame(t(top_feat_data_assay)) %>%
        rownames_to_column("metabolite_id") %>%
        left_join(metabolomics_translate,
                  by="metabolite_id") %>%
        dplyr::select(!metabolite_id) %>%
        column_to_rownames("CHEMICAL_NAME")
      top_feat_data_assay_translated <- as.data.frame(t(top_feat_data_assay_translate))
      
    } else {
      
      top_feat_data_assay_translated <- top_feat_data_assay
      top_feat_data_assay_translate <- NULL
      
    }
    
    colnames(top_feat_data_assay_translated) <- paste(assay,colnames(top_feat_data_assay_translated), sep = "_")
    top_feat_data[[assay]] <- top_feat_data_assay_translated
    rm(top_feat_data_assay_translated, 
       top_feat_data_assay_translate, top_feat_data_assay)
  }
}

top_feat_data_df <- bind_cols(top_feat_data)

spear_significant_convalescent_train_participant_mean <- top_feat_data_df %>%
  rownames_to_column("event_id") %>%
  left_join(clinical_data_train %>% select(event_id, participant_id, event_type) %>% distinct(), by="event_id") %>%
  pivot_longer(cols = -c(event_id, participant_id, event_type), names_to = "feature", values_to = "value") %>%
  group_by(participant_id, feature) %>%
  summarize(mean_value = mean(value, na.rm=T)) %>%
  pivot_wider(id_cols = participant_id, names_from = feature, values_from = mean_value) %>%
  na.omit() %>%
  column_to_rownames("participant_id")


## Combined convalescent -----------


part_id_combined <- intersect(rownames(heme_convalescent_train_participant_mean), 
                              intersect(rownames(androgenic_convalescent_train_participant_mean), 
                                        rownames(spear_significant_convalescent_train_participant_mean)))
combined_convalescent <- cbind(heme_convalescent_train_participant_mean[part_id_combined,], 
                               androgenic_convalescent_train_participant_mean[part_id_combined,], 
                               spear_significant_convalescent_train_participant_mean[part_id_combined,])



# AUROC CV all models -----------------------

# Set number of CV folds for all models
CV_folds = 10

# Set number of bootstrap samples for all models
n_bootstrap_samples = 100

models_folder = "./Models"
dir.create(models_folder)

spear_acute <- cv_scores_lasso_model_min_other(X=spear_acute_train_participant_mean,
                                               model_name = "SPEAR\nacute",
                                               Y=spear_acute_train_encoded[,1],
                                               Y_full=as.factor(spear_acute_train_full[,1]),
                                               n_folds = CV_folds,
                                               bootstrap_sample = n_bootstrap_samples)
saveRDS(spear_acute, file.path(models_folder, "spear_acute.rds"))

clinical_acute_subtyping <- cv_scores_lasso_model_min_other(X=as.matrix(clinical_vars_subtyping_acute_train_ordered),
                                                        model_name = "Clinical",
                                                        Y = spear_acute_train_average_labels_encoded_clinical_subset,
                                                        Y_full = as.factor(spear_acute_train_average_labels_full_clinical_subset),
                                                        n_folds = CV_folds,
                                                        bootstrap_sample = n_bootstrap_samples
                                                      )
saveRDS(clinical_acute_subtyping, file.path(models_folder, "clinical_acute_subtyping.rds"))


# clinical_subtyping_visit1 <- cv_scores_lasso_model_min_other(X=as.matrix(clinical_vars_subtyping_train_visit1),
#                                                       model_name = "Clinical\nSubtyping\nVisit1",
#                                                       Y = clinical_vars_subtyping_train_encoded[,1],
#                                                       Y_full = as.factor(clinical_vars_subtyping_train_full[,1]),
#                                                       n_folds = CV_folds,
#                                                       bootstrap_sample = n_bootstrap_samples
# )

clinical_subtyping_additional <- cv_scores_lasso_model_min_other(X = as.matrix(scale(sapply(clinical_acute_with_additional, as.numeric))),
                                                                 model_name = "Clinical\nadditional",
                                                                 Y = clinical_acute_with_additional_labels_encoded,
                                                                 Y_full = as.factor(clinical_acute_with_additional_labels_full),
                                                                 n_folds = CV_folds,
                                                                 bootstrap_sample = n_bootstrap_samples)
saveRDS(clinical_subtyping_additional, file.path(models_folder, "clinical_subtyping_additional.rds"))


ensemble_acute_subtyping <- cv_scores_lasso_model_min_other(X=as.matrix(ensemble_acute_average_subtyping_train),
                                                            model_name = "SPEAR acute\n+\nClinical",
                                                            Y = spear_acute_train_average_labels_encoded_clinical_subset,
                                                            Y_full = as.factor(spear_acute_train_average_labels_full_clinical_subset),
                                                            n_folds = CV_folds,
                                                            bootstrap_sample = n_bootstrap_samples)
saveRDS(ensemble_acute_subtyping, file.path(models_folder, "ensemble_acute_subtyping.rds"))

spear_convalescent <- cv_scores_lasso_model_min_other(X = as.matrix(spear_convalescent_train_participant_mean),
                                                      model_name = "SPEAR\nconv",
                                                      Y = spear_convalescent_train_encoded[,1],
                                                      Y_full = as.factor(spear_convalescent_train_full[,1]),
                                                      n_folds = CV_folds,
                                                      bootstrap_sample = n_bootstrap_samples)
saveRDS(spear_convalescent, file.path(models_folder, "spear_convalescent.rds"))


ensemble_convalescent_subtyping <- cv_scores_lasso_model_min_other(X = as.matrix(ensemble_conv_average_subtyping_train),
                                                                   model_name = "SPEAR conv\n+\nClinical",
                                                                   Y = spear_convalescent_train_average_labels_encoded_clinical_subset,
                                                                   Y_full = as.factor(spear_convalescent_train_average_labels_full_clinical_subset),
                                                                   n_folds = CV_folds,
                                                                   bootstrap_sample = n_bootstrap_samples)
saveRDS(ensemble_convalescent_subtyping, file.path(models_folder, "ensemble_convalescent_subtyping.rds"))

ensemble_acute_subtyping_additional <- cv_scores_lasso_model_min_other(X = as.matrix(ensemble_acute_average_subtyping_additional_train),
                                                                       model_name = "SPEAR acute\n+\nClinical\nadditional",
                                                                       Y = spear_acute_train_average_labels_encoded_clinical_subset_additional,
                                                                       Y_full = as.factor(spear_acute_train_average_labels_full_clinical_subset_additional),
                                                                       n_folds = CV_folds,
                                                                       bootstrap_sample = n_bootstrap_samples)
saveRDS(ensemble_acute_subtyping_additional, file.path(models_folder, "ensemble_acute_subtyping_additional.rds"))

  
ensemble_conv_subtyping_additional <- cv_scores_lasso_model_min_other(X = as.matrix(ensemble_conv_average_subtyping_additional_train),
                                                                      model_name = "SPEAR conv\n+\nClinical\nadditional",
                                                                      Y = spear_convalescent_train_average_labels_encoded_clinical_subset_additional,
                                                                      Y_full = as.factor(spear_convalescent_train_average_labels_full_clinical_subset_additional),
                                                                      n_folds = CV_folds,
                                                                      bootstrap_sample = n_bootstrap_samples)
saveRDS(ensemble_conv_subtyping_additional, file.path(models_folder, "ensemble_conv_subtyping_additional.rds"))

acute_visit_models = list()
for (visit in c("Visit 1", "Visit 2", "Visit 3", "Visit 4", "Visit 6")){
   acute_visit_models[[visit]][["spear"]] <- cv_scores_lasso_model_min_other(X = as.matrix(spear_acute_visits[[visit]][["matrix"]]),
                                                                             model_name = paste0(visit,"_SPEAR"),
                                                                             Y = spear_acute_visits[[visit]][["labels_encoded"]][,1],
                                                                             Y_full = as.factor(spear_acute_visits[[visit]][["labels_full"]][,1]),
                                                                             n_folds = CV_folds,
                                                                             bootstrap_sample = n_bootstrap_samples)
   acute_visit_models[[visit]][["spear_ensembl_clinical"]] <- cv_scores_lasso_model_min_other(X = as.matrix(spear_acute_visits[[visit]][["ensemble_clinical_matrix"]]),
                                                                                              model_name = paste0(visit,"_SPEARClinical"),
                                                                                              Y = spear_acute_visits[[visit]][["ensemble_clinical_labels_encoded"]],
                                                                                              Y_full = as.factor(spear_acute_visits[[visit]][["ensemble_clinical_labels_full"]]),
                                                                                              n_folds = CV_folds,
                                                                                              bootstrap_sample = n_bootstrap_samples)
   acute_visit_models[[visit]][["spear_ensembl_clinical_additional"]] <- cv_scores_lasso_model_min_other(X = as.matrix(spear_acute_visits[[visit]][["ensemble_clinical_additional_matrix"]]),
                                                                                              model_name = paste0(visit,"_SPEARClinicalAdditional"),
                                                                                              Y = spear_acute_visits[[visit]][["ensemble_clinical_additional_labels_encoded"]],
                                                                                              Y_full = as.factor(spear_acute_visits[[visit]][["ensemble_clinical_additional_labels_full"]]),
                                                                                              n_folds = CV_folds,
                                                                                              bootstrap_sample = n_bootstrap_samples)
}
saveRDS(acute_visit_models, file.path(models_folder, "acute_visit_models.rds"))


conv_visit_models = list()
for (visit in conv_visits){
  conv_visit_models[[visit]][["spear"]] <- cv_scores_lasso_model_min_other(X = as.matrix(spear_conv_visits[[visit]][["matrix"]]),
                                                                            model_name = paste0(visit,"_SPEAR"),
                                                                            Y = spear_conv_visits[[visit]][["labels_encoded"]][,1],
                                                                            Y_full = as.factor(spear_conv_visits[[visit]][["labels_full"]][,1]),
                                                                            n_folds = CV_folds,
                                                                            bootstrap_sample = n_bootstrap_samples)
  conv_visit_models[[visit]][["spear_ensembl_clinical"]] <- cv_scores_lasso_model_min_other(X = as.matrix(spear_conv_visits[[visit]][["ensemble_clinical_matrix"]]),
                                                                                             model_name = paste0(visit,"_SPEARClinical"),
                                                                                             Y = spear_conv_visits[[visit]][["ensemble_clinical_labels_encoded"]],
                                                                                             Y_full = as.factor(spear_conv_visits[[visit]][["ensemble_clinical_labels_full"]]),
                                                                                             n_folds = CV_folds,
                                                                                             bootstrap_sample = n_bootstrap_samples)
  conv_visit_models[[visit]][["spear_ensembl_clinical_additional"]] <- cv_scores_lasso_model_min_other(X = as.matrix(spear_conv_visits[[visit]][["ensemble_clinical_additional_matrix"]]),
                                                                                                        model_name = paste0(visit,"_SPEARClinicalAdditional"),
                                                                                                        Y = spear_conv_visits[[visit]][["ensemble_clinical_additional_labels_encoded"]],
                                                                                                        Y_full = as.factor(spear_conv_visits[[visit]][["ensemble_clinical_additional_labels_full"]]),
                                                                                                        n_folds = CV_folds,
                                                                                                       bootstrap_sample = n_bootstrap_samples)
}
saveRDS(conv_visit_models, file.path(models_folder, "conv_visit_models.rds"))

  
  
# heme_convalescent <- cv_scores_lasso_model_min_other(X = as.matrix(heme_convalescent_train_participant_mean),
#                                                      model_name = "Heme met\nconv",
#                                                      Y = spear_convalescent_train_encoded[rownames(heme_convalescent_train_participant_mean),1],
#                                                      Y_full = as.factor(spear_convalescent_train_full[rownames(heme_convalescent_train_participant_mean),1]),
#                                                      n_folds = CV_folds,
#                                                      bootstrap_sample = n_bootstrap_samples)
# 
# androgenic_convalescent <- cv_scores_lasso_model_min_other(X = as.matrix(androgenic_convalescent_train_participant_mean),
#                                                            model_name = "Androgenic\nconv",
#                                                            Y = spear_convalescent_train_encoded[rownames(androgenic_convalescent_train_participant_mean),1],
#                                                            Y_full = as.factor(spear_convalescent_train_full[rownames(androgenic_convalescent_train_participant_mean),1]),
#                                                            n_folds = CV_folds,
#                                                            bootstrap_sample = n_bootstrap_samples)

spear_significant_convalescent <- cv_scores_lasso_model_min_other(X = as.matrix(spear_significant_convalescent_train_participant_mean),
                                                                  model_name = "SPEAR\nsignificant\nconv",
                                                                  Y = spear_convalescent_train_encoded[rownames(spear_significant_convalescent_train_participant_mean),1],
                                                                  Y_full = as.factor(spear_convalescent_train_full[rownames(spear_significant_convalescent_train_participant_mean),1]),
                                                                  n_folds = CV_folds,
                                                                  bootstrap_sample = n_bootstrap_samples)

saveRDS(spear_significant_convalescent, file.path(models_folder, "spear_significant_convalescent.rds"))

# combined_convalescent_model <- cv_scores_lasso_model_min_other(X = as.matrix(combined_convalescent),
#                                                          model_name = "Combined\nconv",
#                                                          Y = spear_convalescent_train_encoded[rownames(combined_convalescent),1],
#                                                          Y_full = as.factor(spear_convalescent_train_full[rownames(combined_convalescent),1]),
#                                                          n_folds = CV_folds,
#                                                          bootstrap_sample = n_bootstrap_samples)



## Plot CV AUROC averages clinical subtyping and subtyping additional -----
all.res.df <- bind_rows(spear_acute$resdf,
                        clinical_acute_subtyping$resdf, ensemble_acute_subtyping$resdf,
                        clinical_subtyping_additional$resdf,   ensemble_acute_subtyping_additional$resdf,
                        )

all.res.df$Model <- factor(all.res.df$Model, levels = c("SPEAR\nacute",
                                                        "Clinical",   "SPEAR acute\n+\nClinical",
                                                        "Clinical\nadditional", "SPEAR acute\n+\nClinical\nadditional"
                                                        ))
model.cols = c("#08605F",
               "#DBB957", "#728D5B",
               "#DBB957",  "#728D5B"
               )
names(model.cols) = levels(all.res.df$Model)

auroc_table <- all.res.df %>%
  group_by(Model, Seed) %>%
  summarise(AUROC_CV_mean = mean(CV), AUROC_CV_SD = sd(CV)) %>%
  #mutate(Model = str_remove_all(Model,"\n")) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()

stat.test <- auroc_table %>%
  t_test(AUROC_CV_mean ~ Model, 
                  paired = T,
         comparisons = list(c("SPEAR\nacute", "Clinical"),
                            c("SPEAR\nacute", "SPEAR acute\n+\nClinical"),
                            c("Clinical", "Clinical\nadditional"),
                            c("SPEAR\nacute", "Clinical\nadditional"),
                            c("Clinical\nadditional",  "SPEAR acute\n+\nClinical\nadditional"),
                            c("SPEAR\nacute", "SPEAR acute\n+\nClinical\nadditional"),
                            c("SPEAR acute\n+\nClinical", "SPEAR acute\n+\nClinical\nadditional")
                            )) %>%
  add_significance() %>%
  add_xy_position(x = "Model", step.increase = 0.25)

cv.auroc.plot <- ggplot(auroc_table) +
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = AUROC_CV_mean, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = AUROC_CV_mean, fill = Model), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0) +
  scale_fill_manual(values = model.cols) +
  scale_y_continuous(breaks = seq(0.5,0.9,by=0.1), limits=c(0.5,0.83)) + 
  ylab("AUROC (CV)") +
  xlab(NULL) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))


cv.auroc.plot

# Save tables
path_to_save_tables <- "Tables/"
dir.create(path_to_save_tables)
write.csv(auroc_table %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_clinical_subtyping_additional_acute_comparison.csv"), row.names = F, quote = F)
write.csv(all.res.df %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_clinical_subtyping_additional_acute_comparison.csv"), row.names = F, quote =F)

# Save figure:
path_to_save_figs <- "Figs/"
path_to_save_rds <- "FigsRDS/"
#ggsave(plot = cv.auroc.plot, filename = paste0(path_to_save_figs, "CV_AUROC.pdf"), width = 3.5, height = 3)
saveRDS(cv.auroc.plot, file = paste0(path_to_save_rds, "CV_AUROC_spear_clinical_subtyping_additional_acute_comparison.rds"))


## Plot CV AUROC averages clinical acute subtyping only -----
all.res.df <- bind_rows(spear_acute$resdf,
                        clinical_acute_subtyping$resdf, ensemble_acute_subtyping$resdf,
)

all.res.df$Model <- factor(all.res.df$Model, levels = c("SPEAR\nacute",
                                                        "Clinical",   "SPEAR acute\n+\nClinical"))
model.cols = c("#08605F",
               "#DBB957", "#728D5B"
)
names(model.cols) = levels(all.res.df$Model)

auroc_table <- all.res.df %>%
  group_by(Model, Seed) %>%
  summarise(AUROC_CV_mean = mean(CV), AUROC_CV_SD = sd(CV)) %>%
  #mutate(Model = str_remove_all(Model,"\n")) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()

stat.test <- auroc_table %>%
  t_test(AUROC_CV_mean ~ Model, 
         paired = T,
         comparisons = list(c("SPEAR\nacute", "Clinical"),
                            c("SPEAR\nacute", "SPEAR acute\n+\nClinical")
         )) %>%
  add_significance() %>%
  add_xy_position(x = "Model", step.increase = 0.25)

cv.auroc.plot <- ggplot(auroc_table) +
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = AUROC_CV_mean, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = AUROC_CV_mean, fill = Model), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0) +
  scale_fill_manual(values = model.cols) +
  scale_y_continuous(breaks = seq(0.5,0.9,by=0.1), limits=c(0.5,0.83)) + 
  ylab("AUROC (CV)") +
  xlab(NULL) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))


cv.auroc.plot

# Save tables
path_to_save_tables <- "Tables/"
dir.create(path_to_save_tables)

auroc_bootstrap_mean <- auroc_table %>%
  group_by(Model) %>%
  summarise(AUROC_bootstrap_mean = mean(AUROC_CV_mean), AUROC_bootstrap_SD = sd(AUROC_CV_mean)) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()
write.csv(auroc_bootstrap_mean %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_clinical_subtyping_acute_comparison.csv"), row.names = F, quote = F)
write.csv(all.res.df %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_clinical_subtyping_acute_comparison_allseeds.csv"), row.names = F, quote =F)

# Save figure:

#ggsave(plot = cv.auroc.plot, filename = paste0(path_to_save_figs, "CV_AUROC.pdf"), width = 3.5, height = 3)
saveRDS(cv.auroc.plot, file = paste0(path_to_save_rds, "CV_AUROC_spear_clinical_subtyping_acute_comparison.rds"))


## Plot CV AUROC averages clinical acute subtyping additional only -----
all.res.df <- bind_rows(spear_acute$resdf,
                        clinical_subtyping_additional$resdf, ensemble_acute_subtyping_additional$resdf,
)

all.res.df$Model <- factor(all.res.df$Model, levels = c("SPEAR\nacute",
                                                        "Clinical\nadditional",   "SPEAR acute\n+\nClinical\nadditional"))
model.cols = c("#08605F",
               "#DBB957", "#728D5B"
)
names(model.cols) = levels(all.res.df$Model)

auroc_table <- all.res.df %>%
  group_by(Model, Seed) %>%
  summarise(AUROC_CV_mean = mean(CV), AUROC_CV_SD = sd(CV)) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()

stat.test <- auroc_table %>%
  t_test(AUROC_CV_mean ~ Model, 
         paired = T,
         comparisons = list(c("SPEAR\nacute", "Clinical\nadditional"),
                            c("SPEAR\nacute", "SPEAR acute\n+\nClinical\nadditional")
         )) %>%
  add_significance() %>%
  add_xy_position(x = "Model", step.increase = 0.25)

cv.auroc.plot <- ggplot(auroc_table) +
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = AUROC_CV_mean, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = AUROC_CV_mean, fill = Model), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0) +
  scale_fill_manual(values = model.cols) +
  scale_y_continuous(breaks = seq(0.5,0.9,by=0.1), limits=c(0.5,0.83)) + 
  ylab("AUROC (CV)") +
  xlab(NULL) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))


cv.auroc.plot

# Save tables
path_to_save_tables <- "Tables/"
dir.create(path_to_save_tables)

auroc_bootstrap_mean <- auroc_table %>%
  group_by(Model) %>%
  summarise(AUROC_bootstrap_mean = mean(AUROC_CV_mean), AUROC_bootstrap_SD = sd(AUROC_CV_mean)) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()
write.csv(auroc_bootstrap_mean %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_clinical_subtyping_additional_acute_comparison.csv"), row.names = F, quote = F)
write.csv(all.res.df %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_clinical_subtyping_additional_acute_comparison.csv"), row.names = F, quote =F)

# Save figure:

#ggsave(plot = cv.auroc.plot, filename = paste0(path_to_save_figs, "CV_AUROC.pdf"), width = 3.5, height = 3)
saveRDS(cv.auroc.plot, file = paste0(path_to_save_rds, "CV_AUROC_spear_clinical_subtyping_additional_acute_comparison.rds"))


## Plot CV AUROC averages conv clinical subtyping only -----
all.res.df <- bind_rows(
  spear_convalescent$resdf, ensemble_convalescent_subtyping$resdf,
  spear_significant_convalescent$resdf
)
all.res.df$Model <- factor(all.res.df$Model, levels = c(
  "SPEAR\nconv", "SPEAR conv\n+\nClinical", "SPEAR\nsignificant\nconv"
))

model.cols = c(
  "#08605F", "#728D5B","#2D0320"
)
names(model.cols) = levels(all.res.df$Model)

auroc_table <- all.res.df %>%
  group_by(Model, Seed) %>%
  summarise(AUROC_CV_mean = mean(CV), AUROC_CV_SD = sd(CV)) %>%
  #mutate(Model = str_remove_all(Model,"\n")) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()

stat.test <- auroc_table %>%
  t_test(AUROC_CV_mean ~ Model, 
         paired = T,
         comparisons = list(c("SPEAR\nconv", "SPEAR conv\n+\nClinical"),
                            c("SPEAR\nconv", "SPEAR\nsignificant\nconv")
         )) %>%
  add_significance() %>%
  add_xy_position(x = "Model", step.increase = 0.12)

cv.auroc.plot <- ggplot(auroc_table) +
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = AUROC_CV_mean, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = AUROC_CV_mean, fill = Model), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0) +
  scale_fill_manual(values = model.cols) +
  scale_y_continuous(breaks = seq(0.5,0.9,by=0.1), limits=c(0.5,0.83)) + 
  ylab("AUROC (CV)") +
  xlab(NULL) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))


cv.auroc.plot

# Save tables
path_to_save_tables <- "Tables/"
auroc_bootstrap_mean <- auroc_table %>%
  group_by(Model) %>%
  summarise(AUROC_bootstrap_mean = mean(AUROC_CV_mean), AUROC_bootstrap_SD = sd(AUROC_CV_mean)) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()
write.csv(auroc_bootstrap_mean %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_conv_clinical_analyte_sets.csv"), row.names = F, quote = F)
write.csv(auroc_table %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_conv_clinical_analyte_sets_allseeds.csv"), row.names = F, quote = F)

# Save figure:
path_to_save_figs <- "Figs/"
path_to_save_rds <- "FigsRDS/"
#ggsave(plot = cv.auroc.plot, filename = paste0(path_to_save_figs, "CV_AUROC.pdf"), width = 3.5, height = 3)
saveRDS(cv.auroc.plot, file = paste0(path_to_save_rds, "CV_AUROC_spear_conv_clinical_analyte_sets.rds"))

## Plot CV AUROC averages conv clinical subtyping additional only -----
all.res.df <- bind_rows(
  spear_convalescent$resdf, ensemble_conv_subtyping_additional$resdf,
  spear_significant_convalescent$resdf
)
all.res.df$Model <- factor(all.res.df$Model, levels = c(
  "SPEAR\nconv", "SPEAR conv\n+\nClinical\nadditional", "SPEAR\nsignificant\nconv"
))

model.cols = c(
  "#08605F", "#728D5B","#2D0320"
)
names(model.cols) = levels(all.res.df$Model)

auroc_table <- all.res.df %>%
  group_by(Model, Seed) %>%
  summarise(AUROC_CV_mean = mean(CV), AUROC_CV_SD = sd(CV)) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()

stat.test <- auroc_table %>%
  t_test(AUROC_CV_mean ~ Model, 
         paired = T,
         comparisons = list(c("SPEAR\nconv", "SPEAR conv\n+\nClinical\nadditional"),
                            c("SPEAR\nconv", "SPEAR\nsignificant\nconv")
         )) %>%
  add_significance() %>%
  add_xy_position(x = "Model", step.increase = 0.12)

cv.auroc.plot <- ggplot(auroc_table) +
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = AUROC_CV_mean, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = AUROC_CV_mean, fill = Model), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0) +
  scale_fill_manual(values = model.cols) +
  scale_y_continuous(breaks = seq(0.5,0.9,by=0.1), limits=c(0.5,0.83)) + 
  ylab("AUROC (CV)") +
  xlab(NULL) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))


cv.auroc.plot

# Save tables
path_to_save_tables <- "Tables/"
write.csv(auroc_table, file = paste0(path_to_save_tables, "CV_AUROC_spear_conv_clinical_additional_analyte_sets.csv"), row.names = F, quote = F)


# Save figure:
path_to_save_figs <- "Figs/"
path_to_save_rds <- "FigsRDS/"
#ggsave(plot = cv.auroc.plot, filename = paste0(path_to_save_figs, "CV_AUROC.pdf"), width = 3.5, height = 3)
saveRDS(cv.auroc.plot, file = paste0(path_to_save_rds, "CV_AUROC_spear_conv_clinical_additional_analyte_sets.rds"))


## Plot CV SPEAR models and clinical subtyping only per visit -------
all.res.df <- bind_rows(
                        acute_visit_models[["Visit 1"]]$spear$resdf,
                        acute_visit_models[["Visit 1"]]$spear_ensembl_clinical$resdf,
                        acute_visit_models[["Visit 2"]]$spear$resdf,
                        acute_visit_models[["Visit 2"]]$spear_ensembl_clinical$resdf,
                        acute_visit_models[["Visit 3"]]$spear$resdf,
                        acute_visit_models[["Visit 3"]]$spear_ensembl_clinical$resdf,
                        acute_visit_models[["Visit 4"]]$spear$resdf,
                        acute_visit_models[["Visit 4"]]$spear_ensembl_clinical$resdf,
                        acute_visit_models[["Visit 6"]]$spear$resdf,
                        acute_visit_models[["Visit 6"]]$spear_ensembl_clinical$resdf,
                        conv_visit_models[["Visit 7"]]$spear$resdf,
                        conv_visit_models[["Visit 7"]]$spear_ensembl_clinical$resdf,
                        conv_visit_models[["Visit 8"]]$spear$resdf,
                        conv_visit_models[["Visit 8"]]$spear_ensembl_clinical$resdf,
                        conv_visit_models[["Visit 9"]]$spear$resdf,
                        conv_visit_models[["Visit 9"]]$spear_ensembl_clinical$resdf,
                        conv_visit_models[["Visit 10"]]$spear$resdf,
                        conv_visit_models[["Visit 10"]]$spear_ensembl_clinical$resdf
)
all.res.df$Visit <- str_split_i(all.res.df$Model, "_", 1)
all.res.df$Type <- str_split_i(all.res.df$Model, "_", 2)
all.res.df$Visit <- factor(all.res.df$Visit, levels = c("Visit 1",
                                                        "Visit 2",
                                                        "Visit 3",
                                                        "Visit 4",
                                                        "Visit 6",
                                                        "Visit 7",
                                                        "Visit 8",
                                                        "Visit 9",
                                                        "Visit 10"
))
all.res.df$Type <- factor(all.res.df$Type, levels = c("SPEAR","SPEARClinical"))
type.cols = c(
               "#08605F", "#728D5B"
)
names(model.cols) = levels(all.res.df$Type)

auroc_table <- all.res.df %>%
  group_by(Model, Visit, Type, Seed) %>%
  summarise(AUROC_CV_mean = mean(CV), AUROC_CV_SD = sd(CV)) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()

stat.test <- auroc_table %>%
  group_by(Visit) %>%
  t_test(AUROC_CV_mean ~ Type, 
         paired = T
         ) %>%
  adjust_pvalue(method = "hochberg") %>%
  add_significance() %>%
  add_xy_position(x = "Visit", dodge = 0.75) 

cv.auroc.plot <- ggplot(auroc_table) +
  ggbeeswarm::geom_quasirandom(aes(x = Visit, y = AUROC_CV_mean, fill = Type), shape = 21, width = .15, dodge.width = 0.75) +
  geom_boxplot(aes(x = Visit, y = AUROC_CV_mean, fill = Type), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0) +
  scale_fill_manual(values = type.cols, labels = c("SPEAR","SPEAR + Clinical")) +
  scale_y_continuous(breaks = seq(0.5,0.9,by=0.1), limits=c(0.5,0.83)) + 
  labs(x="", y="AUROC (CV)") +
  theme(plot.title = element_text(hjust = .5))


cv.auroc.plot

# Save tables
path_to_save_tables <- "Tables/"
auroc_bootstrap_mean <- auroc_table %>%
  group_by(Model) %>%
  summarise(AUROC_bootstrap_mean = mean(AUROC_CV_mean), AUROC_bootstrap_SD = sd(AUROC_CV_mean)) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()
write.csv(auroc_bootstrap_mean %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_per_visit_clinical_comparison.csv"), row.names = F, quote = F)
write.csv(auroc_table %>% mutate(Model = str_remove_all(Model,"\n")), file = paste0(path_to_save_tables, "CV_AUROC_spear_per_visit_clinical_comparison_bootstraps.csv"), row.names = F, quote = F)


# Save figure:
path_to_save_figs <- "Figs/"
path_to_save_rds <- "FigsRDS/"
#ggsave(plot = cv.auroc.plot, filename = paste0(path_to_save_figs, "CV_AUROC.pdf"), width = 3.5, height = 3)
saveRDS(cv.auroc.plot, file = paste0(path_to_save_rds, "CV_AUROC_spear_per_visit_clinical_comparison.rds"))


## Plot CV SPEAR models with clinical subtyping additional per visit -------
all.res.df <- bind_rows(
  acute_visit_models[["Visit 1"]]$spear$resdf,
  acute_visit_models[["Visit 1"]]$spear_ensembl_clinical_additional$resdf,
  acute_visit_models[["Visit 2"]]$spear$resdf,
  acute_visit_models[["Visit 2"]]$spear_ensembl_clinical_additional$resdf,
  acute_visit_models[["Visit 3"]]$spear$resdf,
  acute_visit_models[["Visit 3"]]$spear_ensembl_clinical_additional$resdf,
  acute_visit_models[["Visit 4"]]$spear$resdf,
  acute_visit_models[["Visit 4"]]$spear_ensembl_clinical_additional$resdf,
  acute_visit_models[["Visit 6"]]$spear$resdf,
  acute_visit_models[["Visit 6"]]$spear_ensembl_clinical_additional$resdf,
  conv_visit_models[["Visit 7"]]$spear$resdf,
  conv_visit_models[["Visit 7"]]$spear_ensembl_clinical_additional$resdf,
  conv_visit_models[["Visit 8"]]$spear$resdf,
  conv_visit_models[["Visit 8"]]$spear_ensembl_clinical_additional$resdf,
  conv_visit_models[["Visit 9"]]$spear$resdf,
  conv_visit_models[["Visit 9"]]$spear_ensembl_clinical_additional$resdf,
  conv_visit_models[["Visit 10"]]$spear$resdf,
  conv_visit_models[["Visit 10"]]$spear_ensembl_clinical_additional$resdf
)
all.res.df$Visit <- str_split_i(all.res.df$Model, "_", 1)
all.res.df$Type <- str_split_i(all.res.df$Model, "_", 2)
all.res.df$Visit <- factor(all.res.df$Visit, levels = c("Visit 1",
                                                        "Visit 2",
                                                        "Visit 3",
                                                        "Visit 4",
                                                        "Visit 6",
                                                        "Visit 7",
                                                        "Visit 8",
                                                        "Visit 9",
                                                        "Visit 10"
))
all.res.df$Type <- factor(all.res.df$Type, levels = c("SPEAR","SPEARClinicalAdditional"))
type.cols = c(
  "#08605F", "#728D5B"
)
names(model.cols) = levels(all.res.df$Type)

auroc_table <- all.res.df %>%
  group_by(Model, Visit, Type, Seed) %>%
  summarise(AUROC_CV_mean = mean(CV), AUROC_CV_SD = sd(CV)) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()

stat.test <- auroc_table %>%
  group_by(Visit) %>%
  t_test(AUROC_CV_mean ~ Type, 
         paired = T
  ) %>%
  adjust_pvalue(method = "hochberg") %>%
  add_significance() %>%
  add_xy_position(x = "Visit", dodge = 0.75) 

cv.auroc.plot <- ggplot(auroc_table) +
  ggbeeswarm::geom_quasirandom(aes(x = Visit, y = AUROC_CV_mean, fill = Type), shape = 21, width = .15, dodge.width = 0.75) +
  geom_boxplot(aes(x = Visit, y = AUROC_CV_mean, fill = Type), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0) +
  scale_fill_manual(values = type.cols, labels = c("SPEAR","SPEAR + \nClinical additional")) +
  scale_y_continuous(breaks = seq(0.5,0.9,by=0.1), limits=c(0.5,0.83)) + 
  labs(x="", y="AUROC (CV)") +
  theme(plot.title = element_text(hjust = .5))

cv.auroc.plot

# Save tables
path_to_save_tables <- "Tables/"
auroc_table <- all.res.df %>%
  mutate(Model = str_remove_all(Model,"\n")) %>%
  group_by(Model) %>%
  summarise(AUROC_CV_mean = mean(CV), AUROC_CV_SD = sd(CV))
write.csv(auroc_table, file = paste0(path_to_save_tables, "CV_AUROC_spear_per_visit_clinical_additional_comparison.csv"), row.names = F, quote = F)


# Save figure:
path_to_save_figs <- "Figs/"
path_to_save_rds <- "FigsRDS/"
#ggsave(plot = cv.auroc.plot, filename = paste0(path_to_save_figs, "CV_AUROC.pdf"), width = 3.5, height = 3)
saveRDS(cv.auroc.plot, file = paste0(path_to_save_rds, "CV_AUROC_spear_per_visit_clinical_additional_comparison.rds"))



#TODO update test data
# all.preds.te <- rbind(spear_list$testscores, clinical_data$testscores, ensemble$testscores)
# 
# CV_AUROC_CURVE <- dplyr::group_by(all.preds.te, Fold, Model) %>%
#   yardstick::roc_curve(truth, MIN) %>% 
#   ggplot2::autoplot()
# CV_AUROC_CURVE <- CV_AUROC_CURVE + 
#   scale_color_manual(values = fold.model.cols) +
#   facet_wrap(vars(Fold)) +
#   guides(color = "none")
# CV_AUROC_CURVE
# saveRDS(CV_AUROC_CURVE, file = paste0(path_to_save_rds, "CV_AUROC_fold_curves_spear_clinical.rds"))





# AUROC Test --------------


# Train a lasso classifier using SPEAR(Phys) Factor 1
# NOTE: need to pass 2-dimensional matrix (just use Factor 1 copied twice)
#       and we can confirm that the second column is zeroed out

# Get scores on test data

all_train_test_data_lasso_AUROC <- function(X_train,
                                            X_test,
                                            Y_train,
                                            Y_test,
                                            Y_train_full,
                                            Y_test_full,
                                            foldid,
                                            model_name="model1") {
  
  
  # Train
  lasso_model <- glmnet::cv.glmnet(x = X_train, y = Y_train, foldid = foldid, family = "multinomial")
  # Evaluate
  lasso.preds.tr = stats::predict(lasso_model, X_train, s = "lambda.min", type="response")[,,1]
  lasso.preds.te = stats::predict(lasso_model, X_test, s = "lambda.min", type="response")[,,1]
  # Convert to df:
  lasso.preds.tr <- as.data.frame(lasso.preds.tr)
  colnames(lasso.preds.tr) <- c("MIN", "OTHER")
  lasso.preds.tr$estimate = c("MIN", "OTHER")[apply(lasso.preds.tr, 1, which.max)]
  lasso.preds.tr$truth = Y_train_full
  lasso.preds.tr <- lasso.preds.tr %>% dplyr::select(truth, estimate, dplyr::everything())
  
  lasso.preds.te <- as.data.frame(lasso.preds.te)
  colnames(lasso.preds.te) <- c("MIN", "OTHER")
  lasso.preds.te$estimate = c("MIN", "OTHER")[apply(lasso.preds.te, 1, which.max)]
  lasso.preds.te$truth = Y_test_full
  lasso.preds.te <- lasso.preds.te %>% dplyr::select(truth, estimate, dplyr::everything())
  
  lasso.preds.tr$truth <- factor(lasso.preds.tr$truth)
  lasso.preds.te$truth <- factor(lasso.preds.te$truth)
  # Assess:
  in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
  cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate
  
  # Add model name
  lasso.preds.tr$Model <- model_name
  lasso.preds.te$Model <- model_name
  
  return(list(train_scores = lasso.preds.tr, test_scores = lasso.preds.te, AUROC_train = in.sample, AUROC_test = cv.sample))
}

# # Spear model
# SPEARphys.test.fs <- SPEARobj.phys$get.factor.scores(data = "test")
# 
# X.tr <- SPEARobj.phys$get.factor.scores()
# X.tr <- cbind(X.tr[,1], X.tr[,1])
# X.te <- SPEARphys.test.fs
# X.te <- cbind(X.te[,1], X.te[,1])
# Y.tr <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", 0, 1)
# Y.te <- ifelse(SPEARobj.phys$data$test$pro_group_labels == "MIN", 0, 1)
# Y.tr.full <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", "MIN", "OTHER")
# Y.te.full <- ifelse(SPEARobj.phys$data$test$pro_group_labels == "MIN", "MIN", "OTHER")
# foldid <- SPEARobj.phys$params$fold.ids
# 
# spear_all_scores <- all_train_test_data_lasso_AUROC(X_train = X.tr,
#                                                     X_test = X.te,
#                                                     Y_train = Y.tr,
#                                                     Y_test = Y.te,
#                                                     Y_train_full = Y.tr.full,
#                                                     Y_test_full = Y.te.full,
#                                                     foldid = foldid,
#                                                     model_name = "SPEAR\nPhysical")
# 
# # Clinical model
# 
# # Order in same order as SPEAR object
# clinical_vars_test_ordered <- clinical_vars_test[rownames(SPEARphys.test.fs),]
# 
# clinical_all_scores <- all_train_test_data_lasso_AUROC(X_train = clinical_vars_ordered,
#                                                        X_test = clinical_vars_test_ordered,
#                                                        Y_train = Y.tr,
#                                                        Y_test = Y.te,
#                                                        Y_train_full = Y.tr.full,
#                                                        Y_test_full = Y.te.full,
#                                                        foldid = foldid,
#                                                        model_name = "Clinical"
# )
# 
# 
# ## Ensemble model
# ensemble_features_test <- cbind(SPEARphys.test.fs, clinical_vars_test_ordered)
# ensemble_all_scores <- all_train_test_data_lasso_AUROC(X_train = ensembl_features,
#                                                        X_test = ensemble_features_test,
#                                                        Y_train = Y.tr,
#                                                        Y_test = Y.te,
#                                                        Y_train_full = Y.tr.full,
#                                                        Y_test_full = Y.te.full,
#                                                        foldid = foldid,
#                                                        model_name = "Ensemble"
# )
# 
# # AUROC:
# all.preds <- bind_rows(spear_all_scores$test_scores, clinical_all_scores$test_scores, ensemble_all_scores$test_scores)
# 
# test.auroc <- all.preds %>% 
#               group_by(Model) %>% 
#               yardstick::roc_curve(truth, MIN) %>%
#               autoplot()
# 
# 
# test.auroc <- test.auroc +
#   ggplot2::ggtitle(("AUROC (Test)"), 
#                    paste0("Spear Physical | AUC=", round(spear_all_scores$AUROC_test, 3), "\n",
#                   "Clinical | AUC=", round(clinical_all_scores$AUROC_test, 3), "\n",
#                    "Ensemble | AUC=", round(ensemble_all_scores$AUROC_test, 3))
#                    ) +
#   ggplot2::xlab("FPR") +
#   ggplot2::ylab("TPR")
# test.auroc
# path_to_save_figs <- "Figs/"
# #ggsave(plot = test.auroc, filename = paste0(path_to_save_figs, "Test_AUROC_Full.pdf"), width = 3, height = 3)
# saveRDS(test.auroc, file = paste0(path_to_save_rds, "Test_AUROC_Full.RDS"))
# 
# # Precision Recall:
# test.precision.recall <- yardstick::pr_curve(data = lasso.preds.te, truth, MIN) %>% 
#   ggplot2::autoplot()
# test.precision.recall <- test.precision.recall +
#   ggplot2::ggtitle(("Precision-Recall (Test)"), "Male + Female") +
#   ggplot2::xlab("Recall") +
#   ggplot2::ylab("Precision")
# path_to_save_figs <- "Figs/"
# #ggsave(plot = test.precision.recall, filename = paste0(path_to_save_figs, "Test_PR_Full.pdf"), width = 3, height = 3)
# saveRDS(test.precision.recall, file = paste0(path_to_save_rds, "Test_PR_Full.RDS"))
# 
