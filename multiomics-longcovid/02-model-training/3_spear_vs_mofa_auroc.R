library(yardstick)
library(MultiAssayExperiment)
library(tidyverse)
select <- dplyr::select
#remotes::install_bitbucket("kleinstein/SPEAR@main")
library(SPEAR)
predict <- stats::predict


# CV Data: AUROC for all Models ----

# Loading MOFA and SPEAR objects...

MOFAobj <- readRDS("/scratch/data-fullcohort/MOFA_factors/MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.rds")
MOFA.fs <- MOFA2::get_factors(MOFAobj)$group1

SPEARobj <- readRDS("/scratch/data-fullcohort/spear_model/results/020624/021624_binomial_SPEAR.rds")
SPEARobj$set.weights(method = "min")
SPEAR.fs <- SPEARobj$get.factor.scores()
SPEAR.cv.fs <- SPEARobj$get.factor.scores(cv = TRUE)
f.ids <- SPEARobj$params$fold.ids
pro_groups <- SPEARobj$data$train$pro_group_labels
pro_groups_encoded <- ifelse(SPEARobj$data$train$pro_group_labels == "MIN", 0, 1)

SPEARobj.phys <- readRDS("/scratch/data-fullcohort/spear_model/results/020624/021924_physical_gaussianNEW_SPEAR.rds")
SPEARobj.phys$set.weights(method = "min")
SPEARphys.fs <- SPEARobj.phys$get.factor.scores()
SPEARphys.cv.fs <- SPEARobj.phys$get.factor.scores(cv = TRUE)
fphys.ids <- SPEARobj.phys$params$fold.ids
physpro_groups <- factor(ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", "MIN", "OTHER"), levels = c("MIN", "OTHER"))
physpro_groups_encoded <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", 0, 1)
tmp <- SPEARobj.phys$data$test@ExperimentList$PPG
tmp[is.na(tmp)] <- 0
SPEARobj.phys$data$test@ExperimentList$PPG <- tmp

data_env_train <- readRDS("/scratch/data-fullcohort/batch-effects/impacc_convalescent_data_env_train_Feb_11_2024.RDS")
clinical_data_train <- data_env_train$clinical_data
spear_object <- readRDS("/scratch/data-fullcohort/spear_model/results/020624/032124_spearphysical_finalscores_allsamples.rds")

spear_convalescent_train <- as.data.frame(rbind(spear_object$Train, spear_object$TrainNoPhysical_Overlap, spear_object$TrainNoPhysical_NoOverlap))
#spear_convalescent_train <- spear_convalescent_train %>% select(Factor1)

spear_convalescent_train_labels <- clinical_data_train %>% 
  select(event_id, pro_2_groups) %>%
  distinct() %>%
  filter(event_id %in% rownames(spear_convalescent_train)) %>%
  mutate(label_encoded = if_else(pro_2_groups == "MIN", 0 , 1)) %>%
  column_to_rownames("event_id")

spear_convalescent_train_labels_encoded <- spear_convalescent_train_labels[rownames(spear_convalescent_train),c("label_encoded")]
spear_convalescent_train_labels_full <- spear_convalescent_train_labels[rownames(spear_convalescent_train),c("pro_2_groups")]


# Male
SPEARobj.male <- readRDS("/scratch/data-fullcohort/spear_model/results/050124/050124_physical_Male_SPEAR.rds")
SPEARobj.male$set.weights(method = "min")
tmp <- SPEARobj.male$data$test@ExperimentList$PPG
tmp[is.na(tmp)] <- 0
SPEARobj.male$data$test@ExperimentList$PPG <- tmp
SPEARobj.male$data$train$pro_simplified <- ifelse(SPEARobj.male$data$train$pro_group_labels == "MIN", "MIN", "OTHER")
SPEARobj.male$data$test$pro_simplified <- ifelse(SPEARobj.male$data$test$pro_group_labels == "MIN", "MIN", "OTHER")

# Female
SPEARobj.female <- readRDS("/scratch/data-fullcohort/spear_model/results/050124/050124_physical_Female_SPEAR.rds")
SPEARobj.female$set.weights(method = "min")
tmp <- SPEARobj.female$data$test@ExperimentList$PPG
tmp[is.na(tmp)] <- 0
SPEARobj.female$data$test@ExperimentList$PPG <- tmp
SPEARobj.female$data$train$pro_simplified <- ifelse(SPEARobj.female$data$train$pro_group_labels == "MIN", "MIN", "OTHER")
SPEARobj.female$data$test$pro_simplified <- ifelse(SPEARobj.female$data$test$pro_group_labels == "MIN", "MIN", "OTHER")


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

# MOFA, SPEAR LC and SPEAR Physical Train CV AUROC scores loop  -------
# Run CV Loop 3 times (MOFA, SPEAR(Pro) and SPEAR(Phys))
# Calculate AUROC for each fold (CV = left out of training)



all.res.df <- data.frame()
all.mofa.preds.te <- data.frame()
for(f in 1:10){
  print(paste0("Fold", f, "..."))
  # Get ids
  train_ids <- f != f.ids
  test_ids <- f == f.ids
  f.ids.tmp <- f.ids[train_ids]
  f.ids.tmp[f.ids.tmp > f] = f.ids.tmp[f.ids.tmp > f] - 1
  # Subset:
  X.tr <- MOFA.fs[train_ids,]
  X.te <- MOFA.fs[test_ids,]
  Y.tr <- pro_groups_encoded[train_ids]
  Y.te <- pro_groups_encoded[test_ids]
  Y.tr.full <- pro_groups[train_ids]
  Y.te.full <- pro_groups[test_ids]
  # Train
  lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = f.ids.tmp, family = "multinomial")
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
  all.mofa.preds.te <- rbind(all.mofa.preds.te, lasso.preds.te)
  # Assess:
  in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
  cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate

  all.res.df <- rbind(all.res.df, data.frame(
    Model = "MOFA\n(Unsupervised)",
    Fold = f,
    Train = in.sample,
    CV = cv.sample
  ))
}

# SPEAR:
all.spear.preds.te <- data.frame()
for(f in 1:10){
  print(paste0("Fold", f, "..."))
  # Get ids
  train_ids <- f != f.ids
  test_ids <- f == f.ids
  f.ids.tmp <- f.ids[train_ids]
  f.ids.tmp[f.ids.tmp > f] = f.ids.tmp[f.ids.tmp > f] - 1
  # Subset:
  X.tr <- SPEAR.cv.fs[train_ids,]
  X.te <- SPEAR.cv.fs[test_ids,]
  Y.tr <- pro_groups_encoded[train_ids]
  Y.te <- pro_groups_encoded[test_ids]
  Y.tr.full <- pro_groups[train_ids]
  Y.te.full <- pro_groups[test_ids]
  # Train
  lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = f.ids.tmp, family = "multinomial")
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
  all.spear.preds.te <- rbind(all.spear.preds.te, lasso.preds.te)
  # Assess:
  in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
  cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate
  
  all.res.df <- rbind(all.res.df, data.frame(
    Model = "SPEAR\nLC",
    Fold = f,
    Train = in.sample,
    CV = cv.sample
  ))
}
# SPEAR Physical:
fig.2.a.all.spearphys.preds.te <- data.frame()
for(f in 1:10){
  print(paste0("Fold", f, "..."))
  # Get ids
  train_ids <- f != fphys.ids
  test_ids <- f == fphys.ids
  f.ids.tmp <- fphys.ids[train_ids]
  f.ids.tmp[f.ids.tmp > f] = f.ids.tmp[f.ids.tmp > f] - 1
  # Subset:
  X.tr <- SPEARphys.cv.fs[train_ids,]
  X.te <- SPEARphys.cv.fs[test_ids,]
  Y.tr <- physpro_groups_encoded[train_ids]
  Y.te <- physpro_groups_encoded[test_ids]
  Y.tr.full <- physpro_groups[train_ids]
  Y.te.full <- physpro_groups[test_ids]
  # Train
  lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = f.ids.tmp, family = "multinomial")
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
  fig.2.a.all.spearphys.preds.te <- rbind(fig.2.a.all.spearphys.preds.te, lasso.preds.te)
  # Assess:
  in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
  cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate
  
  all.res.df <- rbind(all.res.df, data.frame(
    Model = "SPEAR\nPhysical",
    Fold = f,
    Train = in.sample,
    CV = cv.sample
  ))
}
all.res.df$Model <- factor(all.res.df$Model, levels = c("MOFA\n(Unsupervised)", "SPEAR\nLC", "SPEAR\nPhysical"))

# Set model colors:
model.cols = c("#73C6B6", "#C39BD3", "#8E44AD")
names(model.cols) = c("MOFA\n(Unsupervised)", "SPEAR\nLC", "SPEAR\nPhysical")
fold.model.cols = c(rep(model.cols[1], 10),
                    rep(model.cols[2], 10),
                    rep(model.cols[3], 10))
names(fold.model.cols) = c(paste0(1:10, "_", names(model.cols)[1]),
                           paste0(1:10, "_", names(model.cols)[2]),
                           paste0(1:10, "_", names(model.cols)[3]))

cv.auroc.plot <- ggplot(all.res.df) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = CV, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = CV, fill = Model), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  scale_fill_manual(values = model.cols) +
  ylab("AUROC (CV)") +
  ylim(c(0.45, NA)) +
  xlab(NULL) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))

cv.auroc.plot

# Save figure:
path_to_save_figs <- "Figs/"
path_to_save_rds <- "FigsRDS/"
#ggsave(plot = cv.auroc.plot, filename = paste0(path_to_save_figs, "CV_AUROC.pdf"), width = 3.5, height = 3)
saveRDS(cv.auroc.plot, file = paste0(path_to_save_rds, "CV_AUROC_MOFA_SPEARLC_SPEARPhys_nobootstrap.rds"))

# Bootstrapped MOFA, SPEAR LC and SPEAR Physical Train CV AUROC scores loop  -------
mofa_model <- cv_scores_lasso_model_min_other(X=MOFA.fs,
                                              model_name = "MOFA\n(Unsupervised)",
                                              Y=pro_groups_encoded,
                                              Y_full=pro_groups,
                                              n_folds = 10,
                                              bootstrap_sample = 100)
SPEAR_LC <- cv_scores_lasso_model_min_other(X=SPEAR.cv.fs,
                                            model_name = "SPEAR\nLC",
                                            Y=pro_groups_encoded,
                                            Y_full=pro_groups,
                                            n_folds = 10,
                                            bootstrap_sample = 100)
SPEAR_physical <- cv_scores_lasso_model_min_other(X=as.matrix(spear_convalescent_train),
                                                  model_name = "SPEAR\nPhysical",
                                                  Y=spear_convalescent_train_labels_encoded,
                                                  Y_full=as.factor(spear_convalescent_train_labels_full),
                                                  n_folds = 10,
                                                  bootstrap_sample = 100)

path_to_save_models <- "Models/"
saveRDS(mofa_model, paste0(path_to_save_models,"MOFA_model.rds"))
saveRDS(SPEAR_LC, paste0(path_to_save_models,"SPEAR_LC_model.rds"))
saveRDS(SPEAR_physical, paste0(path_to_save_models,"SPEAR_Physical_model.rds"))

mofa_model <- readRDS( paste0(path_to_save_models,"MOFA_model.rds"))
SPEAR_LC <- readRDS(paste0(path_to_save_models,"SPEAR_LC_model.rds"))
SPEAR_physical <- readRDS(paste0(path_to_save_models,"SPEAR_Physical_model.rds"))

all.res.df <- bind_rows(mofa_model$resdf, SPEAR_LC$resdf, SPEAR_physical$resdf)
all.res.df$Model <- factor(all.res.df$Model, levels = c("MOFA\n(Unsupervised)", "SPEAR\nLC", "SPEAR\nPhysical"))


auroc_table <- all.res.df %>%
  group_by(Model, Seed) %>%
  summarise(AUROC_CV_mean = mean(CV), AUROC_CV_SD = sd(CV)) %>%
  mutate(Model = as.factor(Model)) %>%
  ungroup()

stat.test <- auroc_table %>%
  t_test(AUROC_CV_mean ~ Model, 
         paired = T,
         comparisons = list(c("MOFA\n(Unsupervised)", "SPEAR\nLC"),
                            c("SPEAR\nLC", "SPEAR\nPhysical"),
                            c("MOFA\n(Unsupervised)", "SPEAR\nPhysical")
         )) %>%
  add_significance() %>%
  add_xy_position(x = "Model")

cv.auroc.plot <- ggplot(auroc_table) +
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = AUROC_CV_mean, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = AUROC_CV_mean, fill = Model), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0) +
  scale_fill_manual(values = model.cols) +
  scale_y_continuous(breaks = seq(0.6,0.8,by=0.05), limits=c(0.65,0.8)) + 
  ylab("AUROC (CV)") +
  xlab(NULL) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))


cv.auroc.plot

# Save figure:
path_to_save_figs <- "Figs/"
path_to_save_rds <- "FigsRDS/"
#ggsave(plot = cv.auroc.plot, filename = paste0(path_to_save_figs, "CV_AUROC.pdf"), width = 3.5, height = 3)
saveRDS(cv.auroc.plot, file = paste0(path_to_save_rds, "CV_AUROC_withbootstrapping.rds"))


all.mofa.preds.te$Model <- "MOFA\n(Unsupervised)"
all.spear.preds.te$Model <- "SPEAR\nLC"
fig.2.a.all.spearphys.preds.te$Model <- "SPEAR\nPhysical"
all.preds.te <- rbind(all.mofa.preds.te, all.spear.preds.te, fig.2.a.all.spearphys.preds.te)

CV_AUROC_CURVE <- dplyr::group_by(all.preds.te, Fold, Model) %>%
  yardstick::roc_curve(truth, MIN) %>% 
  ggplot2::autoplot()
CV_AUROC_CURVE <- CV_AUROC_CURVE + 
  scale_color_manual(values = fold.model.cols) +
  facet_wrap(vars(Fold)) +
  guides(color = "none")
CV_AUROC_CURVE


# CV DATA: MSE to predict PROMIS ---------


# Run CV Loop 3 times (MOFA, SPEAR(Pro) and SPEAR(Phys))
# Calculate AUROC for each fold (CV = left out of training)

all.res.df <- data.frame()

# SPEAR Models:
other.models <- c("physical", "psycho", "mental", "impact", "dyspnea")
for(cur.mod in other.models){
  print(paste0("Doing ", cur.mod, "..."))
  if(cur.mod != "physical"){
    tmp.SPEARobj <- readRDS(file = paste0("/scratch/data-fullcohort/spear_model/results/020624/021924_", cur.mod, "_gaussian_SPEAR.rds"))
  } else {
    tmp.SPEARobj <- SPEARobj.phys
  }
  tmp.SPEARobj$set.weights(method = "min")
  SPEARgaussian.cv.fs <- tmp.SPEARobj$get.factor.scores(cv = TRUE)
  fphys.ids <- tmp.SPEARobj$params$fold.ids
  names(fphys.ids) <- colnames(tmp.SPEARobj$data$train@ExperimentList$SO)
  gaussian.response <- tmp.SPEARobj$data$train$response
  for(f in 1:max(fphys.ids)){
    print(paste0("Fold", f, "..."))
    # Get ids
    train_ids <- f != fphys.ids
    test_ids <- f == fphys.ids
    f.ids.tmp <- fphys.ids[train_ids]
    f.ids.tmp[f.ids.tmp > f] = f.ids.tmp[f.ids.tmp > f] - 1
    # Subset:
    X.tr <- SPEARgaussian.cv.fs[train_ids,]
    X.te <- SPEARgaussian.cv.fs[test_ids,]
    Y.tr <- gaussian.response[train_ids]
    Y.te <- gaussian.response[test_ids]
    # Train
    lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = f.ids.tmp)
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
    # Assess:
    in.sample <- yardstick::rmse(data = lasso.preds.tr, truth, estimate)$.estimate
    cv.sample <- yardstick::rmse(data = lasso.preds.te, truth, estimate)$.estimate
    
    label = if_else(cur.mod=="psycho","cognitive",cur.mod)
    all.res.df <- rbind(all.res.df, data.frame(
      Model = paste0("SPEAR"),
      Fold = f,
      Train = in.sample,
      CV = cv.sample,
      Scale = paste0("PROMIS ",str_to_title(label), " score")
    ))
    
    
    # Now repeat with MOFA:
    print(paste0("Fold", f, "..."))
    # Subset:
    X.tr <- MOFA.fs[names(train_ids[train_ids == TRUE]),]
    X.te <- MOFA.fs[names(test_ids[test_ids == TRUE]),]
    Y.tr <- gaussian.response[train_ids]
    Y.te <- gaussian.response[test_ids]
    # Train
    lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = f.ids.tmp)
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
    # Assess:
    in.sample <- yardstick::rmse(data = lasso.preds.tr, truth, estimate)$.estimate
    cv.sample <- yardstick::rmse(data = lasso.preds.te, truth, estimate)$.estimate
    
    label = if_else(cur.mod=="psycho","cognitive",cur.mod)
    all.res.df <- rbind(all.res.df, data.frame(
      Model = paste0("MOFA"),
      Fold = f,
      Train = in.sample,
      CV = cv.sample,
      Scale = paste0("PROMIS ",str_to_title(label), " score")
    ))
  }
}

all.res.df$Model <- factor(all.res.df$Model, levels = unique(all.res.df$Model))
all.res.df$Scale <- factor(all.res.df$Scale, levels = c("PROMIS Physical score",
                                                          "PROMIS Cognitive score",
                                                          "PROMIS Mental score",
                                                          "PROMIS Impact score",
                                                          "PROMIS Dyspnea score"))
# Set model colors:
model.cols = rep(c("#8E44AD", "#73C6B6"), length(other.models))
names(model.cols) = unique(all.res.df$Model)


# Make plot:
cv.mse.plot <- ggplot(all.res.df) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = CV, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = CV, fill = Model), outlier.shape=NA,  alpha = .6, color = "black") +
  facet_wrap(~Scale, ncol = 5, scales = "free_y") +
  #geom_vline(xintercept = c(2, 4, 6, 8) + .5, color = "grey") +
  theme_bw() + theme(strip.background = element_rect(fill="white")) +
  scale_fill_manual(values = model.cols) +
  scale_y_continuous(limits = c(0.75,1.25)) +
  ylab("Model Error (RMSE)") +
  #ylim(c(0.5, 1)) +
  xlab(NULL) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))

cv.mse.plot


# Save figure:
path_to_save_figs <- "Figs/"
#ggsave(plot = cv.mse.plot, filename = paste0(path_to_save_figs, "CV_MSE.pdf"), width = 6, height = 3)
saveRDS(cv.mse.plot, file = paste0(path_to_save_rds, "CV_MSE.RDS"))

# CV DATA: ALL Promis models predicting PRO ------

# Run CV Loop 3 times (MOFA, SPEAR(Pro) and SPEAR(Phys))
# Calculate AUROC for each fold (CV = left out of training)

all.res.df <- data.frame()
# SPEAR Models:
other.models <- c("physical", "psycho", "mental", "impact", "dyspnea")
for(cur.mod in other.models){
  print(paste0("Doing ", cur.mod, "..."))
  if(cur.mod != "physical"){
    tmp.SPEARobj <- readRDS(file = paste0("/scratch/data-fullcohort/spear_model/results/020624/021924_", cur.mod, "_gaussian_SPEAR.rds"))
  } else {
    tmp.SPEARobj <- SPEARobj.phys
  }
  tmp.SPEARobj$set.weights(method = "min")
  SPEARgaussian.cv.fs <- tmp.SPEARobj$get.factor.scores(cv = TRUE)
  fphys.ids <- tmp.SPEARobj$params$fold.ids
  names(fphys.ids) <- colnames(tmp.SPEARobj$data$train@ExperimentList$SO)
  pro_groups <- factor(ifelse(tmp.SPEARobj$data$train$pro_group_labels == "MIN", "MIN", "OTHER"), levels = c("MIN", "OTHER"))
  names(pro_groups) <- rownames(SPEARgaussian.cv.fs)
  pro_groups_encoded <- ifelse(tmp.SPEARobj$data$train$pro_group_labels == "MIN", 0, 1)
  names(pro_groups_encoded) <- rownames(SPEARgaussian.cv.fs)
  for(f in 1:max(fphys.ids)){
    print(paste0("Fold", f, "..."))
    # Get ids
    train_ids <- f != fphys.ids
    test_ids <- f == fphys.ids
    f.ids.tmp <- fphys.ids[train_ids]
    f.ids.tmp[f.ids.tmp > f] = f.ids.tmp[f.ids.tmp > f] - 1
    # Subset:
    X.tr <- SPEARgaussian.cv.fs[train_ids,]
    X.te <- SPEARgaussian.cv.fs[test_ids,]
    Y.tr <- pro_groups_encoded[train_ids]
    Y.te <- pro_groups_encoded[test_ids]
    Y.tr.full <- pro_groups[train_ids]
    Y.te.full <- pro_groups[test_ids]
    # Train
    lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = f.ids.tmp, family = "multinomial")
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
    # Assess:
    in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
    cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate
    
    all.res.df <- rbind(all.res.df, data.frame(
      Model = paste0("SPEAR\n(", cur.mod, ")"),
      Fold = f,
      Train = in.sample,
      CV = cv.sample
    ))
    
    
    # Now repeat with MOFA:
    print(paste0("Fold", f, "..."))
    # Get ids
    train_ids <- f != fphys.ids
    test_ids <- f == fphys.ids
    f.ids.tmp <- fphys.ids[train_ids]
    f.ids.tmp[f.ids.tmp > f] = f.ids.tmp[f.ids.tmp > f] - 1
    # Subset:
    X.tr <- MOFA.fs[names(train_ids[train_ids == TRUE]),]
    X.te <- MOFA.fs[names(test_ids[test_ids == TRUE]),]
    Y.tr <- pro_groups_encoded[train_ids]
    Y.te <- pro_groups_encoded[test_ids]
    Y.tr.full <- pro_groups[train_ids]
    Y.te.full <- pro_groups[test_ids]
    # Train
    lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = f.ids.tmp, family = "multinomial")
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
    # Assess:
    in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
    cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate
    
    all.res.df <- rbind(all.res.df, data.frame(
      Model = paste0("MOFA\n(", cur.mod, ")"),
      Fold = f,
      Train = in.sample,
      CV = cv.sample
    ))
  }
}

all.res.df$Model <- factor(all.res.df$Model, levels = unique(all.res.df$Model))

# Set model colors:
model.cols = rep(c("#8E44AD", "#73C6B6"), length(other.models))
names(model.cols) = unique(all.res.df$Model)

# Make plot:
cv.auroc.total.plot <- ggplot(all.res.df) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  geom_hline(yintercept = 1, color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = CV, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = CV, fill = Model ), outlier.shape=NA, alpha = .6, color = "black") +
  geom_vline(xintercept = c(2, 4, 6, 8) + .5, color = "grey") +
  theme_bw() +
  scale_fill_manual(values = model.cols) +
  ylab("AUROC") +
  ylim(c(0.5, 1)) +
  xlab(NULL) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))

cv.auroc.total.plot

# Save figure:
path_to_save_figs <- "Figs/"
#ggsave(plot = cv.auroc.total.plot, filename = paste0(path_to_save_figs, "CV_TOTAL_AUROC.pdf"), width = 6, height = 3)
saveRDS(cv.auroc.total.plot, file = paste0(path_to_save_rds, "CV_TOTAL_AUROC.RDS"))


# TEST DATA: ALL MODELS: AUROC and PR CURVE -----

# Train a lasso classifier using SPEAR(Phys) Factor 1
# NOTE: need to pass 2-dimensional matrix (just use Factor 1 copied twice)
#       and we can confirm that the second column is zeroed out

all.res.df <- data.frame()
# SPEAR Models:
other.models <- c("physical", "psycho", "mental", "impact", "dyspnea")
for(cur.mod in other.models){
  print(paste0("Doing ", cur.mod, "..."))
  if(cur.mod != "physical"){
    tmp.SPEARobj <- readRDS(file = paste0("/scratch/data-fullcohort/spear_model/results/020624/021924_", cur.mod, "_gaussian_SPEAR.rds"))
    tmp <- tmp.SPEARobj$data$test@ExperimentList$PPG
    tmp[is.na(tmp)] <- 0
    tmp.SPEARobj$data$test@ExperimentList$PPG <- tmp
  } else {
    tmp.SPEARobj <- SPEARobj.phys
  }
  tmp.SPEARobj$set.weights(method = "min")
  fphys.ids <- tmp.SPEARobj$params$fold.ids
  names(fphys.ids) <- colnames(tmp.SPEARobj$data$train@ExperimentList$SO)
  pro_groups <- factor(ifelse(tmp.SPEARobj$data$train$pro_group_labels == "MIN", "MIN", "OTHER"), levels = c("MIN", "OTHER"))
  names(pro_groups) <- rownames(SPEARgaussian.cv.fs)
  pro_groups_encoded <- ifelse(tmp.SPEARobj$data$train$pro_group_labels == "MIN", 0, 1)
  names(pro_groups_encoded) <- rownames(SPEARgaussian.cv.fs)
  # Subset:
  X.tr <- tmp.SPEARobj$get.factor.scores(data = "train")
  X.tr <- cbind(X.tr[,1], X.tr[,1])
  X.te <- tmp.SPEARobj$get.factor.scores(data = "test")
  X.te <- cbind(X.te[,1], X.te[,1])
  Y.tr <- ifelse(tmp.SPEARobj$data$train$pro_group_labels == "MIN", 0, 1)
  Y.te <- ifelse(tmp.SPEARobj$data$test$pro_group_labels == "MIN", 0, 1)
  Y.tr.full <- ifelse(tmp.SPEARobj$data$train$pro_group_labels == "MIN", "MIN", "OTHER")
  Y.te.full <- ifelse(tmp.SPEARobj$data$test$pro_group_labels == "MIN", "MIN", "OTHER")
  # Train
  lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = tmp.SPEARobj$params$fold.ids, family = "multinomial")
  # Evaluate
  lasso.preds.tr = stats::predict(lasso_model, X.tr, s = "lambda.min", type="response")[,,1]
  lasso.preds.te = stats::predict(lasso_model, X.te, s = "lambda.min", type="response")[,,1]
  # Convert to df:
  lasso.preds.tr <- as.data.frame(lasso.preds.tr)
  colnames(lasso.preds.tr) <- c("MIN", "OTHER")
  lasso.preds.tr$estimate = c("MIN", "OTHER")[apply(lasso.preds.tr, 1, which.max)]
  lasso.preds.tr$truth = Y.tr.full
  lasso.preds.tr <- lasso.preds.tr %>% dplyr::select(truth, estimate, dplyr::everything())
  lasso.preds.te <- as.data.frame(lasso.preds.te)
  colnames(lasso.preds.te) <- c("MIN", "OTHER")
  lasso.preds.te$estimate = c("MIN", "OTHER")[apply(lasso.preds.te, 1, which.max)]
  lasso.preds.te$truth = Y.te.full
  lasso.preds.te <- lasso.preds.te %>% dplyr::select(truth, estimate, dplyr::everything())
  lasso.preds.tr$truth <- factor(lasso.preds.tr$truth)
  lasso.preds.te$truth <- factor(lasso.preds.te$truth)
  # Assess:
  in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
  cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate
  # AUROC:
  test.auroc <- yardstick::roc_curve(data = lasso.preds.te, truth, MIN) %>% 
    ggplot2::autoplot()
  test.auroc <- test.auroc +
    ggplot2::ggtitle(("AUROC (Test)"), paste0("SPEAR ",str_to_title(cur.mod), " | AUC=", round(cv.sample, 2))) +
    ggplot2::xlab("FPR") +
    ggplot2::ylab("TPR")
  path_to_save_figs <- "Figs/"
  #ggsave(plot = test.auroc, filename = paste0(path_to_save_figs, "Test_AUROC_", cur.mod, ".pdf"), width = 3, height = 3)
  saveRDS(test.auroc, file = paste0(path_to_save_rds, "Test_AUROC_", cur.mod, ".RDS"))
  
  # Precision Recall:
  test.precision.recall <- yardstick::pr_curve(data = lasso.preds.te, truth, MIN) %>% 
    ggplot2::autoplot()
  test.precision.recall <- test.precision.recall +
    ggplot2::ggtitle(("Precision-Recall (Test)"), cur.mod) +
    ggplot2::xlab("Recall") +
    ggplot2::ylab("Precision")
  path_to_save_figs <- "Figs/"
  #ggsave(plot = test.precision.recall, filename = paste0(path_to_save_figs, "Test_PR_", cur.mod, ".pdf"), width = 3, height = 3)
  saveRDS(test.precision.recall, file = paste0(path_to_save_rds, "Test_PR_", cur.mod, ".RDS"))
}

# TEST DATA: AUROC and PR CURVE ------

# Train a lasso classifier using SPEAR(Phys) Factor 1
# NOTE: need to pass 2-dimensional matrix (just use Factor 1 copied twice)
#       and we can confirm that the second column is zeroed out

# Subset:
X.tr <- SPEARobj.phys$get.factor.scores()
X.tr <- cbind(X.tr[,1], X.tr[,1])
X.te <- SPEARobj.phys$get.factor.scores(data = "test")
X.te <- cbind(X.te[,1], X.te[,1])
Y.tr <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", 0, 1)
Y.te <- ifelse(SPEARobj.phys$data$test$pro_group_labels == "MIN", 0, 1)
Y.tr.full <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", "MIN", "OTHER")
Y.te.full <- ifelse(SPEARobj.phys$data$test$pro_group_labels == "MIN", "MIN", "OTHER")
# Train
lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = SPEARobj.phys$params$fold.ids, family = "multinomial")
# Evaluate
lasso.preds.tr = stats::predict(lasso_model, X.tr, s = "lambda.min", type="response")[,,1]
lasso.preds.te = stats::predict(lasso_model, X.te, s = "lambda.min", type="response")[,,1]
# Convert to df:
lasso.preds.tr <- as.data.frame(lasso.preds.tr)
colnames(lasso.preds.tr) <- c("MIN", "OTHER")
lasso.preds.tr$estimate = c("MIN", "OTHER")[apply(lasso.preds.tr, 1, which.max)]
lasso.preds.tr$truth = Y.tr.full
lasso.preds.tr <- lasso.preds.tr %>% dplyr::select(truth, estimate, dplyr::everything())
lasso.preds.te <- as.data.frame(lasso.preds.te)
colnames(lasso.preds.te) <- c("MIN", "OTHER")
lasso.preds.te$estimate = c("MIN", "OTHER")[apply(lasso.preds.te, 1, which.max)]
lasso.preds.te$truth = Y.te.full
lasso.preds.te <- lasso.preds.te %>% dplyr::select(truth, estimate, dplyr::everything())
lasso.preds.tr$truth <- factor(lasso.preds.tr$truth)
lasso.preds.te$truth <- factor(lasso.preds.te$truth)
# Assess:
in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate

# NOTE: Adding cv from fig.2.a.all.spearphys.preds.te, run first chunk above!
cv.preds <- fig.2.a.all.spearphys.preds.te
te.preds <- lasso.preds.te
cv.preds$Group <- "CV"
cv.preds$Fold <- NULL
cv.preds$Model <- NULL
te.preds$Group <- "Test"
all.preds <- rbind(cv.preds, te.preds)

# AUROC:
test.auroc <- yardstick::roc_curve(data = all.preds, truth, MIN) %>% 
  ggplot2::autoplot()
test.auroc <- test.auroc +
  ggplot2::ggtitle(("AUROC (Test)"), paste0("Full | AUC=", round(cv.sample, 2))) +
  ggplot2::xlab("FPR") +
  ggplot2::ylab("TPR")
path_to_save_figs <- "Figs/"
#ggsave(plot = test.auroc, filename = paste0(path_to_save_figs, "Test_AUROC_Full.pdf"), width = 3, height = 3)
saveRDS(test.auroc, file = paste0(path_to_save_rds, "Test_AUROC_Full.RDS"))

# Precision Recall:
test.precision.recall <- yardstick::pr_curve(data = lasso.preds.te, truth, MIN) %>% 
  ggplot2::autoplot()
test.precision.recall <- test.precision.recall +
  ggplot2::ggtitle(("Precision-Recall (Test)"), "Male + Female") +
  ggplot2::xlab("Recall") +
  ggplot2::ylab("Precision")
path_to_save_figs <- "Figs/"
#ggsave(plot = test.precision.recall, filename = paste0(path_to_save_figs, "Test_PR_Full.pdf"), width = 3, height = 3)
saveRDS(test.precision.recall, file = paste0(path_to_save_rds, "Test_PR_Full.RDS"))

# Male vs. Female (for supplemental figures)
male.ids.tr <- colnames(SPEARobj.male$data$train@ExperimentList$SO)
female.ids.tr <- colnames(SPEARobj.female$data$train@ExperimentList$SO)
all.ids.tr <- colnames(SPEARobj.phys$data$train@ExperimentList$SO)
male.ids.te <- colnames(SPEARobj.male$data$test@ExperimentList$SO)
female.ids.te <- colnames(SPEARobj.female$data$test@ExperimentList$SO)
all.ids.te <- colnames(SPEARobj.phys$data$test@ExperimentList$SO)

# Male:
X.tr <- SPEARobj.phys$get.factor.scores()
X.tr <- cbind(X.tr[,1], X.tr[,1])
X.tr <- X.tr[all.ids.tr %in% male.ids.tr,]
X.te <- SPEARobj.phys$get.factor.scores(data = "test")
X.te <- cbind(X.te[,1], X.te[,1])
X.te <- X.te[all.ids.te %in% male.ids.te,]
Y.tr <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", 0, 1)
Y.tr <- Y.tr[all.ids.tr %in% male.ids.tr]
Y.te <- ifelse(SPEARobj.phys$data$test$pro_group_labels == "MIN", 0, 1)
Y.te <- Y.te[all.ids.te %in% male.ids.te]
Y.tr.full <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", "MIN", "OTHER")
Y.tr.full <- Y.tr.full[all.ids.tr %in% male.ids.tr]
Y.te.full <- ifelse(SPEARobj.phys$data$test$pro_group_labels == "MIN", "MIN", "OTHER")
Y.te.full <- Y.te.full[all.ids.te %in% male.ids.te]
# Train
lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = SPEARobj.phys$params$fold.ids[all.ids.tr %in% male.ids.tr], family = "multinomial")
# Evaluate
lasso.preds.tr = stats::predict(lasso_model, X.tr, s = "lambda.min", type="response")[,,1]
lasso.preds.te = stats::predict(lasso_model, X.te, s = "lambda.min", type="response")[,,1]
# Convert to df:
lasso.preds.tr <- as.data.frame(lasso.preds.tr)
colnames(lasso.preds.tr) <- c("MIN", "OTHER")
lasso.preds.tr$estimate = c("MIN", "OTHER")[apply(lasso.preds.tr, 1, which.max)]
lasso.preds.tr$truth = Y.tr.full
lasso.preds.tr <- lasso.preds.tr %>% dplyr::select(truth, estimate, dplyr::everything())
lasso.preds.te <- as.data.frame(lasso.preds.te)
colnames(lasso.preds.te) <- c("MIN", "OTHER")
lasso.preds.te$estimate = c("MIN", "OTHER")[apply(lasso.preds.te, 1, which.max)]
lasso.preds.te$truth = Y.te.full
lasso.preds.te <- lasso.preds.te %>% dplyr::select(truth, estimate, dplyr::everything())
lasso.preds.tr$truth <- factor(lasso.preds.tr$truth)
lasso.preds.te$truth <- factor(lasso.preds.te$truth)
# Assess:
in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate
# AUROC:
test.male.auroc <- yardstick::roc_curve(data = lasso.preds.te, truth, MIN) %>% 
  ggplot2::autoplot()
test.male.auroc <- test.male.auroc +
  ggplot2::ggtitle(("AUROC (Test)"), paste0("Male | AUC=", round(cv.sample, 2))) +
  ggplot2::xlab("FPR") +
  ggplot2::ylab("TPR")
#ggsave(plot = test.male.auroc, filename = paste0(path_to_save_figs, "Test_AUROC_Male.pdf"), width = 3, height = 3)
saveRDS(test.male.auroc, file = paste0(path_to_save_rds, "Test_AUROC_Male.RDS"))

# Precision Recall:
test.male.precision.recall <- yardstick::pr_curve(data = lasso.preds.te, truth, MIN) %>% 
  ggplot2::autoplot()
test.male.precision.recall <- test.male.precision.recall +
  ggplot2::ggtitle(("Precision-Recall (Test)"), "Male") +
  ggplot2::xlab("Recall") +
  ggplot2::ylab("Precision")
#ggsave(plot = test.male.precision.recall, filename = paste0(path_to_save_figs, "Test_PR_Male.pdf"), width = 3, height = 3)
saveRDS(test.male.precision.recall, file = paste0(path_to_save_rds, "Test_PR_Male.RDS"))

# Female:
X.tr <- SPEARobj.phys$get.factor.scores()
X.tr <- cbind(X.tr[,1], X.tr[,1])
X.tr <- X.tr[all.ids.tr %in% female.ids.tr,]
X.te <- SPEARobj.phys$get.factor.scores(data = "test")
X.te <- cbind(X.te[,1], X.te[,1])
X.te <- X.te[all.ids.te %in% female.ids.te,]
Y.tr <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", 0, 1)
Y.tr <- Y.tr[all.ids.tr %in% female.ids.tr]
Y.te <- ifelse(SPEARobj.phys$data$test$pro_group_labels == "MIN", 0, 1)
Y.te <- Y.te[all.ids.te %in% female.ids.te]
Y.tr.full <- ifelse(SPEARobj.phys$data$train$pro_group_labels == "MIN", "MIN", "OTHER")
Y.tr.full <- Y.tr.full[all.ids.tr %in% female.ids.tr]
Y.te.full <- ifelse(SPEARobj.phys$data$test$pro_group_labels == "MIN", "MIN", "OTHER")
Y.te.full <- Y.te.full[all.ids.te %in% female.ids.te]
# Train
lasso_model <- glmnet::cv.glmnet(x = X.tr, y = Y.tr, foldid = SPEARobj.phys$params$fold.ids[all.ids.tr %in% female.ids.tr], family = "multinomial")
# Evaluate
lasso.preds.tr = stats::predict(lasso_model, X.tr, s = "lambda.min", type="response")[,,1]
lasso.preds.te = stats::predict(lasso_model, X.te, s = "lambda.min", type="response")[,,1]
# Convert to df:
lasso.preds.tr <- as.data.frame(lasso.preds.tr)
colnames(lasso.preds.tr) <- c("MIN", "OTHER")
lasso.preds.tr$estimate = c("MIN", "OTHER")[apply(lasso.preds.tr, 1, which.max)]
lasso.preds.tr$truth = Y.tr.full
lasso.preds.tr <- lasso.preds.tr %>% dplyr::select(truth, estimate, dplyr::everything())
lasso.preds.te <- as.data.frame(lasso.preds.te)
colnames(lasso.preds.te) <- c("MIN", "OTHER")
lasso.preds.te$estimate = c("MIN", "OTHER")[apply(lasso.preds.te, 1, which.max)]
lasso.preds.te$truth = Y.te.full
lasso.preds.te <- lasso.preds.te %>% dplyr::select(truth, estimate, dplyr::everything())
lasso.preds.tr$truth <- factor(lasso.preds.tr$truth)
lasso.preds.te$truth <- factor(lasso.preds.te$truth)
# Assess:
in.sample <- yardstick::roc_auc(data = lasso.preds.tr, truth, MIN)$.estimate
cv.sample <- yardstick::roc_auc(data = lasso.preds.te, truth, MIN)$.estimate
# AUROC:
test.female.auroc <- yardstick::roc_curve(data = lasso.preds.te, truth, MIN) %>% 
  ggplot2::autoplot()
test.female.auroc <- test.female.auroc +
  ggplot2::ggtitle(("AUROC (Test)"), paste0("Female | AUC=", round(cv.sample, 2))) +
  ggplot2::xlab("FPR") +
  ggplot2::ylab("TPR")
#ggsave(plot = test.female.auroc, filename = paste0(path_to_save_figs, "Test_AUROC_Female.pdf"), width = 3, height = 3)
saveRDS(test.female.auroc, file = paste0(path_to_save_rds, "Test_AUROC_Female.RDS"))

# Precision Recall:
test.female.precision.recall <- yardstick::pr_curve(data = lasso.preds.te, truth, MIN) %>% 
  ggplot2::autoplot()
test.female.precision.recall <- test.female.precision.recall +
  ggplot2::ggtitle(("Precision-Recall (Test)"), "Female") +
  ggplot2::xlab("Recall") +
  ggplot2::ylab("Precision")
#ggsave(plot = test.female.precision.recall, filename = paste0(path_to_save_figs, "Test_PR_Female.pdf"), width = 3, height = 3)
saveRDS(test.female.precision.recall, file = paste0(path_to_save_rds, "Test_PR_Female.RDS"))


# RMSE to predict SPEAR Factor 1 score -------

# Model function
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


# Load clinical data
data_env_train <- readRDS("/scratch/data-fullcohort/batch-effects/impacc_convalescent_data_env_train_Feb_11_2024.RDS")
clinical_data_train <- data_env_train$clinical_data

# Load analyte data
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

top_feat_data <- list()

for (assay in names(feature_data_train_imputed)) {
  print(assay)
  
  top_feature_data <- list()
  
  top_feat_data_assay <- feature_data_train_imputed[[assay]]
  
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

all_feat_data_df <- bind_cols(top_feat_data)

# Load SPEAR object
spear_object <- readRDS("/scratch/data-fullcohort/spear_model/results/020624/032124_spearphysical_finalscores_allsamples.rds")

spear_convalescent_train <- as.data.frame(rbind(spear_object$Train, spear_object$TrainNoPhysical_Overlap, spear_object$TrainNoPhysical_NoOverlap))

spear_convalescent_train_event_labels <- clinical_data_train %>% 
  select(event_id, pro_2_groups) %>%
  distinct() %>%
  filter(event_id %in% rownames(spear_convalescent_train)) %>%
  mutate(label_encoded = if_else(pro_2_groups == "MIN", 0 , 1)) %>%
  column_to_rownames("event_id")

spear_convalescent_train_encoded <- spear_convalescent_train_event_labels[rownames(spear_convalescent_train),c("label_encoded"), drop=F]
spear_convalescent_train_full <- spear_convalescent_train_event_labels[rownames(spear_convalescent_train),c("pro_2_groups"), drop=F]

set.seed(1)
spear_convalescent_train_fold <- createFolds(as.factor(spear_convalescent_train_full$pro_2_groups), k = 10, list = F, returnTrain = FALSE)
spear_convalescent_train_fold_eventid <- data.frame(fold = spear_convalescent_train_fold)
rownames(spear_convalescent_train_fold_eventid) <- rownames(spear_convalescent_train_full)



## RMSE spear conv per participant --------

SO_mod <- cv_scores_lasso_model_regression(X = as.matrix(all_feat_data_df[rownames(spear_convalescent_train),] %>% select(starts_with("SO_"))),
                                           Y = spear_convalescent_train[,1],
                                           model_name = "SO",
                                           n_folds = 10,
                                           fold_ids = spear_convalescent_train_fold_eventid$fold
)

PPG_mod <- cv_scores_lasso_model_regression(X = as.matrix(all_feat_data_df[rownames(spear_convalescent_train),] %>% select(starts_with("PPG_"))),
                                            Y = spear_convalescent_train[,1],
                                            model_name = "PPG",
                                            n_folds = 10,
                                            fold_ids = spear_convalescent_train_fold_eventid$fold
)

PGX_mod <- cv_scores_lasso_model_regression(X = as.matrix(all_feat_data_df[rownames(spear_convalescent_train),] %>% select(starts_with("PGX_"))),
                                            Y = spear_convalescent_train[,1],
                                            model_name = "PGX",
                                            n_folds = 10,
                                            fold_ids = spear_convalescent_train_fold_eventid$fold
)


PPT_mod <- cv_scores_lasso_model_regression(X = as.matrix(all_feat_data_df[rownames(spear_convalescent_train),] %>% select(starts_with("PPT_"))),
                                            Y = spear_convalescent_train[,1],
                                            model_name = "PPT",
                                            n_folds = 10,
                                            fold_ids = spear_convalescent_train_fold_eventid$fold
)

PMG_mod <- cv_scores_lasso_model_regression(X = as.matrix(all_feat_data_df[rownames(spear_convalescent_train),] %>% select(starts_with("PMG_"))),
                                            Y = spear_convalescent_train[,1],
                                            model_name = "PMG",
                                            n_folds = 10,
                                            fold_ids = spear_convalescent_train_fold_eventid$fold
)

BCT_mod <- cv_scores_lasso_model_regression(X = as.matrix(all_feat_data_df[rownames(spear_convalescent_train),] %>% select(starts_with("BCT_"))),
                                            Y = spear_convalescent_train[,1],
                                            model_name = "BCT",
                                            n_folds = 10,
                                            fold_ids = spear_convalescent_train_fold_eventid$fold
)

ALL_mod <- cv_scores_lasso_model_regression(X = as.matrix(all_feat_data_df[rownames(spear_convalescent_train),]),
                                            Y = spear_convalescent_train[,1],
                                            model_name = "ALL",
                                            n_folds = 10,
                                            fold_ids = spear_convalescent_train_fold_eventid$fold
)


# Plot CV AUROC per fold
all.res.df <- bind_rows(SO_mod$resdf, 
                        PPG_mod$resdf,
                        PPT_mod$resdf,
                        PMG_mod$resdf,
                        PGX_mod$resdf,
                        BCT_mod$resdf,
                        ALL_mod$resdf
)


rmse.plot <- ggplot(all.res.df) +
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey20") +
  ggbeeswarm::geom_quasirandom(aes(x = Model, y = CV, fill = Model), shape = 21, width = .15) +
  geom_boxplot(aes(x = Model, y = CV, fill = Model), outlier.shape = NA, alpha = .6, color = "black") +
  theme_bw() +
  #scale_fill_manual(values = model.cols) +
  labs(x = "", y = "Model Error (RMSE)", title = "Error (RMSE) of models trained on each analyte set to reconstruct SPEAR Factor 1 scores") +
  #ylim(c(0.5, NA)) +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5))

rmse.plot

saveRDS(rmse.plot, file = paste0(path_to_save_rds, "RMSE_TOTAL_AUROC.RDS"))

