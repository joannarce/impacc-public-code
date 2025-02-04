library(MultiAssayExperiment)
library(SPEAR)

# Path to save results to:
path_to_results <- ...
# Date to append to front of saved SPEAR object (e.g. '021924_physical_gaussian_SPEAR.rds')
date_str <- "031324"
# Path to PROMIS Scores:
path_to_promis_scores <- "/scratch/multiomic_factors_PRO_scores_testing/20240229_PRO_Scores_Cleaned.csv"
model_to_predict = "physical"
# Paths to data used to train model:
path_to_tr_data <- "/scratch/data-fullcohort/spear/imputed_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212/"
data_file_tr_intro <- "spear_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_"
metadata_tr_file <- "spear_metadata_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.csv"
data_file_tr_intro <- "spear_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_"
metadata_tr_file <- "spear_metadata_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.csv"
path_to_te_data <- "/scratch/data-fullcohort/spear/imputed_data_MOFA_test_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212/"
data_file_te_intro <- "spear_data_MOFA_test_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_"
metadata_te_file <- "spear_metadata_MOFA_test_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.csv"
data_file_te_intro <- "spear_data_MOFA_test_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_"
metadata_te_file <- "spear_metadata_MOFA_test_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.csv"

# Read Training Data:
data.BCT <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "BCT.csv"))
data.PGX <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "PGX.csv"))
data.PPG <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "PPG.csv"))
data.PMG <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "PMG.csv"))
data.PPT <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "PPT.csv"))
data.SO <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "SO.csv"))
data.meta <- read.csv(path_to_tr_data, metadata_tr_file)

# Preprocessing:
all.data <- list()
all.data[["BCT"]] <- t(as.matrix(data.BCT[2:ncol(data.BCT)]))
all.data[["PGX"]] <- t(as.matrix(data.PGX[2:ncol(data.PGX)]))
all.data[["PPG"]] <- t(as.matrix(data.PPG[2:ncol(data.PPG)]))
all.data[["PMG"]] <- t(as.matrix(data.PMG[2:ncol(data.PMG)]))
all.data[["PPT"]] <- t(as.matrix(data.PPT[2:ncol(data.PPT)]))
all.data[["SO"]] <- t(as.matrix(data.SO[2:ncol(data.SO)]))
rownames(data.meta) <- data.meta$event_id
data.meta$event_id <- NULL
all.data$BCT <- t(apply(all.data$BCT, 1, scale))
all.data$PGX <- t(apply(all.data$PGX, 1, scale))
all.data$PPG <- t(apply(all.data$PPG, 1, scale))
all.data$PMG <- t(apply(all.data$PMG, 1, scale))
all.data$PPT <- t(apply(all.data$PPT, 1, scale))
all.data$SO <- t(apply(all.data$SO, 1, scale))
colnames(all.data[["BCT"]]) <- data.BCT$Features
colnames(all.data[["PGX"]]) <- data.PGX$Features
colnames(all.data[["PPG"]]) <- data.PPG$Features
colnames(all.data[["PMG"]]) <- data.PMG$Features
colnames(all.data[["PPT"]]) <- data.PPT$Features
colnames(all.data[["SO"]]) <- data.SO$Features

# Save train separately:
all.data.tr <- all.data
df.tr <- data.meta

# Test load:
data.BCT <- read.csv(paste0(path_to_te_data, data_file_te_intro, "BCT.csv"))
data.PGX <- read.csv(paste0(path_to_te_data, data_file_te_intro, "PGX.csv"))
data.PPG <- read.csv(paste0(path_to_te_data, data_file_te_intro, "PPG.csv"))
data.PMG <- read.csv(paste0(path_to_te_data, data_file_te_intro, "PMG.csv"))
data.PPT <- read.csv(paste0(path_to_te_data, data_file_te_intro, "PPT.csv"))
data.SO <- read.csv(paste0(path_to_te_data, data_file_te_intro, "SO.csv"))
data.meta <- read.csv(path_to_te_data, metadata_te_file)

# Preprocessing:
all.data <- list()
all.data[["BCT"]] <- t(as.matrix(data.BCT[2:ncol(data.BCT)]))
all.data[["PGX"]] <- t(as.matrix(data.PGX[2:ncol(data.PGX)]))
all.data[["PPG"]] <- t(as.matrix(data.PPG[2:ncol(data.PPG)]))
all.data[["PMG"]] <- t(as.matrix(data.PMG[2:ncol(data.PMG)]))
all.data[["PPT"]] <- t(as.matrix(data.PPT[2:ncol(data.PPT)]))
all.data[["SO"]] <- t(as.matrix(data.SO[2:ncol(data.SO)]))
rownames(data.meta) <- data.meta$event_id
data.meta$event_id <- NULL
all.data$BCT <- t(apply(all.data$BCT, 1, scale))
all.data$PGX <- t(apply(all.data$PGX, 1, scale))
all.data$PPG <- t(apply(all.data$PPG, 1, scale))
all.data$PMG <- t(apply(all.data$PMG, 1, scale))
all.data$PPT <- t(apply(all.data$PPT, 1, scale))
all.data$SO <- t(apply(all.data$SO, 1, scale))
colnames(all.data[["BCT"]]) <- data.BCT$Features
colnames(all.data[["PGX"]]) <- data.PGX$Features
colnames(all.data[["PPG"]]) <- data.PPG$Features
colnames(all.data[["PMG"]]) <- data.PMG$Features
colnames(all.data[["PPT"]]) <- data.PPT$Features
colnames(all.data[["SO"]]) <- data.SO$Features

# Save test separately:
all.data.te <- all.data
df.te <- data.meta

# Promis:
promis <- read.csv(path_to_promis_scores)

promis.df.tr <- data.frame()
for(e.id in rownames(df.tr)){
  promis.row <- promis[promis$event_id == e.id,]
  if(nrow(promis.row) < 1){
    promis.df.tr <- rbind(promis.df.tr, promis.df.tr[nrow(promis.df.tr)-1,])
    promis.df.tr[nrow(promis.df.tr),] <- NA
    promis.df.tr[nrow(promis.df.tr),]$event_id <- e.id
  } else {
    promis.df.tr <- rbind(promis.df.tr, promis.row)
  }
}
promis.df.te <- data.frame()
for(e.id in rownames(df.te)){
  promis.row <- promis[promis$event_id == e.id,]
  if(nrow(promis.row) < 1){
    promis.df.te <- rbind(promis.df.te, promis.df.te[nrow(promis.df.te)-1,])
    promis.df.te[nrow(promis.df.te),] <- NA
    promis.df.te[nrow(promis.df.te),]$event_id <- e.id
  } else {
    promis.df.te <- rbind(promis.df.te, promis.row)
  }
}

# Make SPEAR object:
assays.tr <- all.data.tr
assays.te <- all.data.te
if(model_to_predict == "impact"){
  df.tr$response = as.vector(scale(promis.df.tr$PROMIS_Impact_score))
  df.te$response = as.vector(scale(promis.df.te$PROMIS_Impact_score))
} else if(model_to_predict == "dyspnea"){
  df.tr$response = as.vector(scale(promis.df.tr$PROMIS_Dyspnea_score))
  df.te$response = as.vector(scale(promis.df.te$PROMIS_Dyspnea_score))
} else if(model_to_predict == "mental"){
  df.tr$response = as.vector(scale(promis.df.tr$PROMIS_Mental_score))
  df.te$response = as.vector(scale(promis.df.te$PROMIS_Mental_score))
} else if(model_to_predict == "physical"){
  df.tr$response = as.vector(scale(promis.df.tr$PROMIS_Physical_score))
  df.te$response = as.vector(scale(promis.df.te$PROMIS_Physical_score))
} else if(model_to_predict == "psycho"){
  df.tr$response = as.vector(scale(promis.df.tr$PROMIS_Psycho_score))
  df.te$response = as.vector(scale(promis.df.te$PROMIS_Psycho_score))
} else if(model_to_predict == "combo"){
  df.tr$score.impact = as.vector(scale(promis.df.tr$PROMIS_Impact_score))
  df.tr$score.dyspnea = as.vector(scale(promis.df.tr$PROMIS_Dyspnea_score))
  df.tr$score.mental = as.vector(scale(promis.df.tr$PROMIS_Mental_score))
  df.tr$score.physical = as.vector(scale(promis.df.tr$PROMIS_Physical_score))
  df.tr$score.psycho = as.vector(scale(promis.df.tr$PROMIS_Psycho_score))
  df.te$score.impact = as.vector(scale(promis.df.te$PROMIS_Impact_score))
  df.te$score.dyspnea = as.vector(scale(promis.df.te$PROMIS_Dyspnea_score))
  df.te$score.mental = as.vector(scale(promis.df.te$PROMIS_Mental_score))
  df.te$score.physical = as.vector(scale(promis.df.te$PROMIS_Physical_score))
  df.te$score.psycho = as.vector(scale(promis.df.te$PROMIS_Psycho_score))
}

# Filter for missing train scores:
if(model_to_predict != "combo"){
  not.missing.scores = which(!is.na(df.tr$response))
  assays.tr <- lapply(assays.tr, function(d){
    return(d[,not.missing.scores])
  })
  df.tr <- df.tr[not.missing.scores,]
} else {
  # can't have any missing values in Y...
  not.missing.scores = which(!apply(df.tr, 1, function(row){return(any(is.na(row)))}))
  assays.tr <- lapply(assays.tr, function(d){
    return(d[,not.missing.scores])
  })
  df.tr <- df.tr[not.missing.scores,]
}

# MAE:
mae.tr <- MultiAssayExperiment::MultiAssayExperiment(experiments = assays.tr,
                                                     colData = df.tr)
mae.te <- MultiAssayExperiment::MultiAssayExperiment(experiments = assays.te,
                                                     colData = df.te)

if(model_to_predict != "combo"){
  SPEARobj.gaussian <- SPEAR::new.spear(data = mae.tr,
                                        response = "response",
                                        print.out = 5,
                                        num.folds = 10,
                                        num.factors = 80, # 71, determined from "estimate.num.factors()",  so rounding up
                                        weights.x = c(0, 0.1, 0.5, 1, 2)
  )
  # Add test data:
  SPEARobj.gaussian$add.data(data = mae.te, response = "response", name = "test")
  
} else {
  SPEARobj.gaussian <- SPEAR::new.spear(data = mae.tr,
                                        response = c("score.impact",
                                                     "score.dyspnea",
                                                     "score.mental",
                                                     "score.physical",
                                                     "score.psycho"),
                                        print.out = 5,
                                        num.folds = 10,
                                        num.factors = 80, # 71, determined from "estimate.num.factors()",  so rounding up
                                        weights.x = c(0, 0.1, 0.5, 1, 2)
  )
  # Add test data:
  SPEARobj.gaussian$add.data(data = mae.te, response = c("score.impact",
                                                         "score.dyspnea",
                                                         "score.mental",
                                                         "score.physical",
                                                         "score.psycho"), name = "test")
  
}

# Generate fold ids that balance participant_id across folds:
set.seed(42)
ALL_train_participants_split <- stringr::str_split(rownames(df.tr), pattern = "-")
ALL_train_participants <- sapply(1:length(ALL_train_participants_split), function(i){
  return(paste0(ALL_train_participants_split[[i]][1], "-", ALL_train_participants_split[[i]][2]))
})
df.tmp <- data.frame(
  participant_id = ALL_train_participants,
  pro_group = df.tr$pro_group_labels
)
# Make unique: 306 participants:
df.tmp <- dplyr::distinct(df.tmp)
nfolds = 10 # for 90% train, 10% test
N <- length(df.tmp$pro_group)
df.tmp$foldids = sample(x = rep(1:nfolds, ceiling(N/nfolds)), size = N, replace = F)
# Save for later:
#saveRDS(df.tmp, file = paste0(path_to_results, "021624_10cv_gaussian_foldids.rds"))
# Load fold.ids:
#df.tmp <- readRDS(paste0(path_to_results, "021624_10cv_gaussian_foldids.rds"))
all.fold.ids <- sapply(ALL_train_participants, function(p_id){
  return(df.tmp$foldids[which(df.tmp$participant_id == p_id)])
})
SPEARobj.gaussian$data$train$fold.ids = all.fold.ids

# Train model
SPEARobj.gaussian$train.spear(fold.ids = all.fold.ids)

# Save Object:
saveRDS(SPEARobj.gaussian, file = paste0(path_to_results, paste0(date_str, "_", model_to_predict, "_gaussian_SPEAR.rds")))
