library(MultiAssayExperiment)
#remotes::install_bitbucket("kleinstein/SPEAR@main")
library(SPEAR)

# NOTE: Run after: '1_spear_physical_model_training.R'
# Path to where SPEAR object is saved results to:
path_to_results <- "/scratch/data-fullcohort/spear_model/results/020624/"
# Name of saved SPEAR object in 'path_to_results'
spear_obj_name <- '021924_physical_gaussian_SPEAR.rds'
# Name of file to save factor scores to:
fs_file_name <- "032124_spearphysical_finalscores_allsamples.rds"

# Data:
# Path to PROMIS Scores:
path_to_promis_scores <- "/scratch/multiomic_factors_PRO_scores_testing/20240229_PRO_Scores_Cleaned.csv"
model_to_predict = "physical"
# Paths to data used to train model:
path_to_tr_data <- "/scratch/data-fullcohort/spear/imputed_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212/"
data_file_tr_intro <- "spear_data_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_"
metadata_tr_file <- "spear_metadata_MOFA_train_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.csv"
path_to_te_data <- "/scratch/data-fullcohort/spear/imputed_data_MOFA_test_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212/"
data_file_te_intro <- "spear_data_MOFA_test_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212_"
metadata_te_file <- "spear_metadata_MOFA_test_set_250_factors_SO_PPT_PPG_PMG_PGX_BCT_20240212.csv"
path_to_acute_data <- "/scratch/data-fullcohort/spear/acute_preprocessed_data_20240323/acute_imputed_data.rds"

# Read Training Data:
data.BCT <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "BCT.csv"))
data.PGX <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "PGX.csv"))
data.PPG <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "PPG.csv"))
data.PMG <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "PMG.csv"))
data.PPT <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "PPT.csv"))
data.SO <- read.csv(paste0(path_to_tr_data, data_file_tr_intro, "SO.csv"))
data.meta <- read.csv(paste0(path_to_tr_data, metadata_tr_file))

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
data.meta <- read.csv(paste0(path_to_te_data, metadata_te_file))

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

# Load SPEAR object:
SPEARobj <- readRDS(paste0(path_to_results, spear_obj_name))

# Update NA's in test data with 0 for PPG:
all.data.te$PPG[is.na(all.data.te$PPG)] <- 0
tmp <- SPEARobj$data$test@ExperimentList$PPG
tmp[is.na(tmp)] <- 0
SPEARobj$data$test@ExperimentList$PPG <- tmp

# Set weights to 'min' (best prediction model):
SPEARobj$set.weights(w.x = 0)


# Get factor scores for the following:
fs.list <- list()
# Train (all event_ids used to train the SPEAR model, each with a PROMIS Physical score (not in the object))
# CV (all event_ids used to train the SPEAR model, each with a PROMIS Physical score, using CV score)
# Test (all event_ids with a PROMIS Physical score NOT included in the SPEAR model training)
# "TrainNoPhysical" (open to other names!) - those who would pertain to the training set (80%) but do not have a PROMIS physical score
# TestNoPhysical (same as above, but pertaining to the 20%)
# AcuteTrain (all event_ids corresponding to participants in the train cohort, using acute molecular profiles to get factor score, no physical score obviously)
# AcuteTest (all event_ids corresponding to participants in the test cohort using acute molecular profiles to construct factor score)

# Train:
fs.list[['Train']] <- SPEARobj$get.factor.scores(data = "train", cv = FALSE)

# CV:
fs.list[['CV']] <- SPEARobj$get.factor.scores(data = "train", cv = TRUE)

# Test: (with promis)
fs.te <- SPEARobj$get.factor.scores(data = "test")
test.with.promis <- !is.na(SPEARobj$data$test$response)
fs.list[['Test']] <- fs.te[test.with.promis,]

# Test: (without physical):
fs.list[['TestNoPhysical']] <- fs.te[which(!test.with.promis),]




# Train: (without physical):
# Make SPEAR object:
train.with.no.promis <- is.na(promis.df.tr$PROMIS_Physical_score)
participants.in.train <- unique(sapply(colnames(SPEARobj$data$train@ExperimentList$SO), function(e.id){
  tmp <- stringr::str_split(e.id, "-")
  return(paste0(tmp[[1]][1], "-", tmp[[1]][2]))
}))
promis.df.participants <- sapply(promis.df.tr$event_id, function(e.id){
  tmp <- stringr::str_split(e.id, "-")
  return(paste0(tmp[[1]][1], "-", tmp[[1]][2]))
})
train.id.overlap <- promis.df.participants %in% participants.in.train
train.with.promis.overlap <- train.id.overlap & train.with.no.promis
train.with.promis.nooverlap <- (!train.id.overlap) & train.with.no.promis
assays.tr.nophys.overlap <- lapply(all.data.tr, function(d){
  return(d[,train.with.promis.overlap])
})
assays.tr.nophys.nooverlap <- lapply(all.data.tr, function(d){
  return(d[,train.with.promis.nooverlap])
})
df.tr.nophys.overlap <- df.tr[train.with.promis.overlap,]
df.tr.nophys.nooverlap <- df.tr[train.with.promis.nooverlap,]
df.tr.nophys.overlap$response = NA
df.tr.nophys.nooverlap$response = NA

# MAE:
mae.tr.nophys.overlap <- MultiAssayExperiment::MultiAssayExperiment(experiments = assays.tr.nophys.overlap,
                                                                    colData = df.tr.nophys.overlap)
SPEARobj$add.data(mae.tr.nophys.overlap, name = "train_nophys_overlap")
fs.list[['TrainNoPhysical_Overlap']] <- SPEARobj$get.factor.scores(data = "train_nophys_overlap")
mae.tr.nophys.nooverlap <- MultiAssayExperiment::MultiAssayExperiment(experiments = assays.tr.nophys.nooverlap,
                                                                    colData = df.tr.nophys.nooverlap)
SPEARobj$add.data(mae.tr.nophys.nooverlap, name = "train_nophys_nooverlap")
fs.list[['TrainNoPhysical_NoOverlap']] <- SPEARobj$get.factor.scores(data = "train_nophys_nooverlap")

# Sanity check:
train.ids <- unique(sapply(rownames(fs.list$Train), function(e.id){
  tmp <- stringr::str_split(e.id, "-")
  return(paste0(tmp[[1]][1], "-", tmp[[1]][2]))
}))
train.overlap.ids <- unique(sapply(rownames(fs.list$TrainNoPhysical_Overlap), function(e.id){
  tmp <- stringr::str_split(e.id, "-")
  return(paste0(tmp[[1]][1], "-", tmp[[1]][2]))
}))
train.nooverlap.ids <- sapply(rownames(fs.list$TrainNoPhysical_NoOverlap), function(e.id){
  tmp <- stringr::str_split(e.id, "-")
  return(paste0(tmp[[1]][1], "-", tmp[[1]][2]))
})


# Acute:
acute_data <- readRDS(path_to_acute_data)
all.acute.tr = acute_data$train
all.acute.te = acute_data$test
# Replace NA's in PPG with 0
all.acute.tr$PPG[is.na(all.acute.tr$PPG)] <- 0
all.acute.te$PPG[is.na(all.acute.te$PPG)] <- 0

# Which have a PRO group?
participants.acute.tr <- rownames(all.acute.tr$SO)
participants.acute.te <- rownames(all.acute.te$SO)
participants.acute.tr <- sapply(participants.acute.tr, function(e.id){
  tmp <- stringr::str_split(e.id, "-")
  return(paste0(tmp[[1]][1], "-", tmp[[1]][2]))
})
participants.acute.te <- sapply(participants.acute.te, function(e.id){
  tmp <- stringr::str_split(e.id, "-")
  return(paste0(tmp[[1]][1], "-", tmp[[1]][2]))
})
participants.with.pro.tr <- sapply(rownames(df.tr), function(e.id){
  tmp <- stringr::str_split(e.id, "-")
  return(paste0(tmp[[1]][1], "-", tmp[[1]][2]))
})
participants.with.pro.te <- sapply(rownames(df.te), function(e.id){
  tmp <- stringr::str_split(e.id, "-")
  return(paste0(tmp[[1]][1], "-", tmp[[1]][2]))
})
participants.with.pro <- unique(c(participants.with.pro.tr, participants.with.pro.te))
train.to.keep = participants.acute.tr %in% participants.with.pro
test.to.keep = participants.acute.te %in% participants.with.pro
df.all <- rbind(df.tr, df.te)
df.all$participant_id <- c(participants.with.pro.tr, participants.with.pro.te)
pro.train = sapply(participants.acute.tr[train.to.keep], function(p.id){
  return(df.all$pro_group_labels[which(df.all$participant_id == p.id)[1]])
})
pro.test = sapply(participants.acute.te[test.to.keep], function(p.id){
  return(df.all$pro_group_labels[which(df.all$participant_id == p.id)[1]])
})

# Make data:
all.acute.withPRO <- lapply(1:length(all.acute.tr), function(d){
  d1 <- all.acute.tr[[d]]
  d2 <- all.acute.te[[d]]
  return(rbind(d1[train.to.keep,], d2[test.to.keep,]))
})
names(all.acute.withPRO) <- names(all.acute.tr)
all.acute.order <- list()
for(assay in names(SPEARobj$data$train@ExperimentList)){
  all.acute.order[[assay]] <- t(all.acute.withPRO[[assay]])
}

# MAEs:
mae.acute <- MultiAssayExperiment::MultiAssayExperiment(experiments = all.acute.order,
                                                        colData = data.frame(
                                                          PRO_group = c(pro.train, pro.test)
                                                        ))
SPEARobj$add.data(mae.acute, name = "acute")
fs.list[['Acute']] <- SPEARobj$get.factor.scores(data = "acute")


# Add analyte scores:
feature.df <- SPEARobj$get.analyte.scores(probability.cutoff = 0)
feature.df <- dplyr::filter(feature.df, Factor == "Factor1")
fs.list[['analyte.scores']] <- feature.df

# Save factor scores:
saveRDS(fs.list, paste0(path_to_results, fs_file_name))

