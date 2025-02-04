#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("MultiAssayExperiment")
library(MultiAssayExperiment)
#remotes::install_bitbucket("kleinstein/SPEAR@main")
library(SPEAR)

# NOTE: Run after: '1_spear_physical_model_training.R'
# Path to where SPEAR object is saved results to:
path_to_results <- "/scratch/data-fullcohort/spear/scores_data_acutefullcohort/"
dir.create(path_to_results)
path_to_model <- "/scratch/data-fullcohort/spear_model/results/020624/"
# Name of saved SPEAR object in 'path_to_results'
spear_obj_name <- '021924_physical_gaussian_SPEAR.rds'
# Name of file to save factor scores to:
fs_file_name <- "082024_spearphysical_finalscores_acute_fullcohort_allsamples.rds"

# Data:
# Path to PROMIS Scores:
path_to_promis_scores <- "/scratch/multiomic_factors_PRO_scores_testing/20240229_PRO_Scores_Cleaned.csv"
model_to_predict = "physical"

# Paths to data:
path_to_acute_data <- "/scratch/data-fullcohort/spear/acute_fullcohort_preprocessed_data_20240812/acute_allcohort_imputed_data.rds"


# Load SPEAR object:
SPEARobj <- readRDS(paste0(path_to_model, spear_obj_name))

# Set weights to 'min' (best prediction model):
SPEARobj$set.weights(w.x = 0)


# Get factor scores for the following:
fs.list <- list()


# Acute:
acute_data <- readRDS(path_to_acute_data)
all.acute.tr = acute_data$train
all.acute.te = acute_data$test
# Replace NA's in PPG with 0
all.acute.tr$PPG[is.na(all.acute.tr$PPG)] <- 0
all.acute.te$PPG[is.na(all.acute.te$PPG)] <- 0

# Make data:
all.acute <- lapply(1:length(all.acute.tr), function(d){
  d1 <- all.acute.tr[[d]]
  d2 <- all.acute.te[[d]]
  return(rbind(d1, d2))
})

names(all.acute) <- names(all.acute.tr)
all.acute.order <- list()
for(assay in names(SPEARobj$data$train@ExperimentList)){
  all.acute.order[[assay]] <- t(all.acute[[assay]])
}

# MAEs:
mae.acute <- MultiAssayExperiment::MultiAssayExperiment(experiments = all.acute.order)
SPEARobj$add.data(mae.acute, name = "acute")
fs.list[['Acute_fullcohort']] <- SPEARobj$get.factor.scores(data = "acute")



# Add analyte scores:
feature.df <- SPEARobj$get.analyte.scores(probability.cutoff = 0)
feature.df <- dplyr::filter(feature.df, Factor == "Factor1")
fs.list[['analyte.scores']] <- feature.df

# Save factor scores:
saveRDS(fs.list, paste0(path_to_results, fs_file_name))



