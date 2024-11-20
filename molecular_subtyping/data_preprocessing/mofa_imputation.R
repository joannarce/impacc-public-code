# Dependencies:
require(MOFA2)
predict = stats::predict # otherwise masked by MOFA2::predict

# Define paths:
path_to_codebases <- "/scratch/data-integration/pipelines/combined/"
data_path <- "/scratch/data-integration/pipelines/processedData/"
MOFA_imputation_file_train <- "mofa_imputed_train_jeremy_20220517.RDS"
MOFA_imputation_file_test <- "mofa_imputed_test_jeremy_20220517.RDS"

# Source integration codebase:
source(paste0(path_to_codebases, "codebase_integration.R"))

# Read in the newest data_env:
preprocessed_data <- readRDS(paste0(data_path, "preprocessed_data_ravi_20220513.RDS"))
data_env <- readRDS(paste0(data_path, "data_env_ravi_20220513.RDS"))

# Load the data objects into the environment
for( n in names(data_env) ) {
  assign( n, data_env[[n]] )
}

# Clinical data, add part.event:
a = strsplit(clinical_data$event_type," ")
a1 = sapply(a,function(z) z[[1]])
a2 = sapply(a,function(z) z[[2]])
clinical_data$part.event<-paste(clinical_data$participant_id,a1,a2,sep = ".")

# Designate omics to be used:
omics_used_DR = c(
  "ppt" = "plasma_proteomics_targeted",
  "ppgd" = "plasma_proteomics_global_dda",
  "so" = "serum_olink",
  "pmg" = "plasma_metabolomics_global",
  "nt" = "nasal_transcriptomics",
  "pbmc" = "pbmc_transcriptomics"
)


# Functions to put into codebase:

get.event_id.from.sample_id <- function(s_id, clinical_data){
  dict <- distinct(clinical_data, sample_id, event_id) %>% filter(sample_id==s_id)
  if(nrow(dict) > 0){
    return(dict$event_id[1])
  } else {
    return(NA)
  }
}

get.training_test.from.sample_id <- function(s_id, clinical_data){
  dict <- distinct(clinical_data, sample_id, training_test) %>% filter(sample_id==s_id)
  if(nrow(dict) > 0){
    return(dict$training_test[1])
  } else {
    return(NA)
  }
}

data_prepare_DR = function(omics_used, data_env, clinical_data, rna_keep = 5000){
  # assay data
  datasets.mofa = list()
  #### Identify outlying samples and filter out them.
  for( i in 1:length(omics_used) ) {
    o_name = omics_used[i]
    count_df <- data_env[[omics_used[i]]]
    datasets.mofa[[o_name]] = count_df
  }
  datasets <- datasets.mofa
  train.datasets <- list()
  test.datasets <- list()
  datasets.names <- names(datasets.mofa)
  for(i in 1:length(datasets)){
    dataset.name <- datasets.names[i]
    dataset <- datasets[[i]]
    print(dataset.name)
    ##filter and transformation
    ###filtering transcriptomics
    sds = apply(dataset, 2, sd)
    dataset = dataset[,!is.na(sds)]
    sds = apply(dataset, 2, sd)
    nas = apply(dataset!=0, 2, sum)
    if(grepl("transcriptomics",dataset.name) & !is.null(rna_keep)){
      thr_sd =sort(sds,decreasing = T)[rna_keep]
    }else{
      thr_sd = -Inf
    }
    dataset = dataset[,sds>=thr_sd & nas >= 10]

    # Change sample name to event_id
    train_test = sapply(rownames(dataset), function(z) get.training_test.from.sample_id(z, clinical_data=clinical_data))
    train.dataset = dataset[which(train_test == "training"),]
    test.dataset = dataset[which(train_test == "test"),]
    
    tmp =  sapply(rownames(train.dataset), function(z) get.event_id.from.sample_id(z, clinical_data=clinical_data))
    train.dataset <- train.dataset[!is.na(tmp),]
    rownames(train.dataset) <- tmp
    train.datasets[[i]] <- train.dataset
    
    tmp =  sapply(rownames(test.dataset), function(z) get.event_id.from.sample_id(z, clinical_data=clinical_data))
    test.dataset <- test.dataset[!is.na(tmp),]
    rownames(test.dataset) <- tmp
    test.datasets[[i]] <- test.dataset
  }
  
  # Name the training/testing
  names(train.datasets) <- names(test.datasets) <- names(datasets)
  
  return(list(train = train.datasets, test = test.datasets))
}

generate_full_matrices = function(datasets_combined){
  event_ids <- sort(unique(unlist(sapply(datasets_combined, function(i){return(rownames(i))}))))
  new_datasets_combined <- list()
  for(i in 1:length(datasets_combined)){
    print(i)
    d <- datasets_combined[[i]]
    m <- matrix(NA, nrow = length(event_ids), ncol = ncol(d))
    rownames(m) <- event_ids
    colnames(m) <- colnames(d)
    m <- as.data.frame(m)
    m[rownames(d),] <- d
    new_datasets_combined[[i]] <- m
  }
  names(new_datasets_combined) <- names(datasets_combined)
  return(new_datasets_combined)
}

mofa_run = function(datasets_combined, num_factors = 50,
                    seed = 42, maxiter = 1000, drop_factor_threshold = -1, 
                    convergence_mode = "fast", startELBO = 2, stochastic=F,
                    mofa_file = NULL){
  # Run MOFA+:
  x.mofa <- list()
  for(d in 1:length(datasets_combined)){
    x.mofa[[d]] = t(datasets_combined[[d]])
  }
  names(x.mofa) = names(datasets_combined)
  MOFAobject <- create_mofa(x.mofa)
  # MOFA+ specific parameters:
  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  train_opts <- get_default_training_options(MOFAobject)
  model_opts$num_factors =num_factors
  train_opts$convergence_mode = convergence_mode
  train_opts$startELBO = startELBO
  train_opts$maxiter = maxiter
  train_opts$stochastic = stochastic
  train_opts$seed = seed
  train_opts$drop_factor_threshold =drop_factor_threshold
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  MOFAobject.trained <- run_mofa(MOFAobject,use_basilisk = TRUE, save_data = F)
  return(MOFAobject.trained)
}


# Prepare data for MOFA imputation (will return $train and $test)
prepared_data_DR = data_prepare_DR(omics_used_DR, rna_keep = NULL, clinical_data = clinical_data, data_env = preprocessed_data)

# Make matrices full (i.e. make the event_ids consistent across all omics, and add NA if missing)
prepared_data_DR_full <- list(
  train = generate_full_matrices(prepared_data_DR$train),
  test = generate_full_matrices(prepared_data_DR$test)
)

# IMPUTING TRAINING DATASET VIA MOFA:
cur_data <- prepared_data_DR_full$train
names(cur_data) <- omics_used_DR
for(d in 1:length(cur_data)){
  colnames(cur_data[[d]]) <- paste0(names(cur_data)[d], "_", colnames(cur_data[[d]]))
}
# Train MOFA model:
MOFAobject.trained = mofa_run(datasets_combined = cur_data, num_factors = 200,
                              seed = 2021, maxiter = 1000, drop_factor_threshold = 1e-4, convergence_mode = "fast", stochastic = F,  startELBO = 1)
# Get predicted values from MOFA model:
imputedMOFA <- MOFA2::predict(MOFAobject.trained, factors = "all", add_intercept = T)
# Replace any NA in data with imputedMOFA:
tmp = prepared_data_DR_full$train
names(tmp) <- omics_used_DR
for(i in 1:length(tmp)){
  non_observed <- is.na(tmp[[i]])
  tmp[[i]][non_observed] = t(imputedMOFA[[i]][[1]])[non_observed]
}
# Generate return object and save it:
mofa_example = list()
mofa_example[["datasets"]] =  prepared_data_DR_full$train
mofa_example[["datasets_imputed"]] = tmp
mofa_example[["MOFAobject.trained"]] = MOFAobject.trained 
saveRDS(mofa_example, file = paste0(data_path, MOFA_imputation_file_train))






# IMPUTING TEST DATASET VIA MOFA:
cur_data <- prepared_data_DR_full$test
names(cur_data) <- omics_used_DR
for(d in 1:length(cur_data)){
  colnames(cur_data[[d]]) <- paste0(names(cur_data)[d], "_", colnames(cur_data[[d]]))
}
# Train MOFA model:
MOFAobject.trained = mofa_run(datasets_combined = cur_data, num_factors = 200,
                              seed = 2021, maxiter = 1000, drop_factor_threshold = 1e-4, convergence_mode = "fast", stochastic = F,  startELBO = 1)
# Get predicted values from MOFA model:
imputedMOFA <- MOFA2::predict(MOFAobject.trained, factors = "all", add_intercept = T)
# Replace any NA in data with imputedMOFA:
tmp = prepared_data_DR_full$test
names(tmp) <- omics_used_DR
for(i in 1:length(tmp)){
  non_observed <- is.na(tmp[[i]])
  tmp[[i]][non_observed] = t(imputedMOFA[[i]][[1]])[non_observed]
}
# Generate return object and save it:
mofa_example = list()
mofa_example[["datasets"]] =  prepared_data_DR_full$test
mofa_example[["datasets_imputed"]] = tmp
mofa_example[["MOFAobject.trained"]] = MOFAobject.trained 
saveRDS(mofa_example, file = paste0(data_path, MOFA_imputation_file_test))