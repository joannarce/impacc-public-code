# Here I perform LASSO regression for nasal + CT counts, filtering for DEG genes 
# with a logFC threshold as inputs

## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("janitor")
library("readr")
library("qs")
library("ggrepel")

## Load packages for analysis ======

library("DESeq2")
library("glmnet")
library("edgeR")
library("caret")
library(pROC)


## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")


here()
output_dir <- here("05d_nasal_complete_LASSO_deg_filter")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("04_assemble_bucket_folds")



## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## READ IN FILES ======
cv.folds.nasal <- read_csv(here(input_dir, "nasal_complete_train_folds.csv"))
nasal_counts <- read_csv("../../../data/nasal-transcriptomics/legacy/2022-10-19/nasal-transcriptomics-Counts.csv")
pbmc_counts_rowFeature <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-RowFeature.csv", skip = 1)
rpM_ct_complete <- read_csv("relationship_CT_rpM/rpM_ct_complete.csv")
filter_thresh_list <- qread("04a_nasal_ct_train_extract_deg_thresh/nasal_gene_list_threshold.qs")


## PREPARE DATA ====

protein_coding <- pbmc_counts_rowFeature %>%
  dplyr::filter(gene_biotype == "protein_coding")

# Select protein coding genes only

protein_coding_list <- protein_coding$gene_id
ensembl_to_gene <- setNames(protein_coding$gene_name, protein_coding$gene_id)

set.seed(555)


# Add ct value to metadata
cv.folds.nasal <- rpM_ct_complete %>% 
  select("participant_id", "ct") %>%
  inner_join(., cv.folds.nasal, by = "participant_id")



# INPUTS:
# mod, the output of glmnet, an object with S3 class "glmnet"
# selected_lambda, the selected lambda to obtain the right number of nonzero coefficients

# OUTPUT:
# A dataframe containing the nonzero coefficients, their coefficient, and their fold
# The number of rows is dependent on the penalty
lasso_coef_df <- function(mod, selected_lambda) {
  # Extract the coefficient matrix
  coef_matrix <- coef(mod, s=selected_lambda)
  
  # Filter to non-zero coefficients and remove the intercept
  non_zero_genes <- rownames(coef_matrix)[coef_matrix[, 1] != 0 & rownames(coef_matrix) != "(Intercept)"]
  non_zero_coefs <- coef_matrix[non_zero_genes, 1]
  
  # Create a dataframe with the non-zero coefficients
  data.frame(gene = non_zero_genes, coef = non_zero_coefs)
}

# For storing ROC results
fold.roc <- list()

set.seed(555)
# INPUTS:
# Counts, which is a dataframe with sample_ids as rows and genes/features as observations, or columns
# Meta_folds, which is a dataframe with patient_ids as rows, and metadata variables as columns
# An additional variable in meta_folds is the fold number each patient is a part of
# Fold, which is an integer between 1 and 5
# numFeatures, which is an integer representing the desired number of features
# with non-zero coefficients per L1 normalization

prepare_counts_folds <- function(counts, meta_folds, fold, numFeatures) {
  
  
  
  print(paste("Starting fold...", fold))
  
  # Make sample_id a row
  # Transpose the counts dataframe
  
  counts <- counts %>%
    column_to_rownames(var = "sample_id") %>%
    t(.)
  
  # Only keep protein-coding genes
  counts <- counts[rownames(counts) %in% protein_coding_list, ]
  
  
  # The counts dataframe is now a dataframe with features as rows and sample_ids as columns
  # This is the format DeSeqDataSetFromMatrix is able to read
  
  # Split the metadata into training and test based on the fold
  # 1/5 of the data will be test, 4/5 of the data will be train
  test.fold <- meta_folds[meta_folds$fold == fold,]
  train.folds <- meta_folds[meta_folds$fold != fold,]
  
  
  # Now, set up training data
  
  # Only keep the samples (columns) that are used for training
  counts.train <- counts[,train.folds$sample_id]
  
  # Creates a boolean of the genes/features that pass the QC threshold
  keep.train <- rowSums(counts.train >= 10) >= (0.2*ncol(counts.train))
  
  # Set up the test data
  counts.test <- counts[,test.fold$sample_id]
  
  
  # Create a DESeqDataset from our counts matrix
  dds.train <- DESeqDataSetFromMatrix(
    countData = counts.train[keep.train,], # Only keep the genes/features
    # that pass the QC threshold
    colData = train.folds,
    design = ~1) # We are using DESeq not to gauge differential expression but to 
  # # create a matrix to normalize, so we do not need a formula
  # 
  # # Normalize training data
  dds.train <- estimateSizeFactors(dds.train)
  dds.train <- estimateDispersions(dds.train)
  vsd.train <- varianceStabilizingTransformation(dds.train) %>% 
    assay %>% 
    round(., digits=2)
  
  
  # Now, work with test data
  
  dds.test <- DESeqDataSetFromMatrix(
    countData = counts.test[keep.train,],
    colData = test.fold,
    design = ~1)
  dds.test <- estimateSizeFactors(dds.test)
  dispersionFunction(dds.test) <- dispersionFunction(dds.train)
  vsd.test <- varianceStabilizingTransformation(dds.test, blind=FALSE) %>% 
    assay %>% 
    round(., digits=2)
  
  # Only keep genes that have logFC >= 1 or <= -1 and are DE
  vsd.train <- vsd.train[rownames(vsd.train) %in% filter_thresh_list, ]
  vsd.test <- vsd.test[rownames(vsd.test) %in% filter_thresh_list, ]
  
  # In order to run glmnet, the first argument must be a matrix where rows
  # are observations (samples). We set this up here for both our training
  # and testing.
  
  print("Fitting lasso logistic regression using glmnet...")
  
  x <- t(vsd.train)
  test <- t(vsd.test)
  
  ## Add additional features
  train_meta_subset <- train.folds %>%
    select("sample_id", "admit_age", "ct")
  
  test_meta_subset <- test.fold %>%
    select("sample_id", "admit_age", "ct")
  
  x <- x %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "sample_id") %>% 
    inner_join(train_meta_subset, by = "sample_id") %>% 
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.)
  
  test <- test %>% 
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "sample_id") %>% 
    inner_join(test_meta_subset, by = "sample_id") %>% 
    tibble::column_to_rownames(var = "sample_id") %>%
    as.matrix(.)
  
  ## Set up y data
  
  y_train <- train.folds$group == 1
  y_test <- test.fold$group == 1
  y_train_binary <- ifelse(y_train == TRUE, 1, 0)
  y_test_binary <- ifelse(y_test == TRUE, 1, 0)
  
  # Standardize train and test data
  
  preProcValues <- preProcess(x, method = c("center", "scale"))
  
  trainTransformed <- predict(preProcValues, x)
  testTransformed <- predict(preProcValues, test)
  
  # Set penalty factor to 0 for age and ct value
  penalty.factor <- ifelse(colnames(trainTransformed) %in% c('admit_age', 'ct'), 0, 1)
  
  # Fit lasso logistic regression on the training set ONLY
  # We only use the corresponding trajectory boolean for samples that
  # are not part of the current fold (all train folds only)
  
  mod <- glmnet::glmnet(trainTransformed, 
                        y_train_binary, 
                        family = 'binomial',
                        alpha = 1, 
                        intercept = TRUE,
                        standardize = FALSE,
                        trace.it = 1, 
                        nlambda = 1000,
                        penalty.factor = penalty.factor
  )
  
  # Determine the penalty value required to obtain the desired number
  # of features with non-zero coefficients
  print("Determining the penalty value...")
  print(mod$df)[1]
  lambda_index <- which(mod$df == 7)[1]
  
  print(mod$lambda[lambda_index])
  selected_lambda <- mod$lambda[lambda_index]
  
  pred <- predict(mod, 
                  testTransformed, 
                  type='response',
                  s=selected_lambda)[,1] %>%
    {data.frame(pred=.)} %>%
    tibble::rownames_to_column("sample_id")
  
  
  print("Calculating AUC...")
  
  
  # Build a ROC curve
  # The response is the trajectory group of individuals in the test fold
  fold.roc[[fold]] <<- roc(
    response = y_test_binary, 
    predictor = pred$pred,
    direction = "<",
    plot = TRUE)
  print(sprintf(
    "Fold %d AUC: %.3f", fold, fold.roc[[fold]]$auc))
  
  
  
  # Get the nonzero coefficients and store which split/fold they were from
  print("Extracting non-zero coefficients...")
  lasso_coef_df(mod, selected_lambda) %>%
    dplyr::mutate(fold=fold) ->
    nonzero_coefs
  
  print("Returning results...")
  return(list(mod=mod, coef=nonzero_coefs))
}

# Function to run prepare_counts_folds for multiple folds
# Returns model, nonzero coefficients, and roc curve stats
run_all_folds <- function(counts, meta_folds, numFeatures) {
  results <- list()
  
  for (fold in 1:5) { # 5 fold cv
    print(paste("Processing fold...", fold))
    result <- prepare_counts_folds(counts, meta_folds, fold, numFeatures)
    results[[fold]] <- result
  }
  
  # Get the list of AUCs
  fold.auc <- unlist(lapply(
    fold.roc,
    FUN=function(x) x$auc))
  
  return(list(results = results, fold.auc = fold.auc))
}

## Collect Data =====
n_size <- c(2, 3, 4, 5, 6) # Various gene sizes
for (n in n_size)
{
  print(paste("Processing gene size...", n))
  feature_length =  n + 2
  print(paste("The number of features are: ", feature_length))
  ct_age_fold_results_nasal <- run_all_folds(nasal_counts, cv.folds.nasal, feature_length)
  print(paste("Done with gene size...", n))
  print("Moving on to post processing!")
  
  ## GENE NAME EXTRACTION ====
  
  # Initialize a list to store gene names for each fold
  gene_names_per_fold <- list()
  
  for (i in seq_along(ct_age_fold_results_nasal$results)) {
    # Extract the coef data frame 
    coef_df <- ct_age_fold_results_nasal$results[[i]]$coef
    
    # Extract gene names
    gene_names <- coef_df$gene
    
    # Store the gene names
    gene_names_per_fold[[i]] <- gene_names
  }
  nasal_output <- data.frame(Gene_name = character(), fold = numeric(), stringsAsFactors = FALSE)
  
  # Iterate through each fold and determine the genes/ensembl ids in each fold
  
  for (i in seq_along(gene_names_per_fold))
  {
    ensembl_ids <- gene_names_per_fold[[i]]
    gene_names <- sapply(ensembl_ids, function(id) {
      if (id %in% names(ensembl_to_gene)) {
        return(ensembl_to_gene[[id]])
      } else {
        return(NA)  # Return NA if ID is not found
      }
    })
    
    fold_df <- data.frame(Gene_name = gene_names, fold = i, stringsAsFactors = FALSE)
    valid_indices <- !is.na(gene_names)
    gene_names <- gene_names[valid_indices]
    nasal_output <- rbind(nasal_output, fold_df)
  }
  
  ## EXPORT ====
  print("Exporting...")
  gene_names_df <- do.call(rbind, lapply(gene_names_per_fold, function(x) data.frame(Gene = x)))
  write_csv(gene_names_df, here(output_dir, paste(n, "nasal_gene_names_per_fold.csv", sep = "_")))
  write_csv(nasal_output, here(output_dir, paste(n, "nasal_geneid_output.csv", sep = "_")))
  qsave(ct_age_fold_results_nasal, here(output_dir, paste(n, "nasal_output.qs", sep = "_")))
  qsave(gene_names_per_fold, here(output_dir, paste(n, "nasal_feature_list.qs", sep = "_")))
  print("Done with post-processing!")
}












