# Here I perform glm for pbmc counts to evaluate AUC for gene sets on train data for filtered list

## No ct, for supp analysis

## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("janitor")
library("readr")
library("qs")
library("ggrepel")
library("reshape2")
library("ggpubr")

## Load packages for analysis ======

library("DESeq2")
library("glmnet")
library("edgeR")
library("caret")
library(pROC)


## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")


here()
output_dir <- here("06_train_eval_a_age")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("05c_complete_LASSO_deg_filter")



## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## READ IN FILES ======
cv.folds.pbmc <- read_csv("04_assemble_bucket_folds/pbmc_complete_train_folds.csv")
pbmc_counts <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-Counts.csv")
pbmc_counts_rowFeature <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-RowFeature.csv", skip = 1)
rpM_ct_complete <- read_csv("relationship_CT_rpM/rpM_ct_complete.csv")
two_feature_list <- qread(here(input_dir, "2_pbmc_feature_list.qs"))
three_feature_list <- qread(here(input_dir, "3_pbmc_feature_list.qs"))
four_feature_list <- qread(here(input_dir, "4_pbmc_feature_list.qs"))
five_feature_list <- qread(here(input_dir, "5_pbmc_feature_list.qs"))
six_feature_list <- qread(here(input_dir, "6_pbmc_feature_list.qs"))
eight_feature_list <- qread(here(input_dir, "8_pbmc_feature_list.qs"))
ten_feature_list <- qread(here(input_dir, "10_pbmc_feature_list.qs"))
filter_thresh_list <- qread("04a_pbmc_ct_train_extract_deg_thresh/gene_list_threshold.qs")


## PREPARE DATA ====

protein_coding <- pbmc_counts_rowFeature %>%
  dplyr::filter(gene_biotype == "protein_coding")

# Select protein coding genes only

protein_coding_list <- protein_coding$gene_id
ensembl_to_gene <- setNames(protein_coding$gene_name, protein_coding$gene_id)

set.seed(555)

# Add ct value to metadata
cv.folds.pbmc <- rpM_ct_complete %>% 
  select("participant_id", "ct") %>%
  inner_join(., cv.folds.pbmc, by = "participant_id")


# INPUTS:
# Counts, which is a dataframe with sample_ids as rows and genes/features as observations, or columns
# Meta_folds, which is a dataframe with patient_ids as rows, and metadata variables as columns
# feature_list, which is a nested list containing each set of features to be tested

# OUTPUTS:
# a matrix with:
#nrows = # of feature sets
#ncols = # of bootstraps * number of folds
#values = AUC the feature set for that iteration test
set.seed(555)
determine_3_cv <- function(counts, meta_folds, feature_list, num_folds = 3, num_bootstraps = 50) {
  
  set.seed(555)
  
  # Initialize matrix to store AUCs for each feature set across all bootstraps and folds
  auc_matrix <- matrix(0, nrow = length(feature_list), ncol = num_folds * num_bootstraps)
  # Make sample_id a row
  # Transpose the counts dataframe
  
  counts <- counts %>%
    column_to_rownames(var = "sample_id") %>%
    t(.)
  
  # The counts dataframe is now a dataframe with features as rows and sample_ids as columns
  # This is the format DeSeqDataSetFromMatrix is able to read
  
  # Only keep protein-coding genes
  counts <- counts[rownames(counts) %in% protein_coding_list, ]
  
  
  auc_index <- 1
  NA_counter <- 0
  
  for (bootstrap in seq_len(num_bootstraps)) {
    print(paste("Bootstrap iteration", bootstrap))
    
    # Create folds for cross-validation
    folds <- createFolds(meta_folds$group, k = num_folds, list = TRUE)
    
    fold_auc_results <- numeric(num_folds)
    
    
    for (i in seq_len(num_folds)) {
      
      print(paste("Starting fold...", i))
      
      # Split the metadata into training and test based on the fold
      # 1/3 of the data will be test, 2/3 of the data will be train
      test_indices <- folds[[i]]
      train_indices <- setdiff(seq_len(nrow(meta_folds)), test_indices)
      
      # Split the metadata into training and test based on the fold
      test.fold <- meta_folds[test_indices, ]
      train.folds <- meta_folds[train_indices, ]
      
      
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
      
      print("DESEQ done!")
      # In order to run glm, the first argument must be a matrix where rows
      # are observations (samples). We set this up here for both our training
      # and testing. We also add other info.
      
      
      # Only keep genes that have logFC >= 1 or <= -1 and are DE
      vsd.train <- vsd.train[rownames(vsd.train) %in% filter_thresh_list, ]
      vsd.test <- vsd.test[rownames(vsd.test) %in% filter_thresh_list, ]
      
      print("Preprocessing...")
      
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
      
      # Create df for storing results
      aucs <- numeric(length(feature_list))
      
      
      print("Preprocessing complete!")
      
      # Create matrices
      for (j in seq_along(feature_list)) {
        top_features <- feature_list[[j]]
        print(top_features)
        
        top_features <- top_features[-length(top_features)]
        print(top_features)
        
        # Fit logistic regression on the training set ONLY
        # We only use the corresponding trajectory boolean for samples that
        # are not part of the current fold (all train folds only)
        # Check if top_features indices are within bounds
        if (all(top_features %in% colnames(trainTransformed))) {
          print("Features found, subsetting matrices...")
          train_matrix <- as.matrix(trainTransformed[, top_features, drop = FALSE])
          test_matrix <- as.matrix(testTransformed[, top_features, drop = FALSE])
          
          print("Fitting logistic regression using glm...")
          
          glm_model <- glm(y_train_binary ~ ., 
                           data = as.data.frame(train_matrix), 
                           family = binomial())
          
          pred_logistic <- predict(glm_model, 
                                   newdata = as.data.frame(test_matrix), 
                                   type = "response")
          
          roc.car <- roc(
            response = y_test_binary, 
            predictor = pred_logistic,
            direction = "<",
            plot = TRUE)
          
          auc_val <- as.numeric(auc(roc.car))
          print(auc_val)
          auc_matrix[j, auc_index] <- auc_val
        }
        else {
          # If indices are out of bounds, set NA and break the loop
          print("Some top_features indices are out of bounds.")
          # Identify and print the gene(s) not found in colnames
          missing_genes <- setdiff(top_features, colnames(trainTransformed))
          print(paste("The following genes are not found in colnames:", paste(missing_genes, collapse=", ")))
          auc_matrix[j, auc_index] <- NA
          NA_counter <- NA_counter + 1
          next
        } 
      }
      auc_index <- auc_index + 1
    }
  }
  print(NA_counter)
  return(as.data.frame(auc_matrix))
}

eight <- determine_3_cv(pbmc_counts, cv.folds.pbmc, eight_feature_list)        
five <- determine_3_cv(pbmc_counts, cv.folds.pbmc, five_feature_list)
two <- determine_3_cv(pbmc_counts, cv.folds.pbmc, two_feature_list)
three <- determine_3_cv(pbmc_counts, cv.folds.pbmc, three_feature_list)
ten <- determine_3_cv(pbmc_counts, cv.folds.pbmc, ten_feature_list)
six <- determine_3_cv(pbmc_counts, cv.folds.pbmc, six_feature_list)
four <- determine_3_cv(pbmc_counts, cv.folds.pbmc, four_feature_list)


# OUTPUT ANALYSIS ========

## Process output


# INPUTS:
# feature_list, which is a nested list of features
#the length of the feature_list is five, one for each fold
#unless there are zero genes, in which case it may differ
# AUC_matrix, the output of determine_3_cv
#the number of rows in auc matrix should be equivalent and correspond
#to the length of feature_list
#there are 150 columns, each containing the AUC value for the 50x 3-fold
#cv regimen conducted in determine_3_cv
# n_value, which is the number of genes in each nested list in feature list
# OUTPUTS:
# result_df a dataframe with 5 rows and 155 columns
#each row corresponds to a feature set from feature_list
#n refers to the number of genes
#fold_set refers to the fold number that feature set is from, except that
#identifier, which is a unique identifier for that feature set
#mean_auc, which is the mean auc value for each feature set
#the 150 AUC values
# max_auc_row, which is a dataframe that contains the
#row for each n from results_df with the maximum mean_auc value
create_result_df <- function(feature_list, AUC_matrix, n_value) {
  result_df <- data.frame(
    j = integer(),
    fold_set = integer(),
    identifier = character(),
    features = character(),
    mean_auc = numeric()
  )
  
  # Loop through each fold set
  for (fold_set in seq_along(feature_list)) {
    features <- feature_list[[fold_set]][- (n_value+2)]
    auc_values <- as.numeric(AUC_matrix[fold_set, ])
    
    n_genes <- n_value
    
    # Create a unique identifier for that feature set
    identifier <- paste0(n_value, ".", fold_set)
    
    # Determine the mean of the auc values for that row
    mean_auc <- mean(auc_values)
    
    # Collapse the feature list into string so that it is readable and
    # easily exportable
    features_str <- paste(features, collapse = ", ")
    
    # Create a new row for each feature list
    new_row <- data.frame(
      j = n_genes,
      fold_set = fold_set,
      identifier = identifier,
      features = features_str,
      mean_auc = mean_auc
    )
    
    # Add AUC values from AUC dataframe
    for (i in seq_along(auc_values)) {
      new_row[paste0("AUC_", i)] <- auc_values[i]
    }
    
    result_df <- rbind(result_df, new_row)
  }
  # Extract the max auc (defined by mean) feature set for each n
  max_auc_row <- result_df[which.max(result_df$mean_auc), ]
  return(list(result_df = result_df, max_auc_row = max_auc_row))
}

## Run the above function for the feature lists
results_two <- create_result_df(two_feature_list, two, 2)
results_three <- create_result_df(three_feature_list, three, 3)
results_four <- create_result_df(four_feature_list, four, 4)
results_five <- create_result_df(five_feature_list, five, 5)
results_six <- create_result_df(six_feature_list, six, 6)
results_eight <- create_result_df(eight_feature_list, eight, 8)
results_ten <- create_result_df(ten_feature_list, ten, 10)


## Combine all n complete auc dataframes

auc_results_frame <- do.call("rbind", list(results_two$result_df,
                                           results_three$result_df,
                                           results_four$result_df, results_five$result_df,
                                           results_six$result_df, results_eight$result_df,
                                           results_ten$result_df))

## Extract the best perform feature set for each n, except for n = 0

best_results_df <- do.call("rbind", list(results_two$max_auc_row,
                                         results_three$max_auc_row,
                                         results_four$max_auc_row, results_five$max_auc_row,
                                         results_six$max_auc_row, results_eight$max_auc_row,
                                         results_ten$max_auc_row))

write_csv(auc_results_frame, here(output_dir, "age_genes_only_filter_complete_auc_results.csv"))

write_csv(best_results_df, here(output_dir, "age_genes_only_filter_best_results_df.csv"))

## COMPARE DISTRIBUTIONS =====

# Refactor in order and relabel

new_labels <- c(
  "2" = "Two genes", 
  "3" = "Three genes", 
  "4" = "Four genes", 
  "5" = "Five genes", 
  "6" = "Six genes", 
  "8" = "Eight genes", 
  "10" = "Ten genes")

# Update the factor levels in the dataframe
best_results_df <- best_results_df %>%
  mutate(description = ifelse(j != 0, j, features))

best_results_df$description <- factor(best_results_df$description, 
                                      levels = names(new_labels),
                                      labels = new_labels)


# Establish color scheme

custom_colors <- c("#c89fca", "#e600e6",
                   "#99a8d3", "#c5b88d", "#6fb6d5", "#82b38f", "#a8d8be")


# Pivot so boxplot and violin plot can be made

compare <- best_results_df %>%
  pivot_longer(cols = starts_with("AUC_"), 
               names_to = "AUC_variable", 
               values_to = "AUC_value")
# Plot

pretty_plot <- ggplot(compare, aes(x = description, y = AUC_value, fill = description)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  labs(title = "", x = "Input features", y = "AUC distribution") +
  scale_fill_manual(values = custom_colors) +
  my.theme +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

qsave(pretty_plot, here(output_dir, "age_filtered_violin_genes_only.qs"))


ggsave(here(output_dir, "violin.svg"), pretty_plot, width = 4, height = 4)
