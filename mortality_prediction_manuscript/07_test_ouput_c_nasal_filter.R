# Here I perform glm for nasal only counts to evaluate AUC for gene sets on test data
## These gene sets have been filtered 


## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("readr")
library("qs")
library("ggrepel")
library("reshape2")

## Load packages for analysis ======

library("DESeq2")
library("glmnet")
library("edgeR")
library("caret")
library(pROC)


## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")


here()
output_dir <- here("07_test_ouput_c_nasal_filter")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("06_train_eval_c_nasal")



## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## READ IN FILES ======
train_meta <- read_csv("04_assemble_bucket_folds/nasal_complete_train_folds.csv")
test_meta <- read_csv("04_assemble_bucket_folds/nasal_complete_test_split.csv")
nasal_counts <- read_csv("../../../data/nasal-transcriptomics/legacy/2022-10-19/nasal-transcriptomics-Counts.csv")
pbmc_counts_rowFeature <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-RowFeature.csv", skip = 1)
rpM_ct_complete <- read_csv("relationship_CT_rpM/rpM_ct_complete.csv")
feature_table <- read_csv(here(input_dir, "filter_best_results_df_nasal.csv"))


## MERGE META WITH CT =====

test_meta <- rpM_ct_complete %>% 
  select("participant_id", "ct") %>%
  inner_join(., test_meta, by = "participant_id")

train_meta <- rpM_ct_complete %>% 
  select("participant_id", "ct") %>%
  inner_join(., train_meta, by = "participant_id")

## SET UP COUNTS =====

counts <- nasal_counts

counts <- counts %>%
  column_to_rownames(var = "sample_id") %>%
  t(.)


protein_coding <- pbmc_counts_rowFeature %>%
  dplyr::filter(gene_biotype == "protein_coding")

# Select protein coding genes only

protein_coding_list <- protein_coding$gene_id
ensembl_to_gene <- setNames(protein_coding$gene_name, protein_coding$gene_id)

# The counts dataframe is now a dataframe with features as rows and sample_ids as columns
# This is the format DeSeqDataSetFromMatrix is able to read

# Only keep protein-coding genes
counts <- counts[rownames(counts) %in% protein_coding_list, ]


## EVALUATE FEATURES
set.seed(555)

set.seed(555)

## INPUTS:
# counts is a dataframe containing samples as columns and genes as rows
# train_meta is a dataframe containing the metadata, including the ct value and age
#for the training subset
# test_meta is a dataframe containing the metadata, including the ct value and age
#for the test subset
# features is a list containing the features for testing
determine_performance_gene <- function(counts, train_meta, test_meta, features) {
  
  
  # Set up training counts
  
  # Only keep the samples (columns) that are used for training
  counts.train <- counts[,train_meta$sample_id]
  
  # Creates a boolean of the genes/features that pass the QC threshold
  keep.train <- rowSums(counts.train >= 10) >= (0.2*ncol(counts.train))
  
  # Set up the test data
  counts.test <- counts[,test_meta$sample_id]
  
  
  # Create a DESeqDataset from our counts matrix
  dds.train <- DESeqDataSetFromMatrix(
    countData = counts.train[keep.train,], # Only keep the genes/features
    # that pass the QC threshold
    colData = train_meta,
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
    colData = test_meta,
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
  
  
  print("Preprocessing...")
  
  x <- t(vsd.train)
  test <- t(vsd.test)
  
  ## Add additional features
  train_meta_subset <- train_meta %>%
    select("sample_id", "admit_age", "ct")
  
  test_meta_subset <- test_meta %>%
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
  
  y_train <- train_meta$group == 1
  y_test <- test_meta$group == 1
  y_train_binary <- ifelse(y_train == TRUE, 1, 0)
  y_test_binary <- ifelse(y_test == TRUE, 1, 0)
  
  # Standardize train and test data
  
  preProcValues <- preProcess(x, method = c("center", "scale"))
  
  trainTransformed <- predict(preProcValues, x)
  testTransformed <- predict(preProcValues, test)
  
  
  print("Preprocessing complete!")
  
  # Create matrices
  top_features <- features
  
  # Fit logistic regression on the training set ONLY
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
    
    results_predictions <- data.frame(
      observed = y_test_binary,
      predicted = pred_logistic
    )
    
    roc.car <- roc(
      response = y_test_binary, 
      predictor = pred_logistic,
      direction = "<",
      plot = TRUE)
    
    auc_val <- as.numeric(auc(roc.car))
    roc_plot <- ggroc(roc.car, legacy.axes = TRUE)
    print(auc_val)
    
    ci <- ci.auc(roc.car, method = "bootstrap", boot.n = 5000,
                 boot.stratified = TRUE)
    
    # Youden's index
    youden <- pROC::coords(
      roc.car, x="best", best.method="youden",
      ret=c("threshold","specificity","sensitivity","accuracy","precision"))
    
    # 90% sensitivity
    sens90 <- pROC::coords(
      roc.car, x=0.9, input="sensitivity",
      ret=c("threshold","specificity","sensitivity","accuracy","precision"))
    
    # Confusion matrix
    optimal_threshold <- youden$threshold
    
    # Convert predicted probabilities to binary predictions
    predicted_binary <- ifelse(pred_logistic >= optimal_threshold, 1, 0)
    
    # Create a data frame for binary predictions
    binary_results_predictions <- data.frame(
      observed = y_test_binary,
      predicted = predicted_binary
    )
    
    confusion_matrix <- confusionMatrix(as.factor(predicted_binary), as.factor(y_test_binary),
                                        positive = "1")
    
    
  }
  else {
    # If indices are out of bounds, set NA and break the loop
    print("Some top_features indices are out of bounds.")
    # Identify and print the gene(s) not found in colnames
    missing_genes <- setdiff(top_features, colnames(trainTransformed))
    print(paste("The following genes are not found in colnames:", paste(missing_genes, collapse=", ")))
  } 
  
  
  return(list(model = glm_model, 
              roc = roc.car, 
              auc = auc_val, 
              features = features,
              roc_plot = roc_plot,
              sens90 = sens90,
              youden = youden,
              test_pred = results_predictions,
              confusion = confusion_matrix, 
              ci = ci))
}

## SELECT FEATURES ====

my_features <- feature_table %>%
  filter(j == 3) %>%
  pull(features) %>%
  strsplit(., ",") %>%
  .[[1]] %>%
  trimws()

three_features <- c(my_features[1:3], "admit_age", "ct")

three_genes <- c(my_features[1:3])


# EVALUATE =====

OLAH_integrated <- determine_performance_gene(counts, train_meta, test_meta, c("ENSG00000152463", "admit_age", "ct"))
OLAH_solo <- determine_performance_gene(counts, train_meta, test_meta, c("ENSG00000152463"))

three_eval <- determine_performance_gene(counts, train_meta, test_meta, three_features)
three_genes_eval <- determine_performance_gene(counts, train_meta, test_meta, three_genes)

evaluation_list <- list(three_eval, three_genes_eval, OLAH_integrated)
names(evaluation_list) <- c("Three features + Age + CT", "Three features", "NS OLAH + Age + CT")

## BUILD DATAFRAME OF OUTPUT =======
result_df <-  data.frame(features = character(),
                         auc_values = numeric(),
                         youden_threshold = numeric(),
                         specificity_youden = numeric(),
                         sensitivity_youden = numeric(),
                         accuracy_youden = numeric(),
                         interval_low = numeric(),
                         interval_high = numeric())

for (j in seq_along(evaluation_list)) {
  output <- evaluation_list[[j]]
  name <- names(evaluation_list)[[j]]
  auc <- output[["auc"]]
  youden_threshold <- output$youden$threshold
  specificity <- output$youden$specificity
  sensitivity <- output$youden$sensitivity
  accuracy <- output$youden$accuracy
  i_low <- output$ci[["2.5%"]]
  i_high <- output$ci[["97.5%"]]
  
  new_row <- data.frame(
    features = name,
    auc_values = auc,
    youden_threshold = youden_threshold,
    specificity_youden = specificity,
    sensitivity_youden = sensitivity,
    accuracy_youden = accuracy,
    interval_low = i_low,
    interval_high = i_high
  )
  result_df <- rbind(result_df, new_row)
}

write_csv(result_df, here(output_dir, "nasal_filter_auc_roc_result_df.csv"))

## PLOTTING =======
custom_colors <- c("#94d6e0",
                   "#acaba8",
                   "#639e6c",
                   "#eadca1")

ct_col <- custom_colors[1]
age_col <- custom_colors[2]
col_complete <- custom_colors[3]
col_three <- custom_colors[4]

# Nasal supp roc
main_roc <- c(evaluation_list[1:2])
roc_data_list <- lapply(main_roc, function(x) {
  ggplot_build(x$roc_plot)$data[[1]]
})

# Combine all data into a single dataframe
combined_roc_data <- do.call(rbind, roc_data_list)

combined_roc_data$Source <- rep(names(roc_data_list), sapply(roc_data_list, nrow))

color_mapping <- c(
  "Three features + Age + CT" = col_complete,
  "Three features" = col_three
)

combined_roc_data$Source <- factor(combined_roc_data$Source)

write_csv(combined_roc_data, here(output_dir, "nasal_main_roc_data.csv"))


combined_plot <- ggplot(combined_roc_data, aes(x = x, y = y, color = Source)) +
  geom_line() +
  scale_color_manual(values = color_mapping) +
  labs(x = "(1 - Specificity)", y = "(Sensitivity)", color = "") +
  my.theme +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add border around the plot
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dashed", color = "grey") 

ggsave(here(output_dir, "nasal_combined_plot.svg"), combined_plot, width = 4.5, height = 5)


## Extract and save OLAH data

OLAH_roc <- c(evaluation_list[3])
OLAH_roc_data_list <- lapply(OLAH_roc, function(x) {
  ggplot_build(x$roc_plot)$data[[1]]
})

# Combine all data into a single dataframe
OLAH_roc_data <- do.call(rbind, OLAH_roc_data_list)

OLAH_roc_data$Source <- rep(names(OLAH_roc_data_list), sapply(OLAH_roc_data_list, nrow))


OLAH_roc_data$Source <- factor(OLAH_roc_data$Source)

write_csv(OLAH_roc_data, here(output_dir, "nasal_OLAH_roc_data.csv"))


