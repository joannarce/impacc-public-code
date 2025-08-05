# Here I perform glm for pbmc only counts to evaluate AUC for gene sets on test data
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
output_dir <- here("07_test_output_b_filter")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("06_train_eval_a")



## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## READ IN FILES ======
train_meta <- read_csv("04_assemble_bucket_folds/pbmc_complete_train_folds.csv")
test_meta <- read_csv("04_assemble_bucket_folds/pbmc_complete_test_split.csv")
pbmc_counts <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-Counts.csv")
pbmc_counts_rowFeature <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-RowFeature.csv", skip = 1)
rpM_ct_complete <- read_csv("relationship_CT_rpM/rpM_ct_complete.csv")
feature_table <- read_csv(here(input_dir, "filter_best_results_df.csv"))


## MERGE META WITH CT =====

basic_metadata <- read.csv(here(input_dir, "basic_metadata.csv")) 

test_meta <- rpM_ct_complete %>% 
  select("participant_id", "ct") %>%
  inner_join(., test_meta, by = "participant_id") %>%
  inner_join(., basic_metadata, by = "participant_id")

train_meta <- rpM_ct_complete %>% 
  select("participant_id", "ct") %>%
  inner_join(., train_meta, by = "participant_id") %>%
  inner_join(., basic_metadata, by = "participant_id")

## SET UP COUNTS =====

counts <- pbmc_counts

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
    select("sample_id", "admit_age", "ct", "SOFA", "baseline_lab_lymph", "baseline_lab_crp",
           "resp_status_v1")
  
  test_meta_subset <- test_meta %>%
    select("sample_id", "admit_age", "ct", "SOFA", "baseline_lab_lymph", "baseline_lab_crp",
           "resp_status_v1")
  
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

#sweeney <- c("ENSG00000160883", "ENSG00000112799", "ENSG00000105329",
#            "ENSG00000164821", "ENSG00000156127", "ENSG00000223865",
#            "ct", "admit_age")


# EVALUATE =====
s.OLAH_eval <- determine_performance_gene(counts, train_meta, test_meta, c("ENSG00000152463", "SOFA"))
s.OLAH_integrated <- determine_performance_gene(counts, train_meta, test_meta, c("ENSG00000152463", "admit_age", "ct", "SOFA"))
s.three_eval <- determine_performance_gene(counts, train_meta, test_meta, c(three_features, "SOFA"))
s.three_genes_eval <- determine_performance_gene(counts, train_meta, test_meta, c(three_genes, "SOFA"))

s.s <- determine_performance_gene(counts, train_meta, test_meta, "SOFA")
lymph <- determine_performance_gene(counts, train_meta, test_meta, "baseline_lab_lymph")
crp <- determine_performance_gene(counts, train_meta, test_meta, "baseline_lab_crp")
resp <- determine_performance_gene(counts, train_meta, test_meta, "resp_status_v1")

OLAH_eval <- determine_performance_gene(counts, train_meta, test_meta, c("ENSG00000152463"))
OLAH_integrated <- determine_performance_gene(counts, train_meta, test_meta, c("ENSG00000152463", "admit_age", "ct"))
three_eval <- determine_performance_gene(counts, train_meta, test_meta, three_features)
OLAH_friends <- determine_performance_gene(counts, train_meta, test_meta, c(three_features, "ENSG00000152463"))
age_three_features <- determine_performance_gene(counts, train_meta, test_meta, c(three_genes, "admit_age"))
three_genes_eval <- determine_performance_gene(counts, train_meta, test_meta, three_genes)
sweeney_complete <- determine_performance_gene(counts, train_meta, test_meta, sweeney)
sweeney_genes_subset <- determine_performance_gene(counts, train_meta, test_meta, sweeney[1:6])
age_only <- determine_performance_gene(counts, train_meta, test_meta, c("admit_age"))
ct_only <- determine_performance_gene(counts, train_meta, test_meta, c("ct"))

## De long test

compare_values <- list(
  list(s.s$roc, three_eval$roc),
  list(s.s$roc, three_genes_eval$roc),
  list(s.s$roc, OLAH_integrated$roc),
  list(s.s$roc, OLAH_eval$roc)
)

model_names.s <- list(
  c("sofa", "three_eval"),
  c("sofa", "three_genes_eval"),
  c("sofa", "OLAH_integrated"),
  c("sofa", "OLAH_eval")
)
res_compare_sofa_df <- do.call(rbind, lapply(seq_along(compare_values), function(i) {
  test_res <- roc.test(compare_values[[i]][[1]], compare_values[[i]][[2]], method = "delong")
  data.frame(
    model1 = model_names.s[[i]][1],
    model2 = model_names.s[[i]][2],
    p_value = test_res$p.value
  )
}))

lymph_compare <- list(
  list(lymph$roc, three_eval$roc),
  list(lymph$roc, three_genes_eval$roc),
  list(lymph$roc, OLAH_integrated$roc),
  list(lymph$roc, OLAH_eval$roc)
)


model_names <- list(
  c("lymph", "three_eval"),
  c("lymph", "three_genes_eval"),
  c("lymph", "OLAH_integrated"),
  c("lymph", "OLAH_eval")
)

# Perform roc.test and extract p-values with names
res_compare_lymph_df <- do.call(rbind, lapply(seq_along(lymph_compare), function(i) {
  test_res <- roc.test(lymph_compare[[i]][[1]], lymph_compare[[i]][[2]], method = "delong")
  data.frame(
    model1 = model_names[[i]][1],
    model2 = model_names[[i]][2],
    p_value = test_res$p.value
  )
}))

write.csv(res_compare_lymph_df, here(output_dir, "res_compare_lymph_df.csv"))

## Downstream analysis

evaluation_list <- list(three_eval, three_genes_eval, sweeney_complete, sweeney_genes_subset,
                        age_only, ct_only, OLAH_eval, OLAH_integrated)
names(evaluation_list) <- c("Three features + Age + CT", "Three features", "Sweeney et al. + Age + CT",
                            "Sweeney et al. features only", "Age", "CT", "OLAH", "PBMC OLAH + Age + CT")

OLAH_list <- list(OLAH_eval, OLAH_integrated)
names(OLAH_list) <- c("OLAH", "OLAH + Age + CT")

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

write_csv(result_df, here(output_dir, "filter_auc_roc_result_df.csv"))


## AGE SUBSET EVAL ====== 
three_eval_age_pred_subset_young <- three_eval$test_pred %>%
  rownames_to_column("sample_id") %>%
  merge(., test_meta, by = "sample_id") %>%
  filter(admit_age <= 45)


predicted_binary_young <- ifelse(three_eval_age_pred_subset_young$predicted >= three_eval$youden$threshold, 1, 0)
accuracy_young <- mean(predicted_binary_young == three_eval_age_pred_subset_young$observed)


three_eval_age_pred_subset_old <- three_eval$test_pred %>%
  rownames_to_column("sample_id") %>%
  merge(., test_meta, by = "sample_id") %>%
  filter(admit_age >= 70)

predicted_binary_old <- ifelse(three_eval_age_pred_subset_old$predicted >= three_eval$youden$threshold, 1, 0)

# Calculate accuracy or other metrics
accuracy_old <- mean(predicted_binary_old == three_eval_age_pred_subset_old$observed)


roc_old <- roc(
  response = three_eval_age_pred_subset_old$observed, 
  predictor = three_eval_age_pred_subset_old$predicted,
  direction = "<",
  plot = TRUE
)

auc_old <- as.numeric(auc(roc_old))

youden_old <- pROC::coords(
  roc_old, x="best", best.method="youden",
  ret=c("threshold","specificity","sensitivity","accuracy","precision"))

## PLOTTING =======
custom_colors <- c("#74d6e0",
                   "#ecaba8",
                   "#ff00ff",
                   "#badca1",
                   "#EBB41D",
                   "#a28bc4")

ct_col <- custom_colors[1]
age_col <- custom_colors[2]
col_complete <- custom_colors[3]
col_three <- custom_colors[4]
col_sweeney_complete <- custom_colors[5]
col_sweeney_genes <- custom_colors[6]


# Figure 4b
main_roc <- c(evaluation_list[1:2], evaluation_list[5:6])
roc_data_list <- lapply(main_roc, function(x) {
  ggplot_build(x$roc_plot)$data[[1]]
})

# Combine all data into a single dataframe
combined_roc_data <- do.call(rbind, roc_data_list)

combined_roc_data$Source <- rep(names(roc_data_list), sapply(roc_data_list, nrow))

color_mapping <- c(
  "Three features + Age + CT" = col_complete,
  "Three features" = col_three,
  "CT" = ct_col,
  "Age" = age_col
)

combined_roc_data$Source <- factor(combined_roc_data$Source)

write_csv(combined_roc_data, here(output_dir, "main_roc_data_pbmc.csv"))


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

ggsave(here(output_dir, "combined_plot.svg"), combined_plot, width = 4, height = 3.5)


## supp roc, compare with sweeney


supp_roc <- c(evaluation_list[1:4])
supp_roc_data_list <- lapply(supp_roc, function(x) {
  ggplot_build(x$roc_plot)$data[[1]]
})

# Combine all data into a single dataframe
supp_roc_data <- do.call(rbind, supp_roc_data_list)

supp_roc_data$Source <- rep(names(supp_roc_data_list), sapply(supp_roc_data_list, nrow))

supp_color_mapping <- c(
  "Three features + Age + CT" = col_complete,
  "Three features" = col_three,
  "Sweeney et al. + Age + CT" = col_sweeney_complete,
  "Sweeney et al. features only" = col_sweeney_genes
)

supp_roc_data$Source <- factor(supp_roc_data$Source)

write_csv(supp_roc_data, here(output_dir, "supp_roc_data_pbmc.csv"))


supp_combined_plot <- ggplot(supp_roc_data, aes(x = x, y = y, color = Source)) +
  geom_line() +
  scale_color_manual(values = supp_color_mapping) +
  labs(x = "(1 - Specificity)", y = "(Sensitivity)", color = "") +
  my.theme +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add border around the plot
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dashed", color = "grey") 

ggsave(here(output_dir, "supp_roc_plot.svg"), supp_combined_plot, width = 6, height = 6)


## Figure 4 e ========

OLAH_roc <- OLAH_list[2]
roc_OLAH_list <- lapply(OLAH_roc, function(x) {
  ggplot_build(x$roc_plot)$data[[1]]
})

# Combine all data into a single dataframe
OLAH_roc_data <- do.call(rbind, roc_OLAH_list)

OLAH_roc_data$Source <- rep(names(roc_OLAH_list), sapply(roc_OLAH_list, nrow))

NS_OLAH <- read_csv(here("07_test_ouput_c_nasal_filter", "nasal_OLAH_roc_data.csv"))

OLAH_roc_data <- rbind(OLAH_roc_data, NS_OLAH)

olah_mapping <- c(
  "NS OLAH + Age + CT" = "darkblue",
  "PBMC OLAH + Age + CT" = "#0afa1d"
)

olah_plot <- ggplot(OLAH_roc_data, aes(x = x, y = y, color = Source)) +
  geom_line() +
  scale_color_manual(values = olah_mapping) +
  labs(x = "(1 - Specificity)", y = "(Sensitivity)", color = "") +
  my.theme +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add border around the plot
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dashed", color = "grey") 

ggsave(here(output_dir, "olah_combined_plot.svg"), olah_plot, width = 4, height = 3.5)


## Figure S5A =====
# Extract ROC data from a single ggplot object
roc_data <- ggplot_build(OLAH_friends$roc_plot)$data[[1]]
roc_data$Source <- "Friends"

# Create the plot
friends_plot <- ggplot(roc_data, aes(x = x, y = y, color = "brown")) +
  geom_line() +
  labs(x = "(1 - Specificity)", y = "(Sensitivity)", color = "") +
  my.theme +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dashed", color = "grey")

# Save the plot
ggsave(here(output_dir, "friends_combined_plot.svg"), friends_plot, width = 4, height = 3.5)


write_csv(roc_data, here(output_dir, "supp_olah_friends.csv")) 
