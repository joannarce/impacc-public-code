# Train test on OLINK


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
output_dir <- here("OLINK_01_train_test")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("OLINK_00_setup")



## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## READ IN FILES ======
train_meta <- read.csv("OLINK_00_setup/train_meta.csv")
test_meta <- read.csv("OLINK_00_setup/test_meta.csv")
counts <- read.csv("OLINK_00_setup/olink_counts.csv", row.names = 1)
rpM_ct_complete <- read.csv("relationship_CT_rpM/rpM_ct_complete.csv")



## MERGE META WITH CT =====

test_meta <- rpM_ct_complete %>% 
  select("participant_id", "ct") %>%
  inner_join(., test_meta, by = "participant_id")

train_meta <- rpM_ct_complete %>% 
  select("participant_id", "ct") %>%
  inner_join(., train_meta, by = "participant_id")

## SET UP COUNTS =====

counts <- counts %>%
  column_to_rownames(var = "sample_id") %>%
  t(.)




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
  
  
  # Set up counts
  counts.train <- counts[,train_meta$sample_id]
  counts.test <- counts[,test_meta$sample_id]
  
  
  x <- t(counts.train)
  test <- t(counts.test)
  
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
  
  y_train <- train_meta$trajectory == 1
  y_test <- test_meta$trajectory == 1
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


# EVALUATE =====
IL_6 <- determine_performance_gene(counts, train_meta, test_meta, c("admit_age", "ct", "IL6"))
IL_solo <- determine_performance_gene(counts, train_meta, test_meta, c("IL6"))

integrated.roc <- ggplot_build(IL_6$roc_plot)$data[[1]]
solo.roc <- ggplot_build(IL_solo$roc_plot)$data[[1]]



integrated.roc$Source <- "IL-6 + age + CT, AUC = 0.73"
solo.roc$Source <- "IL-6, AUC = 0.63"


roc <- rbind(integrated.roc, solo.roc)

IL_6.plot <- ggplot(roc, aes(x = x, y = y, color = Source)) +
  geom_line() +
  labs(x = "(1 - Specificity)", y = "(Sensitivity)", color = "") +
  my.theme +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none"# Add border around the plot
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dashed", color = "grey") 

