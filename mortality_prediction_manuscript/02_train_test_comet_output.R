# Here I train and test a glm model using a three gene+age classifier on the COMET data
## The three gene sets have been chosen based on the IMPACC data


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

i_am("COMET.Rproj")


here()
output_dir <- here("02_train_test_comet_output")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

#input_dir <- here("00_train_test")



## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## READ IN FILES ======
meta <- read.csv(here("comet_meta.csv"))
#test_meta <- read.csv(here(input_dir, "test_meta.csv"))

vax <-read.csv(here("COMET_vax_meta.csv"))

meta <- meta %>%
  dplyr::inner_join(., vax, by = "Patient.ID")

counts <- read.csv("counts_preQC.csv")



## SET UP COUNTS =====


counts <- counts %>%
  column_to_rownames(var = "X") %>%
  dplyr::select(-gene_name)


# The counts dataframe is now a dataframe with features as rows and File.names as columns
# This is the format DeSeqDataSetFromMatrix is able to read


## EVALUATE FEATURES


## INPUTS:
# counts is a dataframe containing samples as columns and genes as rows
# train_meta is a dataframe containing the metadata, including the age and trajectory
#for the training subset
# test_meta is a dataframe containing the metadata, including the age and trajectory
#for the test subset
# features is a list containing the features for testing
determine_performance_gene <- function(counts, train_meta, test_meta, features) {
  
  
  # Set up training counts
  
  # Only keep the samples (columns) that are used for training
  counts.train <- counts[,train_meta$File.name]
  
  # Creates a boolean of the genes/features that pass the QC threshold
  keep.train <- rowSums(counts.train >= 10) >= (0.2*ncol(counts.train))
  
  # Set up the test data
  counts.test <- counts[,test_meta$File.name]
  
  
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
    dplyr::select("File.name", "Age.at.Admission")
  
  test_meta_subset <- test_meta %>%
    dplyr::select("File.name", "Age.at.Admission")
  
  x <- x %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "File.name") %>% 
    inner_join(train_meta_subset, by = "File.name") %>% 
    tibble::column_to_rownames(var = "File.name") %>%
    as.matrix(.)
  
  test <- test %>% 
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "File.name") %>% 
    inner_join(test_meta_subset, by = "File.name") %>% 
    tibble::column_to_rownames(var = "File.name") %>%
    as.matrix(.)
  
  ## Set up y data
  
  y_train <- train_meta$Deceased == "Yes"
  y_test <- test_meta$Deceased == "Yes"
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

my_features <- c("Age.at.Admission", "ENSG00000112149", "ENSG00000146122", "ENSG00000129244")
sweeney <- c("ENSG00000160883", "ENSG00000112799", "ENSG00000105329",
             "ENSG00000164821", "ENSG00000156127", "ENSG00000223865",
             "Age.at.Admission")
olah <- c("Age.at.Admission", "ENSG00000152463")
# Perform 5-fold cross-validation with AUC calculation
set.seed(555)

meta <- meta %>%
  mutate(stratify_group = interaction(Deceased, Batch.year, drop = TRUE))
folds <- createFolds(c(meta$stratify_group), k = 5, list = TRUE)

process_fold <- function(fold_index, folds, meta, counts, features) {
  # Split into training and testing folds
  test_indices <- folds[[fold_index]]
  train_indices <- setdiff(seq_len(nrow(meta)), test_indices)
  
  train_fold <- meta[train_indices, ]
  test_fold <- meta[test_indices, ]
  
  # Run the model
  results <- determine_performance_gene(counts, train_fold, test_fold, features)
  
  # Extract test predictions and merge with metadata
  test_results <- results$test_pred
  test_results$Patient_ID <- test_fold$File.name
  test_results$Sex <- test_fold$Sex
  test_results$Age <- test_fold$Age.at.Admission
  test_results$Vax <- test_fold$vax_status
  
  return(list(
    auc = results$auc,
    roc = results$roc,
    roc_plot = results$roc_plot,
    results_df = test_results,
    youdens = results$youden
  ))
}

# Initialize storage
auc_values <- numeric(length(folds))
roc_list <- list()
roc_plot <- list()
youdens_index <- numeric(length(folds))

olah.y <- numeric(length(folds))
olah.auc <- numeric(length(folds))
olah.out <- list()

# Run the function for each fold and feature set
final_results <- data.frame()
#sweeney_results_df <- data.frame()
olah_results_df <- data.frame()

for (i in seq_along(folds)) {
  # Process results for each feature set
  my_features_result <- process_fold(i, folds, meta, counts, my_features)
 # sweeney_result <- process_fold(i, folds, meta, counts, sweeney)
  olah_result <- process_fold(i, folds, meta, counts, olah)
  
  # Store AUC and ROC results
  auc_values[i] <- my_features_result$auc
  roc_list[[i]] <- my_features_result$roc
  roc_plot[[i]] <- my_features_result$roc_plot
  youdens_index[i] <- my_features_result$youdens$threshold
  olah.y[i] <- olah_result$youdens$threshold
  olah.auc[i] <- olah_result$auc
  olah.out[[i]] <- olah_result$roc_plot
  
  # Combine results
  final_results <- rbind(final_results, my_features_result$results_df)
 # sweeney_results_df <- rbind(sweeney_results_df, sweeney_result$results_df)
  olah_results_df <- rbind(olah_results_df, olah_result$results_df)
}

# Rename columns for clarity in all results
col_names <- c("True_Value", "Predicted_Value", "Patient_ID", "Sex", "Age", "Vaccine_Status")
colnames(final_results) <- col_names
#colnames(sweeney_results_df) <- col_names
colnames(olah_results_df) <- col_names

write.csv(x = final_results, here(output_dir, "three_gene_output.csv"))
write.csv(x = olah_results_df, here(output_dir, "olah_output.csv"))

final_results$c <- "Three genes + age classifier"
olah_results_df$c <- "OLAH + age classifier"

plotting.scores <- rbind(final_results, olah_results_df)


three.youden <- mean(youdens_index)
olah.youden <- mean(olah.y)

youden <- data.frame(
  yi = c(three.youden, olah.youden),
  c = c("Three genes + age classifier", "OLAH + age classifier" )
)

# Plot scores =========

plotting.scores$c <- as.factor(plotting.scores$c)

plotting.scores.survival <- plotting.scores %>%
  dplyr::filter(True_Value == 0)

plotting.scores.mortality <- plotting.scores %>%
  dplyr::filter(True_Value == 1)


plot.s <- ggplot(data = na.omit(plotting.scores.survival), aes(x = Vaccine_Status, y = Predicted_Value, fill = Vaccine_Status)) +
  geom_boxplot() +
  facet_wrap(~ c) +
  xlab("Vaccinated status") +
  ylab("Classifier score") +
  ggtitle("Scores for survival trajectory") +
  scale_fill_manual(values = c("purple", "orange")) + 
  my.theme + 
  theme(legend.position = "none") +
  ylim(c(0, 1)) +
  geom_hline(data = youden, aes(yintercept = yi), linetype = "dashed", color = "grey60") 
  

plot.d <- ggplot(data = na.omit(plotting.scores.mortality), aes(x = Vaccine_Status, y = Predicted_Value, fill = Vaccine_Status)) +
  geom_boxplot() +
  facet_grid(~ c) +
  xlab("Vaccinated status") +
  ylab("Classifier score") +
  scale_fill_manual(values = c("purple", "orange")) + 
  my.theme + 
  theme(legend.position = "none") +
  geom_hline(data = youden, aes(yintercept = yi), linetype = "dashed", color = "grey60") +
  ylim(c(0, 1))


ggsave(here("survival_scores.svg"), plot.s, width = 4, height = 4)

ggsave(here("mortality_scores.svg"), plot.d, width = 4.5, height = 4.5)

# For survival
survival_pvals <- na.omit(plotting.scores.survival) %>%
  group_by(c) %>%
  summarise(p_value = wilcox.test(Predicted_Value ~ Vaccine_Status)$p.value)

# Mortality p-values
mortality_pvals <- na.omit(plotting.scores.mortality) %>%
  group_by(c) %>%
  summarise(p_value = wilcox.test(Predicted_Value ~ Vaccine_Status)$p.value)


# Plot ROC for vaccinated group based on classifier ======

vax <- plotting.scores %>%
  dplyr::filter(Vaccine_Status == "Yes")

three.vy <- vax %>%
  dplyr::filter(c == "Three genes + age classifier")

olah.vy <- vax %>%
  dplyr::filter(c == "OLAH + age classifier")

three.roc.vy <- roc(
  response = as.factor(three.vy$True_Value), 
  predictor = as.numeric(three.vy$Predicted_Value),
  direction = "<",
  plot = TRUE)

three.roc.vy.auc <- as.numeric(auc(three.roc.vy))
three.roc.vy.auc.plot <- ggroc(three.roc.vy, legacy.axes = TRUE)


olah.vy <- roc(
  response = as.factor(olah.vy$True_Value), 
  predictor = as.numeric(olah.vy$Predicted_Value),
  direction = "<",
  plot = TRUE)

olah.vy.auc <- as.numeric(auc(olah.vy))
olah.vy.auc.plot <- ggroc(olah.vy, legacy.axes = TRUE)

olah.vy.ci <- ci.auc(olah.vy, method = "bootstrap", boot.n = 5000,
             boot.stratified = TRUE)

three.vy.ci <- ci.auc(three.roc.vy, method = "bootstrap", boot.n = 5000,
                     boot.stratified = TRUE)


three.vy.df <- ggplot_build(three.roc.vy.auc.plot)$data[[1]]
three.vy.df$c <- "3-gene"

olah.vy.df <- ggplot_build(olah.vy.auc.plot)$data[[1]]
olah.vy.df$c <- "OLAH"


vy.roc.df <- rbind(three.vy.df, olah.vy.df)

vy.roc.plot <- ggplot(vy.roc.df, aes(x = x, y = y, color = c)) +
  geom_line() + 
  scale_color_manual(values = c("OLAH" = "darkorange", "3-gene" = "deepskyblue3")) +
  labs(x = "(1 - Specificity)", y = "(Sensitivity)", color = "") +
  my.theme +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add border around the plot
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dashed", color = "grey") 

ggsave(here("vaccinated_roc.svg"), vy.roc.plot, width = 4.4, height = 4.5)



# Combine ROC curves for plotting ======


roc_data_list <- lapply(seq_along(roc_plot), function(i) {
  df <- ggplot_build(roc_plot[[i]])$data[[1]]
  df$Fold <- paste0("Fold_", i)
  df
})

# Combine all fold data
combined_roc_data <- bind_rows(roc_data_list)

# Calculate the average ROC curve
x_seq <- sort(unique(combined_roc_data$x))
interpolated_rocs <- lapply(unique(combined_roc_data$Fold), function(fold) {
  fold_data <- combined_roc_data[combined_roc_data$Fold == fold, ]
  x_vals <- c(0, fold_data$x)
  y_vals <- c(0, fold_data$y)

  step_f <- stepfun(x_vals, c(y_vals, tail(y_vals, 1)))
  data.frame(
    x = x_seq,
    y = step_f(x_seq),
    Fold = fold
  )
})


interpolated_combined <- bind_rows(interpolated_rocs)

# Calculate average
average_roc <- interpolated_combined %>%
  group_by(x) %>%
  summarize(
    y = mean(y),
    Fold = "Average"
  )



final_roc_data <- bind_rows(combined_roc_data, average_roc)



## Now, add the OLAH results =====

roc_olah_list <- lapply(seq_along(olah.out), function(i) {
  roc_plot_i <- olah.out[[i]]
  df <- ggplot_build(roc_plot_i)$data[[1]]
  df$Fold <- paste0("O.Fold_", i)
  
  return(df)
})

# Combine all folds into a single data frame
roc_olah_df <- do.call(rbind, roc_olah_list)


olah_x_seq <- sort(unique(roc_olah_df$x))
interpolated_rocs.o <- lapply(unique(roc_olah_df$Fold), function(fold) {
  fold_data <- roc_olah_df[roc_olah_df$Fold == fold, ]
  x_vals <- c(0, fold_data$x)
  y_vals <- c(0, fold_data$y)
  
  step_f <- stepfun(x_vals, c(y_vals, tail(y_vals, 1)))
  data.frame(
    x = olah_x_seq,
    y = step_f(olah_x_seq),
    Fold = fold
  )
})


interpolated_combined.o <- bind_rows(interpolated_rocs.o)

# Calculate average
average_roc.o <- interpolated_combined.o %>%
  group_by(x) %>%
  summarize(
    y = mean(y),
    Fold = "O.Average"
  )



final_roc_data.o <- bind_rows(roc_olah_df , average_roc.o)


complete_roc <- rbind(final_roc_data.o, final_roc_data)


## Plot =====


color_mapping <- c(
  "Fold_1" = "#f795f7",
  "Fold_2" = "#f795f7",
  "Fold_3" = "#f795f7",
  "Fold_4" = "#f795f7",
  "Fold_5" = "#f795f7",
  "Average" = "#ff00ff",
  "O.Fold_1" = "#92bd96",
  "O.Fold_2" = "#92bd96",
  "O.Fold_3" = "#92bd96",
  "O.Fold_4" = "#92bd96",
  "O.Fold_5" = "#92bd96",
  "O.Average" = "#0afa1d"
)

avgs <- c("O.Average", "Average")

complete_roc$Linewidth <- ifelse(complete_roc$Fold %in% avgs, 1, 0.3)

combined_plot <- ggplot(complete_roc, aes(x = x, y = y, color = Fold)) +
  geom_line(aes(size = Linewidth)) + 
  scale_size_identity() +
  scale_color_manual(values = color_mapping) +
  labs(x = "(1 - Specificity)", y = "(Sensitivity)", color = "") +
  my.theme +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add border around the plot
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dashed", color = "grey") 


qsave(x = complete_roc, file = "rocplot.qs")
write.csv(x = complete_roc, file = "rocplot.csv")
qsave(x = folds, file = "folds.qs")
ggsave(here("roc_combined_plot.svg"), combined_plot, width = 4.4, height = 4.5)



olah_aucs <- c(0.9130435, 0.7130435, 0.8869565, 0.8369565, 0.6145833)
mean(olah_aucs)
sd(olah_aucs)

mean(auc_values)
sd(auc_values)



## Plot main AUC ==============

t.scores <- plotting.scores %>%
  dplyr::filter(c == "Three genes + age classifier")

o.scores <- plotting.scores %>%
  dplyr::filter(c == "OLAH + age classifier")


three.roc <- roc(
  response = as.factor(t.scores$True_Value), 
  predictor = as.numeric(t.scores$Predicted_Value),
  direction = "<",
  plot = TRUE)

o.roc <- roc(
  response = as.factor(o.scores$True_Value), 
  predictor = as.numeric(o.scores$Predicted_Value),
  direction = "<",
  plot = TRUE)


three.roc.auc <- as.numeric(auc(three.roc))
three.roc.auc.plot <- ggroc(three.roc, legacy.axes = TRUE)
three.roc.ci <- ci.auc(three.roc, method = "bootstrap", boot.n = 5000,
                      boot.stratified = TRUE)

o.roc.auc <- as.numeric(auc(o.roc))
o.roc.auc.plot <- ggroc(o.roc, legacy.axes = TRUE)
o.roc.ci <- ci.auc(o.roc, method = "bootstrap", boot.n = 5000,
                       boot.stratified = TRUE)


three.df <- ggplot_build(three.roc.auc.plot)$data[[1]]
three.df$c <- "3-gene"

o.df <- ggplot_build(o.roc.auc.plot)$data[[1]]
o.df$c <- "OLAH"

roc.df <- rbind(three.df, o.df)

roc.plot <- ggplot(roc.df, aes(x = x, y = y, color = c)) +
  geom_line() + 
  scale_color_manual(values = c("OLAH" = "#0afa1d", "3-gene" = "#ff00ff")) +
  labs(x = "(1 - Specificity)", y = "(Sensitivity)", color = "") +
  my.theme +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add border around the plot
  ) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dashed", color = "grey") 

ggsave(here("oof_roc.svg"), roc.plot, width = 4.4, height = 4.5)

