# Here I divide each bucket into five folds for five fold CV


## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("janitor")
library("readr")
library("qs")
library("ggrepel")


## Plotting =====
site_colors <- c("#a6aa5c",
                 "#915cca",
                 "#89b432",
                 "#cd4698",
                 "#4db854",
                 "#d3404a",
                 "#5bb98b",
                 "#c280c2",
                 "#d2a237",
                 "#627dc8",
                 "#ce6f3e",
                 "#46b3d1",
                 "#8a6e30",
                 "#bc6173",
                 "#48793b")

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")


here()
output_dir <- here("04_assemble_bucket_folds")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("01_pbmc_nasal")

## READ IN FILES ======
nasal_only_meta <- read_csv(here(input_dir, "nasal_only_metadata_clin.csv"))
pbmc_only_meta <- read_csv(here(input_dir, "pbmc_only_metadata_clin.csv"))
nasal_complete_meta <- read_csv(here(input_dir, "nasal_complete_metadata_clin.csv"))
pbmc_complete_meta <- read_csv(here(input_dir, "pbmc_complete_metadata_clin.csv"))



## Create boolean based on trajectory group, factorize ====

#inputs: a data frame with information about participants, including their trajectory group
#outputs: the same dataframe with an additional variable; if the value is 1, the corresponding
# trajectory group value is 5. if the value is 0, the corresponding trajectory group value is not 5
# this is a factor with two levels
boolean_trajectory_factor <- function(metadata)
{
  metadata <- metadata %>%
    mutate(group = factor(ifelse(trajectory_group == 5, 1, 0)))
}

nasal_only_meta <- boolean_trajectory_factor(nasal_only_meta)
pbmc_only_meta <- boolean_trajectory_factor(pbmc_only_meta)
nasal_complete_meta <- boolean_trajectory_factor(nasal_complete_meta)
pbmc_complete_meta <- boolean_trajectory_factor(pbmc_complete_meta)


## GENERATE 5 FOLDS ===========



# Split into train and test groups, 70% and 30%

#inputs: a dataframe with information about participants
#a descriptive and distinctive string
#outputs: a dataframe with only 70% of the rows selected, as this is the training dataframe
#the rest of the input dataframe is saved using the descriptive string as a prefix
#this is the testing dataframe
test_train <- function(metadata, file_prefix)
{
  sample <- sample.int(n = nrow(metadata), size = floor(.70*nrow(metadata)), replace = F)
  train <- metadata[sample, ]
  test  <- metadata[-sample, ]
  write_csv(test, here(output_dir, paste0(file_prefix, "_test_split.csv")))
  return(train)
}

# Generate 5-fold cv
#inputs: the train dataframe with information about participants
# a descriptive and distinctive string
#outputs: the same train dataframe with a new variable which indicates which CV
# fold that participant is part of
# The number of severe COVID individuals per fold is saved using the file prefix
fold <- function(metadata, file_prefix)
{
  set.seed(720)
  print(paste0("Generating 5 folds for ", file_prefix))
  # split into train and test
  meta_train <- test_train(metadata, file_prefix)
  print("Finished splitting into train and test...")
  
  # Goal: approximate the ratio of "severe" to "non-severe" COVID-19 case
  
  min.severe <- floor(sum(meta_train$group==1)/5) # minimum number of severe COVID per fold
  print(min.severe)
  while (TRUE) {
    # Generate 5 fold
    cv.folds <- meta_train %>%
      mutate(fold=sample(rep(1:5, length.out=nrow(.))))
    
    # Count number of severe COVID individuals per fold
    cv.folds.table <- cv.folds %>%
      group_by(fold) %>%
      dplyr::count(group)
    
    if (min(cv.folds.table[cv.folds.table$group==1,"n"]) < min.severe) {
      print("At least one fold has too few severe COVID individuals. Regenerating CV folds...")
    } else {
      print("Success generating 5 folds...")
      # Export folds table
      write_csv(cv.folds.table, here(output_dir, paste0(file_prefix, "_cv.folds.table.csv")))
      print(cv.folds.table)
      return(cv.folds)
    }
  }
}

pbmc_complete_cv.folds <- fold(pbmc_complete_meta, "pbmc_complete")
pbmc_only_cv.folds <- fold(pbmc_only_meta, "pbmc_only")
nasal_complete_cv.folds <- fold(nasal_complete_meta, "nasal_complete")
nasal_only_cv.folds <- fold(nasal_only_meta, "nasal_only")

## EXPORT CV FOLDS =====
write_csv(pbmc_complete_cv.folds, here(output_dir, "pbmc_complete_train_folds.csv"))
write_csv(pbmc_only_cv.folds, here(output_dir, "pbmc_only_train_folds.csv"))
write_csv(nasal_only_cv.folds, here(output_dir, "nasal_only_train_folds.csv"))
write_csv(nasal_complete_cv.folds, here(output_dir, "nasal_complete_train_folds.csv"))

## EXPLORE SEVERITY STATUS AND ENROLLMENT SITE WITHIN EACH FOLD =====


# fold vs enrollment_site

#inputs: the train dataframe with information about participants as well as their fold
# a descriptive and distinctive string
#outputs: a visualization of the distribution of enrollment sites across each fold
# the dataframe required to build that visualization is saved using the input string as a prefix
fold_enrollment_site <- function(cv, file_prefix)
{
  # Determine the number of times each enrollment site appears per fold
  fold_enrollment_site <- cv %>%
    group_by(fold) %>%
    dplyr::count(enrollment_site)
  
  fold_enrollment_site$site <- factor(fold_enrollment_site$enrollment_site, levels = unique(fold_enrollment_site$enrollment_site))
  fold_enrollment_site$fold <- factor(fold_enrollment_site$fold, levels = unique(fold_enrollment_site$fold))
  
  #export
  write_csv(fold_enrollment_site, here(output_dir, paste0(file_prefix, "_fold_enrollment_site.csv")))
  #plot distribution
  a <- ggplot(fold_enrollment_site, aes(x = fold, y = n, fill = site)) +
    geom_bar(stat = "identity") +
    labs(title = paste0("Distribution of Enrollment Sites Across Folds ", file_prefix),
         x = "Fold",
         y = "Count",
         fill = "Enrollment Site") +
    scale_fill_manual(values = site_colors) +
    my.theme
  
  return(a)
}

pbmc_complete_enrollment <- fold_enrollment_site(pbmc_complete_cv.folds, "pbmc_ct")
pbmc_only_enrollment <- fold_enrollment_site(pbmc_only_cv.folds, "pbmc_only")
nasal_complete_enrollment <- fold_enrollment_site(nasal_complete_cv.folds, "nasal_ct")
nasal_only_enrollment <- fold_enrollment_site(nasal_only_cv.folds, "nasal_only")




