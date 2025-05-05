# Here, I explore the metadata for the five sample bins

## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("janitor")
library("readr")
library("gridExtra")

## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

bin_colors <- c("#bc6a84",
                "#59aa54",
                "#9762ca",
                "#999a3e",
                "#ce51a0",
                "#45b0a4",
                "#d2485a",
                "#6b85cd",
                "#c95631",
                "#c88744")
nasal_complete_color <- bin_colors[1]
nasal_only_color <- bin_colors[3]
pbmc_complete_color <- bin_colors[5]
pbmc_only_color <- bin_colors[7]
rpM_ct_only_color <- bin_colors[9]

## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")
here()
output_dir <- here("02_bucket_metadata_exploration")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("01_pbmc_nasal")

## READ IN FILES ======
pbmc_complete_metadata <- read_csv(here(input_dir, "pbmc_complete_metadata_clin.csv"))
nasal_complete_metadata <- read_csv(here(input_dir, "nasal_complete_metadata_clin.csv"))
nasal_only_metadata <- read_csv(here(input_dir, "nasal_only_metadata_clin.csv"))
pbmc_only_metadata <- read_csv(here(input_dir, "pbmc_only_metadata_clin.csv"))
rpM_ct_only_metadata <- read_csv(here(input_dir, "rpM_ct_metadata_clin.csv"))

## VISUALIZE METADATA =====

#inputs:
# dataset is a dataframe consisting of the metadata for each participant
# fill_color is a string representing the color that most graphs in the output should be
# output_prefix is a descriptive and distinct string
#outputs:
# combined_plots is a gtable object containing 6 plots based on metadata variables
# ks_test_result is a list containing the output of a Kolmogorov–Smirnov test of 
# the difference in age distributions between trajectory groups 5 and 1
process_and_visualize <- function(dataset, fill_color, output_prefix) {
  # Summarize variables
  p_trajectory_counts <- dataset %>%
    group_by(trajectory_group) %>%
    summarize(participant_count = n_distinct(participant_id), .groups = 'drop')
  
  p_site_counts <- dataset %>%
    group_by(enrollment_site) %>%
    summarize(site_count = n_distinct(participant_id), .groups = 'drop')
  
  # Visualizations of ungrouped variables
  p_hist_age <- ggplot(dataset, aes(x = admit_age)) +
    geom_histogram(binwidth = 5, color = "black", fill = fill_color) +
    labs(title = "Histogram of Admission Age", x = "Age", y = "Frequency") +
    xlim(0, NA) +
    my.theme
  
  p_bar_site <- ggplot(dataset, aes(x = enrollment_site)) +
    geom_bar(color = "black", fill = fill_color) +
    labs(title = "Bar Plot of Enrollment Site", x = "Enrollment Site", y = "Count") +
    my.theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  p_bar_trajectory <- ggplot(dataset, aes(x = factor(trajectory_group))) +
    geom_bar(color = "black", fill = fill_color) +
    labs(title = "Bar Plot of Trajectory Group", x = "Trajectory Group", y = "Count") +
    my.theme
  
  # Visualizations of grouped variables
  p_violin_trajectory <- ggplot(dataset, aes(x = factor(trajectory_group), y = admit_age)) +
    geom_violin(fill = fill_color, color = "black") +
    labs(title = "Violin Plot of Admission Age by Trajectory Group", 
         x = "Trajectory Group", 
         y = "Admission Age") +
    geom_jitter(width = 0.2, alpha = 0.5) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    my.theme
  
  p_violin_age_site <- ggplot(dataset, aes(x = enrollment_site, y = admit_age)) +
    geom_violin(fill = fill_color, color = "black") +
    labs(title = "Violin Plot of Admission Age by Enrollment Site", 
         x = "Enrollment Site", 
         y = "Admission Age") +
    my.theme + 
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  p_violin_trajectory_site <- ggplot(dataset, aes(x = enrollment_site, y = trajectory_group)) +
    geom_violin(fill = fill_color, color = "black") +
    labs(title = "Violin Plot of Trajectory Group by Enrollment Site", 
         x = "Enrollment Site", 
         y = "Trajectory Group") +
    my.theme + 
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Age and trajectory group ks test
  group_1 <- dataset %>% filter(trajectory_group == 1) %>% pull(admit_age)
  group_5 <- dataset %>% filter(trajectory_group == 5) %>% pull(admit_age)
  ks_test_result <- ks.test(group_1, group_5)
  
  # Export summary tables
  write_csv(p_trajectory_counts, here(output_dir, paste0(output_prefix, "_trajectory_counts.csv")))
  write_csv(p_site_counts, here(output_dir, paste0(output_prefix, "_site_counts.csv")))
  
  # Export combined plots
  p_combined_plots <- grid.arrange(p_hist_age, p_bar_site, p_bar_trajectory, p_violin_age_site, 
                                   p_violin_trajectory, p_violin_trajectory_site, ncol = 2)
  ggsave(here(output_dir, paste0(output_prefix, "_combined_plots.png")), plot = p_combined_plots, width = 16, height = 12)
  
  # Export ks test
  qsave(ks_test_result, here(output_dir, paste0(output_prefix, "_ks_test.qs")))
  
  return(list(combined_plots = p_combined_plots, ks_test_result = ks_test_result))
}

pbmc_complete <- process_and_visualize(pbmc_complete_metadata, pbmc_complete_color, "pbmc_complete")
nasal_complete <- process_and_visualize(nasal_complete_metadata, nasal_complete_color, "nasal_complete")
pbmc_only <- process_and_visualize(pbmc_only_metadata, pbmc_only_color, "pbmc_only")
nasal_only <- process_and_visualize(nasal_only_metadata, nasal_only_color, "nasal_only")
rpM_ct_only <- process_and_visualize(rpM_ct_only_metadata, rpM_ct_only_color, "rpM_ct_only")

## EXPLORE CT IN DEPTH ====

violin_ct <- ggplot(rpM_ct_only_metadata, aes(x = factor(trajectory_group), y = ct)) +
  geom_violin(fill = rpM_ct_only_color, color = "black") +
  labs(title = "Violin Plot of CT value by Trajectory Group", 
       x = "Trajectory Group", 
       y = "CT value") +
  geom_jitter(width = 0.2, alpha = 0.5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  my.theme

age_ct <- ggplot(rpM_ct_only_metadata, aes(x = admit_age, y = ct)) +
  geom_point() +
  labs(title = "Scatter Plot of CT vs. Admit Age", 
       x = "Admit Age", 
       y = "CT Value") +
  my.theme

not_5 <- rpM_ct_only_metadata %>% filter(trajectory_group != 5) %>% pull(ct)
group_5 <- rpM_ct_only_metadata %>% filter(trajectory_group == 5) %>% pull(ct)
ks.test(not_5, group_5)
