library(dplyr)
library(tidyr)
remove(list = ls()) #clean-up environment first
# Function to load data
load_data <- function(filepath) {
  read.csv(filepath)
}

# Function to create and save combined data frames per visit_num and metabolite
process_data <- function(data_frame, metabolite, output_dir) {
  print("Starting process_data function...")
  
  unique_visits <- unique(data_frame$visit_num)
  print(paste("Unique visits found:", toString(unique_visits)))
  
  for (visit in unique_visits) {
    print(paste("Processing visit:", visit))
    
    visit_data <- data_frame %>% filter(visit_num == visit)
    print(paste("Data filtered for visit. Number of rows:", nrow(visit_data)))
    
    max_samples <- max(table(visit_data$trajectory_group_new))
    print(paste("Max samples calculated:", max_samples))
    
    group_data_frames <- list()
    
    for (group in unique(visit_data$trajectory_group_new)) {
      print(paste("Processing group:", group))
      
      group_data <- visit_data %>% filter(trajectory_group_new == group)
      print(paste("Data filtered for group. Number of rows:", nrow(group_data)))
      
      transposed <- as.vector(t(group_data[[metabolite]]))
      print("Data transposed.")
      
      pad_length <- max_samples - (length(transposed) %% max_samples)
      transposed_padded <- c(transposed, rep(NA, pad_length))
      
      df <- data.frame(matrix(transposed_padded, ncol = max_samples, byrow = TRUE))
      df <- cbind(Group = group, df)  # Add group identifier as the first column
      print("Data frame created for the group.")
      
      group_data_frames[[as.character(group)]] <- df
    }
    
    combined_df <- do.call(rbind, group_data_frames)
    print("All group data frames combined.")
    
    output_filename <- file.path(output_dir, paste("combined_output_", metabolite, "_visit_", gsub(" ", "_", visit), ".csv", sep = ""))
    
    print(paste("Output filename set:", output_filename))
    
    write.table(combined_df, file = output_filename, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE, na = "")
    print(paste("File written:", output_filename))
  }
}

# Specify the input file, output directory, and metabolite
input_file <- "/Users/boryanapetrova/Dropbox/NaamaLab/projects/IMPACC_study/redone_analysis_post_locked_data/global_for_prism_MTHFR_GG.csv"
output_directory <- "/Users/boryanapetrova/Dropbox/NaamaLab/projects/IMPACC_study/redone_analysis_post_locked_data"
metabolite <- "cysteine"  # Replace with the actual compound name (change any "-" to "." as they get converted when R reads them)

# Load data
data_frame <- load_data(input_file)
print("Data loaded.")

# Process and save data
process_data(data_frame, metabolite, output_directory)
