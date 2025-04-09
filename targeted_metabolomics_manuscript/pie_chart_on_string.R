#install.packages("gridExtra")

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Replace "path_to_your_file.csv" with the actual path to your CSV file
file_path <- "/Users/boryanapetrova/Dropbox/NaamaLab/projects/IMPACC_study/Global_metabo_data/#secure# Re_ Methionine analysis updates 2/visit_counts_gpt_grouped.csv"

# Read the CSV file
data <- read.csv(file_path)

# Assuming the structure is:
# Trajectory (in rows), Group (Global, Targeted), and Visits (in columns)
# We need to reshape the data to a long format
data_long <- data %>%
  gather(key = "Visit", value = "Count", -Trajectory, -Group)

# Calculate proportions
data_long <- data_long %>%
  group_by(Trajectory, Visit) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Modified function to generate a pie chart with adjusted sizes
plot_pie <- function(trajectory, visit) {
  filtered_data <- data_long %>%
    filter(Trajectory == trajectory, Visit == visit)
  
  ggplot(filtered_data, aes(x = "", y = Proportion, fill = Group)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(legend.position = "none", # Remove the legend
          # theme(legend.title = element_text(size = 6), # Adjust legend title size
          # legend.text = element_text(size = 6), # Adjust legend text size
          # legend.key.size = unit(0.5, "cm"), # Adjust legend key size
          plot.title = element_text(size = 8), # Adjust plot title size
          strip.text = element_text(size = 8)) + # Adjust facet strip text size, if used
    labs(title = paste("Trajectory:", trajectory, "Visit:", visit), fill = "Group")
}

# Generate and arrange plots as before
plots <- list()
for (trajectory in unique(data_long$Trajectory)) {
  for (visit in unique(data_long$Visit)) {
    plots[[paste(trajectory, visit, sep = "_")]] <- plot_pie(trajectory, visit)
  }
}

# Arrange plots in a grid with adjusted layout
do.call(gridExtra::grid.arrange, c(plots, ncol = length(unique(data_long$Visit))))