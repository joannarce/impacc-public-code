library(ggplot2)
library(dplyr)
library(tidyr)

# Replace "path_to_your_file.csv" with the actual path to your CSV file
file_path <- "/Users/boryanapetrova/Dropbox/NaamaLab/projects/IMPACC_study/redone_analysis_post_locked_data/visit_counts_gpt_grouped_redone.csv"

# Read the CSV file
data <- read.csv(file_path)

# Transform data from wide to long format, if not already done
data_long <- data %>%
  pivot_longer(cols = -c(Trajectory, Group), names_to = "Visit", values_to = "Count")

# Calculate proportions and total counts
data_long <- data_long %>%
  group_by(Trajectory, Visit) %>%
  mutate(TotalCount = sum(Count),
         Proportion = Count / TotalCount) %>%
  ungroup()

# Check for any potential issues with the data preparation
print(head(data_long))

# Calculate the area for each group based on its proportion
data_long <- data_long %>%
  mutate(Area = (Proportion * TotalCount), # Calculate area as sqrt to adjust for circle area representation
         Size = Area ) # The 'Size' here represents the scaled area for visual representation

# Assuming the data preparation is correct, let's adjust the plotting code
ggplot(data_long, aes(x = Visit, y = Trajectory, color = Group, size = Size)) +
  geom_point(alpha = 0.5) + # Adjust transparency if needed
  scale_size(range = c(1, 20)) + # Adjust this range as needed to fit your visual preference
  theme_minimal() +
  labs(title = "Proportional Area Chart for Visits", x = "Visit", y = "Trajectory", size = "Count") +
  theme(legend.position = "right",
        legend.title.align = 0.5,
        plot.title = element_text(size = 10),
        axis.title = element_text(size =8),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 6))) # Adjust legend appearance

# This adjusts the original plotting code to ensure that size is computed within the aes() function of geom_point