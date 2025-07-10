## Boxplots

## Here, I make boxplots comparing age and CT value with trajectory

library("ggpubr")
library("ggplot2")


metadata <- read_csv("01_pbmc_nasal/pbmc_complete_metadata_clin.csv")
rpM_ct <- read_csv("01_pbmc_nasal/rpM_ct_metadata_clin.csv")

metadata <- metadata %>%
  mutate(traj = ifelse(trajectory_group == 5, 1, 0))

rpM_ct <- rpM_ct %>%
  mutate(traj = ifelse(trajectory_group == 5, 1, 0)) #%>%
# filter(participant_id %in% metadata$participant_id)


## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )


## age plot ======

age <- ggplot(metadata, aes(x = factor(traj), y = admit_age, fill = factor(traj))) +
  geom_boxplot(position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("#6b85cd", "#d2485a")) +
  scale_x_discrete(labels = c("0" = "No mortality", "1" = "Mortality")) + 
  labs(x = "Trajectory", y = "Age", fill = "") + my.theme +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("0", "1")), 
                     bracket.size = 0.3,  # Adjusts the bar thickness
                     label.y = max(metadata$admit_age) * 1.05) +  # Adjusts the height of the bar
  guides(fill = "none")

ggsave("age_boxplot.svg", age, width = 4, height = 4)

## CT value plot

CT_box <- ggplot(rpM_ct, aes(x = factor(traj), y = ct, fill = factor(traj))) +
  geom_boxplot(position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("#6b85cd", "#d2485a")) +
  scale_x_discrete(labels = c("0" = "No mortality", "1" = "Mortality")) + 
  labs(x = "Trajectory", y = "CT value", fill = "") + my.theme +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("0", "1")), 
                     bracket.size = 0.3,  # Adjusts the bar thickness
                     label.y = max(rpM_ct$ct) * 1.05) +  # Adjusts the height of the bar
  guides(fill = "none")

ggsave("CT_boxplot.svg", CT_box, width = 4, height = 4)
