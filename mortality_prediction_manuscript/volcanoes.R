# Here I make figures 2a and 2c

## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("janitor")
library("readr")
library("qs")
library("ggrepel")
library("ggplot2")
library("reshape2")

## Load packages for analysis ======

library("DESeq2")
library("glmnet")
library("edgeR")
library("caret")
library(pROC)

## Load data ======

top_table_for_volcano <- read_csv("03_bucket_DEG/pbmc_ct_age_sex_age_sex_top_table.csv")

## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )


## Plot with cutoffs ====

pretty_plot <- ggplot(top_table_for_volcano, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = diffexpressed), size = 0.8) +
  scale_color_manual(name = "", values = c("grey", "#6b85cd", "#d2485a"),
                     labels = c("Not significant", "Downregulated", "Upregulated"))  +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression(-log[10](P[adj]))) + 
  xlim(c(-2.5, 2.5)) +
  my.theme + 
  theme(legend.position = "none") 

ggsave("volcano_plot.svg", pretty_plot, width = 4.5, height = 5)

## Load nasal data ======
nasal_top_table_for_volcano <- read_csv("03_bucket_DEG/nasal_ct_age_sex_age_sex_top_table.csv")

nasal_pretty_plot <- ggplot(nasal_top_table_for_volcano, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = diffexpressed), size = 0.8, show.legend = FALSE) +
  scale_color_manual(name = "", values = c("grey", "#6b85cd", "#d2485a"),
                     labels = c("Not significant", "Downregulated", "Upregulated")) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  theme(legend.position = NA) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression(-log[10](P[adj]))) + 
  xlim(c(-2.5, 2.5)) +
  my.theme

ggsave("nasal_volcano_plot.svg", nasal_pretty_plot, width = 4.5, height = 5)

# Construct figure 2B ======

combined_top_table <- top_table_for_volcano %>% 
  inner_join(., nasal_top_table_for_volcano, by = "gene_id", suffix = c(".pbmc", ".nasal")) %>%
  mutate(dual_significance = (diffexpressed.nasal != "Not significant" & diffexpressed.pbmc != "Not significant")) %>%
  dplyr::filter(gene_biotype.pbmc == "protein_coding") %>%
  dplyr::mutate(delabel = ifelse(dual_significance, gene_name.pbmc, NA))

set.seed(19)
ggplot(combined_top_table, aes(x = logFC.pbmc, y = logFC.nasal, label = delabel)) +
  geom_point(aes(color = dual_significance), size = 0.8, show.legend = FALSE) +
  scale_color_manual(name = "", values = c("grey", "red")) + 
  geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
  xlab("PBMC logFC") +
  ylab("NS logFC") +
  geom_label_repel(box.padding = 1) +
  my.theme

filtered_top_table <- combined_top_table %>%
  dplyr::filter(diffexpressed.pbmc != "Not significant" | diffexpressed.nasal != "Not significant")


loglog <- ggplot(filtered_top_table, aes(x = logFC.pbmc, y = logFC.nasal, label = delabel)) +
  geom_point(data = subset(filtered_top_table, dual_significance == FALSE),
             aes(color = "grey"), size = 0.8, show.legend = FALSE) +
  geom_point(data = subset(filtered_top_table, dual_significance == TRUE),
             aes(color = "red"), size = 0.8, show.legend = FALSE) +
  scale_color_manual(name = "", values = c("grey", "red")) + 
  geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
  xlab(expression("Blood log"[2]*"FC")) +
  ylab(expression("Nasal log"[2]*"FC")) +
  geom_label_repel(box.padding = 1.5) +
  xlim(c(-2.3, 2.3)) + 
  my.theme


ggsave("2cloglog.svg", loglog, width = 4.5, height = 5)
