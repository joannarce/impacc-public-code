# Here I prepare the OLINK data for analysis

## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("janitor")
library("readr")
library("qs")
library("ggrepel")
library("reshape2")
library("ggpubr")

## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")


here()
output_dir <- here("OLINK_00_setup")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("05d_nasal_complete_LASSO_deg_filter")





## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## READ IN FILES ======
pbmc.meta <- read_csv("01_pbmc_nasal/pbmc_complete_metadata_clin.csv")
olink.counts <- read_csv("../../../data/serum-olink/legacy/2023-01-06/Olink-Counts.csv")
olink.meta <- read_csv("../../../data/serum-olink/legacy/2023-01-06/Olink-Metadata.csv")
olink.rowfeature <- read_csv("../../../data/serum-olink/legacy/2023-01-06/Olink-RowFeature.csv")
clinical_meta <- read_csv("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-individ-locked.csv")
sample_path <- file.path("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-sample-locked.csv")
sample_raw <- read_csv(file = sample_path)
event_raw <-  read_csv("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-event-locked.csv")

## Perform filtering for IMPACC pbmc patients and event date =====
pbmc.patient <- pbmc.meta %>% 
  dplyr::pull(participant_id)

olink.samples <- sample_raw %>%
  dplyr::filter(sample_type == "OLINK") %>%
  dplyr::filter(participant_id %in% pbmc.patient) %>%
  inner_join(., event_raw, by = "event_id") %>%
  dplyr::filter(event_date <= 2) %>%
  filter(event_type != "Escalation 2") %>%
  filter(event_type != "Escalation 1") %>%
  group_by(participant_id.x) %>%
  slice_min(event_date, with_ties = FALSE) %>%
  ungroup() 

olink.meta <- olink.meta %>%
  dplyr::filter(sample_id %in% olink.samples$sample_id) %>%
  dplyr::filter(sample_status == "IMPACC sample assayed and passed QC")

olink.samples <- olink.samples %>%
  dplyr::filter(sample_id %in% olink.meta$sample_id) %>%
  dplyr::select(c(sample_id, participant_id.x)) %>%
  dplyr::rename(participant_id = participant_id.x)

## Filter counts df ==========

olink.counts <- olink.counts %>%
  dplyr::filter(sample_id %in% olink.samples$sample_id)

pbmc.meta <- pbmc.meta %>%
  dplyr::select(-c(sample_id))

olink.complete.metadata <- olink.samples %>%
  dplyr::inner_join(., pbmc.meta, by = "participant_id")


## Train-test split =======
set.seed(19)

olink.complete.metadata <- olink.complete.metadata %>%
  mutate(trajectory = ifelse(trajectory_group == 5, 1, 0))

trainIndex <- createDataPartition(olink.complete.metadata$trajectory, p = .7, 
                                  list = FALSE, 
                                  times = 1)


train.counts <- olink.counts[trainIndex, ]
test.counts <- olink.counts[-trainIndex, ]

train.meta <- olink.complete.metadata[trainIndex, ]
test.meta <- olink.complete.metadata[-trainIndex, ]

write.csv(x = train.meta, here(output_dir, "train_meta.csv"))
write.csv(x = test.meta, here(output_dir, "test_meta.csv"))
write.csv(x = train.counts, here(output_dir, "train_counts.csv"))
write.csv(x = test.counts, here(output_dir, "test_counts.csv"))
write.csv(x = olink.complete.metadata, here(output_dir, "olink_complete_metadata.csv"))
write.csv(x = olink.counts, here(output_dir, "olink_counts.csv"))
