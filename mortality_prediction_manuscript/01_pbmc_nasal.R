## Here I 
# 1) Determine the samples for pbmc and nasal transcriptomics which meet the QC thresholds we have established
# 2) Create a combined dataframe which includes the pbmc + clinical metadata
# 3) Create a combined dataframe which includes the nasal + clinical metadata
# 4) Create a combined dataframe which includes the CT metadata, rpM metadata, values for both, pbmc metadata, clinical metadata
# 5) Create a combined dataframe which includes the CT metadata, rpM metadata, values for both, nasal metadata, clinical metadata
## Load packages for data wrangling, manipulation, saving ======

library("tidyverse")
library("here")
library("qs")
library("janitor")
library("readr")

## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")
here()
output_dir <- here("01_pbmc_nasal")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("relationship_CT_rpM")

## READ IN SAMPLE+EVENT METADATA, CT+RPM INFO, PBMC AND NASAL METADATA ====

sample_path <- file.path("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-sample-locked.csv")
sample_raw <- read_csv(file = sample_path)
event_raw <-  read_csv("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-event-locked.csv")
rpM_ct_complete <- read_csv(here(input_dir, "rpM_ct_complete.csv"))
pbmc_meta <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-Metadata.csv")
nasal_meta <- read_csv("../../../data/nasal-transcriptomics/legacy/2022-10-19/nasal-transcriptomics-Metadata.csv")

## READ IN PATIENT METADATA
clinical_meta <- read_csv("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-individ-locked.csv")
## SELECT COLUMNS OF INTEREST ========

cols_to_select <- c("participant_id", "participant_type", "enrollment_site", "admit_age", "trajectory_group")

clinical_meta <- clinical_meta %>% 
  select(.,all_of(cols_to_select)) %>%
  filter(participant_type == "COVID-19 Positive") %>%
  select(-c("participant_type")) ## exclude controls

## SELECT ONLY SAMPLES THAT HAVE BEEN ASSAYED AND QCd =====

filter_qc_passed <- function(df) {
  df %>%
    filter(sample_status == "IMPACC sample assayed and passed QC")
}

# Filter each data frame using the function
pbmc_meta <- filter_qc_passed(pbmc_meta)
nasal_meta <- filter_qc_passed(nasal_meta)

# Filter to remove controls

pbmc_meta <- pbmc_meta %>%
  filter(Internal_or_IMPACC_Control == "NO")

nasal_meta <- nasal_meta %>%
  filter(!(internal_control == "YES")) %>%
  filter(!(internal_healthy_control == "YES"))

## DETERMINE AVAILABLE SAMPLE_IDS, THEIR PATIENT_IDS, AND EVENT_DATES/TYPES ======

# Inputs: metadata dfs for samples of interest
# Output: metadata dfs for samples of interest with the 
# columns listed in select. Only the rows that have values less than or equal to 2 in 
# event date are selected (only want samples that have been collected before threshold)
# Only want samples collected in hospitals.
process_dataframe <- function(df) {
  df %>%
    inner_join(sample_raw, by = "sample_id") %>%
    inner_join(event_raw, by = "event_id") %>%
    select(sample_id, participant_id, event_id, sample_type, event_type, event_date, samples_collected, event_location) %>%
    filter(event_date <= 2) %>%
    filter(event_location == "Hospital")  %>%  
    filter(event_type != "Escalation 2") %>%
    filter(event_type != "Escalation 1") 
}

# Process the dataframes

pbmc_processed <- process_dataframe(pbmc_meta)
nasal_processed <- process_dataframe(nasal_meta)

## REMOVE DUPLICATE PARTICIPANT ID ====

nasal_processed <- nasal_processed %>%
  filter(event_type == "Visit 1") 

pbmc_processed <- pbmc_processed %>%
  filter((participant_id == "005-0012" | event_type == "Visit 1")) ## extra condition
## for participant 005-0012, as they only have visit 2 data

duplicates <- nasal_processed %>%
  group_by(participant_id) %>%
  filter(n() > 1) %>% 
  ungroup()

nasal_processed <- nasal_processed %>%
  filter(!(sample_id == "0865-00DGJG00-001")) # Visit 6

## MERGE WITH CT AND RPM ====
nasal_complete_metadata <- rpM_ct_complete %>%
  inner_join(., nasal_processed, by = "participant_id") %>%
  select(sample_id.ct, sample_id.rpm, sample_id = sample_id, participant_id,
         event_id = event_id.x, event_date = event_date.x, event_type = event_type.x)

pbmc_complete_metadata <- rpM_ct_complete %>%
  inner_join(., pbmc_processed, by = "participant_id") %>%
  select(sample_id.ct, sample_id.rpm, sample_id = sample_id, participant_id,
         event_id = event_id.x, event_date = event_date.x, event_type = event_type.x)

## MERGE WITH CLINICAL METADATA ====
nasal_complete_metadata_clin <- nasal_complete_metadata %>%
  inner_join(., clinical_meta, by = "participant_id")

pbmc_complete_metadata_clin <- pbmc_complete_metadata %>%
  inner_join(., clinical_meta, by = "participant_id")

pbmc_only_metadata <- pbmc_processed %>%
  inner_join(., clinical_meta, by = "participant_id")

nasal_only_metadata <- nasal_processed %>%
  inner_join(., clinical_meta, by = "participant_id")

rpM_ct_metadata_clin <- rpM_ct_complete %>%
  inner_join(., clinical_meta, by = "participant_id")

## EXPORT ====
write_csv(nasal_complete_metadata_clin, here(output_dir, "nasal_complete_metadata_clin.csv"))
write_csv(pbmc_complete_metadata_clin, here(output_dir, "pbmc_complete_metadata_clin.csv"))
write_csv(nasal_only_metadata, here(output_dir, "nasal_only_metadata_clin.csv"))
write_csv(pbmc_only_metadata, here(output_dir, "pbmc_only_metadata_clin.csv"))
write_csv(rpM_ct_metadata_clin, here(output_dir, "rpM_ct_metadata_clin.csv"))

