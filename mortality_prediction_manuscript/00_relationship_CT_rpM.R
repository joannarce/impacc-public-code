## Here I 
# 1) Determine the samples for nasal rpM and Ct which meet the threshold of event_date <=2 and other filters
# 2) Plot a regression to determine the relationship between rpM and Ct
# 3) Create a combined dataframe which holds the Ct value for each patient, NAs are imputed based on regression

## Load packages for data wrangling, manipulation, saving ======

library("tidyverse")
library("here")
library("qs")
library("janitor")
library("readr")

## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

library(ggplot2)
library(robustbase)

## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")
here()
output_dir <- here("relationship_CT_rpM")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

## READ IN rpM and CT metadata =====
ct_meta <- read_csv("../../../data/nasal-viralload/legacy/2022-01-01/nasal-viralload-Metadata.csv")
rpM_meta <- read_csv("../../../data/nasal-metagenomics/legacy/2022-09-03/nasal-metagenomics-Metadata.csv")

## READ IN SAMPLE METADATA ============

sample_path <- file.path("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-sample-locked.csv")

sample_raw <- read_csv(file = sample_path)

## SELECT SAMPLES OF INTEREST =========

samples_of_interest <- c("NAS SWAB", "BLK RNA-SEQ")


## only select nasal swab and bulk RNA seq (pbmc) samples

sample_merged <- sample_raw %>%
  filter(sample_type %in% samples_of_interest)

## READ IN EVENT METADATA ============

event_raw <-  read_csv("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-event-locked.csv")


sample_merged <- sample_merged %>%
  inner_join(., event_raw, by = "event_id") %>%
  filter(event_date <= 2) %>%
  filter(samples_collected == TRUE) ## only select samples that were collected on or before event date 2

## SELECT ONLY SAMPLES THAT HAVE BEEN ASSAYED AND QCd =====

filter_qc_passed <- function(df) {
  df %>%
    filter(sample_status == "IMPACC sample assayed and passed QC")
}

ct_meta <- filter_qc_passed(ct_meta)
rpM_meta <- filter_qc_passed(rpM_meta)


# Filter to remove controls
rpM_meta <- rpM_meta %>%
  filter(internal_control == "no")

## DETERMINE AVAILABLE SAMPLE_IDS, THEIR PARTICIPANT_IDS, AND EVENT_DATES/TYPES ======

# Inputs: metadata dfs for samples of interest
# Output: metadata dfs for samples of interest with the 
# columns listed in select. Only the rows that have values less than or equal to 2 in 
# event date are selected (only want samples that have been collected before threshold)
# Only want samples collected in hospitals.


process_dataframe <- function(df) {
  df %>%
    inner_join(sample_raw, by = "sample_id") %>%
    inner_join(event_raw, by = "event_id") %>%
    filter(event_date <= 2) %>%
    filter(event_location == "Hospital")  %>%  
    filter(event_type != "Escalation 2") %>%
    filter(event_type != "Escalation 1") %>%
    select(sample_id, participant_id = participant_id.x, event_id, sample_type, event_type, event_date, samples_collected, event_location)
}

rpM_processed <- process_dataframe(rpM_meta)
ct_processed <- process_dataframe(ct_meta)


## COMBINE DATAFRAMES ======
rpM_ct_metadata <- ct_processed %>%
  inner_join(., rpM_processed, by = "event_id", suffix = c(".ct", ".rpm"))

# Check to make sure samples are not being omitted
event_ids_in_ct_not_in_rpM <- setdiff(ct_processed$event_id, rpM_processed$event_id)
event_ids_in_rpM_not_in_ct <- setdiff(rpM_processed$event_id, ct_processed$event_id)
participant_id_in_rpM_not_in_ct <- setdiff(rpM_processed$participant_id, ct_processed$participant_id)
participant_id_in_ct_not_in_rpM <- setdiff(ct_processed$participant_id, rpM_processed$participant_id)

## READ IN SAMPLE DATA =====
ct <- read_csv("../../../data/nasal-viralload/legacy/2022-01-01/nasal-viralload-Counts.csv")
rpM <- read_csv("../../../data/nasal-metagenomics/legacy/2022-09-03/nasal-metagenomics-Statistics.csv")

ct <- ct %>%
  select(sample_id, N1_CT)

rpM <- rpM %>%
  select(sample_id, sc2_rpm)

## Remove samples that have mislabelled event date and that failed nasal transcriptomics qc
rpM <- rpM %>%
  filter(!(sample_id == "0865-00F2MG00-001")) %>% # failed nasal transcriptomics
  filter(!(sample_id == "0865-00DGJG00-001")) # Visit 6


## MERGE SAMPLE DATA WITH METADATA ====
combined_rpM_ct <- rpM_ct_metadata %>%
  inner_join(., rpM, by = c("sample_id.rpm" = "sample_id")) %>%
  inner_join(ct, by = c("sample_id.ct" = "sample_id"))

## PLOT AND REGRESSION ====
ggplot(combined_rpM_ct, aes(x = sc2_rpm, y = N1_CT)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(x = "SC2 RPM", y = "N1 CT", title = "SC2 RPM vs N1 CT") +
  my.theme

# Fit the linear model
model <- lm(N1_CT ~ sc2_rpm, data = combined_rpM_ct)
# Get the model summary
model_summary <- summary(model)
# Extract coefficients and R^2
coefficients <- coef(model)
r_squared <- model_summary$r.squared

combined_rpM_ct <- combined_rpM_ct %>%
  mutate(log_sc2_rpm = log1p(sc2_rpm))  

# Plot with log-transformed RPM
ggplot(combined_rpM_ct, aes(x = log_sc2_rpm, y = N1_CT)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, col = "#7570B3") +
  labs(x = "log(SC2 RPM+1)", y = "N1 CT", title = "Log SC2 RPM vs N1 CT") +
  my.theme

# Explore transformation of data

# Fit the linear model
model_transform <- lm(N1_CT ~ log_sc2_rpm, data = combined_rpM_ct)
# Get the model summary
model_summary_transform <- summary(model_transform)
# Extract coefficients and R^2
coefficients_transform <- coef(model_transform)
r_squared_transform <- model_summary_transform$r.squared

# Explore robust regression on untransformed data

model_robust <- lmrob(N1_CT ~ sc2_rpm, data = combined_rpM_ct)
robust_summary <- summary(model_robust)

# Explore robust regression on transformed data
model_robust_transform <- lmrob(N1_CT ~ log_sc2_rpm, data = combined_rpM_ct)
robust_summary_transform <- summary(model_robust_transform)

g <- ggplot(combined_rpM_ct, aes(x = log_sc2_rpm, y = N1_CT)) +
  geom_point() +
  geom_smooth(method = "lmrob", formula = y ~ x, color = "#6b85cd", se = TRUE) +
  labs(x = "log(SARS-CoV-2 RPM + 1)", y = "Cycle threshold value", title = "SARS-CoV-2 viral load measurement imputation") +
  my.theme

ggplot(combined_rpM_ct, aes(x = log_sc2_rpm, y = N1_CT)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, color = "#7570B3", se = FALSE, linetype = "dashed") + # Linear regression
  geom_smooth(method = "lmrob", formula = y ~ x, color = "indianred2", se = FALSE) + # Robust regression
  labs(x = "Log(SC2 RPM+1)", y = "N1 CT", title = "Robust vs. Linear Regression on Log-Transformed Data") +
  my.theme


## IMPUTATION OF MISSING VALUES ====
rpM_ct_impute <- rpM_processed %>%
  left_join(., ct_processed, by = "event_id", suffix = c(".rpm", ".ct")) %>%
  inner_join(rpM, by = c("sample_id.rpm" = "sample_id")) %>%
  left_join(ct, by = c("sample_id.ct" = "sample_id"))

missing_before <- sum(is.na(rpM_ct_impute$N1_CT))

# Log transform rpm data

rpM_ct_impute <- rpM_ct_impute %>%
  mutate(
    log_sc2_rpm = log1p(sc2_rpm)
  )

rpM_ct_impute <- rpM_ct_impute %>%
  mutate(
    N1_CT_imputed = ifelse(is.na(N1_CT),
                           predict(model_robust_transform, newdata = .),
                           N1_CT)
  )

missing_after <- sum(is.na(rpM_ct_impute$N1_CT_imputed))


## REMOVE DUPLICATE PARTICIPANT ID ====
duplicates <- rpM_ct_impute %>%
  group_by(participant_id.rpm) %>%
  filter(n() > 1) %>% 
  ungroup()

rpM_ct_metadata_imputed <- rpM_ct_impute %>%
  group_by(participant_id.rpm) %>%
  filter(event_type.rpm == "Visit 1") %>%
  ungroup()

# Clean up the dataframe
rpM_ct_metadata_imputed_final <- rpM_ct_metadata_imputed %>%
  select(sample_id.rpm, sample_id.ct, participant_id = participant_id.rpm, event_id,
         event_date = event_date.rpm, event_type = event_type.rpm, ct = N1_CT_imputed)

## EXPORT ====

write_csv(rpM_ct_metadata_imputed_final, here(output_dir, "rpM_ct_complete.csv"))
qsave(robust_summary_transform, here(output_dir, "robust_summary_logtransform.qs"))
write_csv(duplicates, here(output_dir, "duplicate_participants.csv"))
write_csv(sample_merged, here(output_dir, "sample_event_merged.csv"))
