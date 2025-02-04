##### Codebase for the data analysis of convalescent data
##### Authors: Ravi K. Patel, Leying Guan, Jeremy Gygi, Gisela Gabernet
#
#.  This codebase contains functions used for the following:
#. 1) Initial data loading (via load_IMPACC_datasets(...))
#. 2) Misc. functions for color schemes / plotting
#. 3) Misc. functions for data pre-processing

#### Use this function to load all available datasets and return an R environment containing all data objects.
load_IMPACC_datasets <- function(DATA_VERSION, 
                                 PHASES = c(1,2,3,4,5,9,NA), # Phases to include (phase 4a = 4, phase 4b = 5)
                                 ALLOWED_SMPL_STATUS = # Allow samples with the following QC:
                                   c("IMPACC sample assayed and passed QC", 
                                     "IMPACC sample assayed with questionable QC"),
                                 KEEP_COVID19_POS = FALSE, 
                                 KEEP_EMORY_CONTROLS = FALSE, # event_type = "Control donor visit"
                                 KEEP_ESCALATIONS = FALSE, # event_type == "Escalation"
                                 FILTER_BY_CORE_ASSAY_COHORT = FALSE, # Use samples from core assay cohort
                                 FILTER_BY_INTEGRATION_COHORT = FALSE, # Use samples from integration cohort
                                 SELECT_PARTICIPANT_ID = NULL, # Select specific participants
                                 EXCLUDE_SAMPLES = NULL, # Exclude specific samples
                                 EVENT_DATE_UPPER = 999, # Upper date limit, all samples
                                 VISIT_LOWER = 1, # Lower visit limit, all samples
                                 VISIT_UPPER = 999, # Upper visit limit, all samples
                                 SELECT_VISITS = NULL, # Select specific visits to load
                                 USE_LOCKED_CLINICAL = TRUE # Use locked clinical file (from 01/01/2023)
) {
  
  #### Load packages necessary for data loading and processing purpose
  load_packages(packages)
  
  #### Define Data directories: Use the absolute path of the directory containing "current" and "legacy" directories.
  data_dirs <- c(
    bld_cytof_dir = "/data/bld-cytof",
    bld_gwas_dir = "/data/bld-gwas",
    clinical_dir = "/data/clinical",
    ea_cytof_dir = "/data/ea-cytof",
    ea_metagenomics_dir = "/data/ea-metagenomics",
    ea_transcriptomics_dir = "/data/ea-transcriptomics",
    plasma_metabolomics_global_dir = "/data/metabolomics/plasma-metabolomics-global",
    plasma_metabolomics_targeted_dir = "/data/metabolomics/plasma-metabolomics-targeted",
    serum_metabolomics_global_dir = "/data/metabolomics/serum-metabolomics-global",
    nasal_metagenomics_dir = "/data/nasal-metagenomics",
    nasal_transcriptomics_dir = "/data/nasal-transcriptomics",
    nasal_viralload_dir = "/data/nasal-viralload",
    nasal_viralseq_dir = "/data/nasal-viralseq",
    pbmc_transcriptomics_dir = "/data/pbmc-transcriptomics",
    plasma_proteomics_targeted_dir = "/data/proteomics/plasma-proteomics-targeted",
    plasma_proteomics_global_dda_dir = "/data/proteomics/plasma-proteomics-global-DDA",
    plasma_proteomics_global_dia_dir = "/data/proteomics/plasma-proteomics-global-DIA",
    serum_autoantibody_dir = "/data/serum-autoantibody",
    serum_proteomics_global_dir = "/data/proteomics/serum-proteomics-global",
    serum_olink_dir = "/data/serum-olink",
    serum_rbd_abtiters_dir = "/data/serum-rbd-abtiters",
    serum_sarscov2_abtiters_dir = "/data/serum-sarscov2-abtiters"
  )
  
  #### Create an empty environment for storing various data objects
  data_env <- env(data_dirs = data_dirs)
  
  #### Prepare clinical data variables
  data_env <- prepare_clinical_data(data_env = data_env,
                                    DATA_VERSION = DATA_VERSION,
                                    KEEP_COVID19_POS = KEEP_COVID19_POS,
                                    FILTER_BY_CORE_ASSAY_COHORT = FILTER_BY_CORE_ASSAY_COHORT,
                                    FILTER_BY_INTEGRATION_COHORT = FILTER_BY_INTEGRATION_COHORT,
                                    KEEP_ESCALATIONS = KEEP_ESCALATIONS,
                                    KEEP_EMORY_CONTROLS = KEEP_EMORY_CONTROLS,
                                    EVENT_DATE_UPPER = EVENT_DATE_UPPER,
                                    VISIT_UPPER = VISIT_UPPER,
                                    VISIT_LOWER = VISIT_LOWER,
                                    SELECT_VISITS = SELECT_VISITS,
                                    SELECT_PARTICIPANT_ID = SELECT_PARTICIPANT_ID,
                                    EXCLUDE_SAMPLES = EXCLUDE_SAMPLES,
                                    PHASES = PHASES,
                                    USE_LOCKED_CLINICAL = USE_LOCKED_CLINICAL)
  
  #### Prepare info for the assay data files
  data_env <- prepare_assay_files_info(data_env)
  
  #### Fetch appropriate data files and load assay datasets
  data_env <- load_assay_data(data_env, DATA_VERSION, ALLOWED_SMPL_STATUS)
  
  #### Prepare predefined color schemes
  data_env <- prepare_color_schemes(data_env)
  
  #### Store DATA_VERSION and PHASES info in data_env
  data_env[[ "DATA_VERSION" ]] <- DATA_VERSION
  data_env[[ "PHASES" ]] <- PHASES
  data_env[[ "ALLOWED_SMPL_STATUS" ]] <- ALLOWED_SMPL_STATUS
  
  return(data_env)
  
}

### *** FUNCTIONS FOR load_IMPACC_datasets *** ###

# Load packages (needed for load_IMPACC_datasets above)
load_packages <- function(packages) {
  packages <-
    c("tidyverse",
      "ggplot2",
      "pals",
      "stringr",
      "ggpubr",
      "cowplot",
      "RColorBrewer",
      "rlang",
      "gridExtra",
      "corrr",
      "ComplexHeatmap",
      "missForest",
      "doRNG",
      "doParallel",
      # data preprocessing:
      "impute",
      # run PVCA:
      "lme4",
      "pvca", # devtools::install_github("dleelab/pvca"),
      "ggeffects",
      "sva",
      "limma",
      "edgeR"
      )
  
  for (n in seq_along(packages)) {
    suppressMessages(library(packages[n], character.only = TRUE))
  }
}

# Function used to prepare clinical data and attach it to a data_env object
prepare_clinical_data <- function(data_env,
                                  DATA_VERSION,
                                  KEEP_COVID19_POS,
                                  KEEP_ESCALATIONS,
                                  KEEP_EMORY_CONTROLS,
                                  FILTER_BY_CORE_ASSAY_COHORT,
                                  FILTER_BY_INTEGRATION_COHORT,
                                  SELECT_PARTICIPANT_ID,
                                  EXCLUDE_SAMPLES,
                                  EVENT_DATE_UPPER,
                                  VISIT_UPPER,
                                  VISIT_LOWER,
                                  SELECT_VISITS,
                                  PHASES,
                                  USE_LOCKED_CLINICAL) {
  
  if(FILTER_BY_CORE_ASSAY_COHORT & FILTER_BY_INTEGRATION_COHORT) {
    stop("Did you mean to set both FILTER_BY_CORE_ASSAY_COHORT & FILTER_BY_INTEGRATION_COHORT to TRUE? That is unusual! Please check your input.")
  }
  
  # Validating the provided value for the DATA_VERSION
  if(DATA_VERSION != "current") {
    x <- tryCatch({
      as.Date(DATA_VERSION)
    },
    error = function(e) {
      cat(
        "Error: Invalid input (\"",
        DATA_VERSION,
        "\") to DATA_VERSION variable. Either use \"current\" or a valid date in \"YYYY-MM-DD\" format for DATA_VERSION and rerun.\n",
        sep = ""
      )
      knitr::knit_exit()
    })
  }
  
  clinical_dir <- data_env[[ "data_dirs" ]][ "clinical_dir" ]
  
  ### Fetch clinical data files
  additional_regex_tags = ""
  if(USE_LOCKED_CLINICAL) {
    additional_regex_tags = "-locked"
  }
  clinical_sample_file <- fetch_file_names( DATA_VERSION, clinical_dir, paste0("sample", additional_regex_tags, ".csv$"))
  clinical_event_file <- fetch_file_names( DATA_VERSION, clinical_dir, paste0("event", additional_regex_tags, ".csv$"))
  clinical_individ_file <- fetch_file_names( DATA_VERSION, clinical_dir, paste0("individ", additional_regex_tags, ".csv$"))
  # New (convalescent and supplemental data):
  clinical_individ_supp_file <- fetch_file_names( DATA_VERSION, clinical_dir, paste0("individ-supp", additional_regex_tags, ".csv$"))
  
  if (is.null(clinical_sample_file) |
      is.null(clinical_event_file) | is.null(clinical_individ_file)) {
    stop(
      "Error: Either you don't have access to the clinical data or the required clinical data files don't exist for the version of data you have selected. Please make sure that *sample.csv, *event.csv and *.individ.csv files exist in \"",
      clinical_dir,
      "\".\nUSE_LOCKED_CLINICAL = TRUE to use the locked clinical files from /data/clinical/current, otherwise use a specific date for DATA_VERSION.\nExiting...\n",
      sep = ""
    )
  }
  
  message("Using the following clinical data files:\n")
  message("clinical_sample_file <- ", clinical_sample_file, "\n")
  message("clinical_event_file <- ", clinical_event_file, "\n")
  message("clinical_individ_file <- ", clinical_individ_file, "\n")
  message("clinical_individ_supp_file <- ", clinical_individ_supp_file, "\n")
  #message("clinical_participant_status_file <- ", clinical_participant_status_file, "\n")
  #message("clinical_convalescent_file <- ", clinical_convalescent_file, "\n")
  
  ### Load clinical data
  clinical_sample_data <- read_data_files( clinical_sample_file, header = TRUE )
  clinical_event_data <- read_data_files( clinical_event_file, header = TRUE )
  clinical_individ_data <- read_data_files( clinical_individ_file, header = TRUE )
  clinical_individ_supp_data <- read_data_files( clinical_individ_supp_file, header = TRUE )
  
  ### Combine clinical data
  clinical_data <-
    dplyr::inner_join(clinical_sample_data,
               clinical_event_data,
               by = c("event_id" = "event_id")) %>%    # Join the clinical_sample_data and clinical_event_data tables using "event_id" column
    dplyr::mutate(
      participant_id = participant_id.x,
      participant_id.x = NULL,
      participant_id.y = NULL
    ) %>%   # Remove redundant "participant_id" columns
    dplyr::inner_join(clinical_individ_data,
               by = c("participant_id" = "participant_id"))   # Append the clinical_individ_data table
  
  # Add new clinical data by participant_id:
  clinical_data <- dplyr::left_join(clinical_data, clinical_individ_supp_data, by="participant_id")

  
  ############ FILTER OUT SAMPLES VIA VARIABLES/FLAGS ######################
  
  ### Keep data from COVID-19 Positive individuals
  if(KEEP_COVID19_POS) {
    cat("Removing COVID-19 negative patients\n")
    clinical_data <- clinical_data %>% dplyr::filter(grepl("COVID-19 Positive", participant_type))
  }
  
  ### Keep core assay cohort
  if(FILTER_BY_CORE_ASSAY_COHORT) {
    cat("Removing samples not part of core_assay_cohort\n")
    clinical_data <- clinical_data %>% dplyr::filter( clinical_data$core_assay_cohort )
  }
  
  ### Keep integration cohort
  if(FILTER_BY_INTEGRATION_COHORT) {
    if(is.logical(clinical_data$integration_cohort)) {
      cat("Removing samples not part of integration_cohort\n")
      clinical_data <- clinical_data %>% dplyr::filter( clinical_data$integration_cohort )
    } else if(is.character(clinical_data$integration_cohort)) {
      cat("Removing samples not part of integration_cohort\n")
      clinical_data <- clinical_data %>% dplyr::filter( ! is.na(clinical_data$integration_cohort) )
    }
  }
  
  # Filter by provided participant IDs
  if( !is.null(SELECT_PARTICIPANT_ID) ){

    clinical_data <- clinical_data %>%
                      filter(participant_id %in% SELECT_PARTICIPANT_ID)
    
    cat("After filtering for selected participants the clinical data table contains ", 
            length(unique(clinical_data$participant_id)),
            " participants.")
  }
  
  ### Keep samples with event date <= EVENT_DATE_UPPER
  if (!is.null(EVENT_DATE_UPPER)) {
    cat("Removing samples with event_date greater than", EVENT_DATE_UPPER, "\n")
    clinical_data <- clinical_data %>% dplyr::filter( event_date <= EVENT_DATE_UPPER )
  }
  
  ### Keep samples from visit number >= VISIT_LOWER <= VISIT_UPPER
  if (is.null(VISIT_UPPER) && is.null(VISIT_LOWER)) {
    if (is.null(SELECT_VISITS)) {
      allowed_visits <- paste("Visit", 1:999)
    } else {
      allowed_visits <- paste("Visit", SELECT_VISITS)
    }
  } else if (!is.null(VISIT_UPPER) && !is.null(VISIT_LOWER)) {
    cat("Removing samples with visit number (event_type) greater than Visit", VISIT_UPPER, "\n")
    cat("Removing samples with visit number (event_type) smaller than Visit", VISIT_LOWER, "\n")
    allowed_visits <- paste("Visit", VISIT_LOWER:VISIT_UPPER)
  } else {
    stop("Please provide value for both VISIT_LOWER and VISIT_UPPER.")
  }

  ### Filter for Escalation visits
  if (KEEP_ESCALATIONS) {
    escalation_visits <- grepl("Escalation", clinical_data$event_type)
  } else {
    cat("Removing samples with 'Escalation' event_type\n")
    escalation_visits <- rep(FALSE, nrow(clinical_data))
  }

  # Filter for Control donor visits
  if (KEEP_EMORY_CONTROLS) {
    control_visits <- grepl("Control donor visit", clinical_data$event_type)
  } else {
    cat("Removing samples with 'Control donor visit' event_type\n")
    control_visits <- rep(FALSE, nrow(clinical_data))
  }

  # Apply all visit filters
  clinical_data <- clinical_data %>% dplyr::filter( control_visits | escalation_visits | (event_type %in% allowed_visits) )
  
  ### Keep samples from the selected phases
  clinical_data <- clinical_data %>% 
    dplyr::filter( phase %in% c(PHASES) )
  
  ### Stop if no samples are left
  cat(paste0("There are ", nrow(clinical_data), " event entries in the clinical data left after filtering."))
  if (nrow(clinical_data) == 0){
    stop("Error: there are no entries left in the clinical data after filtering.\n Please revise your filtering criteria.")
  }
  
  ############ ADD NEW VARIABLES / FLAGS ######################
  
  ### Generate discretized age based in 5 quantiles.
  admit_age_levels = sort(unique(clinical_data$admit_age))
  admit_age_levels_tiles = ntile(admit_age_levels, 5)
  label_df_tmp = data.frame( admit_age_levels_tiles, admit_age_levels) %>% 
    group_by(admit_age_levels_tiles) %>% 
    summarize(min=min(admit_age_levels), max=max(admit_age_levels)) %>% 
    mutate(labels = paste0("[",min, ",", max, "]") ) 
  admit_age_levels_tiles_labels = label_df_tmp$labels %>% setNames(label_df_tmp$admit_age_levels_tiles)
  admit_age_levels_labels = admit_age_levels_tiles_labels[admit_age_levels_tiles] %>% setNames(admit_age_levels)
  discretized_admit_age_quantile = admit_age_levels_labels[ as.character(clinical_data$admit_age) ]
  
  ### Generate discretized age based on equal width ages.
  discretized_admit_age_equalWidth = cut(clinical_data$admit_age, breaks = 5)
  clinical_data = clinical_data %>% dplyr::mutate(discretized_admit_age_quantile = discretized_admit_age_quantile,
                                                  discretized_admit_age_equalWidth = discretized_admit_age_equalWidth)
  
  ### Add discharge dates column
  # Extract discharge dates, concatenate the dates for patients that have multiple dates.
  discharge_event_data = clinical_event_data[ clinical_event_data$event_type == "Discharge", ] %>% 
    dplyr::group_by(participant_id) %>% 
    dplyr::summarise(event_dates = paste0(sort(event_date), collapse = ",") )
  clinical_data = clinical_data %>% 
    dplyr::mutate(discharge_dates = 
             discharge_event_data$event_dates[ match(clinical_data$participant_id, discharge_event_data$participant_id) ]  
    )
  
  death_event_data = clinical_event_data[ clinical_event_data$event_type == "Death", ]
  clinical_data = clinical_data %>% 
    dplyr::mutate(death_date = death_event_data$event_date[ match(clinical_data$participant_id, death_event_data$participant_id) ]
    )
  
  ############ SAVE AND RETURN CLINICAL OBJECT ######################
  
  data_env[[ "clinical_data" ]] <- clinical_data
  data_env[[ "clinical_sample_file" ]] <- clinical_sample_file # by sample
  data_env[[ "clinical_event_file" ]] <- clinical_event_file # by event
  data_env[[ "clinical_individ_file" ]] <- clinical_individ_file # by participant
  data_env[[ "clinical_individ_supp_file" ]] <- clinical_individ_supp_file # by participant
  #data_env[[ "clinical_participant_status_file" ]] <- clinical_participant_status_file # by participant
  #data_env[[ "clinical_convalescent_file" ]] <- clinical_convalescent_file # by participant

  return(data_env)
}

# Assign elements of a list to variables
assign_env_vars <- function(env) {
  lnames <- names(list)
  if( ! is.null( lnames ) ) {
    for( n in lnames ) {
      assign( n, list[[n]] )
    }
  }
}

# A function to fetch "current" or "legacy" version of the data. The function will return data files that were current on the indicated date.
fetch_file_names <- function(DATA_VERSION, dir, pattern = ".Counts.csv") {
  # Fetch all files within the dir
  files <- list.files(dir, recursive = T)
  
  # Select the files associated with the selection version of the data.
  if (DATA_VERSION == "current") {
    files <- files[grep("current", files)]
    
  } else {
    # Will store the date to use in legacy_date
    legacy_date <- ""
    # Today's date
    today_date <- Sys.Date()
    
    files_legacy <- files[grep("legacy", files)]
    dates <-
      as.Date(stringr::str_match(files_legacy, "/(\\d{4}-\\d{2}-\\d{2})/")[, 2])
    unique_reverse_sorted_dates <-
      sort(unique(dates), decreasing = T)
    
    # If there is no data in legacy folder, use today's date as legacy_date
    if (length(unique_reverse_sorted_dates) == 0) {
      legacy_date = as.Date(NULL)   # Setting as a null date (i.e. length of the date is 0), which will be used later to pick "current" dataset.
      
    } else {
      # Else, find an appropriate version date.
      
      user_version_date <- as.Date(DATA_VERSION)
      # If the user's version date is same as today's date or is more recent than the most recent folder in "legacy" folder, use the today's date as the legacy date.
      if (today_date == user_version_date ||
          unique_reverse_sorted_dates[1] < user_version_date) {
        legacy_date <- as.Date(NULL)   # Setting as a null date (i.e. length of the date is 0), which will be used later to pick "current" dataset.
        
      } else {
        # Else, find the data version from legacy folder that was current on the user's version date.
        for (i in seq_along(unique_reverse_sorted_dates)) {
          if (user_version_date <= unique_reverse_sorted_dates[i]) {
            legacy_date <- unique_reverse_sorted_dates[i]
            
          }
        }
      }
    }
    
    
    if ( length(legacy_date) == 0 ) {
      files <- files[grep("current", files)]
      
    } else {
      files <- files_legacy[dates == legacy_date]
      
    }
    
    
  }
  
  # Select a file that matches the pattern.
  matched_file_name <- files[grepl(pattern, files)]
  # If multiple file exist, and one of them is .csv or .tsv, take that.
  is_csv_tsv <- grepl(".tsv$|.csv$", matched_file_name)
  if (length(matched_file_name) > 1 && sum(is_csv_tsv) > 0) {
    matched_file_name <- matched_file_name[is_csv_tsv][1]
  }
  if (length(matched_file_name) > 1) {
    warning("More than one files match \"",
            pattern  ,
            "\" in \"",
            dir ,
            "\". Skipping this dataset.\n")
    return()
  } else if (length(matched_file_name) == 0) {
    #cat("No files match \"", pattern  ,"\" in \"", dir , "\".\n", sep = "")
    return()
  }
  
  return(file.path(dir, matched_file_name))
}

# Function used to read csv/tsv files (called from load_IMPACC_datasets(...))
read_data_files <- function(file, row.names = NULL, header = FALSE) {
  if(file.access(file, mode = 4) == -1) {
    warning(file, " does not have a read permission.")
    return(data.frame())
  }
  if (grepl(".csv$", file)) {
    x <- tryCatch({
      data <- read.csv(file, row.names = row.names, check.names = F, header=header, comment.char = "#")
    },
    error = function(e) {
      cat("Cannot read ", file,". Please make sure it is a comma-separated file.\nSee the original error message for more detail: ", e$message, "\n")
      return(data.frame())
    })
  } else if (grepl(".tsv$", file)) {
    x <- tryCatch({
      data <- read.table(file, row.names = row.names, check.names = F, header=header)
    },
    error = function(e) {
      cat("Cannot read ", file,". Please make sure it is a tab-separated file.\nSee the original error message for more detail: ", e$message, "\n")
      return(data.frame())
    })
  } else {
    warning("Unrecognized file extension for ", file, ". Currently support only .csv and .tsv files.\n")
    return(data.frame())
  }
  return(data)
}

# Prepare info for the assay data files
prepare_assay_files_info <- function(data_env) {
  ## Individual assay data files are fetched for the selected version of data and stored in a list for easy handling of these large number of data files.
  
  data_dirs <- data_env[[ "data_dirs" ]]
  for( n in names(data_dirs) ) {
    assign( n, data_dirs[ n ] )
  }
  
  ### A list containing, for individual assay files, file-name-pattern, directory, variable name to load the data to and description of the data file/type.
  assay_files_info <- list(
    "plasma_proteomics_targeted" = list(
                            dir = plasma_proteomics_targeted_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "plasma_proteomics_targeted_counts",
                            count_desc = "Count data from targeted proteomics of plasma",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "plasma_proteomics_targeted_rowfeature",
                            rowfeature_desc = "RowFeature data from targeted proteomics of plasma",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "plasma_proteomics_targeted_metadata",
                            metadata_desc = "Metadata from targeted proteomics of plasma",
                            metadata_plate = "Plate.id"
    ),
    "plasma_proteomics_global_dda" = list(
                            dir = plasma_proteomics_global_dda_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "plasma_proteomics_global_dda_counts",
                            count_desc = "Count data from global proteomics (DDA) of plasma",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "plasma_proteomics_global_dda_rowfeature",
                            rowfeature_desc = "RowFeature data from global proteomics (DDA) of plasma",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "plasma_proteomics_global_dda_metadata",
                            metadata_desc = "Metadata from global proteomics (DDA) of plasma",
                            metadata_plate = "plate_num"
    ),
    "plasma_proteomics_global_dia" = list(
                            dir = plasma_proteomics_global_dia_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "plasma_proteomics_global_dia_counts",
                            count_desc = "Count data from global proteomics (DIA) of plasma",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "plasma_proteomics_global_dia_rowfeature",
                            rowfeature_desc = "RowFeature data from global proteomics (DIA) of plasma",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "plasma_proteomics_global_dia_metadata",
                            metadata_desc = "Metadata from global proteomics (DIA) of plasma",
                            metadata_plate = "plate_num"
    ),
    "serum_olink" = list(
                            dir = serum_olink_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "serum_olink_counts",
                            count_desc = "Count data from Olink assay of serum",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "serum_olink_rowfeature",
                            rowfeature_desc = "RowFeature data from Olink assay of serum",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "serum_olink_metadata",
                            metadata_desc = "Metadata from Olink assay of serum",
                            metadata_plate = "plate"
    ),
    "nasal_viralload" = list(
                            dir = nasal_viralload_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "nasal_viralload_counts",
                            count_desc = "Count data for nasal viral load",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "nasal_viralload_rowfeature",
                            rowfeature_desc = "RowFeature data for nasal viral load",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "nasal_viralload_metadata",
                            metadata_desc = "Metadata for nasal viral load",
                            metadata_plate = "plate_num"
    ),
    "serum_rbd_abtiters" = list(
                            dir = serum_rbd_abtiters_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "serum_rbd_abtiters_counts",
                            count_desc = "Count data for RBD ab-titer in serum",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "serum_rbd_abtiters_rowfeature",
                            rowfeature_desc = "RowFeature data for RBD ab-titer in serum",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "serum_rbd_abtiters_metadata",
                            metadata_desc = "Metadata for RBD ab-titer in serum",
                            metadata_plate = "plate_num"
    ),
    "serum_sarscov2_abtiters" = list(
                            dir = serum_sarscov2_abtiters_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "serum_sarscov2_abtiters_counts",
                            count_desc = "Count data for SARS-CoV-2 ab-titer in serum",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "serum_sarscov2_abtiters_rowfeature",
                            rowfeature_desc = "RowFeature data for SARS-CoV-2 ab-titer in serum",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "serum_sarscov2_abtiters_metadata",
                            metadata_desc = "Metadata for SARS-CoV-2 ab-titer in serum",
                            metadata_plate = "plate_num"
    ),
    "nasal_transcriptomics" = list(
                            dir = nasal_transcriptomics_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "nasal_transcriptomics_counts",
                            count_desc = "Count data for nasal transcriptomics",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "nasal_transcriptomics_rowfeature",
                            rowfeature_desc = "RowFeature data for nasal transcriptomics",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "nasal_transcriptomics_metadata",
                            metadata_desc = "Metadata for nasal transcriptomics",
                            metadata_plate = "plate"
    ),
    "plasma_metabolomics_global" = list(
                            dir = plasma_metabolomics_global_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "plasma_metabolomics_global_counts",
                            count_desc = "Count data for global metabolomics of plasma",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "plasma_metabolomics_global_rowfeature",
                            rowfeature_desc = "RowFeature data for global metabolomics of plasma",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "plasma_metabolomics_global_metadata",
                            metadata_desc = "Metadata for global metabolomics of plasma",
                            metadata_plate = "plate"
    ),
    "bld_cytof" = list(
                            dir = bld_cytof_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "bld_cytof_counts",
                            count_desc = "Count data for CyTOF of blood",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "bld_cytof_rowfeature",
                            rowfeature_desc = "RowFeature data for CyTOF of blood",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "bld_cytof_metadata",
                            metadata_desc = "Metadata for CyTOF of blood",
                            metadata_plate = "plate_num"
    ),
    "bld_msi_cytof" = list(
                            dir = bld_cytof_dir,
                            count_pattern = "MeanChannelIntensity-Counts.*$",
                            count_file_path = "",
                            count_var_name = "bld_cytof_msi_counts",
                            count_desc = "MSI count data for CyTOF of blood",
                            rowfeature_pattern = "MeanChannelIntensity-RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "bld_cytof_msi_rowfeature",
                            rowfeature_desc = "RowFeature data for MSI CyTOF of blood",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "bld_cytof_metadata",
                            metadata_desc = "Metadata for CyTOF of blood",
                            metadata_plate = "plate_num"
    ),
    "ea_cytof" = list(
                            dir = ea_cytof_dir,
                            count_pattern = "cytof-Counts.*$",
                            count_file_path = "",
                            count_var_name = "ea_cytof_counts",
                            count_desc = "Count data for CyTOF of EA",
                            rowfeature_pattern = "cytof-RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "ea_cytof_rowfeature",
                            rowfeature_desc = "RowFeature data for CyTOF of EA",
                            metadata_pattern = "cytof-Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "ea_cytof_metadata",
                            metadata_desc = "Metadata for CyTOF of EA",
                            metadata_plate = "plate_num"
    ),
    "ea_transcriptomics" = list(
                            dir = ea_transcriptomics_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "ea_transcriptomics_counts",
                            count_desc = "Count data for EA transcriptomics",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "ea_transcriptomics_rowfeature",
                            rowfeature_desc = "RowFeature data for EA transcriptomics",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "ea_transcriptomics_metadata",
                            metadata_desc = "Metadata for EA transcriptomics",
                            metadata_plate = "plate"
    ),
    "pbmc_transcriptomics" = list(
                            dir = pbmc_transcriptomics_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "pbmc_transcriptomics_counts",
                            count_desc = "Count data for PBMC transcriptomics",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "pbmc_transcriptomics_rowfeature",
                            rowfeature_desc = "RowFeature data for PBMC transcriptomics",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "pbmc_transcriptomics_metadata",
                            metadata_desc = "Metadata for PBMC transcriptomics",
                            metadata_plate = "plate_num"
    ),
    "ea_metagenomics" = list(
                            dir = ea_metagenomics_dir,
                            count_pattern = "BacterialTaxon-Counts_RPM.*$",
                            count_file_path = "",
                            count_var_name = "ea_metagenomics_counts",
                            count_desc = "Count data for EA metagenomics",
                            rowfeature_pattern = "BacterialTaxon-RowFeatures.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "ea_metagenomics_rowfeature",
                            rowfeature_desc = "RowFeature data for EA metagenomics",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "ea_metagenomics_metadata",
                            metadata_desc = "Metadata for EA metagenomics",
                            metadata_plate = "plate_num"
    ),
    "nasal_metagenomics" = list(
                            dir = nasal_metagenomics_dir,
                            count_pattern = "BacterialTaxon-Counts_RPM.*$",
                            count_file_path = "",
                            count_var_name = "nasal_metagenomics_counts",
                            count_desc = "Count data for nasal metagenomics",
                            rowfeature_pattern = "BacterialTaxon-RowFeatures.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "nasal_metagenomics_rowfeature",
                            rowfeature_desc = "RowFeature data for nasal metagenomics",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "nasal_metagenomics_metadata",
                            metadata_desc = "Metadata for nasal metagenomics",
                            metadata_plate = "plate"
    ),
    "serum_autoantibody" = list(
                            dir = serum_autoantibody_dir,
                            count_pattern = "Counts.*$",
                            count_file_path = "",
                            count_var_name = "serum_autoantibody_counts",
                            count_desc = "Count data for serum autoantibody",
                            rowfeature_pattern = "RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "serum_autoantibody_rowfeature",
                            rowfeature_desc = "RowFeature data for serum autoantibody",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "serum_autoantibody_metadata",
                            metadata_desc = "Metadata for serum autoantibody"
    )
  )
  
  data_env[["assay_files_info"]] <- assay_files_info
  return(data_env)
}

# Load assay datasets
load_assay_data <- function(data_env, 
                            DATA_VERSION,
                            ALLOWED_SMPL_STATUS) {
  ### Fetch path of files that match the indicated pattern in the provided directory for each assay file and load the data in the variable names indicated in var_name element of the assay_files_info list. Uses the version of data specified in DATA_VERSION.
  
  assay_files_info <- data_env[[ "assay_files_info" ]]
  clinical_data <- data_env[[ "clinical_data" ]]
  
  for (f in names(assay_files_info)) {
    f_info <- assay_files_info[[f]]
    
    # Loading Metadata table
    if (!is.null(f_info$metadata_pattern)) {
      f_info$metadata_file_path <-
        fetch_file_names(DATA_VERSION, f_info$dir, f_info$metadata_pattern)
      
      if (is.null(f_info$metadata_file_path)) {
        warning(
          "\"",
          f_info$metadata_desc,
          "\" is not accessible. Either you don't have access to this data or the data does not exist for the version of dataset you have selected.\n"
        )
        
      } else {
        ### Load data, keep only those observations for which the clinical data is available and assign the data object to the variable name indicated in "var_name" element.
        cat(
          "Loading ",
          f_info$metadata_desc,
          " from ",
          f_info$metadata_file_path,
          "\n"
        )
        data <-
          read_data_files(f_info$metadata_file_path,
                          row.names = 1,
                          header = TRUE)

        ## Keep those records for which clinical data exists, those that belong to those that have passed/questionable QC.
        data <-
          data %>% dplyr::filter( row.names(.) %in% clinical_data$sample_id &
                            sample_status %in% ALLOWED_SMPL_STATUS )
        
        ## Save the data to the environment.
        data_env[[ f_info$metadata_var_name ]] <- data
        
      }
    }
    
    
    
    # Loading count data
    f_info$count_file_path <-
      fetch_file_names(DATA_VERSION, f_info$dir, f_info$count_pattern)
    
    if (is.null(f_info$count_file_path)) {
      warning(
        "\"",
        f_info$count_desc,
        "\" is not accessible. Either you don't have access to this data or the data does not exist for the version of dataset you have selected.\n"
      )
      
    } else {
      ### Load data, keep only those observations for which the clinical data is available and assign the data object to the variable name indicated in "var_name" element.
      cat("Loading ",
              f_info$count_desc,
              " from ",
              f_info$count_file_path,
              "\n")
      data <-
        read_data_files(f_info$count_file_path,
                        row.names = 1,
                        header = TRUE)

      # In case, if the CSV file has an extra comma at the end of rows, the following line removes the last empty column.
      #data = data %>% select( which(colnames(data) != "") )
      
      
      ## Keep only those records for which clinical data and metadata exist. Clinical data and metadata are filtered using variety of filters (see above).
      data <-
        data %>% dplyr::filter(row.names(.) %in% clinical_data$sample_id &
                          row.names(.) %in% row.names(data_env[[f_info$metadata_var_name]]))
      
      ## Match the order of samples in metadata table with that in the count table.
      data_env[[ f_info$metadata_var_name ]] <-
        data_env[[ f_info$metadata_var_name ]] %>% dplyr::slice( match( row.names(data), row.names(.) ) )
      
      
      ## Save the data to the environment.
      data_env[[ f_info$count_var_name ]] <- data
      
    }
    
    
    
    # Loading rowfeature data
    if (!is.null(f_info$rowfeature_pattern)) {
      f_info$rowfeature_file_path <-
        fetch_file_names(DATA_VERSION, f_info$dir, f_info$rowfeature_pattern)
      
      if (is.null(f_info$rowfeature_file_path)) {
        warning(
          "\"",
          f_info$rowfeature_desc,
          "\" is not accessible. Either you don't have access to this data or the data does not exist for the version of dataset you have selected.\n"
        )
        
      } else {
        ### Load data, keep only those observations for which the clinical data is available and assign the data object to the variable name indicated in "var_name" element.
        cat(
          "Loading ",
          f_info$rowfeature_desc,
          " from ",
          f_info$rowfeature_file_path,
          "\n"
        )
        
        data <-
          read_data_files(f_info$rowfeature_file_path,
                          row.names = 1,
                          header = TRUE)
        
        ## Save the data to the environment.
        data_env[[ f_info$rowfeature_var_name ]] <- data
        
      }
    }

    
    
    assay_files_info[[f]] <- f_info
    
  }
  
  data_env[[ "assay_files_info" ]] <- assay_files_info
  
  return(data_env)
}


################################## Other helper functions #########################################

# Prepare color schemes
prepare_color_schemes <- function(data_env) {
  
  # Load clinical_data table into the local environment
  clinical_data <- data_env[[ "clinical_data" ]]
  
  ### Color for NA values
  data_env[[ "color_na" ]] <- "grey"
  
  ### Color scheme for sex
  data_env[[ "color_sex" ]] <- c("Female" = "#DC2543", "Male" = "#1E94A0")
  
  ### Color scheme for sample_type
  sample_type_levels <- unique( clinical_data$sample_type )
  data_env[[ "color_sample_type" ]] <- get_categorical_color_vector( sample_type_levels )
  
  ### Color scheme for respiratory_status
  respiratory_status_levels <- unique( sort( clinical_data$respiratory_status ) )
  data_env[[ "color_respiratory_status" ]] <- get_ordinal_color_vector( respiratory_status_levels )
  
  ### Color scheme for event_type
  event_type_levels <- unique( sort( clinical_data$event_type ) )
  data_env[[ "color_event_type" ]] <- get_ordinal_color_vector( event_type_levels )
  
  ### Color scheme for race
  race_levels <- unique( sort( clinical_data$race ) )
  data_env[[ "color_race" ]] <- get_categorical_color_vector( race_levels )
  
  
  data_env[[ "color_schemes" ]] <- data.frame(
    object_name = c("color_sex", "color_sample_type", "color_respiratory_status", "color_event_type", "color_race", "color_na"),
    feature_id = c("sex", "sample_type", "respiratory_status", "event_type", "race", "NA"),
    feature_name = c("Sex", "Sample type", "Respiratory status", "Event type", "Race", "NA values")
  )
  
  return(data_env)
}
# Get data summary
# Authors: Ravi K. Patel, Jeremy Gygi
get_data_summary <- function(data_env) { # Print names of data files being used and the object names where the data is stored.
  df <-
    data.frame(
      dataset = "Clinical data",
      file = paste(
        clinical_sample_file,
        clinical_event_file,
        clinical_individ_file,
        sep = "; "
      ),
      object_name = "clinical_data",
      dimensions = paste(dim(data_env$clinical_data), collapse = " x ")
    )
  
  for (f in names(data_env$assay_files_info)) {
    f_info <- data_env$assay_files_info[[f]]
    if (!is.null(f_info$count_file_path)) {
      df <-
        rbind(
          df,
          data.frame(
            dataset = f_info$count_desc,
            file = f_info$count_file_path,
            object_name = f_info$count_var_name,
            dimensions = paste(dim( get(f_info$count_var_name) ), collapse = " x ")
          )
        )
    }
    if (!is.null(f_info$rowfeature_file_path)) {
      df <-
        rbind(
          df,
          data.frame(
            dataset = f_info$rowfeature_desc,
            file = f_info$rowfeature_file_path,
            object_name = f_info$rowfeature_var_name,
            dimensions = paste(dim( get(f_info$rowfeature_var_name) ), collapse = " x ")
          )
        )
    }
    if (!is.null(f_info$metadata_file_path)) {
      df <-
        rbind(
          df,
          data.frame(
            dataset = f_info$metadata_desc,
            file = f_info$metadata_file_path,
            object_name = f_info$metadata_var_name,
            dimensions = paste(dim( get(f_info$metadata_var_name) ), collapse = " x ")
          )
        )
    }
  }
  
  return(df)
}



# Generate and print sample plots (using environmental variables)
generate_sample_plots <- function() {
  
  # Distribution of samples by sex
  bar_sex <- clinical_data %>% count(sex) %>% 
    ggplot( aes(x="Sex", y=n, fill=sex) ) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual( values = color_sex ) + theme_classic() + 
    guides( fill = guide_legend(title="") ) +
    xlab("") + ylab("Frequency")
  
  
  # Distribution of samples by respiratory_status
  bar_respiratory_status <- clinical_data %>% count(respiratory_status) %>% 
    ggplot( aes(x="Respiratory status", y=n, fill = as.factor( respiratory_status ) ) ) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual( values = color_respiratory_status, na.value = color_na ) + theme_classic() + 
    guides( fill = guide_legend(title="Respiratory status") ) +
    xlab("") + ylab("Frequency")
  
  
  # Scatter plot of ITIH4 and SERPINA1 levels
  scatter_ITIH4_vs_SERPINA1 <- plasma_proteomics_targeted_counts %>% ggplot(aes(y = ITIH4, x = SERPINA1 ) ) + 
    geom_point(size = 3, shape = 21, alpha = 0.8) + 
    geom_smooth(method = "lm", alpha = .15) + 
    theme_classic() +
    xlab("SERPINA1 levels") + ylab("ITI45 levels") +
    ggtitle("ITIH4 levels vs. SERPINA1 levels")
  
  
  
  # Jitter plot of ITIH4 levels vs. respiratory_status colored by event_type
  df <- data.frame( 
    ITIH4 = plasma_proteomics_targeted_counts[,"ITIH4"], 
    respiratory_status = clinical_data[ match( row.names(plasma_proteomics_targeted_counts), clinical_data$sample_id ), "respiratory_status" ], 
    event_type = clinical_data[ match( row.names(plasma_proteomics_targeted_counts), clinical_data$sample_id ), "event_type" ] )
  
  jitter_ITIH4_vs_resp <- df %>% ggplot(aes(y = ITIH4, x = as.factor( respiratory_status ), fill = event_type ) ) + 
    geom_point(position = position_jitter(width = 0.25), size = 3, shape = 21, alpha = 0.8) + 
    scale_fill_manual( values = color_event_type ) + 
    theme_classic() +
    labs( fill = "Event type") +
    xlab("respiratory_status") + ylab("ITIH4 levels") +
    ggtitle("ITIH4 levels vs. respiratory_status")
  
  
  
  # Example of imputing missing values
  plasma_proteomics_targeted_counts_imputed = impute::impute.knn(as.matrix(plasma_proteomics_targeted_counts))
  
  na_vals <-
    is.na(
      plasma_proteomics_targeted_counts %>% pivot_longer(cols = colnames(.),  names_to = "features") %>% pull(value)
    )
  
  imputed_vals <-
    plasma_proteomics_targeted_counts_imputed$data %>% data.frame() %>% 
    pivot_longer(cols = colnames(.),  names_to = "features") %>% 
    dplyr::filter(na_vals) %>% pull(value)
  
  unimputed_vals <-
    plasma_proteomics_targeted_counts %>% 
    pivot_longer(cols = colnames(.),  names_to = "features") %>% 
    dplyr::filter(!na_vals) %>% pull(value)
  
  # Comparing the Non-values with the NA values after impution.
  hist_imputed_vals <- rbind(
    data.frame(label = "Non-NA values", vals = unimputed_vals),
    data.frame(label = "NA values", vals = imputed_vals)
  ) %>%
    ggplot(aes(x = vals, fill = label)) +
    geom_density(adjust = 1.5, alpha = .4) +
    scale_fill_manual(values = c("#69b3a2", "#404080")) +
    theme_classic() +
    ggtitle("Distribution of Non-NA values and\n NA values after imputation") +
    labs(fill = "") +
    xlab("Expression values")
  
  
  print(ggarrange(bar_respiratory_status, jitter_ITIH4_vs_resp, scatter_ITIH4_vs_SERPINA1, hist_imputed_vals,
                  widths = c(1,1.5), labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2))
  
  
  # Generate plot:
  print(ggplot(clinical_data) +
          geom_line(aes(x = event_date, y = respiratory_status, group = participant_id),
                    alpha = .1) +
          geom_point(
            aes(
              x = event_date,
              y = respiratory_status,
              group = participant_id,
              fill = as.character(respiratory_status)
            ),
            size = 4,
            shape = 21
          ) +
          # add color scale:
          scale_fill_manual(values = color_respiratory_status, na.value = color_na) +
          ggtitle("Respiratory Status vs. Time") +
          theme_classic() +
          labs(fill = "respiratory_status"))
  
  
  # Generate heatmap:
  melted.df <-
    tidyr::pivot_longer(clinical_data, cols = respiratory_status)
  melted.df <-
    mutate(melted.df, respiratory_status = as.character(value))
  ggplot(melted.df) +
    geom_tile(aes(x = participant_id, y = event_type, fill = respiratory_status)) +
    # add color scale:
    scale_fill_manual(values = color_respiratory_status, na.value = color_na) +
    ggtitle("respiratory_status heatmap") +
    theme_bw() +
    theme(axis.text.x = element_blank())
  
}

# A function to generate a named vector of colors for a categorical variable using the levels_vector (a vector of unique categories) provided by the user.
get_categorical_color_vector <- function(levels_vector) {
  return(pals::alphabet(length(levels_vector)) %>% setNames(levels_vector))
}

# A function to generate a named vector of colors for an ordinal variable using the ordered_levels_vector (a vector of unique categories arranged in an increasing order) provided by the user. The user can provide name of RColorBrewer sequential palette or just assign a value between 1 and 18 to palette_index to choose an RColorBrewer sequential palette at that index.
get_ordinal_color_vector <- function(ordered_levels_vector,
           rcolorbrewer_palette_name = NULL,
           palette_index = 1,
           bidirectional = FALSE) {
    # A list of RColorBrewer sequential palettes
    rcolorbrewer_palettes <-
      brewer.pal.info %>% dplyr::filter(category == "seq") %>% row.names() %>% rev()
    
    # Select a palette that matches user's selection
    if (is.null(rcolorbrewer_palette_name) & palette_index == 1) {
      rcolorbrewer_palette_name <- rcolorbrewer_palettes[1]
    }
    if (palette_index != 1) {
      rcolorbrewer_palette_name <- rcolorbrewer_palettes[palette_index]
    }
    
    # Create a color vector
    palette <-
      colorRampPalette(brewer.pal(9, rcolorbrewer_palette_name))
    color_vector <-
      palette(length(ordered_levels_vector)) %>% setNames(ordered_levels_vector)
    
    # If one of the values is "NA", use grey color for the category.
    # color_vector["NA"] <- "grey"
    return(color_vector)
}

# Generate and return a color-key (legend) displaying the provided color scheme.
display_color_scheme <- function(color_scheme_vector, title = "") {
  p <- data.frame(color_scheme_vector = color_scheme_vector) %>%
    ggplot(aes(
      x = row.names(.),
      y = 1:nrow(.),
      fill = row.names(.)
    )) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = color_scheme_vector) +
    guides(fill = guide_legend(title = title))
  return(get_legend(p))
}

#' @title plot_assay
#' @description Plot a simple heatmap of a matrix (assay). 
#' @author Ravi Patel
#' @param data_env The data environment produced by load_IMPACC_datasets / load_IMPACC_Publication_datasets functions.
#' @param impute_alpha Optional, features with "impute_alpha" or more missing values are discarded within missForest imputation.
#' @param random_seed Optional, random seed used by missForest
#' @param impute_cores Optional, the number of cores to use for imputation.
#' @param assays_to_process Optional, a vector of strings with assay names to be preprocessed. If `NULL`, will preprocess all assays.
#' @param run_pvca logical perform PVCA analysis? Defaults to `TRUE`
#' @return list of preprocessing data matrices and pvca plots: list(preprocessed_data=preprocessed_data, pvca_data = pvca_data)
plot_assay <- function(assay_matrix, row.names = FALSE, col.names = FALSE){
  tidy_mat <- reshape2::melt(as.matrix(assay_matrix))
  g <- ggplot2::ggplot(tidy_mat) +
    ggplot2::geom_tile(ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::scale_fill_gradient2(low = "darkblue", high = "red", na.value = "black") +
    ggplot2::theme_bw() +
    ggplot2::xlab("Analytes") +
    ggplot2::ylab("Samples")
  
  if(!row.names & !col.names){
    g <- g + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  } else if(!row.names){
    g <- g + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  } else if(!col.names){
    g <- g + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  }
  
  return(g)
}



################################ Functions for assay pre-processing ########################

#' @title preprocess_IMPACC_datasets
#' @description preprocess the data using assay specific preprocessing functions.
#' @author Ravi Patel, Jeremy Gygi, Gisela Gabernet
#' @param data_env The data environment produced by load_IMPACC_datasets / load_IMPACC_Publication_datasets functions.
#' @param impute_alpha Optional, features with "impute_alpha" or more missing values are discarded within missForest imputation.
#' @param random_seed Optional, random seed used by missForest
#' @param impute_cores Optional, the number of cores to use for imputation.
#' @param proteomics_missing_fraction_cutoff: maximum allowed fraction of missing values in features. If set to a value (e.g. 0.99), then proteins will be filtered with more than proteomics_missing_fraction_cutoff missing values.
#' @param proteomics_selected_features: Optional, a vector of strings with feature IDs to be selected for downstream analysis. If `NULL`, will instead perform missing fraction cutoff.
#' @param proteomics_batch_correction: perform batch correction on plate and phase using ComBat? Defaults to `TRUE`
#' @param proteomics_combat_matrix_design: design model for ComBat as a string (e.g. ~outcome*event_type), the provided variables need to exist in the clinical file.
#' @param proteomics_combat_matrix_design_variables: design model variables that are used in the design as a vector (e.g. c("outcome","event_type"))
#' @param transcriptomics_batch_correction: perform batch correction on plate and phase using limma removeBatch effects? Defaults to `TRUE`
#' @param transcriptomics_expression_filter: Proportion of samples to apply gene expression feature filtering globally or per outcome group (default 0.05, 5%).
#' @param transcriptomics_expression_filter_grouping: Outcome variable for grouping samples prior to filtering gene features.
#' @param transcriptomics_selected_features: Optional, a vector of strings with feature IDs to be selected for downstream analysis. If `NULL`, will instead perform missing fraction cutoff.
#' @param transcriptomics_limma_matrix_design: design model for limma removeBatcheffect as a string (e.g. ~outcome*event_type), the provided variables need to exist in the clinical file
#' @param transcriptomics_limma_matrix_design_variables: design model variables that are used in the design as a vector (e.g. c("outcome","event_type"))
#' @param transcriptomics_mad_gene_fraction: Proportion of genes to apply gene expression feature filtering globally or per outcome group (default 0.75, 75%).
#' @param metabolomics_selected_features: Optional, a vector of strings with feature IDs to be selected for downstream analysis. If `NULL`, will instead perform IQR feature selection.
#' @param cytof_outlier_samples: Optional, a vector of strings with sample IDs to be removed from the CyTOF data.
#' @param cytof_msi_missing_fraction_cutoff: maximum allowed fraction of missing values in features. If set to a value (e.g. 0.5), then features will be filtered with more than cytof_msi_missing_fraction_cutoff missing values.
#' @param cytof_msi_batch_correction: perform batch correction on plate and phase using limma removeBatch effects? Defaults to `TRUE`.
#' @param cytof_msi_limma_matrix_design: design model for limma removeBatcheffect as a string (e.g. ~outcome*event_type), the provided variables need to exist in the clinical file.
#' @param cytof_msi_limma_matrix_design_variables: design model variables that are used in the design as a vector (e.g. c("outcome","event_type"))
#' @param serum_rbd_abtiters_batch_correction: perform batch correction on plate and phase using limma removeBatch effects? Defaults to `TRUE`
#' @param assays_to_process Optional, a vector of strings with assay names to be preprocessed. If `NULL`, will preprocess all assays.
#' @param run_pvca logical perform PVCA analysis? Defaults to `TRUE`
#' @param additional_pvca_variables character vector of additional variables to include in pvca analysis, they need to exist in the clinical metadata annotation.
#' @return list of preprocessing data matrices and pvca plots: list(preprocessed_data=preprocessed_data, pvca_data = pvca_data)
#' @examples \dontrun{ preprocessed_data <- preprocess_data(data_env) }
# previously called preprocess_data
preprocess_IMPACC_datasets <- function(data_env, 
                                        impute_alpha=0.3, 
                                        random_seed=2021,  
                                        impute_cores=8, 
                                        proteomics_missing_fraction_cutoff = 0.99,
                                        proteomics_selected_features = NULL,
                                        proteomics_batch_correction = FALSE,
                                        proteomics_combat_matrix_design = NULL,
                                        proteomics_combat_matrix_design_variables = NULL,
                                        transcriptomics_batch_correction = FALSE,
                                        transcriptomics_expression_filter_grouping = NULL,
                                        transcriptomics_expression_filter = 0.05,
                                        transcriptomics_selected_features = NULL,
                                        transcriptomics_limma_matrix_design = NULL,
                                        transcriptomics_limma_matrix_design_variables = NULL,
                                        transcriptomics_mad_gene_fraction = 0.75,
                                        metabolomics_selected_features = NULL,
                                        cytof_outlier_samples = NULL,
                                        cytof_msi_missing_fraction_cutoff = 0.5,
                                        cytof_msi_selected_features = NULL,
                                        cytof_msi_batch_correction = FALSE,
                                        cytof_msi_limma_matrix_design = NULL,
                                        cytof_msi_limma_matrix_design_variables = NULL,
                                        serum_rbd_abtiters_batch_correction = FALSE,
                                        preprocess_default_impute = TRUE,
                                        assays_to_process = NULL,
                                        run_pvca = TRUE,
                                        additional_pvca_variables = NULL
                                       ) {
  # List of preprocessing functions
  preprocess_funcs = list(
    "plasma_proteomics_targeted" = preprocess_proteomics,
    "plasma_proteomics_global_dda" = preprocess_proteomics,
    "serum_olink" = preprocess_olink,
    "nasal_viralload" = preprocess_nasal_viralload,
    "serum_rbd_abtiters" = preprocess_serum_rbd_abtiters,
    "serum_sarscov2_abtiters" = preprocess_serum_sarscov2_abtiters,
    "nasal_transcriptomics" = preprocess_transcriptomics,
    "plasma_metabolomics_global"= preprocess_metabolomics,
    "bld_cytof" = preprocess_bld_cytof_child,
    "bld_msi_cytof" = preprocess_bld_cytof_msi,
    "ea_transcriptomics" = preprocess_transcriptomics,
    "pbmc_transcriptomics" = preprocess_transcriptomics,
    "ea_metagenomics" = preprocess_metagenomics,
    "nasal_metagenomics" = preprocess_metagenomics
  )

  # Load the missingness/imputation-related parameters into data_env for an easy access in the downstream functions.
  data_env$impute_alpha = impute_alpha
  data_env$random_seed = random_seed
  data_env$impute_cores = impute_cores
  data_env$proteomics_missing_fraction_cutoff = proteomics_missing_fraction_cutoff
  data_env$proteomics_combat_matrix_design = proteomics_combat_matrix_design
  data_env$proteomics_combat_matrix_design_variables = proteomics_combat_matrix_design_variables
  data_env$proteomics_selected_features = proteomics_selected_features
  data_env$proteomics_batch_correction = proteomics_batch_correction
  data_env$transcriptomics_limma_matrix_design = transcriptomics_limma_matrix_design
  data_env$transcriptomics_expression_filter = transcriptomics_expression_filter
  data_env$transcriptomics_expression_filter_grouping = transcriptomics_expression_filter_grouping
  data_env$transcriptomics_selected_features = transcriptomics_selected_features
  data_env$additional_pvca_variables = additional_pvca_variables
  data_env$transcriptomics_batch_correction = transcriptomics_batch_correction
  data_env$transcriptomics_mad_gene_fraction = transcriptomics_mad_gene_fraction
  data_env$transcriptomics_limma_matrix_design = transcriptomics_limma_matrix_design
  data_env$transcriptomics_limma_matrix_design_variables = transcriptomics_limma_matrix_design_variables
  data_env$metabolomics_selected_features = metabolomics_selected_features
  data_env$cytof_outlier_samples = cytof_outlier_samples
  data_env$cytof_msi_missing_fraction_cutoff = cytof_msi_missing_fraction_cutoff
  data_env$cytof_msi_selected_features = cytof_msi_selected_features
  data_env$cytof_msi_batch_correction = cytof_msi_batch_correction
  data_env$cytof_msi_limma_matrix_design = cytof_msi_limma_matrix_design
  data_env$cytof_msi_limma_matrix_design_variables = cytof_msi_limma_matrix_design_variables
  data_env$serum_rbd_abtiters_batch_correction = serum_rbd_abtiters_batch_correction
  data_env$preprocess_default_impute = preprocess_default_impute
  
  # Validate parameters
  # ----------------------
  # Preprocess all assays if assays_to_preprocess is NULL:
  if(is.null(assays_to_process)) {
    assays_to_process = names(data_env$assay_files_info)
  }
  
  # check that transcriptomics expression filter group is in clinical table
  if (!is.null(transcriptomics_expression_filter_grouping)) {
    if (!(transcriptomics_expression_filter_grouping %in% colnames(data_env$clinical_data))) {
      stop("Provided group ", transcriptomics_expression_filter_grouping, " is not present in the clinical data. Please choose another grouping or perform gene filtering without grouping.")
    }
  }
  
  # check that all the provided variables for pvca are in clinical table
  if (!is.null(additional_pvca_variables)) {
    if (!all(additional_pvca_variables %in% colnames(data_env$clinical_data))) {
      stop("Not all provided variables as additional_pvca_variables are in the clinical data table: ",
           setdiff(additional_pvca_variables, colnames(data_env$clinical_data)))
    }
  }

  # Initialize result tables:
  pvca_data = data.frame()
  preprocessed_data = list()
  pvca_plots = list()
  mad_plots = list()

  # Loop through each assay in assays_to_process:
  for(assay in assays_to_process ) {
    # Preprocess the data for the assays that have preprocessing functions defined.
    if( ! is.null(preprocess_funcs[[assay]]) ) {
      cat("Preprocessing", assay, "\n")
      # Preprocess the data
      preprocess_out = preprocess_funcs[[assay]](assay, data_env, run_pvca)
      preprocessed_data[[ assay ]] = preprocess_out$data
      if(!is.null(preprocess_out$pvca)) {
        pvca_plots[[assay]] = preprocess_out$pvca
      }
      # If PVCA was run, store in data.frame:
      if(!is.null(preprocess_out$pvca)) {
        p = preprocess_out$pvca$p1
        p$data$assay = paste0(assay, "_final")
        pvca_data = bind_rows(pvca_data, p$data)
        if(!is.null(preprocess_out$pvca$p0)) {
          p = preprocess_out$pvca$p0
          p$data$assay = paste0(assay, "_beforeBatchCorr")
          pvca_data = bind_rows(pvca_data, p$data)
        }
      }
      if(!is.null(preprocess_out$mad_plots)){
        mad_plots[[assay]] = preprocess_out$mad_plots
      }
    }
  }
  return(list(preprocessed_data=preprocessed_data, pvca_data = pvca_data, 
              pvca_plots = pvca_plots, mad_plots = mad_plots))
}

#' Pre-process Transcriptomics
#'
#' @param assay: assay name in data_env (e.g. "plasma_proteomics_targeted").
#' @description filter gene biotypes, remove genes not expressed in more than certain percentage of samples globally or according to outcome groups (default 5%),
#' perform VOOM normalization, and batch effect correction according to the phase and plate number. Finally scale data and center to mean of 0 and standard deviation of 1.
#' @param data_env: Data env object loaded with Load_IMPACC_datasets function.
#' @param run_pvca: logical whether to run pvca
#' @param allowed_gene_biotypes: Allowed gene biotypes in transcriptomics matrix. E.g. c("protein_coding")
#' @param transcriptomics_batch_correction: Logical whether to perform batch correction for plate and phase
#' @param transcriptomics_expression_filter: Proportion of samples to apply gene expression feature filtering globally or per outcome group (default 0.05, 5%).
#' @param transcriptomics_expression_filter_grouping: Outcome variable for grouping samples prior to filtering gene features.
#' @param transcriptomics_selected_features: Optional, a vector of strings with feature IDs to be selected for downstream analysis. If `NULL`, will instead perform missing fraction cutoff.
#' @param transcriptomics_limma_matrix_design: design model for limma removeBatcheffect as a string (e.g. ~outcome*event_type), the provided variables need to exist in the clinical file
#' @param transcriptomics_limma_matrix_design_variables: design model variables that are used in the design as a vector (e.g. c("outcome","event_type"))
#' @param transcriptomics_mad_gene_fraction: Proportion of genes to apply gene expression feature filtering globally or per outcome group (default 0.75, 75%).
#' @returns named list (data, plot) where data is the preprocessed matrix and
#' plot the PVCA plot on the preprocessed matrix
preprocess_transcriptomics = function(assay,
                                      data_env,
                                      run_pvca,
                                      allowed_gene_biotypes = c("protein_coding"),
                                      transcriptomics_batch_correction = data_env$transcriptomics_batch_correction,
                                      transcriptomics_expression_filter = data_env$transcriptomics_expression_filter,
                                      transcriptomics_expression_filter_grouping = data_env$transcriptomics_expression_filter_grouping,
                                      transcriptomics_selected_features = data_env$transcriptomics_selected_features,
                                      transcriptomics_mad_gene_fraction = data_env$transcriptomics_mad_gene_fraction,
                                      transcriptomics_limma_matrix_design = data_env$transcriptomics_limma_matrix_design,
                                      transcriptomics_limma_matrix_design_variables = data_env$transcriptomics_limma_matrix_design_variables
                                      ) {
  
  # Make list to store PVCA results (if run_pvca == TRUE)
  if(run_pvca){
    pvca_res <- list()
  } else {
    pvca_res <- NULL
  }
  
  # Get counts, rowfeature, and metadata:
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  rowfeature = data_env[[ data_env$assay_files_info[[assay]]$rowfeature_var_name ]]
  metadata = data_env[[ data_env$assay_files_info[[assay]]$metadata_var_name ]]
  clinical_data <- data_env$clinical_data
  
  # Generate genes_pc table (only protein coding genes)
  genes_pc <- rowfeature %>%
    # get rownames as an explicit column
    tibble::rownames_to_column(var = "gene_id") %>%
    # remove row names
    tibble::remove_rownames() %>%
    # keep only gene ids with valid hgnc symbols
    tidyr::drop_na(gene_name) %>%
    # keep only protein coding genes
    dplyr::filter(gene_biotype %in% allowed_gene_biotypes) %>%
    # remove redundant genes if any and make unique
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    # add gene ids as row names
    tibble::column_to_rownames(var = "gene_id")
  
  # get counts data only for pc genes and transpose:
  counts1 <- counts[,rownames(genes_pc)] %>% t()

  # create DGElist object
  dge_list <- edgeR::DGEList(counts = counts1, genes = genes_pc)
  dge_list <- edgeR::calcNormFactors(dge_list)

  if (is.null(transcriptomics_selected_features)) {
    # expressed gene (CPM >=1) should be present in at least 5% of samples (globally or per outcome group)
    cut.filter <- transcriptomics_expression_filter
    if (!is.null(transcriptomics_expression_filter_grouping)) {
      # filter for genes expressed (CPM>=1) in more than 5% of samples any group
      group = data_env$clinical_data[[transcriptomics_expression_filter_grouping]][ match( rownames(counts), data_env$clinical_data$sample_id ) ]
      agg_res = stats::aggregate(t(round(edgeR::cpm(dge_list$counts)) >= 1), by=list(group=group), mean)
      keepRows = colSums(agg_res[,-1] > cut.filter) > 0
    } else {
      # filter for genes expressed (CPM>1) in at least 5% of samples globally
      agg_res = colMeans(t(round(edgeR::cpm(dge_list$counts)) >= 1))
      keepRows = agg_res > cut.filter
    }
    curDGE <- dge_list[keepRows,]
    curDGE <- edgeR::calcNormFactors(curDGE)
    
    counts1_cpm <- edgeR::cpm(curDGE$counts)
    
    ### Select top P fraction of genes with highest MAD
    P <- transcriptomics_mad_gene_fraction
    
    # Using data in log space
    mad_vals = apply(log1p(counts1_cpm), 1, mad)
    mad_cutoff = min(tail(sort(mad_vals), P*length(mad_vals)))
    df_mad <- data.frame(
      gene_id = rownames(counts1_cpm), 
      mean_exp = rowMeans(log1p(counts1_cpm)),
      mad_val = mad_vals)
    p1 = ggplot(df_mad, aes(y=mad_val)) +
      geom_histogram(bins = 1000) +
      coord_flip() + ggtitle(assay) + 
      geom_hline(yintercept = mad_cutoff, color="red")
    
    p2 = ggplot(df_mad, aes(mean_exp, mad_val)) + 
      #geom_point(size=0.5) + 
      geom_hex(bins=100) +
      geom_smooth() + ggtitle(assay) +
      geom_hline(yintercept = mad_cutoff, color="red") + 
      scale_fill_viridis_c(option = "A")
    
    mad_plot <- ggarrange(p1,p2)
    print(mad_plot)
  } else {
    curDGE <- edgeR::calcNormFactors(dge_list)
  }

  # Voom normalization, conversion to log2CPM values
  voomCounts_1 <- limma::voom(curDGE, design = NULL, plot = FALSE, save.plot = FALSE)
  
  # Harmonize variables for batch correction
  if (assay=="nasal_transcriptomics"){
    metadata$median_cv_coverage <- metadata$`median_cv_coverage_(QCMetrics)`
  }
  metadata$plate <- metadata[, data_env$assay_files_info[[assay]]$metadata_plate]
  
  # Metadata needed for batch correction and PVCA
  metadata_batch_corr <- metadata %>% 
    dplyr::mutate(phase_plate = paste0(phase, "_", plate)) %>%
    dplyr::mutate_at(vars(phase_plate, phase, plate), list(factor)) %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::filter(sample_id %in% colnames(counts1)) %>%
    dplyr::arrange(match(sample_id, colnames(counts1))) %>%
    tidyr::drop_na(phase_plate, median_cv_coverage)
  
  if(transcriptomics_batch_correction){
    
    if (is.null(transcriptomics_limma_matrix_design) && is.null(transcriptomics_limma_matrix_design_variables)) {
      # If no design is specified batch correction without design Matrix is performed.
      cat("\nWarning: performing batch correction without design matrix.\n")
    } else if (!is.null(transcriptomics_limma_matrix_design) && !is.null(transcriptomics_limma_matrix_design_variables)) {
      cat("\nWarning: performing batch correction with design matrix. Samples with NA values in : ", 
          transcriptomics_limma_matrix_design_variables, "will be dropped.\n")
      
      metadata_batch_corr <- metadata_batch_corr %>% 
        #Remove the phase annotated at clinical data, as we are using the annotations in the assay metadata level
        dplyr::left_join(clinical_data %>% 
                          dplyr::select(-phase) %>% 
                          dplyr::select(all_of(append(transcriptomics_limma_matrix_design_variables, c("sample_id")))),
                          by="sample_id") %>%
        dplyr::arrange(match(sample_id, colnames(counts1))) %>%
        tidyr::drop_na(all_of(transcriptomics_limma_matrix_design_variables))
      
      designMat <- model.matrix(eval(parse(text=transcriptomics_limma_matrix_design)), 
                                data = metadata_batch_corr)
    } else {
      stop("\nError: Both `transcriptomics_limma_matrix_design` and `transcriptomics_limma_matrix_design_variables` need to be provided, or none of them.")
    }
    
    
    # Prepare variables for the batch-correction
    phase_plate <- metadata_batch_corr$phase_plate
    medCV <- metadata_batch_corr$median_cv_coverage
  }
  
  # Run PVCA on pre-batch corrected transcriptomics:
  if(run_pvca){
    p0 = run_pvca(assay, data_env, 
                  preprocessed_matrix = t(voomCounts_1$E), 
                  plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                  core_site_metadata_column = "Core", 
                  plot_title = paste0(assay, "- before batch correction") )
    print(p0)
    pvca_res$p0 <- p0
  }

  curDGE_noNA <- curDGE[, metadata_batch_corr$sample_id]
  if (!(all(rownames(curDGE_noNA$samples) == metadata_batch_corr$sample_id))){
    stop("Error: make sure samples are ordered consistently between counts data and metadata. \n")
  } else {
    cat("Samples are ordered consistently between counts data and metadata.\n")
  }
  
  # compute voom counts
  voomCounts_batch_corr_input <- limma::voom(curDGE_noNA, design = NULL, 
                                             plot = FALSE, save.plot = FALSE)
  
  if(transcriptomics_batch_correction){
    if (!is.null(transcriptomics_limma_matrix_design)) {
      # Batch effect correction with design matrix
      exp_batch_corr_output <- limma::removeBatchEffect(voomCounts_batch_corr_input,
                                                        batch = phase_plate,
                                                        covariates = medCV,
                                                        design = designMat)
    } else {
      # Batch effect correction without design matrix
      exp_batch_corr_output <- limma::removeBatchEffect(voomCounts_batch_corr_input,
                                                        batch = phase_plate,
                                                        covariates = medCV)
    }
    
    # Run PVCA on pre-batch corrected transcriptomics:
    if(run_pvca){
      p1 = run_pvca(assay, data_env, 
                    preprocessed_matrix = t(exp_batch_corr_output), 
                    plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                    core_site_metadata_column = "Core", 
                    plot_title = paste0(assay, " - after batch correction") )
      print(p1)
      pvca_res$p1 <- p1
    }
    
  }

  if(transcriptomics_batch_correction){
    # Processed data matrix
    processed_data <- t(exp_batch_corr_output)
  } else {
    processed_data <- t(voomCounts_batch_corr_input$E)
  }
  
  # Final feature selection
  if (is.null(transcriptomics_selected_features)) {
    # Filter genes based on the MAD cutoff from above.
    mad_genes <- df_mad %>%
    filter(mad_val >= mad_cutoff) %>%
    pull(gene_id)
  
    cat("MAD-based filtering retained ", round(length(mad_genes)/dim(df_mad)[1]*100, 2), "% genes.\n" )

    # Processed data matrix filtered by MAD threshold
    processed_data <- processed_data[, mad_genes ]
  } else {
    # Processed data matrix filtered by user-provided gene list
    processed_data <- processed_data[, transcriptomics_selected_features]
  }
  
  # Scale the processed (normalized and log-transformed) data
  # to 0 mean 1 SD. It is specific to MOFA, for DEG this is not desired.
  processed_data <- apply(processed_data, 2, scale) %>% 
                          data.frame(check.names = F) %>% 
                          `rownames<-`(row.names(processed_data))
  if (is.null(transcriptomics_selected_features)){
    return(list(data=processed_data, pvca=pvca_res, mad_plots=mad_plot))
  } else {
    return(list(data=processed_data, pvca=pvca_res))
  }
  
}


#' Pre-process Proteomics
#' 
#' @param assay: assay name in data_env (e.g. "plasma_proteomics_targeted").
#' @param data_env: Data env object loaded with Load_IMPACC_datasets function.
#' @param run_pvca: logical whether to run pvca
#' @param proteomics_missing_fraction_cutoff: fraction of missing values allowed per feature
#' @param proteomics_selected_features: Optional, a vector of strings with feature IDs to be selected for downstream analysis. If `NULL`, will instead perform missing fraction cutoff.
#' @param proteomics_batch_correction: logical whether to perform batch correction for plate and phase
#' @param proteomics_combat_matrix_design: design matrix for batch correction
#' @param proteomics_combat_matrix_design_variables: variables to be used in the design matrix
#' @returns named list (data, plot) where data is the preprocessed matrix and 
#' plot the PVCA plot on the preprocessed matrix

preprocess_proteomics <- function(assay, 
                                  data_env, 
                                  run_pvca = TRUE,
                                  proteomics_missing_fraction_cutoff = data_env$proteomics_missing_fraction_cutoff,
                                  proteomics_selected_features = data_env$proteomics_selected_features,
                                  proteomics_batch_correction = data_env$proteomics_batch_correction,
                                  proteomics_combat_matrix_design = data_env$proteomics_combat_matrix_design,
                                  proteomics_combat_matrix_design_variables = data_env$proteomics_combat_matrix_design_variables
                                  ){
  
  # Extract counts + clinical_data:
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  metadata = data_env[[ data_env$assay_files_info[[assay]]$metadata_var_name ]]
  
  if (is.null(proteomics_selected_features)) {

    # Remove features that have more than "proteomics_missing_fraction_cutoff" missing values
    discard_features = apply( is.na( counts ), 2, mean) > proteomics_missing_fraction_cutoff

    if( sum(discard_features) > 0 ) {
      warning(
      "Discarding the following ",
      sum(discard_features),
      " features of ",
      assay,
      ", because they have > ",
      proteomics_missing_fraction_cutoff * 100,
      "% missing values:\n"
    )
    warning(paste(colnames(counts)[discard_features], collapse = ","), "\n")
    counts <- counts[,!discard_features]
    }

  } else {
    counts <- counts[, proteomics_selected_features]
    # Checking if any features are left with all NAs
    discard_features = apply( is.na( counts ), 2, mean) == 1
    print(paste0("After selecting the provided features, there are ", sum(discard_features), " features with all NAs:"))
    print(discard_features[discard_features])
  }

  counts <- 2^counts #unlog2 the counts
  counts_copy <- counts #copy. This is the one we will change and return
  
  NormalisationFactor <- median(rowSums(counts, na.rm = TRUE)) #calculate row(sample) sums. Take the median of all sums.
  
  #run for loop over each over the rows (samples)
  for(i in seq_len(nrow(counts))){
    #get sample i, apply normalisation factor on it. Change sample i in the normalized df.
    counts_copy[i,] <- counts[i,] * (NormalisationFactor / sum(counts[i,], na.rm = TRUE))
  }

  #Impute with half min value per Protein
  for(i in seq_len(ncol(counts_copy))){       # for-loop over columns
    halfmin <- min(counts_copy[,i], na.rm = TRUE)/2
    if (is.infinite(halfmin)) { halfmin <- 0 }
    counts_copy[,i][is.na(counts_copy[,i])] <- halfmin
  }
  
  # Log-transform data
  counts_copy <- log1p(counts_copy)

  metadata$plate <- metadata[, data_env$assay_files_info[[assay]]$metadata_plate]
  metadata_batch_corr <- metadata %>% 
                          dplyr::mutate(phase_plate = paste0(phase, "_", plate)) %>%
                          dplyr::mutate_at(vars(phase_plate, phase, plate), list(factor)) %>%
                          tibble::rownames_to_column("sample_id") %>%
                          dplyr::filter(sample_id %in% row.names(counts_copy)) %>%
                          dplyr::arrange(match(sample_id, row.names(counts_copy))) %>%
                          tidyr::drop_na(phase_plate)

  if (proteomics_batch_correction) {
    if (is.null(proteomics_combat_matrix_design) && is.null(proteomics_combat_matrix_design_variables)) {
      # If no design is specified batch correction without design Matrix is performed.
      cat("\nWarning: performing batch correction without design matrix.\n")
    } else if (!is.null(proteomics_combat_matrix_design) && !is.null(proteomics_combat_matrix_design_variables)) {
      cat("\nWarning: performing batch correction with design matrix. Samples with NA values in : ", 
          proteomics_combat_matrix_design_variables, "will be dropped.\n")
      
      metadata_batch_corr <- metadata_batch_corr %>% 
        #Remove the phase annotated at clinical data, as we are using the annotations in the assay metadata level
        dplyr::left_join(clinical_data %>% 
                          dplyr::select(-phase) %>% 
                          dplyr::select(all_of(append(proteomics_combat_matrix_design_variables, c("sample_id")))),
                        by="sample_id") %>%
        dplyr::arrange(match(sample_id, row.names(counts_copy))) %>%
        tidyr::drop_na(all_of(proteomics_combat_matrix_design_variables))
      
      designMat <- model.matrix(eval(parse(text=proteomics_combat_matrix_design)), 
                                data = metadata_batch_corr)
    } else {
      stop("\nError: Both `proteomics_combat_matrix_design` and `proteomics_combat_matrix_design_variables` need to be provided, or none of them.")
    }
  }

  # Prepare variables for the batch-correction
  phase_plate <- metadata_batch_corr$phase_plate
  
  # Run PVCA on pre-batch corrected matrix:
  if(run_pvca){
    pvca_res <- list()
    p0 = run_pvca(assay, data_env, 
                  preprocessed_matrix = counts_copy, 
                  plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                  core_site_metadata_column = NULL, 
                  plot_title = paste0(assay, "- before batch correction") )
    
    print(p0)
    pvca_res$p0 <- p0
  } else {
    pvca_res <- NULL
  }

  if (proteomics_batch_correction){
    if (!is.null(proteomics_combat_matrix_design)) {
      # Batch effect correction with design matrix
      exp_batch_corr_counts <- sva::ComBat(t(counts_copy),
                                            batch = phase_plate,
                                            mod = designMat)
    } else {
      # Batch effect correction without design matrix
      exp_batch_corr_counts <- sva::ComBat(t(counts_copy),
                                            batch = phase_plate)
    }

    # Run PVCA after batch correction
    if(run_pvca){
        p1 = run_pvca(assay, 
                      data_env, 
                      preprocessed_matrix = t(exp_batch_corr_counts), 
                      plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                      core_site_metadata_column = NULL, 
                      plot_title = paste0(assay, "- after batch correction") )
        print(p1)
        pvca_res$p1 <- p1
    } else {
      pvca_res = NULL
    }
    counts_copy <- t(exp_batch_corr_counts)
  }

  # Scale data
  counts_return <- counts_copy %>% 
    # Scale by columns (proteins):
    apply(2, scale) %>% 
    # coerce into data.frame:
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(counts_copy))
  
  return(list(data=counts_return, pvca=pvca_res))
}

#' Pre-process Metabolomics
#' 
#' @param assay: assay name in data_env (e.g. "plasma_metabolomics_global").
#' @param data_env: Data env object loaded with Load_IMPACC_datasets function.
#' @param metabolomics_selected_features: Optional, a vector of strings with feature IDs to be selected for downstream analysis. If `NULL`, will instead perform IQR feature selection.
#' @param run_pvca: logical whether to run pvca
#' @returns named list (data, plot) where data is the preprocessed matrix and 
#' plot the PVCA plot on the preprocessed matrix
preprocess_metabolomics <- function(assay,
                                    data_env,
                                    metabolomics_selected_features = data_env$metabolomics_selected_features,
                                    run_pvca = TRUE
                                    ){
  
  # Get counts and row_feature from data_env:
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  rowfeature = data_env[[ data_env$assay_files_info[[assay]]$rowfeature_var_name ]]
  
  # Copy counts and rowfeature:
  counts_return <- counts #copy. This is the one we will change and return
  rowfeature_return <- rowfeature #copy. This is the one we will change and return
  
  # Impute the NA values by half-min value for metabolite:
  # Loop over columns (metabolites):
  for(i in 1:ncol(counts_return)) {
    halfmin <- min(counts_return[,i], na.rm = TRUE)/2
    counts_return[,i][is.na(counts_return[,i])] <- halfmin
  }
  
  if (is.null(metabolomics_selected_features)) {
    # Removing non-Xenobiotic annotated features with IQR == 0
    iqrs = apply(counts_return, 2, stats::IQR)
    if(ncol(counts_return) == nrow(rowfeature_return)) {
      ft_to_remove = iqrs == 0 & 
        ! is.na(rowfeature_return$SUPER_PATHWAY) &
        rowfeature_return$SUPER_PATHWAY != "Xenobiotics"
      counts_return = counts_return[, ! ft_to_remove]
      rowfeature_return = rowfeature_return[! ft_to_remove, ]
    }
    warning("IQR-based filtering retained ", round(ncol(counts_return)/ncol(counts)*100,2), "% features." )

  } else {
    counts_return <- counts_return[, metabolomics_selected_features]
  }
  
  # Perform pareto-scaling (divide by sqrt sd):
  pareto_scaling <- function(x) {(x - mean(x, na.rm = T))/sqrt(sd(x, na.rm = T))}
  
  counts_return <- counts_return %>%
    apply(2, pareto_scaling) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(counts_return))
  
  # Run PVCA analysis:
  if(run_pvca){
    p1 = run_pvca(assay, 
                  data_env, 
                  preprocessed_matrix = counts_return, 
                  plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                  core_site_metadata_column = NULL, 
                  plot_title = paste0(assay) )
    print(p1)
    pvca_res = list(p1=p1)
  } else {
    pvca_res = NULL
  }

  return(list(data=counts_return, pvca=pvca_res))
}

#' Pre-process Olink
#' 
#' @param assay: assay name in data_env (e.g. "serum_olink").
#' @param data_env: Data env object loaded with Load_IMPACC_datasets function.
#' @param run_pvca: logical whether to run pvca
#' @returns preprocessed counts 
preprocess_olink <- function(assay,
                             data_env,
                             run_pvca = TRUE
                             ) {
  # get counts from data_env for OLINK:
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]

  # Remove samples that have all NA values
  counts_processed <- counts[rowMeans(is.na(counts)) < 1, ] %>%
    as.matrix()
  
  # Run PVCA
  if(run_pvca){
    #PVCA cannot handle missing values
    imputed_counts <- impute::impute.knn(t(counts_processed))$data
    p1 = run_pvca(assay, data_env, 
                    preprocessed_matrix = t(imputed_counts), 
                    plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                    core_site_metadata_column = "comment", 
                    plot_title = paste0(assay))
    print(p1)
    pvca_res = list(p1=p1)
  } else {
    pvca_res = NULL
  }
  
  # Scaling the data (log-transformation is not performed since there are negative values)
  counts_processed <- apply(counts_processed, 2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(counts_processed))
  return(list(data=counts_processed, pvca=pvca_res))
}


# Cytof
#' @title preprocesseed blood cytof children populations
#' @description remove debris, red blood cells and undefined populations and normalize on total number of events per sample 
#' @author Slim Fourati
#' @author Brian Lee
#' @author Jingjing Qi
#' @author Ravi Patel
#' @param assay: assay name in data_env (e.g. "bld_cytof" or "ea_cytof")
#' @param data_env: the data environment generated using load_IMPACC_datasets function
#' @param run_pvca: logical whether to run pvca
preprocess_bld_cytof_child <- function(assay,
                                 data_env,
                                 run_pvca = TRUE,
                                 cytof_outlier_samples = data_env$cytof_outlier_samples
                                 ) {

  if (!is.null(cytof_outlier_samples)){
    outlier_samples <- cytof_outlier_samples
  } else {
    outlier_samples <- c()
  }
  
  # Load counts data:
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  output_counts <-  counts %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::filter(!sample_id %in% outlier_samples) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Tier1_Undefined`,
                  -`Tier2_undefined_Undefined`) %>%
    tibble::column_to_rownames(var = "sample_id")
  output_counts <- output_counts/rowSums(output_counts)
  
  return(list(data=output_counts, pvca=pvca_res))
}

# Cytof MSI
#' @title preprocesseed blood cytof MSI data
#' @description remove debris, red blood cells and undefined populations and normalize on total number of events per sample 
#' @author Slim Fourati
#' @author Gisela Gabernet
#' @param assay: assay name in data_env (e.g. "bld_msi_cytof")
#' @param data_env: the data environment generated using load_IMPACC_datasets function
#' @param run_pvca: logical whether to run pvca
#' @param cytof_outlier_samples: vector of sample ids to be removed from the analysis
#' @param cytof_msi_batch_correction: logical whether to perform batch correction for plate and phase
#' @param cytof_msi_missing_fraction_cutoff: fraction of missing values allowed per feature
#' @param cytof_msi_selected_features: Optional, a vector of strings with feature IDs to be selected for downstream analysis. If `NULL`, will instead perform missing fraction cutoff.
#' @param cytof_msi_limma_matrix_design: design matrix for batch correction
#' @param cytof_msi_limma_matrix_design_variables: variables to be used in the design matrix
preprocess_bld_cytof_msi <- function(assay,
                                 data_env,
                                 run_pvca = TRUE,
                                 cytof_outlier_samples = data_env$cytof_outlier_samples,
                                 cytof_msi_batch_correction = data_env$cytof_msi_batch_correction,
                                 cytof_msi_missing_fraction_cutoff = data_env$cytof_msi_missing_fraction_cutoff,
                                 cytof_msi_selected_features = data_env$cytof_msi_selected_features,
                                 cytof_msi_limma_matrix_design = data_env$cytof_msi_limma_matrix_design,
                                 cytof_msi_limma_matrix_design_variables = data_env$cytof_msi_limma_matrix_design_variables
                                 ) {

  if (!is.null(cytof_outlier_samples)){
    outlier_samples <- cytof_outlier_samples
  } else {
    outlier_samples <- c()
  }
  
  # Load counts data:
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  output_counts <-  counts %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::filter(!sample_id %in% outlier_samples) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -contains(match="Platelets"),
                  -contains(match="RBC"),
                  -contains(match="Tier1_Undefined"),
                  -contains(match="Tier2_undefined_Undefined"),
                  -contains(match="193Ir_DNA"),
                  -contains(match="191Ir_DNA"),
                  -contains(match="192Os"),
                  -contains(match="beadDist"),
                  -contains(match="SampleID")
                  ) %>%
    tibble::column_to_rownames(var = "sample_id")

  #Filter out features with more than specified missing values,
  if (is.null(cytof_msi_selected_features)) {
    filter_counts <- output_counts[, colMeans(is.na(output_counts)) < cytof_msi_missing_fraction_cutoff]
  } else {
    filter_counts <- output_counts[, cytof_msi_selected_features]
  }

  #Filter out samples with more than 80% missing values
  filter_counts <- filter_counts[rowMeans(is.na(filter_counts)) < 0.8, ]
  filter_counts <- as.matrix(filter_counts)

  # Scale the data
  scaled_counts <- scale(filter_counts)
  
  # Run PVCA
  if(run_pvca){
    #PVCA cannot handle missing values
    imputed_counts <- t(impute::impute.knn(t(scaled_counts))$data)

    pvca_res = list()
    p0 = run_pvca(assay, data_env,
                  preprocessed_matrix = imputed_counts,
                  plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                  core_site_metadata_column = "Core",
                  plot_title = paste0(assay))
    print(p0)
    pvca_res$p0 = p0
  } else {
    pvca_res = NULL
  }

  if(cytof_msi_batch_correction){

    # Load metadata
    metadata = data_env[[ data_env$assay_files_info[[assay]]$metadata_var_name ]]
    metadata$plate <- metadata[, data_env$assay_files_info[[assay]]$metadata_plate]
    metadata_batch_corr <- metadata %>% 
                          dplyr::mutate(phase_plate = paste0(phase, "_", plate)) %>%
                          dplyr::mutate_at(vars(phase_plate, phase, plate), list(factor)) %>%
                          tibble::rownames_to_column("sample_id") %>%
                          dplyr::filter(sample_id %in% row.names(imputed_counts)) %>%
                          dplyr::arrange(match(sample_id, row.names(imputed_counts))) %>%
                          tidyr::drop_na(phase_plate)
    
    # Load clinical data
    clinical_data <- data_env$clinical_data
    
    # Generate model matrix for batch correction
    if (is.null(cytof_msi_limma_matrix_design) && is.null(cytof_msi_limma_matrix_design_variables)) {
      # If no design is specified batch correction without design Matrix is performed.
      cat("\nWarning: performing batch correction without design matrix.\n")
    } else if (!is.null(cytof_msi_limma_matrix_design) && !is.null(cytof_msi_limma_matrix_design_variables)) {
      cat("\nWarning: performing batch correction with design matrix. Samples with NA values in : ", 
          cytof_msi_limma_matrix_design_variables, "will be dropped.\n")
      
      metadata_batch_corr <- metadata_batch_corr %>% 
        #Remove the phase annotated at clinical data, as we are using the annotations in the assay metadata level
        dplyr::left_join(clinical_data %>% 
                           dplyr::select(-phase) %>% 
                           dplyr::select(all_of(append(cytof_msi_limma_matrix_design_variables, c("sample_id")))),
                         by="sample_id") %>%
        dplyr::arrange(match(sample_id, rownames(scaled_counts))) %>%
        tidyr::drop_na(all_of(cytof_msi_limma_matrix_design_variables))
      
      designMat <- model.matrix(eval(parse(text=cytof_msi_limma_matrix_design)), 
                                data = metadata_batch_corr)
    } else {
      stop("\nError: Both `cytof_msi_limma_matrix_design` and `cytof_msi_limma_matrix_design_variables` need to be provided, or none of them.")
    }
    
    # Prepare variables for the batch-correction
    phase_plate <- metadata_batch_corr$phase_plate

    if (!is.null(cytof_msi_limma_matrix_design)) {
      # Batch effect correction with design matrix
      exp_batch_corr_output <- limma::removeBatchEffect(t(scaled_counts),
                                                        batch = phase_plate,
                                                        design = designMat)
    } else {
      # Batch effect correction without design matrix
      exp_batch_corr_output <- limma::removeBatchEffect(t(scaled_counts),
                                                        batch = phase_plate)
    }
    
    # Run PVCA batch corrected data
    if(run_pvca){
      #PVCA cannot handle missing data
      imputed_counts <- impute::impute.knn(exp_batch_corr_output)$data
      p1 = run_pvca(assay, data_env, 
                    preprocessed_matrix = t(imputed_counts), 
                    plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                    core_site_metadata_column = "Core", 
                    plot_title = paste0(assay, " - after batch correction") )
      print(p1)
      pvca_res$p1 <- p1
    }

    return_counts = t(exp_batch_corr_output)
    
  } else {
    return_counts = scaled_counts
  }

  
  return(list(data=as.data.frame(return_counts), pvca=pvca_res))
}


#' Pre-process default function
#' 
#' @param assay: assay name in data_env (e.g. "").
#' @param data_env: Data env object loaded with Load_IMPACC_datasets function.
#' @param impute: whether to impute missing values
#' @param run_pvca: logical whether to run pvca
#' @returns preprocessed counts 
preprocess_default <- function(assay,
                               data_env,
                               impute=data_env$preprocess_default_impute,
                               run_pvca=FALSE
                               ) {
  # Load in counts from data_env for assay:
  #counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]

  if (impute) {
    # Imputation using missForest and save to data_env
    data_env <- missForest_imputation(assay, 
                                      data_env, 
                                      alpha = data_env$impute_alpha, 
                                      random_seed = data_env$random_seed, 
                                      cores = data_env$impute_cores)
  } else {
    message("Skipping imputation.")
  }

  # Return the imputed count data; the imputed count data frame is stored in data_env by missForest_imputation function
  output_counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]

  # Scaling the data (log-transformation is not performed since there could be negative values for some datasets.
  # If log-transformation is needed for a specific assay, create a new function for that assay and log-transform only for that assay)
  output_counts <- apply(output_counts, 2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(output_counts))

  if(run_pvca){
      pvca_res <- list()
      #PVCA cannot handle missing data
      imputed_counts <- impute::impute.knn(t(output_counts))$data
      p1 = run_pvca(assay, data_env, 
                    preprocessed_matrix = t(imputed_counts), 
                    plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                    core_site_metadata_column = "Core", 
                    plot_title = paste0(assay) )
      print(p1)
      pvca_res$p1 <- p1
  } else {
    pvca_res <- NULL
  }
  
  return( list(data=output_counts, pvca = pvca_res) )
}

# The following functions call preprocess_default:
# nasal_viralload
preprocess_nasal_viralload <- function(assay,
                                       data_env,
                                       run_pvca){
  preprocess_results = preprocess_default(assay, data_env, run_pvca)
  return(preprocess_results)
}

# serum_rbd_abtiters
preprocess_serum_rbd_abtiters <- function(assay,
                                          data_env,
                                          run_pvca,
                                          serum_rbd_abtiters_batch_correction = data_env$serum_rbd_abtiters_batch_correction){

  # Return the imputed count data; the imputed count data frame is stored in data_env by missForest_imputation function
  output_counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]

  # Log transform counts
  output_counts <- log2(output_counts)

  if(run_pvca){
      pvca_res <- list()
      #PVCA cannot handle missing data
      imputed_counts <- impute::impute.knn(t(output_counts))$data
      p0 = run_pvca(assay, data_env, 
                    preprocessed_matrix = t(imputed_counts), 
                    plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                    core_site_metadata_column = "Core", 
                    plot_title = paste0(assay) )
      print(p0)
      pvca_res$p0 <- p0
  } else {
    pvca_res <- NULL
  }

  if(serum_rbd_abtiters_batch_correction){

    # Load metadata
    metadata = data_env[[ data_env$assay_files_info[[assay]]$metadata_var_name ]]
    metadata$plate <- metadata[, data_env$assay_files_info[[assay]]$metadata_plate]
    metadata_batch_corr <- metadata %>% 
                          dplyr::mutate(phase_plate = paste0(phase, "_", plate)) %>%
                          dplyr::mutate_at(vars(phase_plate, phase, plate), list(factor)) %>%
                          tibble::rownames_to_column("sample_id") %>%
                          dplyr::filter(sample_id %in% row.names(output_counts)) %>%
                          dplyr::arrange(match(sample_id, row.names(output_counts))) %>%
                          tidyr::drop_na(phase_plate)
      
    # Prepare variables for the batch-correction
    phase_plate <- metadata_batch_corr$phase_plate

    # Batch effect correction without design matrix
    exp_batch_corr_output <- limma::removeBatchEffect(t(output_counts),
                                                      batch = phase_plate)
    
    # Run PVCA batch corrected data
    if(run_pvca){
      #PVCA cannot handle missing data
      imputed_counts <- impute::impute.knn(exp_batch_corr_output)$data
      p1 = run_pvca(assay, data_env, 
                    preprocessed_matrix = t(imputed_counts), 
                    plate_metadata_column = data_env$assay_files_info[[assay]]$metadata_plate,
                    core_site_metadata_column = "Core", 
                    plot_title = paste0(assay, " - after batch correction") )
      print(p1)
      pvca_res$p1 <- p1
    }

    return_counts = t(exp_batch_corr_output)
    
  } else {
    return_counts = output_counts
  }
  
  return( list(data=return_counts, pvca = pvca_res) )
}

# serum_sarscov2_abtiters
preprocess_serum_sarscov2_abtiters <- function(assay,
                                               data_env,
                                               run_pvca){
  preprocess_results = preprocess_default(assay, data_env, run_pvca, impute=dat_env$preprocess_default_impute)
  return(preprocess_results)
}

# metagenomics
preprocess_metagenomics <- function(assay,
                                    data_env,
                                    run_pvca){
  preprocess_results = preprocess_default(assay, data_env, run_pvca)
  return(preprocess_results)
}


############ OTHER NEEDED FUNCTIONS #########

#' Perform PVCA analysis
#' 
#' @param assay Assay name for which to perform PVCA. Choose from...
#' @param data_env Data env object loaded with Load_IMPACC_datasets
#' @param preprocessed_matrix 
#' @param plate_metadata_column Name of metadata column containing plate info
#' @param core_site_metadata_column Name of metadata column where core facility site is.
#' @param plot_title
#' @returns PVCA plot
#' @examples
#' plot <- run_pvca()
#' plot
run_pvca <- function(assay, 
                    data_env, 
                    preprocessed_matrix, 
                    plate_metadata_column,
                    additional_pvca_variables = data_env$additional_pvca_variables,
                    core_site_metadata_column = NULL, 
                    plot_title) {

  # Set VarCorr function to be from lme4 package:
  VarCorr <- lme4::VarCorr
  
  # Load counts, metadata, and plate information from data_env:
  assay_counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  assay_metadata = data_env[[ data_env$assay_files_info[[assay]]$metadata_var_name ]]
  assay_metadata$plate = assay_metadata[, plate_metadata_column]
  
  # Specify core site column if provided
  if(!is.null(core_site_metadata_column)) {
    assay_metadata$selected_core_site_column = assay_metadata[,core_site_metadata_column]
  }
  clinical_data <- data_env$clinical_data
  
  # specify which assay is going to be used as input to PVCA:
  pvca_input_matrix <- t(preprocessed_matrix)
  if (sum(is.na(pvca_input_matrix)) > 0) {
    stop("pvca input matrix should not have missing values")
  }
  pvca_input_raw <- t(assay_counts)
  pvca_input_metadata <- assay_metadata
  
  # select variable to be included in PVCA
  pvca_input_phenodata <- clinical_data %>%
    dplyr::select(any_of(append(additional_pvca_variables, 
                                c("sample_id", "event_type", "respiratory_status",
                                  "sex", "event_location",
                                  "discretized_admit_age_quantile", "ethnicity", 
                                  "participant_id", "race", "admit_month",
                                  "trajectory_group", "enrollment_site", "bmi"))))
  
  # step 4 append meta data to pvca phenodata
  #   core site: where are the sample processed (this chunk of code needs to be modified for each core assay)
  if(!is.null(core_site_metadata_column)) {
    pvca_input_metadata <- pvca_input_metadata %>%
      dplyr::mutate(core_site = gsub(pattern = "_.+", replacement = "", selected_core_site_column))
    # Print the core_site table
    print(table(pvca_input_metadata$core_site))
  }

  #   plate = interaction term of phase and plate
  pvca_input_phenodata <- pvca_input_metadata %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::mutate(plate     = interaction(phase, plate, drop = TRUE)) %>%
    dplyr::select( which(colnames(.) %in% c("sample_id", "plate", "phase", "core_site")) ) %>%
    merge(x     = pvca_input_phenodata,
          by    = "sample_id",
          all.x = TRUE)
  
  # append rownames
  pvca_input_phenodata <- pvca_input_phenodata %>%
    tibble::column_to_rownames(var = "sample_id")
  pvca_input_phenodata <- pvca_input_phenodata[colnames(pvca_input_matrix), , drop = FALSE]
  
 
  # Keep those phenodata columns for which there are more than one levels
  columns_to_keep = apply(pvca_input_phenodata, 2, function(x){length(na.omit(unique(x))) > 1} )
  pvca_input_phenodata = pvca_input_phenodata[, columns_to_keep]
  
  # run PVCA
  message("Running PVCA")
  fit <- pvca::PVCA(counts    = pvca_input_matrix,
              meta      = pvca_input_phenodata,
              inter     = FALSE,
              threshold = 0.6)
  message("Done running PVCA")
  
  # plot barplot with PVCA results
  pvca_barplot_data <- data.frame(explained = as.vector(fit),
                                  effect    = names(fit)) %>%
    dplyr::arrange(explained) %>%
    dplyr::mutate(effect = factor(effect, levels = effect))
  
  pvca_barplot <- ggplot2::ggplot(data    = pvca_barplot_data,
                        mapping = aes(x = explained, y = effect)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(aes(label = signif(explained, digits = 3)),
              nudge_y   = 0.01,
              size      = 3) +
    ggplot2::labs(x = NULL, y = "Proportion of the variance explained") + 
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                    vjust = 1,
                                    hjust = 1)) +
    ggplot2::ggtitle(plot_title)
    return(pvca_barplot)
}


#' Perform Missforest imputation
#' 
#' @author Leying Guan
#' @param assay_info Assay name for which to perform PVCA. Choose from...
#' @param data_env Data env object loaded with Load_IMPACC_datasets
#' @param alpha parameters fo 
#' @param plate_metadata_column Name of metadata column containing plate info
#' @param core_site_metadata_column Name of metadata column where core facility site is.
#' @param plot_title
#' @returns return data environment with imputed data
#' @examples
#' plot <- run_pvca()
#' plot
missForest_imputation <- function(assay_info,
                                  data_env,
                                  alpha = 0.2,
                                  random_seed = 2021,
                                  cores = 8 ){
  # Get assay_files_info
  assay_info <- data_env$assay_files_info[[assay]]

  # Load count, rowfeature, and count_desc from assay_info:
  count_df <- data_env[[ assay_info$count_var_name ]]
  rowfeature_df <- data_env[[ assay_info$rowfeature_var_name ]]
  count_desc <- assay_info$count_desc

  message(
    "Imputing the missing values (if any) in ",
    count_desc,
    " ...\n"
  )

  # Discard features that have "alpha" or more missing values.
  impute_df <- count_df
  discard_features <- apply(is.na(impute_df), 2, mean) >= alpha
  if( sum(discard_features) > 0 ) {
    warning(
      "Discarding the following ",
      sum(discard_features),
      " features of ",
      count_desc,
      ", because they have >= ",
      alpha * 100,
      "% missing values:\n"
    )
    warning(paste(colnames(impute_df)[discard_features], collapse = ","), "\n")
    impute_df <- impute_df[,!discard_features]
    rowfeature_df <- rowfeature_df %>% dplyr::filter(!discard_features)
  }
  
  # Perform the missing data imputation
  set.seed(random_seed)
  cores = min( cores, (ncol(count_df)-sum(discard_features)) )
  doParallel::registerDoParallel(cores = cores) # set based on number of CPU cores
  doRNG::registerDoRNG(seed = random_seed)
  missF_out <- missForest::missForest(impute_df,  verbose = FALSE, parallelize = "variables")
  count_df_wImputed <- missF_out$ximp
  
  data_env[[ assay_info$count_var_name ]] <- count_df_wImputed
  data_env[[ assay_info$rowfeature_var_name ]] <- rowfeature_df

  return(data_env)
}

#############################################
#### Longitudinal association functions
#############################################

pairwise_lme_modeling <- function(data, variable = "trajectory_group", fixedKnots = F, sex_age = T,
                                  knots = c(1,4,7,14,21), time_var = "event_date") {
  for(i in seq_along(unique(data$name))){
    ## retrieve comparisons
    active_data <- data[data$name == unique(data$name)[i], ]
    pairwise_comparisons <- combn(levels(active_data[[variable]]), 2)
    if(sex_age){
      compout <- matrix(NA, ncol(pairwise_comparisons), nrow = 4)
      rownames(compout) = c("p.slope", "p.intercept", "p.intercept.sex", "p.intercept.age.quantile")
    }else{
      compout <- matrix(NA, ncol(pairwise_comparisons), nrow = 2)
      rownames(compout) = c("p.slope", "p.intercept")
    }
    colnames(compout) <- paste0(pairwise_comparisons[1, ], "v", pairwise_comparisons[2, ])
    for(j in 1:ncol(pairwise_comparisons)){
      data_tmp = active_data[active_data[[variable]] == pairwise_comparisons[1, j] |
                               active_data[[variable]] == pairwise_comparisons[2, j],]
      if(!fixedKnots){
        if(sex_age){
          formula_use= formula(paste0("value ~ ",time_var," * ",variable, "+sex +discretized_admit_age_quantile"))
        }else{
          formula_use= formula(paste0("value ~ ",time_var," * ",variable))
        }
      }else{
        if(sex_age){
          formula_use <- formula(paste0("value ~ splines::bs(",time_var,", knots=c(",paste(knots, collapse = ',') ,")) * ",variable,"+sex +discretized_admit_age_quantile"))
        } else {
          formula_use <- formula(paste0("value ~ splines::bs(",time_var,", knots=c(",paste(knots, collapse = ',') ,")) * ",variable))
        }
      }
      lmfit= lme(fixed = formula_use, random = ~1|enrollment_site/participant_id, data = data_tmp)
      if(sex_age){
        compout[1,j] = (anova( lmfit, type = "marginal")$"p-value")[6]
        compout[2,j] = (anova( lmfit, type = "marginal")$"p-value")[3]
        compout[3,j] = (anova( lmfit, type = "marginal")$"p-value")[4] #sex
        compout[4,j] = (anova( lmfit, type = "marginal")$"p-value")[5] #age_quant
      }else{
        compout[1,j] = (anova( lmfit, type = "marginal")$"p-value")[4]
        compout[2,j] = (anova( lmfit, type = "marginal")$"p-value")[3]
      }
    }
    if(length(unique(data$name)) > 1){
      if(i ==1){
        compout_final <- list(compout)
        names(compout_final)[i] <- unique(data$name)[i]
      } else {
        compout_final <- c(compout_final, list(compout))
        names(compout_final)[i] <- unique(data$name)[i]
      }
    } else {
      compout_final <- compout
    } 
  }
  return(compout_final)
}

pairwise_mgcv_modeling <- function(data, age_sex = T, spline_type='cr', time_var = "event_date", variable = "trajectory_group"){
  for(i in seq_along(unique(data$name))){
    active_data <- data[data$name == unique(data$name)[i], ]
    levels_encode = sort(unique(active_data[[variable]]))
    pairwise_comparisons <- combn( levels_encode, 2)
    if(age_sex){
      compout = matrix(0, ncol = ncol( pairwise_comparisons), nrow = 4)
    }else{
      compout = matrix(0, ncol = ncol( pairwise_comparisons), nrow = 2)
    }
    for(j in 1:ncol(pairwise_comparisons)){
      tmp_data = active_data[active_data[[variable]]==pairwise_comparisons[1,j]|active_data[[variable]]==pairwise_comparisons[2,j],]
      tmp_data$participant_id = as.factor(tmp_data$participant_id )
      tmp_data$enrollment_site = as.factor(tmp_data$enrollment_site )
      tmp_data[[variable]] = ordered(as.factor(as.character(tmp_data[[variable]])), levels = c(as.character(pairwise_comparisons[1,j]), as.character(pairwise_comparisons[2,j])))
      if(age_sex){
        formula_use = formula(paste0("value~ s(",time_var,", bs = '",spline_type,"')+
              s(",time_var,", bs = '",spline_type,"', by =", variable, ")+", variable,
              "+ sex + discretized_admit_age_quantile"))
      }else{
        formula_use = formula(paste0("value~ s(",time_var,", bs = '",spline_type,"')+
              s(",time_var,", bs = '",spline_type,"', by =", variable, ")+",variable))
      }
      fit <- try(gamm4::gamm4(formula_use, data = tmp_data, random = ~(1|enrollment_site/participant_id)))
      a1 = anova(fit$gam)
      compout[2,j] = a1$pTerms.pv[1] # parametric terms pval of endpoint
      compout[1,j] = a1$s.table[2,4] # significance of smooth terms of endpoint
      if(age_sex){
        compout[3,j] = a1$pTerms.pv[2] #sex
        compout[4,j] = a1$pTerms.pv[3] #discretized_admit_age_quantile
      }
    }
    pairwise_comparisons = paste0(pairwise_comparisons[1, ], "v", pairwise_comparisons[2, ])
    colnames(compout) = pairwise_comparisons
    if(age_sex){
      rownames(compout) = c("p.slope", "p.intercept", "p.intercept.sex", "p.intercept.age.quantile")
    }else{
      rownames(compout) = c("p.slope", "p.intercept")
    }
    if(length(unique(data$name)) > 1){
      if(i==1){
        compout_final <- list(compout)
        names(compout_final)[i] <- unique(data$name)[i]
      } else {
        compout_final <- c(compout_final, list(compout))
        names(compout_final)[i] <- unique(data$name)[i]
      }
    } else {
      compout_final <- compout
    } 
  }
  return(compout_final)
}

pairwise_mgcv_modeling_nosex <- function(data, age = T, spline_type='cr', variable = "trajectory_group"){
  for(i in seq_along(unique(data$name))){
    active_data <- data[data$name == unique(data$name)[i], ]
    levels_encode = sort(unique(active_data[[variable]]))
    pairwise_comparisons <- combn( levels_encode, 2)
    if(age){
      compout = matrix(0, ncol = ncol( pairwise_comparisons), nrow = 3)
    }else{
      compout = matrix(0, ncol = ncol( pairwise_comparisons), nrow = 2)
    }
    for(j in 1:ncol(pairwise_comparisons)){
      tmp_data = active_data[active_data[[variable]]==pairwise_comparisons[1,j]|active_data[[variable]]==pairwise_comparisons[2,j],]
      tmp_data$participant_id = as.factor(tmp_data$participant_id )
      tmp_data$enrollment_site = as.factor(tmp_data$enrollment_site )
      tmp_data[[variable]] = ordered(as.factor(as.character(tmp_data[[variable]])), levels = c(as.character(pairwise_comparisons[1,j]), as.character(pairwise_comparisons[2,j])))
      if(age){
        formula_use = formula(paste0("value~ s(event_date, bs = '",spline_type,"')+
              s(event_date, bs = '",spline_type,"', by =", variable, ")+", variable,
              "+ discretized_admit_age_quantile"))
      }else{
        formula_use = formula(paste0("value~ s(event_date, bs = '",spline_type,"')+
              s(event_date, bs = '",spline_type,"', by =", variable, ")+",variable))
      }
      fit <- try(gamm4::gamm4(formula_use, data = tmp_data, random = ~(1|enrollment_site/participant_id)))
      a1 = anova(fit$gam)
      compout[2,j] = a1$pTerms.pv[1] # parametric terms pval of endpoint
      compout[1,j] = a1$s.table[2,4] # significance of smooth terms of endpoint
      if(age){
        compout[3,j] = a1$pTerms.pv[2] #discretized_admit_age_quantile
      }
    }
    pairwise_comparisons = paste0(pairwise_comparisons[1, ], "v", pairwise_comparisons[2, ])
    colnames(compout) = pairwise_comparisons
    if(age){
      rownames(compout) = c("p.slope", "p.intercept", "p.intercept.age.quantile")
    }else{
      rownames(compout) = c("p.slope", "p.intercept")
    }
    if(length(unique(data$name)) > 1){
      if(i==1){
        compout_final <- list(compout)
        names(compout_final)[i] <- unique(data$name)[i]
      } else {
        compout_final <- c(compout_final, list(compout))
        names(compout_final)[i] <- unique(data$name)[i]
      }
    } else {
      compout_final <- compout
    } 
  }
  return(compout_final)
}

pairwise_mgcv_modeling_batch <- function(data, age_sex = T, spline_type='cr', variable = "trajectory_group"){
  for(i in seq_along(unique(data$name))){
    active_data <- data[data$name == unique(data$name)[i], ]
    levels_encode = sort(unique(active_data[[variable]]))
    pairwise_comparisons <- combn( levels_encode, 2)
    if(age_sex){
      compout = matrix(0, ncol = ncol( pairwise_comparisons), nrow = 4)
    }else{
      compout = matrix(0, ncol = ncol( pairwise_comparisons), nrow = 2)
    }
    for(j in 1:ncol(pairwise_comparisons)){
      tmp_data = active_data[active_data[[variable]]==pairwise_comparisons[1,j]|active_data[[variable]]==pairwise_comparisons[2,j],]
      tmp_data$participant_id = as.factor(tmp_data$participant_id )
      tmp_data$enrollment_site = as.factor(tmp_data$enrollment_site )
      tmp_data[[variable]] = ordered(as.factor(as.character(tmp_data[[variable]])), levels = c(as.character(pairwise_comparisons[1,j]), as.character(pairwise_comparisons[2,j])))
      if(age_sex){
        formula_use = formula(paste0("value~ s(event_date, bs = '",spline_type,"')+
              s(event_date, bs = '",spline_type,"', by =", variable, ")+", variable,
              "+ sex + discretized_admit_age_quantile"))
      }else{
        formula_use = formula(paste0("value~ s(event_date, bs = '",spline_type,"')+
              s(event_date, bs = '",spline_type,"', by =", variable, ")+",variable))
      }
      fit <- try(gamm4::gamm4(formula_use, data = tmp_data, random = ~(1|enrollment_site/participant_id) + (1|phase/participant_id)))
      a1 = anova(fit$gam)
      compout[2,j] = a1$pTerms.pv[1] # parametric terms pval of endpoint
      compout[1,j] = a1$s.table[2,4] # significance of smooth terms of endpoint
      if(age_sex){
        compout[3,j] = a1$pTerms.pv[2] #sex
        compout[4,j] = a1$pTerms.pv[3] #discretized_admit_age_quantile
      }
    }
    pairwise_comparisons = paste0(pairwise_comparisons[1, ], "v", pairwise_comparisons[2, ])
    colnames(compout) = pairwise_comparisons
    if(age_sex){
      rownames(compout) = c("p.slope", "p.intercept", "p.intercept.sex", "p.intercept.age.quantile")
    }else{
      rownames(compout) = c("p.slope", "p.intercept")
    }
    if(length(unique(data$name)) > 1){
      if(i==1){
        compout_final <- list(compout)
        names(compout_final)[i] <- unique(data$name)[i]
      } else {
        compout_final <- c(compout_final, list(compout))
        names(compout_final)[i] <- unique(data$name)[i]
      }
    } else {
      compout_final <- compout
    } 
  }
  return(compout_final)
}

mgcv_global = function(inputDF, age_sex = T, spline_type = "cr", 
                       endpoints = "trajectory_group", time_var = "event_date"){
  features_names = sort(unique(inputDF$name))
  if(age_sex){
    res_table =   data.frame(matrix(NA, ncol = 8, nrow = length(features_names)))
    colnames(res_table) = c("p.slope", "p.intercept","p.intercept.sex", "p.intercept.age.quantile",
                            "adjp.slope", "adjp.intercept", "adjp.intercept.sex", 
                            "adjp.intercept.age.quantile")
  } else {
    res_table =   data.frame(matrix(NA, ncol = 4, nrow = length(features_names)))
    colnames(res_table) = c("p.slope", "p.intercept",
                            "adjp.slope", "adjp.intercept")
  }
  
  for(i in 1:length(features_names)){
    message(i)
    feature_name = features_names[i]
    tmp_input = inputDF[inputDF$name==feature_name,]
    tmp_input[[endpoints]] =ordered(as.factor(tmp_input[[endpoints]]))
    tmp_input$participant_id = as.factor(tmp_input$participant_id)
    if(age_sex){
      formula_use = formula(paste0("value~ s(",time_var,", bs = '",spline_type,"')+
              s(",time_var,", bs = 'cr', by =", endpoints, ")+", endpoints,
                                   "+ sex + discretized_admit_age_quantile"))
    }else{
      formula_use = formula(paste0("value~ s(",time_var,", bs = '",spline_type,"')+
              s(",time_var,", bs = 'cr', by =", endpoints, ")+",endpoints))
    }
    fit <- try(gamm4::gamm4(formula_use, data = tmp_input, random = ~(1|enrollment_site/participant_id)))
    if(class(fit)!="try-error"){
      a1 =anova(fit$gam)
      aa1 = sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2 = sum(a1$s.table[-c(1), 2])
      res_table[i,1]=pchisq(aa1,aa2, lower.tail = F)
      res_table[i,2] =a1$pTerms.pv[1]
      if(age_sex){
        res_table[i,3] =a1$pTerms.pv[2] ## sex
        res_table[i,4] =a1$pTerms.pv[3] ## discretized_admit_age_quantile
      }
    }
  }
  if (nrow(res_table) == 1) {
    res_table$adjp.intercept = res_table$p.intercept
    res_table$adjp.slope = res_table$p.slope
  } else {
    res_table$adjp.intercept = qvalue::qvalue(res_table$p.intercept, fdr.level = 0.05, pi0 = 1)$qvalues
    res_table$adjp.slope = qvalue::qvalue(res_table$p.slope, fdr.level = 0.05, pi0 = 1)$qvalues
  }

  if(age_sex){
      if (nrow(res_table) == 1) {
        res_table$adjp.intercept.sex = res_table$p.intercept.sex
        res_table$adjp.intercept.age.quantile = res_table$p.intercept.age.quantile
      } else {
        res_table$adjp.intercept.sex = qvalue::qvalue(res_table$p.intercept.sex, fdr.level = 0.05, pi0 = 1)$qvalues
        res_table$adjp.intercept.age.quantile = qvalue::qvalue(res_table$p.intercept.age.quantile, fdr.level = 0.05, pi0 = 1)$qvalues
      }
  }
  rownames(res_table) = features_names
  return(res_table)
}

mgcv_global_nosex = function(inputDF, age = T, spline_type = "cr", endpoints = "trajectory_group"){
  features_names = sort(unique(inputDF$name))
  if(age){
    res_table =   data.frame(matrix(NA, ncol = 6, nrow = length(features_names)))
    colnames(res_table) = c("p.slope", "p.intercept","p.intercept.age.quantile",
                            "adjp.slope", "adjp.intercept",  
                            "adjp.intercept.age.quantile")
  } else {
    res_table =   data.frame(matrix(NA, ncol = 4, nrow = length(features_names)))
    colnames(res_table) = c("p.slope", "p.intercept",
                            "adjp.slope", "adjp.intercept")
  }
  
  for(i in 1:length(features_names)){
    message(i)
    feature_name = features_names[i]
    tmp_input = inputDF[inputDF$name==feature_name,]
    tmp_input[[endpoints]] =ordered(as.factor(tmp_input[[endpoints]]))
    tmp_input$participant_id = as.factor(tmp_input$participant_id)
    if(age){
      formula_use = formula(paste0("value~ s(event_date, bs = '",spline_type,"')+
              s(event_date, bs = '",spline_type,"', by =", endpoints, ")+", endpoints,
                                   "+ discretized_admit_age_quantile"))
    }else{
      formula_use = formula(paste0("value~ s(event_date, bs = '",spline_type,"')+
              s(event_date, bs = '",spline_type,"', by =", endpoints, ")+",endpoints))
    }
    fit <- try(gamm4::gamm4(formula_use, data = tmp_input, random = ~(1|enrollment_site/participant_id)))
    if(class(fit)!="try-error"){
      a1 =anova(fit$gam)
      aa1 = sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2 = sum(a1$s.table[-c(1), 2])
      res_table[i,1]=pchisq(aa1,aa2, lower.tail = F)
      res_table[i,2] =a1$pTerms.pv[1]
      if(age){
        res_table[i,3] =a1$pTerms.pv[2] ## discretized_admit_age_quantile
      }
    }
  }
  if (nrow(res_table) == 1) {
    res_table$adjp.intercept = res_table$p.intercept
    res_table$adjp.slope = res_table$p.slope
  } else {
    res_table$adjp.intercept = qvalue::qvalue(res_table$p.intercept, fdr.level = 0.05, pi0 = 1)$qvalues
    res_table$adjp.slope = qvalue::qvalue(res_table$p.slope, fdr.level = 0.05, pi0 = 1)$qvalues
  }

  if(age){
      if (nrow(res_table) == 1) {
        res_table$adjp.intercept.age.quantile = res_table$p.intercept.age.quantile
      } else {
        res_table$adjp.intercept.age.quantile = qvalue::qvalue(res_table$p.intercept.age.quantile, fdr.level = 0.05, pi0 = 1)$qvalues
      }
  }
  rownames(res_table) = features_names
  return(res_table)
}

mgcv_global_batch = function(inputDF, age_sex = T, spline_type = "cr", endpoints = "trajectory_group"){
  features_names = sort(unique(inputDF$name))
  if(age_sex){
    res_table =   data.frame(matrix(NA, ncol = 8, nrow = length(features_names)))
    colnames(res_table) = c("p.slope", "p.intercept","p.intercept.sex", "p.intercept.age.quantile",
                            "adjp.slope", "adjp.intercept", "adjp.intercept.sex", 
                            "adjp.intercept.age.quantile")
  } else {
    res_table =   data.frame(matrix(NA, ncol = 4, nrow = length(features_names)))
    colnames(res_table) = c("p.slope", "p.intercept",
                            "adjp.slope", "adjp.intercept")
  }
  
  for(i in 1:length(features_names)){
    message(i)
    feature_name = features_names[i]
    tmp_input = inputDF[inputDF$name==feature_name,]
    tmp_input[[endpoints]] =ordered(as.factor(tmp_input[[endpoints]]))
    tmp_input$participant_id = as.factor(tmp_input$participant_id)
    if(age_sex){
      formula_use = formula(paste0("value~ s(event_date, bs = '",spline_type,"')+
              s(event_date, bs = 'cr', by =", endpoints, ")+", endpoints,
                                   "+ sex + discretized_admit_age_quantile"))
    }else{
      formula_use = formula(paste0("value~ s(event_date, bs = '",spline_type,"')+
              s(event_date, bs = 'cr', by =", endpoints, ")+",endpoints))
    }
    fit <- try(gamm4::gamm4(formula_use, data = tmp_input, random = ~(1|enrollment_site/participant_id) + (1|phase/participant_id)))
    if(class(fit)!="try-error"){
      a1 =anova(fit$gam)
      aa1 = sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2 = sum(a1$s.table[-c(1), 2])
      res_table[i,1]=pchisq(aa1,aa2, lower.tail = F)
      res_table[i,2] =a1$pTerms.pv[1]
      if(age_sex){
        res_table[i,3] =a1$pTerms.pv[2] ## sex
        res_table[i,4] =a1$pTerms.pv[3] ## discretized_admit_age_quantile
      }
    }
  }
  if (nrow(res_table) == 1) {
    res_table$adjp.intercept = res_table$p.intercept
    res_table$adjp.slope = res_table$p.slope
  } else {
    res_table$adjp.intercept = qvalue::qvalue(res_table$p.intercept, fdr.level = 0.05, pi0 = 1)$qvalues
    res_table$adjp.slope = qvalue::qvalue(res_table$p.slope, fdr.level = 0.05, pi0 = 1)$qvalues
  }

  if(age_sex){
      if (nrow(res_table) == 1) {
        res_table$adjp.intercept.sex = res_table$p.intercept.sex
        res_table$adjp.intercept.age.quantile = res_table$p.intercept.age.quantile
      } else {
        res_table$adjp.intercept.sex = qvalue::qvalue(res_table$p.intercept.sex, fdr.level = 0.05, pi0 = 1)$qvalues
        res_table$adjp.intercept.age.quantile = qvalue::qvalue(res_table$p.intercept.age.quantile, fdr.level = 0.05, pi0 = 1)$qvalues
      }
  }
  rownames(res_table) = features_names
  return(res_table)
}
 
#' @title Longitudinal model loop
#' @description Model loop: For a given input dataframe, creates a df of p.slope, p.intercept, their pairwise counterparts
#' and the respective adjusted p.values for all of these. It will loop over all the provided features / factors, and adjust for all the features / factors tested.
#' @author Cole Maguire
#' @author Ravi Patel
#' @author Leying Guan
#' @author Gisela Gabernet
#' @param inputDF: dataframe in long format, with mandatory columns: name (feature or factor name to be modeled), value (feature value), 
#' sample (sample ID), enrollment_site, event_date, participant_id, endpoint column (provide name to the endpoint parameter), 
#' and sex and admit_age if adjustment desired (age_sex param set to TRUE)
#' @param model_type: model type to apply. Available "lme", "smoothSpline", "fixedKnots".
#' @param age_sex: whether to adjust for age and sex. The `sex` and `admit_age` columns need to be provided in the dataframe if True.
#' @param knots: knots for the fixed spline model.
#' @param spline_type: spline type for the fixed spline model.
#' @param endpoint: endpoint column name.
#' @param time_var: variable to use as time values for association. Default is event_date (days after admission).
#' @param old_p_corrections: whether to use the old p-value correction method.
model_loop <- function(inputDF, model_type = "lme", spline_type = "cr", age_sex =TRUE, knots=c(1,2),
                       endpoint = "pro_2_groups", time_var = "event_date", old_p_corrections = FALSE){
  if (!model_type %in% c("lme", "smoothSpline", "fixedKnots")) {
    stop("model_type must be one of 'lme', 'smoothSpline', or 'fixedKnots'")
  }
  if (!spline_type %in% c("cr", "ts")) {
    stop("spline_type must be one of 'cr' or 'ts'")
  }
  inputDF$participant_id = as.factor( inputDF$participant_id )
  if(model_type == "smoothSpline"){
    #inputDF$event_date_transformed = inputDF$event_date
    lmDF <- mgcv_global(inputDF, endpoints = endpoint, spline_type = spline_type, 
                        age_sex = age_sex, time_var = time_var)
    lmDF$name <- rownames(lmDF)
    ## Because of the pairwise function requires there to be a column for "name"
    ## the rownames are copied to a new column here, ultimately this column is removed
  } else {
    lmDF <- inputDF %>%
      group_by(name) %>% 
      do(p.slope = {
        if(model_type == "lme"){
          if(age_sex){
            fit <- try(lme(fixed = as.formula(paste0("value ~ ", time_var, " * ", endpoint,
                            "+sex+discretized_admit_age_quantile")), data = .,
                           random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ ", time_var, " * ", endpoint)), data = .,
                           random =  ~1|enrollment_site/participant_id))
          }
        } else if (model_type == "fixedKnots"){
          if(age_sex){
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(",time_var,", 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                             "+sex+discretized_admit_age_quantile")), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(",time_var,", 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint)), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          }
        }
        if(class(fit)[1]!="try-error"){
            print(anova(fit, type = "marginal"))
            idx <- nrow(anova(fit, type = "marginal"))
            result = (anova(fit, type = "marginal")$"p-value")[idx]
        } else {
          result = NA
        }
      }, 
      p.intercept = {
        if(model_type == "lme"){
          if(age_sex){
            fit <- try(lme(fixed = as.formula(paste0("value ~ ",time_var," * ", endpoint,
                                                     "+sex+discretized_admit_age_quantile")), data = .,
                           random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ ",time_var," * ", endpoint)), data = .,
                           random =  ~1|enrollment_site/participant_id))
          }
        } else if (model_type == "fixedKnots"){
          if(age_sex){
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(",time_var,", 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                                                     "+sex+discretized_admit_age_quantile")), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(",time_var,", 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint)), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          }
          
        }
        if(class(fit)[1]!="try-error"){
          result = (anova(fit, type = "marginal")$"p-value")[3] #endpoint is position 3
        }else{
          result = NA
        }
      },
      p.intercept.sex = {
        if(age_sex){
          if(model_type == "lme"){
              fit <- try(lme(fixed = as.formula(paste0("value ~ ",time_var," * ", endpoint,
                                                       "+sex+discretized_admit_age_quantile")), data = .,
                             random =  ~1|enrollment_site/participant_id))
          } else if (model_type == "fixedKnots"){
              fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(",time_var,", 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                                                       "+sex+discretized_admit_age_quantile")), 
                             data = .,random =  ~1|enrollment_site/participant_id))
          }
          if(class(fit)[1]!="try-error"){
            result = (anova(fit, type = "marginal")$"p-value")[4] ## sex is position 4
          } else {
            result = NA
          }
        }
      },
      p.intercept.age.quantile = {
        if(age_sex){
          if(model_type == "lme"){
            fit <- try(lme(fixed = as.formula(paste0("value ~ ",time_var," * ", endpoint,
                                                     "+sex+discretized_admit_age_quantile")), data = .,
                           random =  ~1|enrollment_site/participant_id))
          } else if (model_type == "fixedKnots"){
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(",time_var,", 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                                                     "+sex+discretized_admit_age_quantile")), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          }
          if(class(fit)[1]!="try-error"){
            result = (anova(fit, type = "marginal")$"p-value")[5] ## discretized_admit_age_quantile is position 4
          } else {
            result = NA
          }
        }
      }
      ) %>% ungroup() %>%
      mutate(p.slope = unlist(p.slope),
             p.intercept  = unlist(p.intercept),
             p.intercept.sex  = unlist(p.intercept.sex),
             p.intercept.age.quantile  = unlist(p.intercept.age.quantile),
      ) %>% as.data.frame()
    lmDF$adjp.intercept=qvalue::qvalue(lmDF$p.intercept, fdr.level = 0.05, pi0 = 1)$qvalues
    lmDF$adjp.slope=qvalue::qvalue(lmDF$p.slope, fdr.level = 0.05, pi0 = 1)$qvalues
    lmDF$adjp.intercept.sex=qvalue::qvalue(lmDF$p.intercept.sex, fdr.level = 0.05, pi0 = 1)$qvalues
    lmDF$adjp.intercept.age.quantile=qvalue::qvalue(lmDF$p.intercept.age.quantile, fdr.level = 0.05, pi0 = 1)$qvalues
  }
  
  message("\n Initial calculation of p.slope and p.intercept complete. Beginning pairwise comparisons.")
  #### Pairwise comparisons ####
  ## Loop the pairwise comparisons ##
  for(i in seq_along(rownames(lmDF))){
    if(age_sex){
      tmp.out = matrix(NA,ncol = ncol(combn(unique(inputDF[[endpoint]]), 2)), nrow =4)
    }else{tmp.out = matrix(NA,ncol = ncol(combn(unique(inputDF[[endpoint]]), 2)), nrow =2)}

    try(
      if(model_type == "lme"){
        tmp.out <- pairwise_lme_modeling(inputDF[inputDF$name == lmDF$name[i],],
                                         variable = endpoint, fixedKnots = F, time_var=time_var)
      } else if(model_type == "fixedKnots"){
        tmp.out <- pairwise_lme_modeling(inputDF[inputDF$name == lmDF$name[i],],
                                         variable = endpoint, fixedKnots = T, knots = knots, time_var=time_var)
      } else if (model_type == "smoothSpline"){
        tmp.out <- pairwise_mgcv_modeling(inputDF[inputDF$name == lmDF$name[i],],
                                          variable = endpoint, spline_type = spline_type, time_var = time_var)
      })
    if(age_sex){
      if(i==1){
        values.out <- c(tmp.out[1,], tmp.out[2,], tmp.out[3,], tmp.out[4,])
        names(values.out) <- c(paste0(rownames(tmp.out)[1],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[2],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[3],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[4],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])))
      } else {
        values.out <- rbind(values.out, c(tmp.out[1,], tmp.out[2,], tmp.out[3,], tmp.out[4,]))
      }
    }else{
      if(i==1){
        values.out <- c(tmp.out[1,], tmp.out[2,])
        names(values.out) <- c(paste0(rownames(tmp.out)[1],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[2],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])))
      } else {
        values.out <- rbind(values.out, c(tmp.out[1,], tmp.out[2,]))
      }
    }

  }

  if(is.vector(values.out)) {
    tmp = matrix(values.out, nrow = 1)
    colnames(tmp) = names(values.out)
    values.out = tmp
  }

  ## Change to use qvalue like in tidyr loop
  if(old_p_corrections){
    p.slope.adj <- apply(values.out[,grepl("p.slope", colnames(values.out)), drop=FALSE], 2, function(x){
      qvalue::qvalue(x, fdr.level = 0.05, pi0 = 1)$qvalues
    })

    p.intercept.adj <- apply(values.out[,grepl("p.intercept", colnames(values.out)), drop=FALSE], 2, function(x){
      qvalue::qvalue(x, fdr.level = 0.05, pi0 = 1)$qvalues
    })
  } else {
    p.slope.adj <- as.data.frame(values.out[,grepl("p.slope", colnames(values.out)), drop=FALSE])
    if(nrow(p.slope.adj) > 1) {
      p.slope.adj <- qvalue::qvalue(as.vector(p.slope.adj[[1]]), fdr.level = 0.05, pi0 = 1)$qvalues
    }
    p.intercept.adj <- as.data.frame(values.out[,grepl("p.intercept", colnames(values.out)), drop=FALSE])
    if(nrow(p.intercept.adj) > 1) {
      p.intercept.adj <- qvalue::qvalue(as.vector(p.intercept.adj[[1]]), fdr.level = 0.05, pi0 = 1)$qvalues
    }
  }
  #colnames(p.slope.adj) <- gsub("_", ".adj_", colnames(p.slope.adj))
  #colnames(p.intercept.adj) <- gsub("_", ".adj_", colnames(p.intercept.adj))
  
  p.slope.adj <- data.frame(p_slope_adj = p.slope.adj)
  p.intercept.adj <- data.frame(p_intercept_adj = p.intercept.adj)

  lmDF_pairwise <- as.data.frame(cbind(lmDF, values.out, p.slope.adj, p.intercept.adj))
  rownames(lmDF_pairwise) <- lmDF_pairwise$name
  lmDF_pairwise <- lmDF_pairwise[,colnames(lmDF_pairwise) != "name"]
  return(lmDF_pairwise)
}

#' @title Longitudinal model loop with correction only by age
#' @description Model loop: For a given input dataframe, creates a df of p.slope, p.intercept, their pairwise counterparts
#' and the respective adjusted p.values for all of these. It will loop over all the provided features / factors, and adjust for all the features / factors tested.
#' @author Cole Maguire
#' @author Ravi Patel
#' @author Leying Guan
#' @author Gisela Gabernet
#' @param inputDF: dataframe in long format, with mandatory columns: name (feature or factor name to be modeled), value (feature value), 
#' sample (sample ID), enrollment_site, event_date, participant_id, endpoint column (provide name to the endpoint parameter), 
#' and sex and admit_age if adjustment desired (age_sex param set to TRUE)
#' @param model_type: model type to apply. Available "lme", "smoothSpline", "fixedKnots".
#' @param age: whether to adjust for age. The `admit_age` columns need to be provided in the dataframe if True.
#' @param knots: knots for the fixed spline model.
#' @param spline_type: spline type for the fixed spline model.
#' @param endpoint: endpoint column name.
#' @param old_p_corrections: whether to use the old p-value correction method.
model_loop_nosex <- function(inputDF, model_type = "lme", spline_type = "cr", age =TRUE, knots=c(1,2),
                       endpoint = "pro_2_groups", old_p_corrections = FALSE){
  if (!model_type %in% c("lme", "smoothSpline", "fixedKnots")) {
    stop("model_type must be one of 'lme', 'smoothSpline', or 'fixedKnots'")
  }
  if (!spline_type %in% c("cr", "ts")) {
    stop("spline_type must be one of 'cr' or 'ts'")
  }
  inputDF$participant_id = as.factor( inputDF$participant_id )
  if(model_type == "smoothSpline"){
    inputDF$event_date_transformed = inputDF$event_date
    lmDF <- mgcv_global_nosex(inputDF, endpoints = endpoint, spline_type = spline_type, age = age)
    lmDF$name <- rownames(lmDF)
    ## Because of the pairwise function requires there to be a column for "name"
    ## the rownames are copied to a new column here, ultimately this column is removed
  } else {
    lmDF <- inputDF %>%
      group_by(name) %>% 
      do(p.slope = {
        if(model_type == "lme"){
          if(age){
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint,
                                                     "+discretized_admit_age_quantile")), data = .,
                           random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint)), data = .,
                           random =  ~1|enrollment_site/participant_id))
          }
        } else if (model_type == "fixedKnots"){
          if(age){
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                                                     "+discretized_admit_age_quantile")), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint)), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          }
        }
        if(class(fit)[1]!="try-error"){
          #print(anova(fit, type = "marginal"))
          idx <- nrow(anova(fit, type = "marginal"))
          result = (anova(fit, type = "marginal")$"p-value")[idx]
        } else {
          result = NA
        }
      }, 
      p.intercept = {
        if(model_type == "lme"){
          if(age){
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint,
                                                     "+discretized_admit_age_quantile")), data = .,
                           random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint)), data = .,
                           random =  ~1|enrollment_site/participant_id))
          }
        } else if (model_type == "fixedKnots"){
          if(age){
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                                                     "+discretized_admit_age_quantile")), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint)), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          }
          
        }
        if(class(fit)[1]!="try-error"){
          result = (anova(fit, type = "marginal")$"p-value")[3] #endpoint is position 3
        }else{
          result = NA
        }
      },
      p.intercept.age.quantile = {
        if(age){
          if(model_type == "lme"){
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint,
                                                     "+discretized_admit_age_quantile")), data = .,
                           random =  ~1|enrollment_site/participant_id))
          } else if (model_type == "fixedKnots"){
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                                                     "+discretized_admit_age_quantile")), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          }
          if(class(fit)[1]!="try-error"){
            result = (anova(fit, type = "marginal")$"p-value")[4] ## discretized_admit_age_quantile is position 4
          } else {
            result = NA
          }
        }
      }
      ) %>% ungroup() %>%
      mutate(p.slope = unlist(p.slope),
             p.intercept  = unlist(p.intercept),
             p.intercept.age.quantile  = unlist(p.intercept.age.quantile),
      ) %>% as.data.frame()
    lmDF$adjp.intercept=qvalue::qvalue(lmDF$p.intercept, fdr.level = 0.05, pi0 = 1)$qvalues
    lmDF$adjp.slope=qvalue::qvalue(lmDF$p.slope, fdr.level = 0.05, pi0 = 1)$qvalues
    lmDF$adjp.intercept.age.quantile=qvalue::qvalue(lmDF$p.intercept.age.quantile, fdr.level = 0.05, pi0 = 1)$qvalues
  }
  
  message("\n Initial calculation of p.slope and p.intercept complete. Beginning pairwise comparisons.")
  #### Pairwise comparisons ####
  ## Loop the pairwise comparisons ##
  for(i in seq_along(rownames(lmDF))){
    if(age){
      tmp.out = matrix(NA,ncol = ncol(combn(unique(inputDF[[endpoint]]), 2)), nrow =3)
    }else{tmp.out = matrix(NA,ncol = ncol(combn(unique(inputDF[[endpoint]]), 2)), nrow =2)}
    
    try(
      if(model_type == "lme"){
        tmp.out <- pairwise_lme_modeling(inputDF[inputDF$name == lmDF$name[i],],
                                         variable = endpoint, fixedKnots = F)
      } else if(model_type == "fixedKnots"){
        tmp.out <- pairwise_lme_modeling(inputDF[inputDF$name == lmDF$name[i],],
                                         variable = endpoint, fixedKnots = T, knots = knots)
      } else if (model_type == "smoothSpline"){
        tmp.out <- pairwise_mgcv_modeling_nosex(inputDF[inputDF$name == lmDF$name[i],],
                                          variable = endpoint, spline_type = spline_type)
      })
    if(age){
      if(i==1){
        values.out <- c(tmp.out[1,], tmp.out[2,], tmp.out[3,])
        names(values.out) <- c(paste0(rownames(tmp.out)[1],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[2],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[3],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])))
      } else {
        values.out <- rbind(values.out, c(tmp.out[1,], tmp.out[2,], tmp.out[3,]))
      }
    }else{
      if(i==1){
        values.out <- c(tmp.out[1,], tmp.out[2,])
        names(values.out) <- c(paste0(rownames(tmp.out)[1],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[2],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])))
      } else {
        values.out <- rbind(values.out, c(tmp.out[1,], tmp.out[2,]))
      }
    }
    
  }
  
  if(is.vector(values.out)) {
    tmp = matrix(values.out, nrow = 1)
    colnames(tmp) = names(values.out)
    values.out = tmp
  }
  
  ## Change to use qvalue like in tidyr loop
  if(old_p_corrections){
    p.slope.adj <- apply(values.out[,grepl("p.slope", colnames(values.out)), drop=FALSE], 2, function(x){
      qvalue::qvalue(x, fdr.level = 0.05, pi0 = 1)$qvalues
    })
    
    p.intercept.adj <- apply(values.out[,grepl("p.intercept", colnames(values.out)), drop=FALSE], 2, function(x){
      qvalue::qvalue(x, fdr.level = 0.05, pi0 = 1)$qvalues
    })
  } else {
    p.slope.adj <- as.data.frame(values.out[,grepl("p.slope", colnames(values.out)), drop=FALSE])
    if(nrow(p.slope.adj) > 1) {
      p.slope.adj <- qvalue::qvalue(as.vector(p.slope.adj[[1]]), fdr.level = 0.05, pi0 = 1)$qvalues
    }
    p.intercept.adj <- as.data.frame(values.out[,grepl("p.intercept", colnames(values.out)), drop=FALSE])
    if(nrow(p.intercept.adj) > 1) {
      p.intercept.adj <- qvalue::qvalue(as.vector(p.intercept.adj[[1]]), fdr.level = 0.05, pi0 = 1)$qvalues
    }
  }
  #colnames(p.slope.adj) <- gsub("_", ".adj_", colnames(p.slope.adj))
  #colnames(p.intercept.adj) <- gsub("_", ".adj_", colnames(p.intercept.adj))
  
  p.slope.adj <- data.frame(p_slope_adj = p.slope.adj)
  p.intercept.adj <- data.frame(p_intercept_adj = p.intercept.adj)
  
  lmDF_pairwise <- as.data.frame(cbind(lmDF, values.out, p.slope.adj, p.intercept.adj))
  rownames(lmDF_pairwise) <- lmDF_pairwise$name
  lmDF_pairwise <- lmDF_pairwise[,colnames(lmDF_pairwise) != "name"]
  return(lmDF_pairwise)
}

#' @title Longitudinal model loop with batch as random effect
#' @description Model loop: For a given input dataframe, creates a df of p.slope, p.intercept, their pairwise counterparts
#' and the respective adjusted p.values for all of these. It will loop over all the provided features / factors, and adjust for all the features / factors tested.
#' @author Cole Maguire
#' @author Ravi Patel
#' @author Leying Guan
#' @author Gisela Gabernet
#' @param inputDF: dataframe in long format, with mandatory columns: name (feature or factor name to be modeled), value (feature value), 
#' sample (sample ID), enrollment_site, event_date, participant_id, endpoint column (provide name to the endpoint parameter), 
#' and sex and admit_age if adjustment desired (age_sex param set to TRUE)
#' @param model_type: model type to apply. Available "lme", "smoothSpline", "fixedKnots".
#' @param age_sex: whether to adjust for age and sex. The `sex` and `admit_age` columns need to be provided in the dataframe if True.
#' @param knots: knots for the fixed spline model.
#' @param spline_type: spline type for the fixed spline model.
#' @param endpoint: endpoint column name.
#' @param old_p_corrections: whether to use the old p-value correction method.
model_loop_batch <- function(inputDF, model_type = "lme", spline_type = "cr", age_sex =TRUE, knots=c(1,2),
                       endpoint = "pro_2_groups", old_p_corrections = FALSE){
  if (!model_type %in% c("lme", "smoothSpline", "fixedKnots")) {
    stop("model_type must be one of 'lme', 'smoothSpline', or 'fixedKnots'")
  }
  if (!spline_type %in% c("cr", "ts")) {
    stop("spline_type must be one of 'cr' or 'ts'")
  }
  inputDF$participant_id = as.factor( inputDF$participant_id )
  if(model_type == "smoothSpline"){
    inputDF$event_date_transformed = inputDF$event_date
    lmDF <- mgcv_global_batch(inputDF, endpoints = endpoint, spline_type = spline_type, age_sex = age_sex)
    lmDF$name <- rownames(lmDF)
    ## Because of the pairwise function requires there to be a column for "name"
    ## the rownames are copied to a new column here, ultimately this column is removed
  } else {
    stop("No other model supported for batch effect in random effects.")
  }
  
  message("\n Initial calculation of p.slope and p.intercept complete. Beginning pairwise comparisons.")
  #### Pairwise comparisons ####
  ## Loop the pairwise comparisons ##
  for(i in seq_along(rownames(lmDF))){
    if(age_sex){
      tmp.out = matrix(NA,ncol = ncol(combn(unique(inputDF[[endpoint]]), 2)), nrow =4)
    }else{tmp.out = matrix(NA,ncol = ncol(combn(unique(inputDF[[endpoint]]), 2)), nrow =2)}

    try(
      if (model_type == "smoothSpline"){
        tmp.out <- pairwise_mgcv_modeling_batch(inputDF[inputDF$name == lmDF$name[i],],
                                          variable = endpoint, spline_type = spline_type)
      } else {
        stop("No other model supported for batch effect in random effects.")
      })
    if(age_sex){
      if(i==1){
        values.out <- c(tmp.out[1,], tmp.out[2,], tmp.out[3,], tmp.out[4,])
        names(values.out) <- c(paste0(rownames(tmp.out)[1],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[2],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[3],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[4],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])))
      } else {
        values.out <- rbind(values.out, c(tmp.out[1,], tmp.out[2,], tmp.out[3,], tmp.out[4,]))
      }
    }else{
      if(i==1){
        values.out <- c(tmp.out[1,], tmp.out[2,])
        names(values.out) <- c(paste0(rownames(tmp.out)[1],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                               paste0(rownames(tmp.out)[2],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])))
      } else {
        values.out <- rbind(values.out, c(tmp.out[1,], tmp.out[2,]))
      }
    }

  }

  if(is.vector(values.out)) {
    tmp = matrix(values.out, nrow = 1)
    colnames(tmp) = names(values.out)
    values.out = tmp
  }

  ## Change to use qvalue like in tidyr loop
  if(old_p_corrections){
    p.slope.adj <- apply(values.out[,grepl("p.slope", colnames(values.out)), drop=FALSE], 2, function(x){
      qvalue::qvalue(x, fdr.level = 0.05, pi0 = 1)$qvalues
    })

    p.intercept.adj <- apply(values.out[,grepl("p.intercept", colnames(values.out)), drop=FALSE], 2, function(x){
      qvalue::qvalue(x, fdr.level = 0.05, pi0 = 1)$qvalues
    })
  } else {
    p.slope.adj <- as.data.frame(values.out[,grepl("p.slope", colnames(values.out)), drop=FALSE])
    if(nrow(p.slope.adj) > 1) {
      p.slope.adj <- qvalue::qvalue(as.vector(p.slope.adj[[1]]), fdr.level = 0.05, pi0 = 1)$qvalues
    }
    p.intercept.adj <- as.data.frame(values.out[,grepl("p.intercept", colnames(values.out)), drop=FALSE])
    if(nrow(p.intercept.adj) > 1) {
      p.intercept.adj <- qvalue::qvalue(as.vector(p.intercept.adj[[1]]), fdr.level = 0.05, pi0 = 1)$qvalues
    }
  }
  #colnames(p.slope.adj) <- gsub("_", ".adj_", colnames(p.slope.adj))
  #colnames(p.intercept.adj) <- gsub("_", ".adj_", colnames(p.intercept.adj))
  
  p.slope.adj <- data.frame(p_slope_adj = p.slope.adj)
  p.intercept.adj <- data.frame(p_intercept_adj = p.intercept.adj)

  lmDF_pairwise <- as.data.frame(cbind(lmDF, values.out, p.slope.adj, p.intercept.adj))
  rownames(lmDF_pairwise) <- lmDF_pairwise$name
  lmDF_pairwise <- lmDF_pairwise[,colnames(lmDF_pairwise) != "name"]
  return(lmDF_pairwise)
}

#' @title Plot model
#' @description Plot individual feature / factor values longitudinally with provided model. It supports only one feature at a time.
#' @author Cole Maguire
#' @author Ravi Patel
#' @author Leying Guan
#' @author Gisela Gabernet
#' @param plotDF: dataframe in long format, with mandatory columns: name (feature or factor name to be modeled - filtered to contain only one feature), value (feature value), 
#' sample (sample ID), enrollment_site, event_date, participant_id, endpoint column (provide name to the endpoint parameter), 
#' and sex and admit_age if adjustment desired (age_sex param set to TRUE).
#' @param model_loop: model loop output dataframe.
#' @param endpoint: endpoint column name.
#' @param model_type: model type to apply, select the same as the one used in the model loop function. Available: "lme", "fixedKnots", "smoothSpline".
#' @param spline_type: spline type for the spline model.
#' @param knots: knots for the fixed spline model.
#' @param age_sex: whether to correct for age and sex.
#' @param time_var: Variable to use as time variable (defaults to event_date).s
#' @param facet: whether to facet the plot by endpoint.
#' @param p_adjust: whether to adjust the p.values for multiple comparisons.
#' @param signif_markers: whether to add significance markers to the plot.
#' @param custom_signif_markers: custom significance markers to use.
#' @param remove_NS: whether to remove non-significant markers.
#' @param remove_NS_threshold: threshold for removing non-significant markers.
#' @param fit: provide previous fit.
#' @param font: font to use for the plot.
#' @param title_size: title size.
#' @param bar_height: bar height.
#' @param colors: colors to use for the plot.
#' @param title: plot title.
#' @param xlabel: x-axis label.
#' @param ylabel: y-axis label.
#' @param custom_theme_graph: custom theme for the plot.
#' @param custom_theme_bars: custom theme for the bars.
#' @param knot_lines: whether to add knot lines.
#' @param group_trendline: whether to add group trendline.
#' @param group_trendline_color: color for the group trendline.
#' @param group_trendline_line_width: line width for the group trendline.
#' @param group_trendline_dropout: whether to drop out the group trendline.
#' @param individual_trendlines: whether to add individual trendlines.
#' @param individual_trendline_size: size for the individual trendlines.
#' @param individual_trendline_alpha: alpha for the individual trendlines.
#' @param individual_trendlines_color: color for the individual trendlines.
#' @param individual_points: whether to add individual points.
#' @param individual_points_size: size for the individual points.
#' @param individual_points_alpha: alpha for the individual points.
#' @param individual_paths: whether to add individual paths.
#' @param individual_paths_size: size for the individual paths.
#' @param individual_paths_alpha: alpha for the individual paths.
#' @param sex_age_intercept_bars: whether to add intercept bars.
#' @param sex_age_intercept_x: x position for the intercept bars.
#' @param sex_bar_seperator_dist: distance for sex bar separator.
#' @param bar_size: bar size.
#' @param facet_space: space for the facet.
#' @param intercept_hist: intercept histogram.
#' @param sex_age_intercept_bars_colors: colors for the intercept bars.
#' @param CI: whether to add confidence intervals.
#' @param CI_alpha: alpha for the confidence intervals.
#' @param CI_interval: interval for the confidence intervals.
#' @param y_axis_reverse: whether to reverse the y-axis.
#' @param y_limits: y-axis limits.
#' @param x_limits: x-axis limits.
#' @param return_multi_obj: whether to return the plot object.
#' @param p_value_text_size: p-value text size.
plot_model <- function(plotDF, model_loop, endpoint = "trajectory_group", model_type = "lme",
                       spline_type = "cr", knots = c(1, 4, 7, 14, 21), age_sex = TRUE, age_only = FALSE, time_var = "event_date",
                       facet = TRUE, p_adjust = NA, signif_markers =TRUE, custom_signif_markers, 
                       remove_NS = FALSE, remove_NS_threshold = 0.05, fit, font, title_size =10, bar_height = 5,
                       colors, title, xlabel, ylabel, custom_theme_graph, custom_theme_bars, knot_lines = TRUE,
                       group_trendline = TRUE, group_trendline_color = "black", group_trendline_line_width = 1.5,
                       group_trendline_dropout = TRUE,
                       individual_trendlines = TRUE, individual_trendline_size = 0.6, 
                       individual_trendline_alpha = 0.9, individual_trendlines_color = "#5F5F5F",
                       individual_points = TRUE, individual_points_size = 1, individual_points_alpha = 0.8,
                       individual_paths = TRUE, individual_paths_size = 1, individual_paths_alpha = 0.4,
                       sex_age_intercept_bars = FALSE, sex_age_intercept_x = -3, sex_bar_seperator_dist = 2,
                       bar_size = 10, facet_space =1, intercept_hist = NULL, sex_age_intercept_bars_colors = NULL,
                       CI = FALSE, CI_alpha = 0.3, CI_interval = 0.95, y_axis_reverse = FALSE, y_limits = NULL, x_limits = NULL,
                       return_multi_obj = FALSE, p_value_text_size = 3.88){
  #### DATA CHECKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (!model_type %in% c("lme", "fixedKnots", "smoothSpline")){
    stop("Please choose a valid model type: lme, fixedKnots, or smoothSpline.")
  }
  if(length(unique(plotDF$name)) > 1){
    stop("Too many different variables in the name column! Please reduce to the one you want to plot.")
  } else if (!(unique(plotDF$name) %in% rownames(model_loop))){
    stop("Your provided model_loop does not contain the variable for plotting in plotDF.")
  }
  if(!(is.factor(plotDF[[endpoint]]))){
    stop("The endpoint you have selected is not a factor. Please convert it to a factor.")
  }
  if(age_sex){
    if(!is.factor(plotDF$sex)){
      stop("Please convert sex to a factor.")
    } else if (!is.factor(plotDF$discretized_admit_age_quantile)){
      stop("Please convert discretized_admit_age_quantile to a factor.")
    }
  }
  if(model_type == "smoothSpline"){
    if(!is.factor(plotDF$participant_id)){
      stop("Please convert participant_id to a factor.")
    # } else if (!is.ordered(plotDF[[endpoint]])){
    #   stop("Please convert your choosen endpoint to an ordered factor.")
    } else if (!is.factor(plotDF$enrollment_site)){
      stop("Please convert your choosen enrollment_site to a factor.")
    }
    if (!spline_type %in% c("cr", "ts")){
      stop("Please choose a valid spline type: cr or ts.")
    }
  }
  ## Labels to replace trajectory_group as appropriate, updated 1/25/22
  trajectory_labels <- c("TG1", "TG2",
                         "TG3",
                         "TG4", "TG5")
  # trajectory_labels <- c("1-Brief length \nof stay", "2-Intermediate length \nof stay",
  #                        "3-Intermediate stay with\ndischarge limitations",
  #                        "4-Prolonged \nhospitalization", "5-Fatal")
  names(trajectory_labels) <- c("1", "2", "3", "4", "5")
  
  ## Labels to replace PRO min groups
  pro_min_labels <- c("MIN","LC")
  names(pro_min_labels) <- c("MIN", "OTHER")
  
  ## Labels to replace PRO groups
  pro_labels <- c("MIN", "COG", "PHY", "MLT")
  names(pro_labels) <- c("MIN", "COG", "PHY", "MLT")
  
  # Time variable
  plotDF$time_var <- plotDF[[time_var]]
  
  #### SIGNIFICANCE BARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(!is.na(p_adjust)){
    comps_loop <- model_loop[rownames(model_loop) == unique(plotDF$name),-1:-8]
    comps_loop <- pivot_longer(comps_loop, cols = 1:ncol(comps_loop), names_sep = "_", names_to = c("test", "comp"))
    if(signif_markers){
      if(missing(custom_signif_markers)){
        signif_mapping <- c("***"=0.001, "**"=0.01, "*"=0.05)
      } else {
        signif_mapping <- custom_signif_markers
      }
      ## This loop replaces p-values according to the signif_mapping named vector
      for(j in 1:length(signif_mapping)){
        if(j==length(signif_mapping)){
          comps_loop$value[!(comps_loop$value %in% names(signif_mapping))] <-
            ifelse(as.numeric(comps_loop$value[!(comps_loop$value %in% names(signif_mapping))])<=signif_mapping[j], names(signif_mapping)[j], "NS")
        } else {
          comps_loop$value[!(comps_loop$value %in% names(signif_mapping))] <-
            ifelse(as.numeric(comps_loop$value[!(comps_loop$value %in%                                              
                                                   names(signif_mapping))])<=signif_mapping[j], names(signif_mapping)[j],
                   comps_loop$value[!(comps_loop$value %in% names(signif_mapping))])
        }
      }
    }
    if(p_adjust){
      adjusted <- "adj.pval"
      if(!(signif_markers)){
        comps <- as.data.frame(paste(signif(comps_loop$value[comps_loop$test == "p.slope.adj"],
                                            digits = 2), "/", signif(comps_loop$value[comps_loop$test == "p.intercept.adj"], digits = 2)))
      } else {
        comps <- as.data.frame(paste(comps_loop$value[comps_loop$test == "p.slope.adj"],
                                     "/", comps_loop$value[comps_loop$test == "p.intercept.adj"]))
      }
    } else {
      adjusted <- "non-adjusted.pval"
      if(!(signif_markers)){
        comps <- as.data.frame(paste(signif(comps_loop$value[comps_loop$test == "p.slope"],
                                            digits = 2), "/", signif(comps_loop$value[comps_loop$test == "p.intercept"], digits = 2)))
      } else {
        comps <- as.data.frame(paste(comps_loop$value[comps_loop$test == "p.slope"],
                                     "/", comps_loop$value[comps_loop$test == "p.intercept"]))
      }
    }
    colnames(comps) <- "p.value"
    rownames(comps) <- unique(comps_loop$comp)
    if(remove_NS){
      if(signif_markers){
        comps <- comps %>% filter(p.value != "NS / NS")
      } else {
        tmp_split <- as.data.frame(str_split(comps[,1], " / "))
        comps <- comps %>% filter(as.numeric(tmp_split[1,]) <= remove_NS_threshold |
                                    as.numeric(tmp_split[2,]) <= remove_NS_threshold)
      }
    }
    comps$group1 <- unlist(strsplit(rownames(comps), "v"))[c(TRUE, FALSE)]
    comps$group2 <- unlist(strsplit(rownames(comps), "v"))[c(FALSE, TRUE)]
    if(nrow(comps) > 0){
      comps$y_pos <- seq(from = 2, to = nrow(comps)*2, by = 2)
      
      p_bar <-  ggplot(data = plotDF, mapping = aes_string(x = endpoint, y = "value")) +
        geom_point(mapping = aes_string(color = endpoint), size = 0, alpha = 0) +
        suppressWarnings(geom_signif(data=comps, inherit.aes =F, vjust = -0.2, textsize = p_value_text_size,
                                     aes(xmin = group1, xmax = group2, annotations = p.value, y_position = y_pos),
                                     manual= TRUE)) + coord_cartesian(ylim = c(2, (nrow(comps)*2+2))) +
        theme_void() + theme(legend.position = "none") 
    } else {
      p_bar <- ggplot() + theme(panel.background = element_rect(fill = 'white', colour = 'white'))
    }
    if(missing(title)){
      if(model_type == "smoothSpline"){
        p_bar <- p_bar + ggtitle(paste0(unique(plotDF$name), " ", "shape/average ", adjusted)) # used to be slope/intercept
      } else {
        p_bar <- p_bar + ggtitle(paste0(unique(plotDF$name), " ", "slope/intercept ", adjusted))
      }
    } else {
      p_bar <- p_bar + ggtitle(title)
    }
    
    if(!(missing(font))){
      p_bar <- p_bar + theme(text=element_text(size=title_size,  family= font))
    }
  } 
  #### LINEAR MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(model_type == "lme"){
    if(missing(fit)){
      if(age_sex){
        if(age_only){
          fit <- lmer(formula(paste0("value ~ time_var * ", endpoint, " + discretized_admit_age_quantile + (1|enrollment_site/participant_id)")), 
                      data = plotDF)
        } else {
          fit <- lmer(formula(paste0("value ~ time_var * ", endpoint, " + sex + discretized_admit_age_quantile + (1|enrollment_site/participant_id)")), 
                      data = plotDF)
        }

      } else {
        fit <- lmer(formula(paste0("value ~ time_var * ", endpoint, " + (1|enrollment_site/participant_id)")), 
                    data = plotDF)
      }
    }
    plotDF$yhat <- stats::predict(fit, newdata = plotDF)
    ### Fixed Effect Intercept Matrix List
    if(age_sex){
      if(age_only){
        stop("LME model with age only is not yet supported.")
      } else {
        intercept_list <- NULL
        fixed_effects <- fixef(fit)
        for(u in levels(plotDF[[endpoint]])){
          if(sum(fixed_effects[grepl(paste0("trajectory_group", u), names(fixed_effects))]) == 0){
            active_effects <- 0
          } else {
            active_effects <- fixed_effects[grepl(paste0("trajectory_group", u), names(fixed_effects))]
            active_effects <- active_effects[!(grepl(":", names(active_effects)))]
          }
          group1 <- c(active_effects + fixed_effects[1] + fixed_effects[grepl("sex", names(fixed_effects))],
                      active_effects + fixed_effects[1] +
                        fixed_effects[grepl("sex", names(fixed_effects))] +
                        fixed_effects[grepl("discretized_admit_age_quantile", names(fixed_effects))])
          group2 <- c(active_effects + fixed_effects[1],
                      active_effects + fixed_effects[1] +
                        fixed_effects[grepl("discretized_admit_age_quantile", names(fixed_effects))])
          active_intercepts <- as.data.frame(rbind(group1, group2))
          rownames(active_intercepts)[1] <- names(fixed_effects)[grepl("sex", names(fixed_effects))]
          rownames(active_intercepts) <- gsub("sex", "", rownames(active_intercepts))
          rownames(active_intercepts)[2] <- levels(fit@frame$sex)[!(grepl(rownames(active_intercepts)[1], levels(fit@frame$sex)))]
          ann_text <- data.frame(time_var = sex_age_intercept_x,
                                 value = unlist(c(active_intercepts[1,], active_intercepts[2,])), lab = "Text",
                                 trajectory_group = factor(u, levels = c("1","2","3", "4", "5")),
                                 sex_age_group = c(paste0(rownames(active_intercepts)[1], colnames(active_intercepts)),
                                                   paste0(rownames(active_intercepts)[2], colnames(active_intercepts))))
          ann_text$sex_age_group <- gsub("\\(Intercept\\)", paste0("discretized_admit_age_quantile",
                                                                   levels(plotDF$discretized_admit_age_quantile)[1]), ann_text$sex_age_group)
          ann_text$sex_age_group <- gsub("trajectory_group.", paste0("discretized_admit_age_quantile",
                                                                     levels(plotDF$discretized_admit_age_quantile)[1]), ann_text$sex_age_group)
          ann_text$time_var <- ifelse(grepl(levels(plotDF$sex)[1], ann_text$sex_age_group),sex_age_intercept_x+sex_bar_seperator_dist, sex_age_intercept_x)
          intercept_list <- c(intercept_list, list(ann_text))
          names(intercept_list)[length(intercept_list)] <- paste0("trajectory_group", u)
        }
      }
    }
    pred <- ggpredict(fit, c("time_var", endpoint), type = "fixed",
                      ci.lvl = CI_interval)
    colnames(pred) <- gsub("x", "time_var", gsub("group", endpoint, colnames(pred)))
    
    p1 <- ggplot(data = plotDF, mapping = aes(x = time_var, y = value)) 
    if(individual_paths){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = value), color ='gray', alpha = individual_paths_alpha)
    }
    if(individual_trendlines){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = yhat),
                           alpha = individual_trendline_alpha,
                           linewidth = individual_trendline_size,
                           color = individual_trendlines_color)
    }
    if(individual_points){
      p1 <- p1 + geom_point(mapping = aes_string(color = endpoint), alpha = individual_points_alpha,
                            size = individual_points_size)
    }
    if(group_trendline == TRUE ){
      if(group_trendline_dropout){
        pred <- pred %>%
          left_join(plotDF %>%
                      mutate({{endpoint}} := as.character(.data[[endpoint]])) %>%
                      group_by(.data[[endpoint]]) %>% dplyr::summarize(max = max(time_var)), by = c({{endpoint}})) %>%
          filter(time_var <= max)
      }
      if(group_trendline_color == "endpoint"){
        p1 <- p1 + geom_line(inherit.aes = F, data = pred, mapping = aes_string(group = endpoint, y = "predicted", x="time_var",
                                           color = endpoint), linewidth = group_trendline_line_width)
      } else {
        p1 <- p1 + geom_line(inherit.aes = F, data = pred, mapping = aes_string(group = endpoint, y = "predicted", x="time_var"),
                             color = group_trendline_color, linewidth = group_trendline_line_width)
      }
    }
    
    if(CI){
      if(age_sex){
        p1 <- p1 + geom_ribbon(inherit.aes = F, data = pred,
                               aes_string(x="time_var", ymin = "conf.low", ymax = "conf.high", fill = endpoint), alpha = CI_alpha)
      } else {
        p1 <- p1 + geom_ribbon(inherit.aes = F, data = pred,
                               aes_string(x="time_var", ymin = "conf.low", ymax = "conf.high", fill = endpoint), alpha = CI_alpha)
      }
      
    }
    ##### Small intercept bars for the different age/sex groups
    if(sex_age_intercept_bars){
      if(age_sex){
        p1 <- p1 + coord_cartesian(xlim=c(0,max(plotDF$time_var)), clip = "off")
        if(is.null(sex_age_intercept_bars_colors)){
          color_sets <- c("#10A7E4", "#50BCE8", "#8ACBE6", "#B8D8E5", "#D2E1E7", "#F10816", "#F32E3A", "#EC767D", "#ECABAF", "#EDD0D2")
        } else {
          color_sets <- sex_age_intercept_bars_colors
        }
        for(w in 1:length(intercept_list)){
          for(y in 1:nrow(intercept_list[[w]])){
            p1 <- p1 + geom_text(data = intercept_list[[w]][y,], label = "-", size = bar_size, color = color_sets[y])
          }
        }
      }
    }
    
  } else if (model_type == "smoothSpline"){
    #### SMOOTH SPLINE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(age_sex){
      if(age_only) {
        formula_use <- formula(paste0("value~ s(time_var, bs = '",spline_type,"')+
              s(time_var, bs = '",spline_type,"', by =", endpoint, ")+", endpoint,
                                      "+ discretized_admit_age_quantile"))
      } else {
        formula_use <- formula(paste0("value~ s(time_var, bs = '",spline_type,"')+
              s(time_var, bs = '",spline_type,"', by =", endpoint, ")+", endpoint,
                                      "+ sex + discretized_admit_age_quantile"))
      }
    } else {
      formula_use <- formula(paste0("value~ 
          s(time_var, bs = '",spline_type,"')+s(time_var, bs = '",spline_type,"', by =", endpoint, ")+",endpoint))
    }
    fit <- gamm4::gamm4(formula_use, data = plotDF, random = ~(1|enrollment_site/participant_id))
    #print(formula_use)
    plotDF$yhat <- stats::predict(fit$mer)
    
    pred <- ggpredict(fit, c("time_var", endpoint), type = "fixed",
                      ci.lvl = CI_interval)
    colnames(pred) <- gsub("x", "time_var", gsub("group", endpoint, colnames(pred)))
    
    p1 <-ggplot(data = plotDF, mapping = aes(x = time_var, y = value)) 
    if(individual_paths){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = value), color ='gray', alpha = individual_paths_alpha)
    }
    if(individual_trendlines){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = yhat),
                           alpha = individual_trendline_alpha,
                           linewidth = individual_trendline_size,
                           color = individual_trendlines_color)
    }
    if(individual_points){
      p1 <- p1 + geom_point(mapping = aes_string(color = endpoint), alpha = individual_points_alpha,
                            size = individual_points_size)
    }
    if(group_trendline == TRUE ){
      if(group_trendline_dropout){
        pred <- pred %>%
          left_join(plotDF %>%
                      mutate({{endpoint}} := as.character(.data[[endpoint]])) %>%
                      group_by(.data[[endpoint]]) %>% dplyr::summarize(max = max(time_var)), by = c({{endpoint}})) %>%
          filter(time_var <= max)
      }
      
      if(group_trendline_color == "endpoint"){
        p1 <- p1 + geom_line(inherit.aes = F, data = pred, mapping = aes_string(group = endpoint, y = "predicted", x="time_var",
                                             color = endpoint), linewidth = group_trendline_line_width)
      } else {
        p1 <- p1 + geom_line(inherit.aes = F, data = pred, mapping = aes_string(group = endpoint, y = "predicted", x="time_var"),
                             color = group_trendline_color, linewidth = group_trendline_line_width)
      }
    }
    
    if(knot_lines){
      p1 <- p1 + geom_vline(xintercept = knots, linetype = 2)
    }
    if(CI){
      p1 <- p1 + geom_ribbon(inherit.aes = F, data = pred,
                             aes_string(x="time_var", ymin = "conf.low", ymax = "conf.high", fill = endpoint), alpha = CI_alpha)
    }
  } else if(model_type =="fixedKnots") {
    ################ FIXED KNOTS ################
    if(missing(fit)){
      if(age_sex){
        if(age_only) {
          fit_formula <- formula(paste0("value  ~ splines::bs(time_var, knots = knots) * ",
                                        endpoint, "+ discretized_admit_age_quantile + (1|enrollment_site/participant_id)"))
          fit <- lmer(fit_formula, data = plotDF)
          plotDFmod <- plotDF
          plotDFmod$sex <- levels(plotDFmod$sex)[1]
          plotDFmod$discretized_admit_age_quantile <- levels(plotDFmod$discretized_admit_age_quantile)[1]
        } else {
          fit_formula <- formula(paste0("value  ~ splines::bs(time_var, knots = knots) * ",
                                        endpoint, "+ sex + discretized_admit_age_quantile + (1|enrollment_site/participant_id)"))
          fit <- lmer(fit_formula, data = plotDF)
          plotDFmod <- plotDF
          plotDFmod$sex <- levels(plotDFmod$sex)[1]
          plotDFmod$discretized_admit_age_quantile <- levels(plotDFmod$discretized_admit_age_quantile)[1]
        }
      } else {
        fit_formula <- formula(paste0("value  ~ splines::bs(time_var, knots = knots) * ",
                                      endpoint, "+ (1|enrollment_site/participant_id)"))
        fit <- lmer(fit_formula, data = plotDF)
        plotDFmod <- plotDF
      }
    }
    plotDF$yhat = NULL; plotDF$yhat_wo_re = NULL
    plotDF$yhat <- stats::predict(fit, newdata = plotDF)
    plotDF$yhat_wo_re <- stats::predict(fit, re.form=NA, newdata=plotDFmod)
    if(CI){
      message("Confidence intervals will be added for fixedKnots soon.")
      #   if(age_sex){ ## In Development and Testing Still
      #     # fit_formula2 <- formula(paste0("value  ~ bs(time_var, knots = knots) * ",
      #     #                               endpoint, "+ (1|enrollment_site/participant_id)"))
      #     fit_formula2 <- formula(paste0("value  ~ bs(time_var) * ",
      #                                   endpoint, "+ (1|enrollment_site/participant_id)"))
      #     fit2 <- lmer(fit_formula2, data = plotDF)
      #     pred <- ggpredict(fit2, terms = c("time_var [all]", endpoint), type = "fixed",
      #                       ci.lvl = CI_interval)
      #     colnames(pred) <- gsub("x", "time_var", gsub("group", endpoint, colnames(pred)))
      #   }else{
      #     #plotDF <- as.data.frame(plotDF)
      #     pred <- ggeffect(fit, terms = c("time_var [all]", endpoint), type = "fixed",
      #                       ci.lvl = CI_interval)
      #     colnames(pred) <- gsub("x", "time_var", gsub("group", endpoint, colnames(pred)))
      #   }
    }
    p1 <-ggplot(data = plotDF, mapping = aes(x = time_var, y = value)) 
    
    if(individual_paths){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = value), color ='gray', alpha = individual_paths_alpha)
    }
    if(individual_trendlines){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = yhat),
                           alpha = individual_trendline_alpha,
                           linewidth = individual_trendline_size,
                           color = individual_trendlines_color)
    }
    if(individual_points){
      p1 <- p1 + geom_point(mapping = aes_string(color = endpoint), alpha = individual_points_alpha,
                            size = individual_points_size)
    }
    if(group_trendline == TRUE){
      if(group_trendline_dropout){
        pred <- pred %>%
          left_join(plotDF %>%
                      mutate({{endpoint}} := as.character(.data[[endpoint]])) %>%
                      group_by(.data[[endpoint]]) %>% dplyr::summarize(max = max(time_var)), by = c({{endpoint}})) %>%
          filter(time_var <= max)
      }
      if(group_trendline_color == "endpoint"){
        p1 <- p1 + geom_line(mapping = aes_string(group = endpoint, y = "yhat_wo_re", color = endpoint),
                             alpha = 1, linewidth = group_trendline_line_width)
      } else {
        p1 <- p1 + geom_line(mapping = aes_string(group = endpoint, y = "yhat_wo_re"),
                             alpha = 1, color = group_trendline_color, linewidth = group_trendline_line_width)
      }
    }
    if(knot_lines){
      p1 <- p1 + geom_vline(xintercept = knots, linetype = 2)
    }
    # if(CI){ ## In Development and Testing Still
    #   p1 <- p1 + geom_ribbon(inherit.aes = F, data = pred,
    #                          aes_string(x="time_var", ymin = "conf.low", ymax = "conf.high", fill = endpoint), alpha = .2)
    # }
  }
  ### Applies to all modeling p1s
  if(endpoint == "trajectory_group"){
    if(facet){
      p1 <- p1 + facet_wrap(facets = ~trajectory_group, nrow = 1,
                            labeller = labeller(trajectory_group = trajectory_labels)) + 
        theme_bw() + 
        labs(y = unique(plotDF$name)) + 
        theme(legend.pos = "none", strip.text.x = element_text(size = 7.5)) 
    } else {
      p1 <- p1 +  
        theme_bw() + 
        labs(y = unique(plotDF$name)) #+ 
        #theme(legend.pos = "none", strip.text.x = element_text(size = 7.5)) 
    }
  } else if (endpoint == "pro_group_labels") {
    if(facet){
      p1 <- p1 + facet_wrap(facets = ~factor(pro_group_labels, levels=c("MIN", "COG", "PHY", "MLT")), nrow = 1) + 
        theme_bw() + 
        labs(y = unique(plotDF$name)) + 
        theme(legend.pos = "none") 
    } else{
      p1 <- p1 + 
        theme_bw() + 
        labs(y = unique(plotDF$name)) + 
        theme(legend.pos = "none") 
    }
  } else if (endpoint == "pro_2_groups") {
    if(facet){
      p1 <- p1 + facet_wrap(facets = ~pro_2_groups, labeller=labeller(pro_2_groups = pro_min_labels), nrow = 1) + 
        theme_bw() + 
        labs(y = unique(plotDF$name)) + 
        theme(legend.pos = "none") 
    } else{
      p1 <- p1 + 
        theme_bw() + 
        labs(y = unique(plotDF$name)) + 
        theme(legend.pos = "none") 
    }
  } else {
    if(facet){
      p1 <- p1 + facet_wrap(facets = ~get(endpoint), nrow = 1) + 
        theme_bw() + 
        labs(y = unique(plotDF$name)) + 
        theme(legend.pos = "none") 
    } else{
      p1 <- p1 + 
        theme_bw() + 
        labs(y = unique(plotDF$name)) + 
        theme(legend.pos = "none") 
    }
  }
  if(!(missing(font))){
    p1 <- p1 + theme(text=element_text(family=font))
  }
  if(!(missing(facet_space))){
    p1 <- p1 + theme(panel.spacing = unit(facet_space, "lines"))
  }
  if(!(missing(colors))){
    p1 <- p1 + scale_color_manual(values=colors)
    p1 <- p1 + scale_fill_manual(values=colors)
  } else if (endpoint == "trajectory_group"){
    p1 <- p1 + scale_color_manual(values=c("#639A21", "#39828C", "#6371AD", "#BD7D31", "#9C3418") %>% setNames(1:5))
    p1 <- p1 + scale_fill_manual(values=c("#639A21", "#39828C", "#6371AD", "#BD7D31", "#9C3418") %>% setNames(1:5))
  }
  if(!(missing(xlabel))){
    p1 <- p1 + xlab(xlabel)
  } else {
    p1 <- p1 + xlab("Days from Admission")
  }
  if(!(missing(ylabel))){
    p1 <- p1 + ylab(ylabel)
  }
  if(!(missing(custom_theme_graph))){
    Theme_graph <-paste0('theme(', custom_theme_graph,')')
    p1 <- p1 + eval(parse(text=as.character(Theme_graph)))
  }
  if(!(missing(custom_theme_bars))){
    Theme_bars <-paste0('theme(', custom_theme_bars,')')
    p_bar <- p_bar + eval(parse(text=as.character(Theme_bars)))
  }
  if(y_axis_reverse){
    p1 <- p1 + scale_y_reverse()
  }
 
  if(!is.null(y_limits)){
    if(!is.null(x_limits)){
      p1 <- p1 + coord_cartesian(ylim = y_limits, xlim = x_limits)
    } else{
      p1 <- p1 + coord_cartesian(ylim = y_limits)
    }
  } else {
    if(!is.null(x_limits)){
      p1 <- p1 + coord_cartesian(xlim = x_limits)
    }
  }
  if(return_multi_obj){
    return(list(p1, p_bar))
  } else {
    if(!is.na(p_adjust)){
      return(egg::ggarrange(p_bar, p1, heights = c(bar_height, 10)))
    } else {
      return(p1)
    }
  }
 
}

mixed_pairwise = function(my.formula0, my.formula1, data_use){
  endpoints0 = sort(unique(data_use$endpoints))
  pair_names = c()
  res_table = c()
  endpoint_coef = c()
  
  for(i in 1:(length(endpoints0)-1)){
    for(j in (i+1):(length(endpoints0))){
      pair_names = c(pair_names, paste0(endpoints0[i],"|", endpoints0[j]))
      data_tmp = data_use[data_use$endpoints == endpoints0[i] | data_use$endpoints == endpoints0[j],]
      av = anova(lme4::lmer(my.formula0,data = data_tmp),
                 lme4::lmer(my.formula1,data = data_tmp))
      lmer_1 = lme4::lmer(my.formula1,data = data_tmp)
      #print(av)
      #print(lme4::lmer(my.formula1, data = data_tmp))
      res_table = c( res_table, av$`Pr(>Chisq)`[2])
      endpoint_coef = c(endpoint_coef, lmer_1@beta[2])
    }
  }
  
  names(res_table) = pair_names
  names(endpoint_coef) = pair_names
  return(list("res_table" = res_table, "coef_endpoint" = endpoint_coef))
}

violin_plot_visits <- function(data, colors, limits, ncol=5) {
  p <- ggplot(data, aes(x=endpoints, y=value, fill=endpoints)) +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
    facet_wrap(~Factor, drop = T, scales = "free_x", ncol = ncol) +
    theme_bw() + theme(legend.position = "none")
  if(!(missing(limits))){
    p <- p + scale_y_continuous(limits = limits)
  }
  if(!(missing(colors))){
    p <- p + scale_color_manual(values=colors)
    p <- p + scale_fill_manual(values=colors)
  }
  return(p)
}
