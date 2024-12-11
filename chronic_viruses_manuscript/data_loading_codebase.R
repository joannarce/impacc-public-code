##### Codebase for the loading IMPACC data into a data environmental object
##### Authors: Ravi K. Patel, Leying Guan, Jeremy Gygi, Gisela Gabernet
#
#.  This codebase contains functions used for the following:
#.  Initial data loading (via load_IMPACC_datasets(...))
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
      "doParallel"
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