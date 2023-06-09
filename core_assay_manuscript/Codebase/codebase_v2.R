##### Codebase for the Data analysis template
##### Authors: Ravi K. Patel, Leying Guan, Jeremy Gygi


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
    "ComplexHeatmap")


# A function to fetch "current" or "legacy" version of the data. The function will return data files that were current on the indicated date.
fetch_file_names <- function(dir, pattern = ".Counts.csv") {
  # Fetch all files within the dir
  files <- list.files(dir, recursive = T)
  
  # Select a file that matches the pattern.
  matched_file_name <- files[grepl(pattern, files)]
  # If multiple file exist, and one of them is .csv or .tsv, take that.
  is_csv_tsv <- grepl(".csv$", matched_file_name)
  if (length(matched_file_name) > 1 & sum(is_csv_tsv) == 1) {
    matched_file_name <- matched_file_name[is_csv_tsv][1]
  }
  if (length(matched_file_name) > 1) {
    stop("More than one files match \"",
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


# A function to generate a named vector of colors for a categorical variable using the levels_vector (a vector of unique categories) provided by the user.
get_categorical_color_vector <- function(levels_vector) {
  return(pals::alphabet(length(levels_vector)) %>% setNames(levels_vector))
}


# A function to generate a named vector of colors for an ordinal variable using the ordered_levels_vector (a vector of unique categories arranged in an increasing order) provided by the user. The user can provide name of RColorBrewer sequential palette or just assign a value between 1 and 18 to palette_index to choose an RColorBrewer sequential palette at that index.
get_ordinal_color_vector <-
  function(ordered_levels_vector,
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


read_data_files <- function(file, row.names = NULL, header = FALSE) {
  if(file.access(file, mode = 4) == -1) {
    warning(file, " does not have a read permission.")
    return(data.frame())
  }
  if (grepl(".csv$", file)) {
    x <- tryCatch({
      data <- read.csv(file, row.names = row.names, check.names = F, header=header, comment.char = "#", fileEncoding="UTF-8-BOM")
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


### Assign elements of a list to variables
assign_env_vars <- function(env) {
  lnames <- names(list)
  if( ! is.null( lnames ) ) {
    for( n in lnames ) {
      assign( n, list[[n]] )
    }
  }
}


#### Prepare clinical data
prepare_clinical_data <- function(data_env, DATA_VERSION, KEEP_COVID19_POS, FILTER_BY_CORE_ASSAY_COHORT, EVENT_DATE_UPPER, VISIT_UPPER, PHASES) {
  
  clinical_dir <- data_env[[ "data_dirs" ]][ "clinical_dir" ]
  
  ### Fetch clinical data files
  clinical_sample_file <- fetch_file_names( clinical_dir, "sample.csv$")
  clinical_event_file <- fetch_file_names( clinical_dir, "event.csv$")
  clinical_individ_file <- fetch_file_names( clinical_dir, "individ.csv$")
  clinical_individ_supp_file <- fetch_file_names( clinical_dir, "individ-supp.csv$")
  
  if (is.null(clinical_sample_file) |
      is.null(clinical_event_file) | is.null(clinical_individ_file)) {
    stop(
      "Error: Either you don't have access to the clinical data or the required clinical data files don't exist for the version of data you have selected. Please make sure that *sample.csv, *event.csv and *.individ.csv files exist in \"",
      clinical_dir,
      "\".\nExiting...\n",
      sep = ""
    )
  }
  
  message("Using the following clinical data files:\n")
  message("clinical_sample_file <- ", clinical_sample_file, "\n")
  message("clinical_event_file <- ", clinical_event_file, "\n")
  message("clinical_individ_file <- ", clinical_individ_file, "\n")
  message("clinical_individ_supp_file <- ", clinical_individ_supp_file, "\n")
  
  ### Load clinical data
  clinical_sample_data <- read_data_files( clinical_sample_file, header = TRUE )
  clinical_event_data <- read_data_files( clinical_event_file, header = TRUE )
  clinical_individ_data <- read_data_files( clinical_individ_file, header = TRUE )
  clinical_individ_supp_data <- read_data_files( clinical_individ_supp_file, header = TRUE )
  
  ### Combine clinical data
  clinical_data <-
    inner_join(clinical_sample_data,
               clinical_event_data,
               by = c("event_id" = "event_id")) %>%    # Join the clinical_sample_data and clinical_event_data tables using "event_id" column
    mutate(
      participant_id = participant_id.x,
      participant_id.x = NULL,
      participant_id.y = NULL
    ) %>%   # Remove redundant "participant_id" columns
    inner_join(clinical_individ_data,
               by = c("participant_id" = "participant_id"))   # Append the clinical_individ_data table
  
  ### Add individ-suppl file.
  if(! any(duplicated(clinical_individ_supp_data$participant_id))) {
   clinical_data <- left_join(clinical_data, clinical_individ_supp_data, by="participant_id")
  }
  
  
  ### Keep data from COVID-19 Positive individuals
  if(KEEP_COVID19_POS) {
    cat("Removing COVID-19 negative patients\n")
    clinical_data <- clinical_data %>% dplyr::filter(grepl("COVID-19 Positive", participant_type))
  }
  
  if(FALSE) {
  ### Keep core assay cohort
  if(FILTER_BY_CORE_ASSAY_COHORT) {
    cat("Removing samples not part of core_assay_cohort\n")
    clinical_data <- clinical_data %>% dplyr::filter( clinical_data$core_assay_cohort )
  }
  }
  
  ### Keep samples with event date <= EVENT_DATE_UPPER
  if( ! is.null(EVENT_DATE_UPPER) ) {
    cat("Removing samples with event_date greater than", EVENT_DATE_UPPER, "\n")
    clinical_data <- clinical_data %>% dplyr::filter( event_date <= EVENT_DATE_UPPER )
  }
  
  ### Keep samples from visit number <= VISIT_UPPER
  if( ! is.null(VISIT_UPPER) ) {
    cat("Removing samples with visit number (event_type) greater than Visit", VISIT_UPPER, "\n")
    allowed_visits = paste("Visit", 1:VISIT_UPPER)
    clinical_data <- clinical_data %>% dplyr::filter( grepl("Escalation|Control donor visit", event_type) | (event_type %in% allowed_visits) )
  }
  
  ### Keep samples from the selected phases
  clinical_data <- clinical_data %>% dplyr::filter(
    ! is.na(phase) &
    phase %in% c(PHASES) 
    )
  
  
  ### Generate discretized age based in 5 quantiles.
  admit_age_levels = sort(unique(clinical_data$admit_age))
  admit_age_levels_tiles = ntile(admit_age_levels, 5)
  label_df_tmp = data.frame( admit_age_levels_tiles, admit_age_levels) %>% 
    group_by(admit_age_levels_tiles) %>% 
    summarize(min=min(admit_age_levels), max=max(admit_age_levels)) %>% 
    mutate(labels = paste0("[",min, ",", max, "]") ) 
  admit_age_levels_tiles_labels = label_df_tmp$labels %>% setNames(label_df_tmp$admit_age_levels_tiles)
  admit_age_levels_labels = admit_age_levels_tiles_labels[admit_age_levels_tiles] %>% setNames(admit_age_levels)
  discretized_admit_age_quantile_old = admit_age_levels_labels[ as.character(clinical_data$admit_age) ]

  age = clinical_data$admit_age
  discretized_admit_age_quantile = case_when( age >= 19 & age <= 35 ~ "[19,35]",
                                              age >= 36 & age <= 50 ~ "[36,50]",
                                              age >= 51 & age <= 65 ~ "[51,65]",
                                              age >= 66 & age <= 80 ~ "[66,80]",
                                              age >= 81 & age <= 95 ~ "[81,95]")
  
  ### Generate discretized age based on equal width ages.
  discretized_admit_age_equalWidth = cut(clinical_data$admit_age, breaks = 5)
  
  clinical_data = clinical_data %>% mutate( discretized_admit_age_quantile = discretized_admit_age_quantile, discretized_admit_age_equalWidth = discretized_admit_age_equalWidth )
  
  
  ### Add discharge dates column
  # Extract discharge dates, concatenate the dates for patients that have multiple dates.
  discharge_event_data = clinical_event_data[ clinical_event_data$event_type == "Discharge", ] %>% 
    group_by(participant_id) %>% 
    summarise(event_dates = paste0(sort(event_date), collapse = ",") )
  clinical_data = clinical_data %>% 
    mutate(discharge_dates = 
             discharge_event_data$event_dates[ match(clinical_data$participant_id, discharge_event_data$participant_id) ]  
    )
  
  death_event_data = clinical_event_data[ clinical_event_data$event_type == "Death", ]
  clinical_data = clinical_data %>% 
    mutate(death_date = death_event_data$event_date[ match(clinical_data$participant_id, death_event_data$participant_id) ]
    )
  

  data_env[[ "clinical_data" ]] <- clinical_data
  data_env[[ "clinical_sample_file" ]] <- clinical_sample_file
  data_env[[ "clinical_event_file" ]] <- clinical_event_file
  data_env[[ "clinical_individ_file" ]] <- clinical_individ_file
  
  return(data_env)
}


#### Prepare info for the assay data files
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
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata from targeted proteomics of plasma"
    ),
    "plasma_proteomics_global_dda" = list(
                            dir = plasma_proteomics_global_dda_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata from global proteomics (DDA) of plasma"
    ),
    "plasma_proteomics_global_dia" = list(
                            dir = plasma_proteomics_global_dia_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata from global proteomics (DIA) of plasma"
    ),
    "serum_olink" = list(
                            dir = serum_olink_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata from Olink assay of serum"
    ),
    "nasal_viralload" = list(
                            dir = nasal_viralload_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata for nasal viral load"
    ),
    "serum_rbd_abtiters" = list(
                            dir = serum_rbd_abtiters_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata for RBD ab-titer in serum"
    ),
    "serum_sarscov2_abtiters" = list(
                            dir = serum_sarscov2_abtiters_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata for SARS-CoV-2 ab-titer in serum"
    ),
    "nasal_transcriptomics" = list(
                            dir = nasal_transcriptomics_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata for nasal transcriptomics"
    ),
    "plasma_metabolomics_global" = list(
                            dir = plasma_metabolomics_global_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata for global metabolomics of plasma"
    ),
    "bld_cytof" = list(
                            dir = bld_cytof_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata for CyTOF of blood"
    ),
    "ea_cytof" = list(
                            dir = ea_cytof_dir,
                            count_pattern = "cytof-Counts\\..*$",
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
                            metadata_desc = "Metadata for CyTOF of EA"
    ),
    "ea_transcriptomics" = list(
                            dir = ea_transcriptomics_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata for EA transcriptomics"
    ),
    "pbmc_transcriptomics" = list(
                            dir = pbmc_transcriptomics_dir,
                            count_pattern = "Counts\\..*$",
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
                            metadata_desc = "Metadata for PBMC transcriptomics"
    ),
    "ea_metagenomics" = list(
                            dir = ea_metagenomics_dir,
                            count_pattern = "BacterialTaxon-Counts_RPM.*$",
                            count_file_path = "",
                            count_var_name = "ea_metagenomics_counts",
                            count_desc = "Count data for EA metagenomics",
                            rowfeature_pattern = "BacterialTaxon-RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "ea_metagenomics_rowfeature",
                            rowfeature_desc = "RowFeature data for EA metagenomics",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "ea_metagenomics_metadata",
                            metadata_desc = "Metadata for EA metagenomics"
    ),
    "nasal_metagenomics" = list(
                            dir = nasal_metagenomics_dir,
                            count_pattern = "BacterialTaxon-Counts_RPM.*$",
                            count_file_path = "",
                            count_var_name = "nasal_metagenomics_counts",
                            count_desc = "Count data for nasal metagenomics",
                            rowfeature_pattern = "BacterialTaxon-RowFeature.*$",
                            rowfeature_file_path = "",
                            rowfeature_var_name = "nasal_metagenomics_rowfeature",
                            rowfeature_desc = "RowFeature data for nasal metagenomics",
                            metadata_pattern = "Metadata.*$",
                            metadata_file_path = "",
                            metadata_var_name = "nasal_metagenomics_metadata",
                            metadata_desc = "Metadata for nasal metagenomics"
    ),
    "serum_autoantibody" = list(
                            dir = serum_autoantibody_dir,
                            count_pattern = "Counts\\..*$",
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


#### Load assay datasets
load_assay_data <- function(data_env, DATA_VERSION, ALLOWED_SMPL_STATUS) {
  ### Fetch path of files that match the indicated pattern in the provided directory for each assay file and load the data in the variable names indicated in var_name element of the assay_files_info list. Uses the version of data specified in DATA_VERSION.
  
  assay_files_info <- data_env[[ "assay_files_info" ]]
  clinical_data <- data_env[[ "clinical_data" ]]
  
  for (f in names(assay_files_info)) {
    f_info <- assay_files_info[[f]]
    
    
    
    # Loading Metadata table
    if (!is.null(f_info$metadata_pattern)) {
      f_info$metadata_file_path <-
        fetch_file_names(f_info$dir, f_info$metadata_pattern)
      
      if (is.null(f_info$metadata_file_path)) {
        warning(
          "\"",
          f_info$metadata_desc,
          "\" is not accessible. Either you don't have access to this data or the data does not exist for the version of dataset you have selected.\n"
        )
        
      } else {
        ### Load data, keep only those observations for which the clinical data is available and assign the data object to the variable name indicated in "var_name" element.
        message(
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
      fetch_file_names(f_info$dir, f_info$count_pattern)
    
    if (is.null(f_info$count_file_path)) {
      warning(
        "\"",
        f_info$count_desc,
        "\" is not accessible. Either you don't have access to this data or the data does not exist for the version of dataset you have selected.\n"
      )
      
    } else {
      ### Load data, keep only those observations for which the clinical data is available and assign the data object to the variable name indicated in "var_name" element.
      message("Loading ",
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
        fetch_file_names(f_info$dir, f_info$rowfeature_pattern)
      
      if (is.null(f_info$rowfeature_file_path)) {
        warning(
          "\"",
          f_info$rowfeature_desc,
          "\" is not accessible. Either you don't have access to this data or the data does not exist for the version of dataset you have selected.\n"
        )
        
      } else {
        ### Load data, keep only those observations for which the clinical data is available and assign the data object to the variable name indicated in "var_name" element.
        message(
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


#### Prepare color schemes
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


#### Generate and print sample plots
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




#### Use this function to load all available datasets and return an R environment containing all data objects.
load_IMPACC_datasets <- function(DATA_VERSION, 
                                 KEEP_COVID19_POS, 
                                 PHASES = 1:100, 
                                 ALLOWED_SMPL_STATUS = 
                                   c("IMPACC sample assayed and passed QC", 
                                     "IMPACC sample assayed with questionable QC"),
                                 FILTER_BY_CORE_ASSAY_COHORT = TRUE,
                                 EVENT_DATE_UPPER = 42,
                                 VISIT_UPPER = 6,
                                 ASSAY_NAMES = NULL,
                                 data_base_dir = "/impacc-study/data_deposition/processed-data/"
                                 ) {
  
  if(! dir.exists(data_base_dir)) {
    cat(data_base_dir, "does not exist. Importing the base directory of computable matrices from Codebase/config.txt\n")
    config = read.table("../../Codebase/config.txt", header = T)
    data_base_dir = config$value[ config$key == "data_base_dir" ]
    if(! dir.exists(data_base_dir)) {
      stop("Error: the data base_dir specified in Codebase/config.txt does not exist. If you have not downloaded the computable matrices from ImmPort (ID:SDY1760) yet, please do so and refer to the usage instruction in the base README.md. Aborting...")
    }
  }
  
  #### Load packages necessary for data loading and processing purpose
  load_packages(c(packages, packages_imputation))
  
  
  
  #### Define Data directories: Use the absolute path of the directory containing "current" and "legacy" directories.
  data_dirs <- c(
      bld_cytof_dir = paste0(data_base_dir, "bld-cytof"),
      bld_gwas_dir = paste0(data_base_dir, "bld-gwas"),
      clinical_dir = paste0(data_base_dir, "clinical"),
      ea_cytof_dir = paste0(data_base_dir, "ea-cytof"),
      ea_metagenomics_dir = paste0(data_base_dir, "ea-metagenomics"),
      ea_transcriptomics_dir = paste0(data_base_dir, "ea-transcriptomics"),
      plasma_metabolomics_global_dir = paste0(data_base_dir, "metabolomics/plasma-metabolomics-global"),
      plasma_metabolomics_targeted_dir = paste0(data_base_dir, "metabolomics/plasma-metabolomics-targeted"),
      serum_metabolomics_global_dir = paste0(data_base_dir, "metabolomics/serum-metabolomics-global"),
      nasal_metagenomics_dir = paste0(data_base_dir, "nasal-metagenomics"),
      nasal_transcriptomics_dir = paste0(data_base_dir, "nasal-transcriptomics"),
      nasal_viralload_dir = paste0(data_base_dir, "nasal-viralload"),
      nasal_viralseq_dir = paste0(data_base_dir, "nasal-viralseq"),
      pbmc_transcriptomics_dir = paste0(data_base_dir, "pbmc-transcriptomics"),
      plasma_proteomics_targeted_dir = paste0(data_base_dir, "proteomics/plasma-proteomics-targeted"),
      plasma_proteomics_global_dda_dir = paste0(data_base_dir, "proteomics/plasma-proteomics-global-DDA"),
      plasma_proteomics_global_dia_dir = paste0(data_base_dir, "proteomics/plasma-proteomics-global-DIA"),
      serum_autoantibody_dir = paste0(data_base_dir, "serum-autoantibody"),
      serum_proteomics_global_dir = paste0(data_base_dir, "proteomics/serum-proteomics-global"),
      serum_olink_dir = paste0(data_base_dir, "serum-olink"),
      serum_rbd_abtiters_dir = paste0(data_base_dir, "serum-rbd-abtiters"),
      serum_sarscov2_abtiters_dir = paste0(data_base_dir, "serum-sarscov2-abtiters")
  )
  
  if(! is.null(ASSAY_NAMES)) {
    ASSAY_DIRS = c(paste0(ASSAY_NAMES, "_dir"), "clinical_dir")
    data_dirs[ ! names(data_dirs) %in% ASSAY_DIRS] = ""
    #data_dirs = data_dirs[ASSAY_DIRS]
  }
  
  #### Create an empty environment for storing various data objects
  data_env <- env(data_dirs = data_dirs)
  
  #### Prepare clinical data variables
  data_env <- prepare_clinical_data(data_env, DATA_VERSION, KEEP_COVID19_POS, FILTER_BY_CORE_ASSAY_COHORT, EVENT_DATE_UPPER, VISIT_UPPER, PHASES)

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



################################## Functions specific for data-integration group #########################################

packages_imputation <- 
  c("missForest", "doRNG", "doParallel")

load_packages <- function(packages) {
  for (n in 1:length(packages)) {
    suppressMessages(library(packages[n], character.only = TRUE))
  }
}


missForest_imputation <- function( assay_info, data_env, alpha = 0.2, random_seed = 2021, cores = 8 ) {
  
  count_df <- data_env[[ assay_info$count_var_name ]]
  rowfeature_df <- data_env[[ assay_info$rowfeature_var_name ]]
  count_desc <- assay_info$count_desc

  message(
    "Imputing the missing values (if any) in ",
    count_desc,
    " ...\n"
  )

  #### Imputation code adapted from Leying Guan
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
  registerDoParallel(cores = cores) # set based on number of CPU cores
  registerDoRNG(seed = random_seed)
  missF_out <- missForest(impute_df,  verbose = FALSE, parallelize = "variables")
  count_df_wImputed <- missF_out$ximp

  #t = sapply(1:ncol(X0[[1]]), function(x) {
  #  imp = median(X[[1]][is.na(X0[[1]][, x]) , x])
  #  nimp = median(X[[1]][!is.na(X0[[1]][, x]) , x])
  #  return(c(imp = imp, noimp = nimp))
  #})
  #t = t[,!is.na(t[1,])]
  #data.frame(t(t)) %>% 
  #  rownames_to_column("rownames") %>% 
  #  pivot_longer(names_to = "class",
  #               values_to = "median",
  #               cols = -rownames) %>% 
  #  ggplot(aes(x = rownames, y = median, fill = class)) + 
  #  geom_bar(stat = "identity", position = "dodge")
  
  data_env[[ assay_info$count_var_name ]] <- count_df_wImputed
  data_env[[ assay_info$rowfeature_var_name ]] <- rowfeature_df
}


preprocess_plasma_proteomics_targeted <- function( count_df ) {
  count_df_processed <- apply( count_df, 2, function(x){
    x[is.na(x)] <- min(na.omit(x))/2
    return(x)
  }) %>% data.frame()
  
  count_df_processed <- log( count_df_processed, 2 )
  return(count_df_processed)
}

preprocess_serum_olink_obsolete <- function(serum_olink_counts, removeOutliers = TRUE) {
  # list of outlier samples to remove
  flag <- NULL
  if (removeOutliers) {
    flag <- c("0865-002K3G00-002", "0865-002NVJ00-002", "0865-0041HH00-002", "0865-0062BG00-002")
  }
  serum_olink_counts_processed <- serum_olink_counts[rowMeans(is.na(serum_olink_counts)) < 1 &
                                                       !(rownames(serum_olink_counts) %in% flag), ] %>%
    as.matrix() %>%
    impute::impute.knn() %>%
    .$data %>%
    return(value = .)
}


load_IMPACC_datasets_vDataIntegration <- function(DATA_VERSION, 
                                                  KEEP_COVID19_POS, 
                                                  PHASES = 1:100, 
                                                  alpha = 0.2, 
                                                  random_seed = 2021,
                                                  ALLOWED_SMPL_STATUS = 
                                                    c("IMPACC sample assayed and passed QC", 
                                                      "IMPACC sample assayed with questionable QC"),
                                                  impute_cores = 8,
                                                  FILTER_BY_CORE_ASSAY_COHORT = TRUE,
                                                  EVENT_DATE_UPPER = 42,
                                                  VISIT_UPPER = 6
                                                  ) {
  
  #### Load packages necessary for data loading and processing purpose
  load_packages(c(packages, packages_imputation))
  
  #### Load all available IMPACC datasets associated with DATA_VERSION.
  data_env <- load_IMPACC_datasets(DATA_VERSION, KEEP_COVID19_POS, PHASES, ALLOWED_SMPL_STATUS, FILTER_BY_CORE_ASSAY_COHORT, EVENT_DATE_UPPER, VISIT_UPPER)
  
  #### Get assay data info
  assay_files_info <- data_env[[ "assay_files_info" ]]
  
  #### Impute each omic dataset using the standard approach
  for( n in names(assay_files_info)) {
    nrow_count_df <- nrow(data_env[[ assay_files_info[[n]]$count_var_name ]])
    if (!is.null(nrow_count_df)) {
      if (nrow_count_df > 0) {
        missForest_imputation(assay_files_info[[n]], data_env, alpha, random_seed, cores = impute_cores)
      }
    }
  }
  
  return(data_env)
}





load_IMPACC_Publication_datasets <- function(DATA_VERSION, 
                                             KEEP_COVID19_POS, 
                                             PHASES = 1:100, 
                                             ALLOWED_SMPL_STATUS = 
                                               c("IMPACC sample assayed and passed QC", 
                                                 "IMPACC sample assayed with questionable QC"),
                                             FILTER_BY_CORE_ASSAY_COHORT = TRUE,
                                             EVENT_DATE_UPPER = 42,
                                             VISIT_UPPER = 6,
                                             PUBLICATION_ID = 1234) {

  #### Load packages necessary for data loading and processing purpose
  load_packages(c(packages, packages_imputation))
  
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
  data_env <- prepare_clinical_data(data_env, DATA_VERSION, KEEP_COVID19_POS, FILTER_BY_CORE_ASSAY_COHORT, EVENT_DATE_UPPER, VISIT_UPPER, PHASES)
  
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

  
  
  #### Identify publication specific samples and remove the rest.
  
  path_to_metadata <- "/scratch/data-integration/pipelines/PublicationFiltering/example_project_df.csv"
  publication_id <- PUBLICATION_ID
  clinical_data <- data_env$clinical_data
  
  ### TODO:
  # Fix site ids to match clinical data!!
  
  ##############################
  #  Load publication_id info: #
  ##############################
  publication.metadata.df <- read_csv(file = path_to_metadata)
  if(!publication_id %in% publication.metadata.df$publication_id){
    stop("ERROR: publication_id provided was not found in the file. Double check that it is the correct publication_id.")
  }
  row.data <- publication.metadata.df[which(publication.metadata.df$publication_id == publication_id),]
  publication_name <- row.data$publication_name
  ###############
  #  Filtering: #
  ###############
  
  # Site ID Dictionary:
  temp <- clinical_data
  temp$site_id <- sapply(temp$participant_id, function(id){
    res <- stringr::str_split(id, "-")
    if(length(res[[1]]) == 1){
      return(id)
    } else {
      return(res[[1]][1])
    }
  })
  id.dict <- dplyr::arrange(
    dplyr::filter(
      dplyr::distinct(temp, enrollment_site, site_id), 
      !grepl("C", site_id)),
    site_id)
  
  # enrollment_site:
  site_cols <- which(grepl(pattern = "site_", x = colnames(publication.metadata.df)))
  site_vals <- unlist(row.data[,site_cols])
  site_vals <- site_vals[site_vals == 1]
  names(site_vals) <- gsub("site_", "", names(site_vals))
  available_sites <- sapply(names(site_vals), function(site_id){
    return(id.dict$enrollment_site[which(id.dict$site_id == site_id)])
  })
  data_env$clinical_data <- dplyr::filter(temp, enrollment_site %in% available_sites)
  
  
  #### Remove samples from each dataset that are not part of clinical_data
  #### Get assay data info
  assay_files_info <- data_env[[ "assay_files_info" ]]
  for( n in names(assay_files_info)) {
    count_var_name = assay_files_info[[n]]$count_var_name
    metadata_var_name = assay_files_info[[n]]$metadata_var_name
    if( ! is.null( data_env[[ count_var_name ]] )  ) {
      data_env[[ count_var_name ]] = 
        data_env[[ count_var_name ]] %>% 
        dplyr::filter( rownames(.) %in% data_env$clinical_data$sample_id )
      data_env[[ metadata_var_name ]] = 
        data_env[[ metadata_var_name ]] %>% 
        dplyr::filter( rownames(.) %in% data_env$clinical_data$sample_id )
    }
  }
  
  cat("Filtered for publication_id ", publication_id, ": ", publication_name, "...\n", sep = "")
  
  return(data_env)
  
}








################################## Functions specific for WGCNA #########################################

packages_WGCNA <- 
  c("WGCNA")

plot_soft_threshold_selection <- function(sft) {
  # Plot the results:
  par(mfrow = c(1, 2))
  
  cex1 = 0.9
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n",
    main = paste("Scale independence"),
    ylim = c(0, 0.9)
  )
  
  powers = sft$fitIndices$Power
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  
  
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.90, col = "red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity")
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  
}


tune_soft_threshold_WGCNA <-
  function(data_df,
           networkType = "signed",
           corFnc = "bicor",
           powers = c(c(1:10), seq(from = 12, to = 30, by = 2)), ...) {
    
  #### Load packages necessary for data loading and processing purpose
  load_packages(c(packages, packages_WGCNA))
  
  # Call the network topology analysis function
  sft = pickSoftThreshold(
    data_df,
    powerVector = powers,
    verbose = 0,
    networkType = networkType,
    corFnc = corFnc,
    ...
  )
  plot_soft_threshold_selection(sft)
  cat("WGCNA's suggestion for the power value: ", sft$powerEstimate, "\n")
  return(sft)
}



generate_WGCNA_modules <-
  function(data_df,
           networkType = "signed",
           power = NULL,
           minModuleSize = NULL,
           reassignThreshold = 1e-6,
           mergeCutHeight = 0.15,
           minKMEtoStay = 0.3,
           minCoreKME = 0.5,
           corType = "bicor",
           maxPOutliers = 0.1,
           assay_alias = NULL,   ### Append this assay_alias to the module names as prefix.
           includ_grey_in_plots = TRUE,
           plot_heatmap = TRUE,
           ...) {
    
    #### Load packages necessary for data loading and processing purpose
    load_packages(c(packages, packages_WGCNA))
    
    if( is.null( c(power, minModuleSize) ) ) {
      cat("ERROR:: Provide values for both power and minModuleSize parameters\n")
      return()
    }
    
    ## Network construction 
    net = blockwiseModules(
      data_df,
      power = power,
      networkType = networkType,
      TOMType = "unsigned",
      minModuleSize = minModuleSize,
      reassignThreshold = reassignThreshold,
      mergeCutHeight = mergeCutHeight,
      numericLabels = TRUE,
      pamRespectsDendro = FALSE,
      minKMEtoStay = minKMEtoStay,
      minCoreKME = minCoreKME,
      corType = corType,
      verbose = 3,
      maxPOutliers = maxPOutliers,
      ...
    )
    
    # Convert labels to colors for plotting
    mergedColors = labels2colors(net$colors)
    module_labels = paste0("mod",net$colors) %>% setNames(names(net$colors))
    module_color_key = labels2colors(unique(net$colors)) %>% setNames(paste0("mod",unique(net$colors)))
    # Plot the dendrogram and the module colors underneath
    plotDendroAndColors(
      net$dendrograms[[1]],
      mergedColors[net$blockGenes[[1]]],
      "Module colors",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05
    )
    
    if( ! is.null(assay_alias) ) {
      module_labels = paste0(assay_alias, "_", module_labels)
      module_color_key = labels2colors(unique(net$colors)) %>% setNames(paste0(assay_alias, "_","mod",unique(net$colors)))
    }
    
    # Recalculate MEs with module labels
    MEs = moduleEigengenes(data_df, module_labels)$eigengenes
    colnames(MEs) = gsub("^ME","",colnames(MEs))
    #MEs = orderMEs(MEs)
    
    module_membership = data.frame( feature = colnames(data_df), module = module_labels )
    
    pp <- NULL
    pp2 <- NULL
    
    if(plot_heatmap) {
      trajectory_group = clinical_data %>% dplyr::slice( match( row.names(data_df), clinical_data$sample_id ) ) %>% pull(trajectory_group)
      
      non_grey_data = data_df[,! grepl("mod0", module_membership$module)]
      data_for_plot = non_grey_data
      module_for_plot = module_membership$module[ ! grepl("mod0", module_membership$module) ]
      if(includ_grey_in_plots) {
        data_for_plot = data_df
        module_for_plot = module_membership$module
      }
      data_for_plot = apply(data_for_plot, 2, scale) %>% `rownames<-`(rownames(data_for_plot))
      p <-
        Heatmap(
          data_for_plot, 
          column_title_rot = 90,
          column_title_gp = gpar(fontsize=10),
          column_split = module_for_plot,
          show_column_names = FALSE,
          show_row_names = FALSE,
          left_annotation = rowAnnotation( traj_grp = trajectory_group, col=list(traj_grp = brewer.greens(5) %>% setNames(1:5) ) ),
          top_annotation = columnAnnotation( modules = module_for_plot, col=list(modules = module_color_key ) )
        )
      draw(p, padding = unit(c(2, 2, 15, 5), "mm"))
      #print(p)
      
      p2 <- 
        Heatmap(
          t(data_for_plot),
          column_title_rot = 90,
          column_title_gp = gpar(fontsize=10),
          row_names_gp = gpar(fontsize=6),
          column_split = trajectory_group,
          show_column_names = FALSE,
          top_annotation = columnAnnotation( traj_grp = trajectory_group, col=list(traj_grp = brewer.greens(5) %>% setNames(1:5) ) ),
          left_annotation = rowAnnotation( modules = module_for_plot, col=list(modules = module_color_key ) )
        )
      
      draw(p2, padding = unit(c(2, 2, 15, 5), "mm"))
      #print(p)
      
      
      
      pp <-
        Heatmap(
          data_for_plot, 
          column_title_rot = 90,
          column_title_gp = gpar(fontsize=10),
          column_split = module_for_plot,
          row_split = trajectory_group,
          show_column_names = FALSE,
          show_row_names = FALSE,
          name = "Z-score",
          left_annotation = rowAnnotation( TG = trajectory_group, 
                                           col = list(TG = c(
                                             "1" = "#639A21",
                                             "2" = "#39828C",
                                             "3" = "#6371AD",
                                             "4" = "#BD7D31",
                                             "5" = "#9C3418")),
                                           show_annotation_name = FALSE),
          top_annotation = columnAnnotation( modules = module_for_plot, col=list(modules = module_color_key ), show_legend = FALSE, show_annotation_name = FALSE)   )
      #draw(p, padding = unit(c(2, 2, 20, 5), "mm"))
      #print(p)
      
      pp2 <- 
        Heatmap(
          t(data_for_plot),
          column_title_rot = 90,
          column_title_gp = gpar(fontsize=10),
          row_names_gp = gpar(fontsize=6),
          column_split = trajectory_group,
          show_column_names = FALSE,
          name = "Z-score",
          top_annotation = columnAnnotation( trajectory = trajectory_group, 
                                             col = list(trajectory = c("1" = "#639A21",
                                                                       "2" = "#39828C",
                                                                       "3" = "#6371AD",
                                                                       "4" = "#BD7D31",
                                                                       "5" = "#9C3418"))),
          left_annotation = rowAnnotation( modules = module_for_plot, col=list(modules = module_color_key ) )
        )
      
      
      
      
      
    }
    
    colnames(MEs) = gsub("^ME", "", colnames(MEs))
    
    print("Note: mod0 = 'grey' module")
    
    return( 
      list( 
        module_membership = module_membership, 
        MEs = MEs, 
        WGCNAPlot1 = pp,
        WGCNAPlot2 = pp2
      ) 
    )
  }



how_to_use_WGCNA_functions <- function() {
  cat(quote(`
WGCNA analysis is a two-step process:
1: Tune the soft-threshold (power)
2: Construct a co-expression network

* Use the following function to help tune the power parameter. The following function will generate two plots, use the left plot to choose a lowest power value for which the R-square of Scale-free topology model reaches close to 0.9 (the red line).
        sft_tuned <-
          tune_soft_threshold_WGCNA(
                 data_df,                  # Preprocessed data in data.frame class
                 networkType = "signed",   # Indicate the type of network to construct
                 corFnc = "bicor",         # Correlation function    
                 ...)                      # Additional parameters that the user may want to pass on to pickSoftThreshold function.

* The following function constructs a co-expression network based on the input parameters (which users should tweak depending on the dataset being used) and returns a list with two elements: 1) a data frame indicating mappings between modules and features, and 2) another data frame (MEs) containing module-eigengenes. This function also generates a heatmap of feature expression across samples, with the features grouped by modules (grey module (mod0) is excluded).
        out_list <- 
          generate_WGCNA_modules(
                   data_df,                    # Preprocessed data in data.frame class
                   networkType = "signed",     # Indicate the type of network to construct
                   power = NULL,               # Provide a power value chosen from the first step, this parameter is mandatory.
                   minModuleSize = NULL,       # Provide a value for the minimum number of features a module can have; use a larger value (> 30) for datasets with thousands of features, and a smaller value otherwise.
                   corType = "bicor",          # Correlation function
                   maxPOutliers = 0.1,         # It is safer to use a smaller value for this parameter when corType = bicor is used.
                   assay_alias = NULL,         # Append this assay_alias to the module names as prefix.
                   reassignThreshold = 1e-6,   # This and rest of the parameters use default values, but users may want to tweak them based on the assay type and data.
                   mergeCutHeight = 0.15,
                   minKMEtoStay = 0.3,
                   minCoreKME = 0.5,
                   ...)                        # Additional parameters that the user may want to pass on to blockwiseModules function.
            `))
}


###cox model###############################################################################
#' prepare_right_censored_clinical_data (modified LG, 06/16/2021)
#' @description create right censored time-to-event data using user defined cutoff. 
#'              return a dataframe extend the clinical_data with columns:
#'              outcomeD14: recoded respiratory_status_day14 to factor "discharge" ~ 2, "moderate" ~3, 4, "severe" ~ 5, 6, "death" ~ 7.
#'              status: right-censored-endpoint time-to-event outcome status "censored" or "event"
#'              day_of_death: numeric number of days from admission to death
#'              day_event_or_censored: numeric number of days from admission to death or censored
#'              censored_event_date: right censored event_date in the clinical_data
#'              event_date_from_sympt, day_of_death_from_sympt, day_event_or_censored_from_sympt, censorred_event_date_from_sympt.  
#'               
#' @param clinical_data clinical_data dataframe generated using the data_analysis_template
#' @param censored_cutoff user defined numeric number of days from admission, 
#'                        it is used to create the right cencored time-to-event status
#' @param visit_filter optioinal integer used to filter for subject with minimal number of visits

prepare_right_censored_clinical_data <- function(clinical_data, censored_cutoff, visit_filter = 0){
  clinical_data %>%
    distinct(event_id, event_date, participant_id, enrollment_site ,
             sex, admit_age, race, symptom_date,
             death, respiratory_status_day14, event_type, event_location,          
             samples_collected, respiratory_status) %>%
    group_by(participant_id) %>% 
    # option to filter subjects with minimal `visit_filter` number of visits 
    dplyr::filter(max(row_number()) > visit_filter) %>%  
    mutate( 
      # code for outcome
      outcomeD14 = recode(respiratory_status_day14,
                          `2` = "discharge", 
                          `3` = "moderate", 
                          `4` = "moderate",
                          `5` = "severe",
                          `6` = "severe",
                          `7` = "death"),
      outcomeD14 = factor(outcomeD14,
                          levels = c("discharge", "moderate", "severe", "death")),
      day_of_death          = ifelse(test = death, 
                                     yes  = max(event_date), 
                                     no   = NA),
      
      # outcome right censored from admission
      status                = ifelse(test = death & (max(event_date) <= (censored_cutoff)), 
                                     yes  = "event", 
                                     no   = "censored"),    
      # right censored outcome time
      day_event_or_censored = ifelse(test = status == "censored", 
                                     yes  = min(max(event_date), censored_cutoff), 
                                     no   = day_of_death),
      # right censored time of each observation 
      censored_event_date   = pmin(event_date, day_event_or_censored),
      # same censoring using symptom onset as reference,
      symptom_date            = ifelse(test = symptom_date %in% -999, yes = NA, no = symptom_date),
      event_date_from_sympt   = event_date - symptom_date,
      day_of_death_from_sympt = ifelse(test = death, yes = max(event_date_from_sympt), no = NA),
      day_event_or_censored_from_sympt = ifelse(test = status == "censored", 
                                                yes  = min(max(event_date_from_sympt), 
                                                           (censored_cutoff - symptom_date)), day_of_death_from_sympt),
      censorred_event_date_from_sympt = pmin(event_date_from_sympt, day_event_or_censored_from_sympt),
      # use sort by subject length of time in study (for plot)
      n                     = max(event_date)) %>%
    # dropped NA symptom date
    dplyr::filter(event_date > 0, !is.na(event_date)) %>% 
    ungroup() %>%
    # arrange subject in order by  length of time in study (for plot)
    arrange(n) %>%
    mutate(participant_id = factor(participant_id, levels = unique(participant_id)))
}

#' display visit timecourse in days by subject, paneled in data without censoring and right censored 
#' @param right_censored_data
#' @param censored_cutoff
#' @param ... argument passed to ggplot to adjut font size for example

plot_right_censored_subject_visits <- function(right_censored_data, censored_cutoff, ...){
  plotTitle <- "no censoring"
  plotSubTitle <- paste(length(unique(right_censored_data$participant_id)), 
                        "participant /",  length(unique(right_censored_data$participant_id[right_censored_data$death])), 
                        "death")
  
  plotNotCensored <- right_censored_data %>%
    group_by(participant_id) %>%
    mutate(event_circle             = ifelse(test = death, 
                                             yes  = max(event_date), 
                                             no   = NA),
           respiratory_status_day14 = factor(respiratory_status_day14)) %>%
    ggplot(mapping = aes(x = event_date, y = participant_id, group = participant_id)) +
    geom_line() +
    geom_point(mapping = aes(x = event_circle), color = "red", shape = 4, size = 3) +
    geom_point(mapping = aes(x = event_date, shape = event_type), size = 2) +
    scale_shape_manual(values = c(21:25, 7:14) ) +
    labs(title    = plotTitle, 
         subtitle = plotSubTitle) +
    theme_classic(base_size = 10) +
    theme(title           = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(shape = 16)))
  
  # censored 
  plotTitle <- paste("censored at day", censored_cutoff, "from admission")
  plotSubTitle <- paste(length(unique(right_censored_data$participant_id)), 
                        "participants /", length(unique(right_censored_data$participant_id[right_censored_data$status == "event"])), "censored event")
  
  plotCensored <- right_censored_data %>%
    group_by(participant_id) %>%
    mutate(event_circle             = ifelse(test = status == "event", 
                                             yes  = max(censored_event_date), 
                                             no   = NA),
           respiratory_status_day14 = factor(respiratory_status_day14)) %>%
    dplyr::filter(event_date <= day_event_or_censored) %>%
    ggplot(mapping = aes(x = censored_event_date, y = participant_id, group = participant_id)) +
    geom_line() +
    geom_point(mapping = aes(x = event_circle), color = "red", shape = 4, size = 3) +
    geom_point(mapping = aes(x = event_date, shape = event_type), size = 2) +
    scale_shape_manual(values = c(21:25, 7:14)) +
    labs(title = plotTitle, subtitle = plotSubTitle) +
    theme_classic() +
    theme(title           = element_text(size=10, face='bold', hjust = 0.5),
          legend.position = "bottom")+
    guides(color = guide_legend(override.aes = list(shape = 16)))
  
  ggarrange(plotNotCensored, plotCensored, common.legend = TRUE, nrow = 1, legend = "bottom")
}


#' right join right censored clinical data to assay data (modified by Leying on 06/16/2021)
#' @description right join the right censored time-to-event data to assay data
#'              returns dataframe of long form of data 
#' example  
#'        participant_id  event_id        death     status    status_endpoint   censored_event_date     event_start    marker
#'        001-0002	      001-0002-1	    FALSE	    0	              0              12	                     0	            xxx    
#'        001-0002	      001-0002-2	    FALSE     0	              0              15	                     12	            xxx    
#'        003-0014	      003-0014-1	    TRUE      0	              1              25	                     0	            xxx   
#'        003-0014	      003-0014-2	    TRUE      0	              1              29	                     25	            xxx   
#'        003-0014	      003-0014-3	    TRUE      1	              1              32	                     29	            xxx  
#' status: binary form event-to-outcome status 1: death, 0: censored at a visit observasion, there the 'censored_event_date'
#' status_endpoint: binary form event-to-outcome status 1: death, 0: censored refer to the censored cutoff 
#' censored_event_date: days from admission that a visit that the assay data was observed
#' event_start: days from admission that the  privious timepoint to censored_event_date
#'             
#' @param assay_dataframe an imputated dataframe of assay data, with subject 'sample_id' as rowname and features as colnames
#'                        example:      
#'                                            feature_1   feature_2   feature_3  ....
#'                               sample_id_1     xxx         xxx         xxx     ....   
#'                               sample_id_2     xxx         xxx         xxx     .... 
#'                               sample_id_3     xxx         xxx         xxx     ....    
#' @param clinical_data
#' @param right_censored_data

prepare_right_censored_clinical_event_assay <- function(assay_dataframe, 
                                                        clinical_data,
                                                        right_censored_data){
  assay_dataframe %>%
    rownames_to_column(var = "sample_id") %>%
    merge(y  = clinical_data, by = "sample_id") %>%
    merge(y  = right_censored_data %>%
            # LG, add  day_event_or_censored
            distinct(participant_id, event_id, event_date, status, 
                     censored_event_date,  day_event_or_censored),
          by = c("participant_id", "event_id", "event_date")) %>%
    group_by(participant_id) %>% 
    dplyr::filter(event_date <= censored_event_date) %>% 
    dplyr::filter(max(row_number()) > visit_filter) %>%
    # sort by observation time, create long form data 
    arrange(censored_event_date) %>%
    mutate(event_start = lag(censored_event_date, default = 0),
           status_endpoint = ifelse(test = (status == "event"),
                                    yes  = 1, 
                                    no   = 0),
           status      = ifelse(test = (status == "event") & (censored_event_date == max(censored_event_date)),
                                yes  = 1, 
                                no   = 0)
    ) %>%
    arrange(participant_id, censored_event_date)%>%
    ungroup()
}

#' extract coxme results (LG, 06/16/2021)
#' x: fitted coxme object
#'
coxme_extract = function (x, rcoef = FALSE, digits = 2) {
  coxme_results =list()
  beta <- x$coefficients
  nvar <- length(beta)
  nfrail <- nrow(x$var) - nvar
    se <- sqrt(diag(x$var)[nfrail + 1:nvar])
    tmp <- cbind(beta, exp(beta), se, round(beta/se, 2), 
                 signif(1 - pchisq((beta/se)^2, 1), 2))
  dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", 
                                        "se(coef)", "z", "p"))
  coxme_results$coefficients = tmp
  coxme_results
}

##########visit based random effect model##############
##ordinal mixed effect model testing for ordinal effects for a given module
mixed_ordinal = function(my.formula0, my.formula1, data_use){
  av = anova(ordinal::clmm(formula(my.formula0), data = data_use),
             ordinal::clmm(formula(my.formula1), data = data_use))
  res = rep(0, 2)
  names(res) = c("AIC", "pval")
  res[1] = av$AIC[2]
  res[2] = (av$`Pr(>Chisq)`)[2]
  return(res)
}

##mixed effect pairwise comparison
mixed_pairwise = function(my.formula0, my.formula1, data_use){
  endpoints0 = sort(unique(data_use$endpoints))
  pair_names = c()
  res_table = c()
  for(i in 1:(length(endpoints0)-1)){
    for(j in (i+1):(length(endpoints0))){
      pair_names = c(pair_names, paste0(endpoints0[i],"|", endpoints0[j]))
      data_tmp = data_use[data_use$endpoints == endpoints0[i] | data_use$endpoints == endpoints0[j],]
      av = anova(lme4::lmer(my.formula0,data = data_tmp),
                 lme4::lmer(my.formula1,data = data_tmp))
      res_table = c( res_table, av$`Pr(>Chisq)`[2])
    }
  }
  names(res_table) = pair_names
  return(res_table)
}


violin_plot_module_baseline = function(data_use, res_table){
  data_plot <- data_use[,colnames(data_use)%in%rownames(res_table)]
  data_plot$endpoints = data_use$endpoints
  data_plot<- data_plot%>%gather(key = 'module', value = 'expression', -endpoints)
  #data_plot <- data_use%>%gather(key = 'module', value = 'expression', -endpoints, -control, -sites)
  data_plot$expression = as.numeric(data_plot$expression)
  
  # Adding p-values:
  data_plot$p.val <- sapply(data_plot$module, function(mod){
    return(res_table[which(rownames(res_table) == mod), 2])
  })
  # Adding q-values:
  data_plot$q.val <- sapply(data_plot$module, function(mod){
    return(res_table[which(rownames(res_table) == mod), 3])
  })
  
  res <- ggplot() + 
    geom_violin(data = data_plot, aes(x=endpoints, y=expression, fill = endpoints)) +
    geom_text(data = distinct(data_plot, module, p.val), aes(x = 3.5, y = .2, label = paste0("p.val - ", round(p.val, 4)))) +
    geom_text(data = distinct(data_plot, module, q.val), aes(x = 3.5, y = .15, label = paste0("q.val - ", round(q.val, 4)))) +
    facet_wrap(~module) +
    theme_bw()
  return(res)
}

#Maintainer: Leqi Xu
#' modules: NULL=plot all; 
#' label: q.signif or  q.val
association_plot = function(data_use,  res_table_pairwise, res_table,
                            modules = NULL, ncol = 4, 
                            label = "q.signif",
                            colors = c("#639A21", "#39828C",  "#6371AD","#BD7D31", "#9C3418")){
  
  pairwise_res_table_formatted= pairwise_format(pairwise_table = res_table_pairwise)
  ordinal_res_table_formatted= ordinal_format(ordinal_table = res_table)
  ordinal_res_table_formatted = ordinal_res_table_formatted[,c(4,6,7)]
  data_plot <- data_use[,colnames(data_use)%in%rownames(res_table)]
  if(is.null(dim(data_plot))){
    data_plot = data.frame(data_plot)
    colnames(data_plot)=rownames(res_table)
  }
  data_plot$endpoints = data_use$endpoints
  data_plot<- data_plot%>%gather(key = 'module', value = 'expression', -endpoints)
  #data_plot <- data_use%>%gather(key = 'module', value = 'expression', -endpoints, -control, -sites)
  data_plot$expression = as.numeric(data_plot$expression)
  # Adding p-values:
  data_plot$p.val <- sapply(data_plot$module, function(mod){
    return(res_table[which(rownames(res_table) == mod), 2])
  })
  # Adding q-values:
  data_plot$q.val <- sapply(data_plot$module, function(mod){
    return(res_table[which(rownames(res_table) == mod), 3])
  })
  
  module_pvalue <- ordinal_res_table_formatted[,2]
  names(module_pvalue) <- ordinal_res_table_formatted$module
  module_psig <- ordinal_res_table_formatted[,3]
  names(module_psig) <-ordinal_res_table_formatted$module
  if(!is.null(modules)){
    data_plot_01 <- data_plot[which(data_plot$module %in% modules),]
    res_table_pairwise_all01 <- pairwise_res_table_formatted[which(pairwise_res_table_formatted$module%in%modules),]
    module_pvalue = module_pvalue[modules]
    module_psig = module_psig[modules]
  }else{
    data_plot_01 <- data_plot
    res_table_pairwise_all01 <- pairwise_res_table_formatted
  }
  data_plot_01$endpoints = factor(paste0("TG",data_plot_01$endpoints),levels = c("TG1","TG2","TG3","TG4","TG5"))
  
  mod01_box0 <- ggboxplot(data_plot_01, x = "endpoints", y = "expression",
                            color = "endpoints", palette =colors,
                            add = "jitter") +
  stat_pvalue_manual(data = res_table_pairwise_all01[res_table_pairwise_all01$q.signif!="ns",], 
                         y.position = quantile(data_plot_01$expression,0.99), step.increase = 0.05, step.group.by = "module", bracket.size = 0.2,  remove.bracket  = F,
                         label = "q.signif") +
  facet_wrap(~ module,labeller = labeller(module = module_pvalue), ncol = min(ncol, length(unique(data_plot_01$module))))+
  theme_bw()
  return(mod01_box0)
}

pairwise_format <- function(pairwise_table){
  pairwise_table12 <- data.frame(module = rownames(pairwise_table),group1="1",group2="2",q.val = pairwise_table$`1|2`)
  pairwise_table23 <- data.frame(module = rownames(pairwise_table),group1="2",group2="3",q.val = pairwise_table$`2|3`)
  pairwise_table34 <- data.frame(module = rownames(pairwise_table),group1="3",group2="4",q.val = pairwise_table$`3|4`)
  pairwise_table45 <- data.frame(module = rownames(pairwise_table),group1="4",group2="5",q.val = pairwise_table$`4|5`)
  pairwise_table13 <- data.frame(module = rownames(pairwise_table),group1="1",group2="3",q.val = pairwise_table$`1|3`)
  pairwise_table24 <- data.frame(module = rownames(pairwise_table),group1="2",group2="4",q.val = pairwise_table$`2|4`)
  pairwise_table35 <- data.frame(module = rownames(pairwise_table),group1="3",group2="5",q.val = pairwise_table$`3|5`)
  pairwise_table14 <- data.frame(module = rownames(pairwise_table),group1="1",group2="4",q.val = pairwise_table$`1|4`)
  pairwise_table25 <- data.frame(module = rownames(pairwise_table),group1="2",group2="5",q.val = pairwise_table$`2|5`)
  pairwise_table15 <- data.frame(module = rownames(pairwise_table),group1="1",group2="5",q.val = pairwise_table$`1|5`)
  pairwise_table_all <- rbind(pairwise_table12,pairwise_table23,pairwise_table34,pairwise_table45,pairwise_table13,pairwise_table24,pairwise_table35,pairwise_table14,pairwise_table25,pairwise_table15)
  pairwise_table_all$q.val <- signif(pairwise_table_all$q.val,3)
  pairwise_table_all$q.signif <- cut(pairwise_table_all$q.val, 
                                     breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                     labels = c("****", "***", "**", "*", "ns"))
  pairwise_table_all$group1 = paste0("TG",pairwise_table_all$group1)
  pairwise_table_all$group2 = paste0("TG",pairwise_table_all$group2)
  return(pairwise_table_all)
}

ordinal_format <- function(ordinal_table){
  ordinal_table$module <- rownames(ordinal_table)
  ordinal_table$qval <- signif(ordinal_table$qval,3)
  ordinal_table$q.signif <- cut(ordinal_table$qval, 
                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                labels = c("****", "***", "**", "*", "ns"))
  ordinal_table$modp <- paste0(ordinal_table$module, ": ", ordinal_table$qval)
  ordinal_table$modsig <- paste0(ordinal_table$module, ": ", ordinal_table$q.signif) 
  
  return(ordinal_table)
}

##########trajectory analysis
visit_summary_plot = function(clinical_data_use){
  plotDF <- clinical_data_use %>%
    dplyr::filter(grepl(pattern = "Visit", event_type)) %>%
    add_column(flag = 1) %>%
    dplyr::select(participant_id, event_type, flag) %>%
    distinct() %>%
    pivot_wider(names_from = event_type, values_from = flag) %>%
    dplyr::select(-participant_id) %>%
    group_by_all() %>%
    summarize(n = n()) %>%
    arrange(desc(n)) %>%
    rowid_to_column() %>%
    pivot_longer(cols = -c(rowid, n), names_to = "event_type") %>%
    dplyr::filter(!is.na(value))
  
  ggplot(data = plotDF,
         mapping = aes(x = event_type, y = rowid)) +
    geom_point(mapping = aes(size = n)) +
    geom_line(mapping = aes(group = rowid)) +
    scale_y_continuous(breaks = 1:max(plotDF$rowid), 
                       labels = distinct(dplyr::select(plotDF,rowid, n))$n) +
    labs(y = "Number of participants") +
    theme_bw() +
    theme(axis.ticks.y = element_blank(),
          axis.text.x  = element_text(angle = 45, hjust = 1))
}


mgcv_global = function(inputDF, age_sex = T, endpoints = "trajectory_group"){
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
    print(i)
    feature_name = features_names[i]
    tmp_input = inputDF[inputDF$name==feature_name,]
    tmp_input[[endpoints]] =ordered(as.factor(tmp_input[[endpoints]]))
    tmp_input$participant_id = as.factor(tmp_input$participant_id)
    if(age_sex){
      formula_use = formula(paste0("value~ s(event_date, bs = 'cr')+
              s(event_date, bs = 'cr', by =", endpoints, ")+", endpoints,
                                   "+ sex + discretized_admit_age_quantile"))
    }else{
      formula_use = formula(paste0("value~ s(event_date, bs = 'cr')+
              s(event_date, bs = 'cr', by =", endpoints, ")+",endpoints))
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
  res_table$adjp.intercept= qvalue::qvalue(res_table$p.intercept, fdr.level = 0.05, pi0 = 1)$qvalues
  res_table$adjp.slope= qvalue::qvalue(res_table$p.slope, fdr.level = 0.05, pi0 = 1)$qvalues
  if(age_sex){
    res_table$adjp.intercept.sex= qvalue::qvalue(res_table$p.intercept.sex, fdr.level = 0.05, pi0 = 1)$qvalues
    res_table$adjp.intercept.age.quantile= qvalue::qvalue(res_table$p.intercept.age.quantile, fdr.level = 0.05, pi0 = 1)$qvalues
  }
  rownames(res_table) = features_names
  return(res_table)
}

pairwise_lme_modeling <- function(data, variable = "trajectory_group", fixedKnots = F, sex_age = T,
                                  knots = c(1,4,7,14,21)) {
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
          formula_use= formula(paste0("value ~ event_date * ",variable, "+sex +discretized_admit_age_quantile"))
        }else{
          formula_use= formula(paste0("value ~ event_date * ",variable))
        }
      }else{
        if(sex_age){
          formula_use <- formula(paste0("value ~ splines::bs(event_date, knots=c(",paste(knots, collapse = ',') ,")) * ",variable,"+sex +discretized_admit_age_quantile"))
        } else {
          formula_use <- formula(paste0("value ~ splines::bs(event_date, knots=c(",paste(knots, collapse = ',') ,")) * ",variable))
        }
      }
      lmfit= lme(fixed = formula_use, random = ~1|enrollment_site/participant_id, data = data_tmp)
      if(sex_age){
        compout[1,j] = (anova( lmfit)$"p-value")[6]
        compout[2,j] = (anova( lmfit)$"p-value")[3]
        compout[3,j] = (anova( lmfit)$"p-value")[4] #sex
        compout[4,j] = (anova( lmfit)$"p-value")[5] #age_quant
      }else{
        compout[1,j] = (anova( lmfit)$"p-value")[4]
        compout[2,j] = (anova( lmfit)$"p-value")[3]
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

pairwise_mgcv_modeling <- function(data, age_sex = T, variable = "trajectory_group"){
  ## collect pairwise comparisons for outcomeD14 and prep output vector/df
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
        formula_use = formula(paste0("value~ s(event_date, bs = 'cr')+
              s(event_date, bs = 'cr', by =", variable, ")+", variable,
              "+ sex + discretized_admit_age_quantile"))
      }else{
        formula_use = formula(paste0("value~ s(event_date, bs = 'cr')+
              s(event_date, bs = 'cr', by =", variable, ")+",variable))
      }
      fit <- try(gamm4::gamm4(formula_use, data = tmp_data, random = ~(1|enrollment_site/participant_id)))
      a1 = anova(fit$gam)
      compout[2,j] = a1$pTerms.pv[1]
      compout[1,j] = a1$s.table[2,4]
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

## For a given input, creates a df of p.slope, p.intercept, their pairwise 
## counterparts and the respective adjusted p.values for all of these
model_loop <- function(inputDF, modelType = "lme", age_sex =TRUE, knots = c(1, 4, 7, 14, 21),
                       endpoint = "trajectory_group", old_p_corrections = FALSE){
  inputDF$participant_id = as.factor(inputDF$participant_id )
  if(modelType == "smoothSpline"){
    inputDF$event_date_transformed = inputDF$event_date
    lmDF <- mgcv_global(inputDF, endpoint, age_sex = age_sex)
    lmDF$name <- rownames(lmDF)
    ## Because of the pairwise function requires there to be a column for "name"
    ## the rownames are copied to a new column here, ultimately this column is removed
  } else {
    lmDF <- inputDF %>%
      group_by(name) %>% 
      do(p.slope = {
        if(modelType == "lme"){
          if(age_sex){
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint,
                            "+sex+discretized_admit_age_quantile")), data = .,
                           random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint)), data = .,
                           random =  ~1|enrollment_site/participant_id))
          }
        } else if (modelType == "fixedKnots"){
          if(age_sex){
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                             "+sex+discretized_admit_age_quantile")), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint)), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          }
        }
        if(class(fit)[1]!="try-error"){
            idx <- nrow(anova(fit))
            result = (anova(fit)$"p-value")[idx]
        } else {
          result = NA
        }
      }, 
      p.intercept = {
        if(modelType == "lme"){
          if(age_sex){
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint,
                                                     "+sex+discretized_admit_age_quantile")), data = .,
                           random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint)), data = .,
                           random =  ~1|enrollment_site/participant_id))
          }
        } else if (modelType == "fixedKnots"){
          if(age_sex){
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                                                     "+sex+discretized_admit_age_quantile")), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          } else{
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint)), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          }
          
        }
        if(class(fit)[1]!="try-error"){
          result = (anova(fit)$"p-value")[3] #trajectory_group is position 3
        }else{
          result = NA
        }
      },
      p.intercept.sex = {
        if(age_sex){
          if(modelType == "lme"){
              fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint,
                                                       "+sex+discretized_admit_age_quantile")), data = .,
                             random =  ~1|enrollment_site/participant_id))
          } else if (modelType == "fixedKnots"){
              fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                                                       "+sex+discretized_admit_age_quantile")), 
                             data = .,random =  ~1|enrollment_site/participant_id))
          }
          if(class(fit)[1]!="try-error"){
            result = (anova(fit)$"p-value")[4] ## sex is position 4
          } else {
            result = NA
          }
        }
      },
      p.intercept.age.quantile = {
        if(age_sex){
          if(modelType == "lme"){
            fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint,
                                                     "+sex+discretized_admit_age_quantile")), data = .,
                           random =  ~1|enrollment_site/participant_id))
          } else if (modelType == "fixedKnots"){
            fit <- try(lme(fixed = as.formula(paste0("value ~ splines::bs(event_date, 
                             knots =c(",paste(knots, collapse = ',') ,")) * ", endpoint,
                                                     "+sex+discretized_admit_age_quantile")), 
                           data = .,random =  ~1|enrollment_site/participant_id))
          }
          if(class(fit)[1]!="try-error"){
            result = (anova(fit)$"p-value")[5] ## discretized_admit_age_quantile is position 4
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
      if(modelType == "lme"){
        tmp.out <- pairwise_lme_modeling(inputDF[inputDF$name == lmDF$name[i],],
                                         variable = endpoint, fixedKnots = F)
      } else if(modelType == "fixedKnots"){
        tmp.out <- pairwise_lme_modeling(inputDF[inputDF$name == lmDF$name[i],],
                                         variable = endpoint, fixedKnots = T, knots = knots)
      } else if (modelType == "smoothSpline"){
        tmp.out <- pairwise_mgcv_modeling(inputDF[inputDF$name == lmDF$name[i],],
                                          variable = endpoint)
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
    p.slope.adj <- qvalue::qvalue(as.vector(p.slope.adj), fdr.level = 0.05, pi0 = 1)$qvalues
    p.intercept.adj <- as.data.frame(values.out[,grepl("p.intercept", colnames(values.out)), drop=FALSE])
    p.intercept.adj <- qvalue::qvalue(as.vector(p.intercept.adj), fdr.level = 0.05, pi0 = 1)$qvalues
  }
  colnames(p.slope.adj) <- gsub("_", ".adj_", colnames(p.slope.adj))
  colnames(p.intercept.adj) <- gsub("_", ".adj_", colnames(p.intercept.adj))
  
  lmDF_pairwise <- as.data.frame(cbind(lmDF, values.out, p.slope.adj, p.intercept.adj))
  rownames(lmDF_pairwise) <- lmDF_pairwise$name
  lmDF_pairwise <- lmDF_pairwise[,colnames(lmDF_pairwise) != "name"]
  return(lmDF_pairwise)
}

## plotDF takes an IMPACC inputDF data format. model_loop takes the output from the model_loop function
## for p-values and adjusted p-values. modelType can be "lme", "fixedKnots", or "smoothSpline", if "fixedKnots"
## provide the days at which knots should be place (defaults to 1,4,7,14, and 21). 
## p_adjust = T will use adj. p-values from model_loop. signif_markers = T will switch
## the comparison bars to use a * based system. Defaults to p<0.5 = *, p<0.01 = **, and p<0.001 = ***
## NS (no significance). This system can be overrode by providing a named list of cutoffs from
## least to greatest to custom_signif_markers (e.g. custom_signif_markers = c("*" = 0.05, "." = 0.1, " "=1) )
## These provided markers will be used for the first marker that comparison qualifies for
## by being less than or equal to the cutoff (hence why it is important for it to be least to greatest).
## remove_NS = T will remove any bars for which NS is detected for both p.slope and p.intercept.
## remove_NS_threshold allows you to change the criteria for which comparisons should be removed.
## fit allows users to provide an alternate model they generated (will conflict with group trendline if used)
## font changes font broadly, title_size is the title font size, bar_height is the
## relative proportion of the significance bars to the actual plots 
## (the plots are at height 10, so if this is set to 10 then sign bars will be half the image)
## custom_theme_graph allows for custom ggtheme parameters to be provided for the plots (bottom half of image)
## custom_theme_bars allows for custom ggtheme parameters to be provided for the significance bar image (top half)
## knot_lines will add dotted lines at the knots by default for the fixedKnot and smoothSpline model.
## group_trendline adds that models overall group model (ie without random effects).
## There are associated arguments for its color and size. Same for the next two arguments:
## individual_trendlines, individual_points, and individual_paths. Alpha is a 
## parameter for transparency.

plot_model <- function(plotDF, model_loop, endpoint = "trajectory_group", modelType = "lme",
                       knots = c(1, 4, 7, 14, 21), age_sex = TRUE, facet = TRUE,
                       p_adjust = NA, signif_markers =TRUE, custom_signif_markers, remove_NS = FALSE,
                       remove_NS_threshold = 0.05, fit, font, title_size =10, bar_height = 5,
                       colors, title, xlabel, ylabel, custom_theme_graph, custom_theme_bars, knot_lines = TRUE,
                       group_trendline = TRUE, group_trendline_color = "black", group_trendline_line_width = 1.5,
                       group_trendline_dropout = TRUE,
                       individual_trendlines = TRUE, individual_trendline_size = 0.6, individual_trendline_alpha = 0.9, individual_trendlines_color = "#5F5F5F",
                       individual_points = TRUE, individual_points_size = 1, individual_points_alpha = 0.8,
                       individual_paths = TRUE, individual_paths_size = 1, individual_paths_alpha = 0.4,
                       sex_age_intercept_bars = FALSE, sex_age_intercept_x = -3, sex_bar_seperator_dist = 2,
                       bar_size = 10, facet_space =1, intercept_hist = NULL, sex_age_intercept_bars_colors = NULL,
                       CI = FALSE, CI_alpha = 0.3, CI_interval = 0.95, y_axis_reverse = FALSE, y_limits = NULL, x_limits = NULL,
                       return_multi_obj = FALSE, p_value_text_size = 3.88){
  #### DATA CHECKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  if(modelType == "smoothSpline"){
    if(!is.factor(plotDF$participant_id)){
      stop("Please convert participant_id to a factor.")
    # } else if (!is.ordered(plotDF[[endpoint]])){
    #   stop("Please convert your choosen endpoint to an ordered factor.")
    } else if (!is.factor(plotDF$enrollment_site)){
      stop("Please convert your choosen enrollment_site to a factor.")
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
        comps <- comps %>% dplyr::filter(p.value != "NS / NS")
      } else {
        tmp_split <- as.data.frame(str_split(comps[,1], " / "))
        comps <- comps %>% dplyr::filter(as.numeric(tmp_split[1,]) <= remove_NS_threshold |
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
      if(modelType == "smoothSpline"){
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
  if(modelType == "lme"){
    if(missing(fit)){
      if(age_sex){
        fit <- lmer(formula(paste0("value ~ event_date * ", endpoint, " + sex + discretized_admit_age_quantile + (1|enrollment_site/participant_id)")), 
                    data = plotDF)
      } else {
        fit <- lmer(formula(paste0("value ~ event_date * ", endpoint, " + (1|enrollment_site/participant_id)")), 
                    data = plotDF)
      }
    }
    plotDF$yhat <- stats::predict(fit, newdata = plotDF)
    ### Fixed Effect Intercept Matrix List
    if(age_sex){
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
        ann_text <- data.frame(event_date = sex_age_intercept_x,
                               value = unlist(c(active_intercepts[1,], active_intercepts[2,])), lab = "Text",
                               trajectory_group = factor(u, levels = c("1","2","3", "4", "5")),
                               sex_age_group = c(paste0(rownames(active_intercepts)[1], colnames(active_intercepts)),
                                                 paste0(rownames(active_intercepts)[2], colnames(active_intercepts))))
        ann_text$sex_age_group <- gsub("\\(Intercept\\)", paste0("discretized_admit_age_quantile",
                                                                 levels(plotDF$discretized_admit_age_quantile)[1]), ann_text$sex_age_group)
        ann_text$sex_age_group <- gsub("trajectory_group.", paste0("discretized_admit_age_quantile",
                                                                   levels(plotDF$discretized_admit_age_quantile)[1]), ann_text$sex_age_group)
        ann_text$event_date <- ifelse(grepl(levels(plotDF$sex)[1], ann_text$sex_age_group),sex_age_intercept_x+sex_bar_seperator_dist, sex_age_intercept_x)
        intercept_list <- c(intercept_list, list(ann_text))
        names(intercept_list)[length(intercept_list)] <- paste0("trajectory_group", u)
      }
    }
    pred <- ggpredict(fit, c("event_date", endpoint), type = "fixed",
                      ci.lvl = CI_interval)
    colnames(pred) <- gsub("x", "event_date", gsub("group", endpoint, colnames(pred)))
    
    p1 <- ggplot(data = plotDF, mapping = aes(x = event_date, y = value)) 
    if(individual_paths){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = value), color ='gray', alpha = individual_paths_alpha)
    }
    if(individual_trendlines){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = yhat),
                           alpha = individual_trendline_alpha,
                           size = individual_trendline_size,
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
                      group_by(.data[[endpoint]]) %>% summarize(max = max(event_date)), by = c({{endpoint}})) %>%
          filter(event_date <= max)
      }
      
      if(group_trendline_color == "endpoint"){
        p1 <- p1 + geom_line(inherit.aes = F, data = pred, mapping = aes_string(group = endpoint, y = "predicted", x="event_date",
                                           color = endpoint), size = group_trendline_line_width)
      } else {
        p1 <- p1 + geom_line(inherit.aes = F, data = pred, mapping = aes_string(group = endpoint, y = "predicted", x="event_date"),
                             color = group_trendline_color, size = group_trendline_line_width)
      }
    }
    
    
    if(CI){
      if(age_sex){
        p1 <- p1 + geom_ribbon(inherit.aes = F, data = pred,
                               aes_string(x="event_date", ymin = "conf.low", ymax = "conf.high", fill = endpoint), alpha = CI_alpha)
      } else {
        p1 <- p1 + geom_ribbon(inherit.aes = F, data = pred,
                               aes_string(x="event_date", ymin = "conf.low", ymax = "conf.high", fill = endpoint), alpha = CI_alpha)
      }
      
    }
    ##### Small intercept bars for the different age/sex groups
    if(sex_age_intercept_bars){
      if(age_sex){
        p1 <- p1 + coord_cartesian(xlim=c(0,max(plotDF$event_date)), clip = "off")
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
    
  } else if (modelType == "smoothSpline"){
    #### SMOOTH SPLINE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(age_sex){
      formula_use <- formula(paste0("value~ s(event_date, bs = 'cr')+
              s(event_date, bs = 'cr', by =", endpoint, ")+", endpoint,
                                    "+ sex + discretized_admit_age_quantile"))
    } else {
      formula_use <- formula(paste0("value~ 
          s(event_date, bs = 'cr')+s(event_date, bs = 'cr', by =", endpoint, ")+",endpoint))
    }
    fit <- gamm4::gamm4(formula_use, data = plotDF, random = ~(1|enrollment_site/participant_id))
    plotDF$yhat <- predict(fit$mer)
    
    pred <- ggpredict(fit, c("event_date", endpoint), type = "fixed",
                      ci.lvl = CI_interval)
    colnames(pred) <- gsub("x", "event_date", gsub("group", endpoint, colnames(pred)))
    
    p1 <-ggplot(data = plotDF, mapping = aes(x = event_date, y = value)) 
    if(individual_paths){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = value), color ='gray', alpha = individual_paths_alpha)
    }
    if(individual_trendlines){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = yhat),
                           alpha = individual_trendline_alpha,
                           size = individual_trendline_size,
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
                      group_by(.data[[endpoint]]) %>% summarize(max = max(event_date)), by = c({{endpoint}})) %>%
          filter(event_date <= max)
      }
      
      if(group_trendline_color == "endpoint"){
        p1 <- p1 + geom_line(inherit.aes = F, data = pred, mapping = aes_string(group = endpoint, y = "predicted", x="event_date",
                                             color = endpoint), size = group_trendline_line_width)
      } else {
        p1 <- p1 + geom_line(inherit.aes = F, data = pred, mapping = aes_string(group = endpoint, y = "predicted", x="event_date"),
                             color = group_trendline_color, size = group_trendline_line_width)
      }
    }
    
    if(knot_lines){
      p1 <- p1 + geom_vline(xintercept = knots, linetype = 2)
    }
    if(CI){
      p1 <- p1 + geom_ribbon(inherit.aes = F, data = pred,
                             aes_string(x="event_date", ymin = "conf.low", ymax = "conf.high", fill = endpoint), alpha = CI_alpha)
    }
  } else if(modelType =="fixedKnots") {
    ################ FIXED KNOTS ################
    if(missing(fit)){
      if(age_sex){
        fit_formula <- formula(paste0("value  ~ splines::bs(event_date, knots = knots) * ",
                                      endpoint, "+ sex + discretized_admit_age_quantile + (1|enrollment_site/participant_id)"))
        fit <- lmer(fit_formula, data = plotDF)
        plotDFmod <- plotDF
        plotDFmod$sex <- levels(plotDFmod$sex)[1]
        plotDFmod$discretized_admit_age_quantile <- levels(plotDFmod$discretized_admit_age_quantile)[1]
      } else {
        fit_formula <- formula(paste0("value  ~ splines::bs(event_date, knots = knots) * ",
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
      #     # fit_formula2 <- formula(paste0("value  ~ bs(event_date, knots = knots) * ",
      #     #                               endpoint, "+ (1|enrollment_site/participant_id)"))
      #     fit_formula2 <- formula(paste0("value  ~ bs(event_date) * ",
      #                                   endpoint, "+ (1|enrollment_site/participant_id)"))
      #     fit2 <- lmer(fit_formula2, data = plotDF)
      #     pred <- ggpredict(fit2, terms = c("event_date [all]", endpoint), type = "fixed",
      #                       ci.lvl = CI_interval)
      #     colnames(pred) <- gsub("x", "event_date", gsub("group", endpoint, colnames(pred)))
      #   }else{
      #     #plotDF <- as.data.frame(plotDF)
      #     pred <- ggeffect(fit, terms = c("event_date [all]", endpoint), type = "fixed",
      #                       ci.lvl = CI_interval)
      #     colnames(pred) <- gsub("x", "event_date", gsub("group", endpoint, colnames(pred)))
      #   }
    }
    p1 <-ggplot(data = plotDF, mapping = aes(x = event_date, y = value)) 
    
    if(individual_paths){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = value), color ='gray', alpha = individual_paths_alpha)
    }
    if(individual_trendlines){
      p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = yhat),
                           alpha = individual_trendline_alpha,
                           size = individual_trendline_size,
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
                      group_by(.data[[endpoint]]) %>% summarize(max = max(event_date)), by = c({{endpoint}})) %>%
          filter(event_date <= max)
      }
      
      if(group_trendline_color == "endpoint"){
        p1 <- p1 + geom_line(mapping = aes_string(group = endpoint, y = "yhat_wo_re", color = endpoint),
                             alpha = 1, size = group_trendline_line_width)
      } else {
        p1 <- p1 + geom_line(mapping = aes_string(group = endpoint, y = "yhat_wo_re"),
                             alpha = 1, color = group_trendline_color, size = group_trendline_line_width)
      }
    }
    if(knot_lines){
      p1 <- p1 + geom_vline(xintercept = knots, linetype = 2)
    }
    # if(CI){ ## In Development and Testing Still
    #   p1 <- p1 + geom_ribbon(inherit.aes = F, data = pred,
    #                          aes_string(x="event_date", ymin = "conf.low", ymax = "conf.high", fill = endpoint), alpha = .2)
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
    p1 <- p1 + scale_color_manual(values=c("#639A21", "#39828C", "#6371AD", "#BD7D31", "#9C3418"))
    p1 <- p1 + scale_fill_manual(values=c("#639A21", "#39828C", "#6371AD", "#BD7D31", "#9C3418"))
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

### Simplified plot_model function for non-pairwise compairsons and only linear modeling
longitudinal_plot <- function(plotDF, endpoint = "trajectory_group", model_formula,
                              colors, title, xlabel, ylabel, font,
                              group_trendline = T, group_trendline_color = "blue", group_trendline_line_width = 1.5,
                              individual_points = T, individual_points_size = 1, individual_points_alpha = 0.8,
                              individual_paths = T, individual_paths_size = 1, individual_paths_alpha = 0.4){
  
  ##### Modeling General Trendlines
  if(missing(model_formula)){
    formula_lmer <- formula(paste0("value ~ event_date * ", endpoint, " + sex + discretized_admit_age_quantile + (1|enrollment_site/participant_id)"))
    
  } else {
    formula_lmer <- formula(paste0(model_formula))
  }
  plotDFmod <- plotDF
  plotDFmod$sex <- "Female"
  plotDFmod$discretized_admit_age_quantile <- unique(plotDFmod$discretized_admit_age_quantile)[1]
  fit <- lmer(formula_lmer, data = plotDF)
  plotDF$yhat <- stats::predict(fit, newdata = plotDF)
  plotDF$yhat_wo_re <- stats::predict(fit, re.form=NA, newdata=plotDFmod)
  ################ Plotting ################ 
  p1 <- ggplot(data = plotDF, mapping = aes(x = event_date, y = value)) 
  if(individual_paths){
    p1 <- p1 + geom_line(mapping = aes(group = participant_id, y = value), color ='gray', alpha = individual_paths_alpha)
  }
  if(individual_points){
    p1 <- p1 + geom_point(mapping = aes_string(color = endpoint), alpha = individual_points_alpha,
                          size = individual_points_size)
  }
  if(group_trendline){
    p1 <- p1 + geom_line(mapping = aes_string(group = endpoint, y = "yhat_wo_re"),
                         color = group_trendline_color, size = group_trendline_line_width)
  }
  p1 <- p1 + #facet_grid(facets = get(facet_rows)~get(endpoint)) + 
    theme_bw() + 
    labs(y = unique(plotDF$name)) + 
    theme(legend.pos = "none") 
  if(missing(title)){
    p1 <- p1 + ggtitle(paste0(unique(plotDF$name), " by ", endpoint))
  } else {
    p1 <- p1 + ggtitle(title)
  }
  if(!(missing(font))){
    p1 <- p1 + theme(text=element_text(family=font))
  }
  if(!(missing(colors))){
    p1 <- p1 + scale_color_manual(values=colors)
  } else if (endpoint == "trajectory_group"){
    p1 <- p1 + scale_color_manual(values=c("#639A21", "#39828C", "#BD7D31", "#6371AD", "#9C3418"))
  }
  if(!(missing(xlabel))){
    p1 <- p1 + xlab(xlabel)
  }
  if(!(missing(ylabel))){
    p1 <- p1 + ylab(ylabel)
  }
  return(p1)
}

##code for visualize the results
plot_result_list <- function(result_list, 
                             what_to_plot
){
  result_list <- result_list[result_list[[what_to_plot]] <= 0.05, ]
  
  pright <- result_list%>%
    group_by(type)%>%
    summarise(total = n())%>%
    ggplot(aes(y = total, x = type))+
    geom_bar( stat = "identity", alpha = 0.5)+
    geom_text(aes(label = total), alpha = 0.9)+
    labs(y = "", x = "")+
    coord_flip()+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "none")
  
  size_tit <- paste0("-log10(", what_to_plot, ")") 
  pmain <- result_list%>%
    ggplot(aes(name, type))+
    geom_point(aes( size = -log10(get(what_to_plot))), shape = 21)+
    geom_line(aes(group = name), alpha = 0.6)+
    labs(size = size_tit, x = "", y = "")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom",
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  ptop <- result_list%>%
    group_by(name)%>%
    summarise(total = n())%>%
    ggplot()+
    geom_bar(aes(x = name, y = total), stat = "identity", position = "stack", alpha = 0.5)+
    labs(y = "", x = "")+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "none")
  
  return((ptop + plot_spacer() + pmain + pright) +  plot_layout(widths = c(4.3, 1), height = c(1, 2)))
  
}
