# About 

This repo houses code for the analyses in the manuscript "Minimalistic Transcriptomic Signatures Permit Accurate Early Prediction of COVID-19 Mortality." The repository includes scripts for analyzing data and developing classifiers on the IMPACC cohort. It includes scripts for analyzing data from the EARLI cohort. It includes scripts for analyzing data from the COMET cohort for subsequent validation. Full analysis details are available in the publication Methods section. 

# Files

## Datasets

Found in `Data.zip`

- `strandedness_metadata.csv`: Contains the strandedness information for the COMET data.

- `unique_COMET_patients.csv`: Contains the unique COMET (non-IMPACC) participants from the COMET data. 

- `metadata_COMET_validation_simplified.csv`: Metadata for the COMET cohort.

- `genecounts`: Genecounts for the COMET data from each batch.

- `data_files`: Metadata and genecounts for the EARLI data. 

- `top_DE_EARLI.csv`: EARLI DE output.

- `pbmc_ct_age_sex_top_table.xlsx`: IMPACC PBMC DE output.

## Scripts

### IMPACC SCRIPTS

These can generally be run in the order of numbers. They require data from the IMPACC server as detailed in the methods section of the manuscript. 

- `00_OLINK.R` prepares the OLINK proteomic data for analysis. It performes sample collection and QC and a train-test split. This should be run after the 01 scripts. 

- `00_relationship_CT_rpM`: Prepares the cycle threshold data and performs imputation as detailed in the methods section of the manuscript. 

- `01_OLINK.R`: Trains and tests a model using IL-6 to predict trajectory outcomes. 

- `01_pbmc_nasal.R`: Prepares the metadata for the pbmc and nasal swab samples that pass the requirements for inclusion in the analysis. 

- `02_bucket_metadata_exploration.R`: Not relevant for any of the figures or analysis in the manuscript. Provides a general framework for exploring the metadata compiled in `01_pbmc_nasal.R`.

- `03_bucket_DEG.R`: Performes DE analysis on the pbmc and nasal data. Key outputs include top tables with the DE results. 

- `04_assemble_bucket_folds.R`: Splits the nasal and pbmc data into a train-test split. Further divides the training data into balanced folds for CV analysis. 

- `04a_nasal_ct_train_extract_deg_thresh.R`: Determines nasal genes which should be considered as inputs in the LASSO model. Based on the training data only. 

- `04a_pbmc_ct_train_extract_deg_thresh.R`: Determines pbmc genes which should be considered as inputs in the LASSO model. Based on the training data only. 

- `05c_complete_LASSO_deg_filter.R`: Takes the output of `04a_pbmc_ct_train_extract_deg_thresh.R`as input as well as the pbmc metadata. Performs LASSO regression and 5-fold cross validation to determine the gene sets of various lengths. 

- `05d_nasal_complete_LASSO_deg_filter`: Takes the output of `04a_nasal_ct_train_extract_deg_thresh.R`as input as well as the nasal metadata. Performs LASSO regression and 5-fold cross validation to determine the gene sets of various lengths.

- `06_train_eval_a`: Evaluates the pbmc gene sets with additional features CT value and age on the training data using repeated random partitioning. 

- `06_train_eval_a_age`: Evaluates the pbmc gene sets with additional feature, age, on the training data using repeated random partitioning. 

- `06_train_eval_b_no_age_ct.R`: Evaluates the pbmc gene sets with no additional features on the training data using repeated random partitioning. 

- `06_train_eval_c_nasal.R`: Evaluates the nasal gene sets with additional features CT value and age on the training data using repeated random partitioning. 

- `07_test_output_b_filter.R`: Evaluates the selected pbmc gene set with age and CT on the test data. Compares with Sweeney et. al. Evaluates the OLAH + age + CT classifier. 

- `07_test_output_c_nasal_filter.R`: Evaluates the selected nasal gene set with age and CT on the test data. Evaluates the OLAH + age + CT classifier. 

- `boxplots.R`: Takes in the pbmc and nasal metadata and creates boxplots comparing the distribution of age and CT value between the two trajectory groups. 

- `pathways.R`: Takes in the DE results from `03_bucket_DEG.R`and performs pathway analysis. Generates pathway plots. 

- `volcanoes.R`: Takes in the DE results from `03_bucket_DEG.R`and generates volcano plots and the log-log plot between pbmc and nasal DE results.

### EARLI SCRIPTS
These scripts compare the IMPACC results with data from a Sepsis cohort. 

- `01_DE_EARLI.R`: Performs QC and DE analysis on the EARLI data. Outputs `top_DE_EARLI.csv`

- `fig3a_combined_logfc_earli_pbmc.R`: Creates the log-log plot and the venn diagram plot. Takes in `top_DE_EARLI.csv` and the pbmc DE top table from `03_bucket_DEG.R` as inputs. 

- `de_overlap.R`: Creates the Sepsis-PBMC pathway overlap plot from the DE data.

### COMET SCRIPTS
These scripts perform DE analysis on the COMET data. They evaluate the classifiers on the COMET data as well. 

- `COMET_DE_unstrand.R`: Analyzes the COMET data. Performs DE analysis and creates a train-test split. Additionally, creates boxplots and pathway plots.

- `02_train_test_comet_output.R`: Evaluates the classifiers on the COMET data. Performs out-of-fold ROC analysis. Performs vaccinated stratified analysis.






