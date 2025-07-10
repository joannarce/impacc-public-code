# Here, I perform DE analysis on nasal+CT train samples to extract only
# Genes that are significantly DE with logFC >= 1 or <= -1 
# For inputs to the LASSO classiier

## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("janitor")
library("readr")
library("qs")
library("ggrepel")


## Load packages for differential gene expression and pathway analysis ====

library(limma)
library(fgsea)

## Set up plotting specs ====

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")
here()
output_dir <- here("04a_nasal_ct_train_extract_deg_thresh")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("04_assemble_bucket_folds")

## READ IN METADATA ======
nasal_complete_metadata_train <- read_csv(here(input_dir, "nasal_complete_train_folds.csv"))
clinical_meta <- read_csv("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-individ-locked.csv")

# Metadata pre-processing (clinical)

clinical_meta <- clinical_meta %>%
  select(participant_id, sex)

clinical_meta$sex <- as.factor(clinical_meta$sex)

## READ IN OMICS ======
nasal_counts <- read_csv("../../../data/nasal-transcriptomics/legacy/2022-10-19/nasal-transcriptomics-Counts.csv")
pbmc_rowFeature <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-RowFeature.csv", skip = 1)


## FUNCTION FOR DEG extraction ====

#inputs:
# counts is a dataframe with gene ensembl ids, or features, as columns and sample_ids as a single column
# metadata is a dataframe with information about each participant, including the sample_id

#outputs:
# a list of the genes passing the threshold
# a csv of the DEG results

de_and_list <- function(counts, metadata) {
  counts <- counts %>%
    column_to_rownames(var = "sample_id")
  
  metadata <- metadata %>%
    left_join(., clinical_meta, by = "participant_id")
  
  metadata$sex <- as.factor(metadata$sex)
  
  # Quality control. Keep samples with >= 5,000 genes
  # Keep genes with >= 10 counts in >= 20% of samples
  counts <- counts[rowSums(counts > 0) >= 5000, ]
  counts <- as.data.frame(t(counts))
  keep <- rowSums(counts >= 10) >= 0.2*ncol(counts)
  counts <- counts[keep, ]
  
  # Select samples for which metadata exists and vice versa
  metadata <- metadata %>%
    subset((sample_id %in% colnames(counts)))
  
  counts <- counts[,metadata$sample_id]
  stopifnot(colnames(counts)==metadata$sample_id)
  
  # Make trajectory group a factor
  metadata <- metadata %>%
    mutate(group = factor(ifelse(trajectory_group == 5, 1, 0)))
  metadata$group <- as.factor(metadata$group)
  ## Make design matrix, controlling for age and sex
  design <- model.matrix(~group + admit_age + sex,
                         data = metadata)
  
  ## Normalize with voom, perform DE analysis
  vwts <- voom(counts, 
               design = design,
               normalize.method = "quantile",
               plot = T)
  
  vfit <- lmFit(vwts)
  vfit <- eBayes(vfit)
  
  ## Extract DE genes with topTable
  
  top_DE <- topTable(vfit, coef = "group1", sort.by = "none",
                     number = Inf, p.value = 1)
  
  top_DE$diffexpressed <- "Not significant"
  top_DE$diffexpressed[top_DE$logFC > 0 & top_DE$adj.P.Val < 0.05] <- "Significantly upregulated"
  top_DE$diffexpressed[top_DE$logFC < 0 & top_DE$adj.P.Val < 0.05] <- "Significantly downregulated"
  
  top_DE$threshold[top_DE$logFC <= -1 & top_DE$adj.P.Val < 0.05] <- "Thresh"
  top_DE$threshold[top_DE$logFC >= 1 & top_DE$adj.P.Val < 0.05] <- "Thresh"
  
  ## Get feature info
  
  analyze_top <- top_DE %>%
    rownames_to_column(var = "gene_id") %>%
    inner_join(pbmc_rowFeature, by = "gene_id")
  
  sig_DE <-  topTable(vfit, coef = "group1", sort.by = "none",
                      number = Inf, p.value = 0.05)
  
  write_csv(analyze_top, here(output_dir, "de_results_train_only_nasal_CT.csv"))
  
  ## Save threshold genes
  
  geneListThreshold <- analyze_top %>%
    filter(threshold == "Thresh") %>%
    pull("gene_id") %>%
    as.list()
  
  qsave(geneListThreshold, here(output_dir, "nasal_gene_list_threshold.qs"))
  
  
}

de_and_list(nasal_counts, nasal_complete_metadata_train)
