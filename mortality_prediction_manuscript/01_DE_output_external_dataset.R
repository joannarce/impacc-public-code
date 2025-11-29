## Here I conduct DE analysis on the external dataset linked by the reviewer

## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("janitor")
library("readr")
library("qs")
library("ggrepel")
library("reshape2")
library(readr)
library(readxl)

## Load packages for analysis ======

library("DESeq2")
library("limma")
library("biomaRt")
library("edgeR")
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

i_am("Additional_dataset.Rproj")


here()
output_dir <- here("01_DE_output")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("data_files")

## READ IN DATA ======

data <- read_tsv(here(input_dir, "counts.tsv"))
annot <- read_tsv(here(input_dir, "annot.tsv"))
meta <- read_excel(here(input_dir, "meta.xlsx"))
dict <- read_excel(here(input_dir, "dict.xlsx"))
library(tidyverse)

# 1. List all .txt.gz count files from the extracted tar directory
count_files <- list.files(here(input_dir, "GSE217948_RAW"), pattern = "\\.txt\\.gz$", full.names = TRUE)

# 2. Read each compressed count file and format it
library(purrr)
counts_list <- lapply(count_files, function(f) {
  dat <- read.delim(gzfile(f), header = TRUE, stringsAsFactors = FALSE)
  stopifnot(ncol(dat) >= 2)
  dat <- dat[, 1:2]  # assume 1st column = gene, 2nd column = count
  colnames(dat) <- c("Gene", tools::file_path_sans_ext(basename(f)))
  return(dat)
})

# 3. Merge all data frames on Gene
counts_combined <- purrr::reduce(counts_list, full_join, by = "Gene")

# 4. Clean up: set rownames, remove Gene column
rownames(counts_combined) <- counts_combined$Gene
counts_combined <- counts_combined[, -1]

# 5. Clean up:: column names

colnames(counts_combined) <- sub(".*_(.*?)\\.counts.*", "\\1", colnames(counts_combined))

## METADATA_CLEANING ====

meta <- meta %>%
  dplyr::filter(Status == "COVID-19") %>%
  dplyr::select(c(Patient_ID, Gender, Age, Clinical_outcome))

dict_clean <- dict %>%
  janitor::clean_names() %>%        # converts "Patient Identifier" -> "patient_identifier"
  filter(!is.na(sample_identifier) & !is.na(geo_accession)) %>%
  dplyr::rename(
    Patient_ID = patient_identifier,
    GEO_accession = geo_accession
  ) %>%
  distinct(Patient_ID, .keep_all = T)


# Step 2: Join GEO accession into meta
meta_annotated <- meta %>%
  left_join(dict_clean, by = "Patient_ID")


filtered_meta <- meta_annotated %>%
  filter(Clinical_outcome %in% c("Discharged alive", "Dead"),
         !is.na(Age),
         !is.na(Gender))

# Counts
# Quality control. Keep samples with >= 5,000 genes
# Keep genes with >= 10 counts in >= 20% of samples
counts <- counts_combined %>% t()
counts <- counts[rowSums(counts > 0) >= 5000, ]
counts <- as.data.frame(t(counts))
keep <- rowSums(counts >= 10) >= 0.2*ncol(counts)
counts <- counts[keep, ]


 annot$GeneID <- as.character(annot$GeneID)
 # Ensure rownames are characters
 rownames(counts) <- as.character(rownames(counts))
 
 # Create a mapping: Ensembl ID -> Gene Symbol
 geneid_to_symbol <- setNames(annot$Symbol, annot$EnsemblGeneID)
 
 # Map rownames from Ensembl to symbols
 new_rownames <- geneid_to_symbol[rownames(counts)]
 
 # Remove duplicated gene symbols
 counts <- counts[!duplicated(new_rownames) & !is.na(new_rownames), ]
 
 # Set new rownames (gene symbols)
 rownames(counts) <- new_rownames[!duplicated(new_rownames) & !is.na(new_rownames)]
 
 #rownames(counts) <- new_rownames[!duplicated(new_rownames)]

c <- colnames(counts)
# Select samples for which metadata exists and vice versa
metadata <- filtered_meta %>%
  subset((sample_identifier %in% c)) %>%
  distinct(Patient_ID, .keep_all = T)

counts <- counts[,metadata$sample_identifier]
stopifnot(colnames(counts)==metadata$sample_identifier)

write_csv(metadata, here(output_dir, "cleaned_metadata.csv"))
write_csv(counts, here(output_dir, "cleaned_counts.csv"))

## DE ANALYSIS =====

metadata$Gender <- as.factor(metadata$Gender)
metadata$Clinical_outcome <- factor(metadata$Clinical_outcome,
                                    levels = rev(c("Dead", "Discharged alive")))


metadata$Age <- as.numeric(metadata$Age)
## Make design matrix, controlling for age and Gender
design <- model.matrix(~Clinical_outcome + Age + Gender,
                       data = metadata)

## Normalize with voom, perform DE analysis
vwts <- voom(counts, 
             design = design,
             normalize.method = "quantile",
             plot = T)

vfit <- lmFit(vwts)
vfit <- eBayes(vfit)

## Extract DE genes with topTable

top_DE <- topTable(vfit, coef = "Clinical_outcomeDead", sort.by = "none",
                   number = Inf, p.value = 1)
top_DE$name <- rownames(top_DE)




## Gene expression analysis ======

reactome <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = "C2", 
                             subcategory = "CP:REACTOME")
reactome.list <- split(x = reactome$gene_symbol, f = reactome$gs_description)

gene.ranks <- -log10(top_DE$P.Value) * sign(top_DE$logFC)
names(gene.ranks) <- top_DE$name
gene.ranks <- sort(gene.ranks, decreasing = TRUE)


set.seed(19) ## Reproducibility
reactome.gsea <- fgseaMultilevel(
  pathways = reactome.list,
  stats = gene.ranks,
  minSize = 10,
  maxSize = 500,
  nproc = 1
)


