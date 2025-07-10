# Here, I perform DE and pathway analysis on the four sample bins with genetic information

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

bin_colors <- c("#bc6a84",
                "#59aa54",
                "#9762ca",
                "#999a3e",
                "#ce51a0",
                "#45b0a4",
                "#d2485a",
                "#6b85cd",
                "#c95631",
                "#c88744")
nasal_complete_up <- bin_colors[1]
nasal_only_up <- bin_colors[3]
pbmc_complete_up <- bin_colors[5]
pbmc_only_up <- bin_colors[7]
nasal_complete_down <- bin_colors[2]
nasal_only_down <- bin_colors[4]
pbmc_complete_down <- bin_colors[6]
pbmc_only_down <- bin_colors[8]

## Set up here and output/input directories =====

i_am("IMPACC_071824.Rproj")
here()
output_dir <- here("03_bucket_DEG")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
  
} else{
  print("Output directory already exists")
}

input_dir <- here("01_pbmc_nasal")

## READ IN METADATA ======
pbmc_complete_metadata <- read_csv(here(input_dir, "pbmc_complete_metadata_clin.csv"))
nasal_complete_metadata <- read_csv(here(input_dir, "nasal_complete_metadata_clin.csv"))
nasal_only_metadata <- read_csv(here(input_dir, "nasal_only_metadata_clin.csv"))
pbmc_only_metadata <- read_csv(here(input_dir, "pbmc_only_metadata_clin.csv"))
clinical_meta <- read_csv("../../../data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-individ-locked.csv")

# uncomment if running DE analysis with sex as a factor
clinical_meta <- clinical_meta %>%
  select(participant_id, sex)

clinical_meta$sex <- as.factor(clinical_meta$sex)

## READ IN OMICS ======
pbmc_counts <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-Counts.csv")
pbmc_rowFeature <- read_csv("../../../data/pbmc-transcriptomics/legacy/2023-07-19/pbmc-transcriptomics-RowFeature.csv", skip = 1)
nasal_counts <- read_csv("../../../data/nasal-transcriptomics/legacy/2022-10-19/nasal-transcriptomics-Counts.csv")


## FUNCTION FOR DE AND PATHWAY ANALYSIS ====

#inputs:
# counts is a dataframe with gene ensembl ids, or features, as columns and sample_ids as a single column
# metadata is a dataframe with information about each participant, including the sample_id
# up_color is a string with the color for the fill for dots representing significantly upregulated genes in the volcano plot
# down color is a string with the color for the fill for dots representing significantly upregulated genes in the volcano plot
#outputs: a list
# sig_de is a dataframe with information about the significant DE genes
# volcano is a volcano plot based on the result of the DE analysis
# hallmark is a dataframe consisting of information about each of the 50 hallmark pathways,
#including their p value, NES, and leading edge genes
# annotated_top is a dataframe with information about every gene in the DE analysis,
#including whether or not it was significant, its logFC value, its gene name, its chr location,
#and its p_adj value
de_and_pathway <- function(counts, metadata, up_color, down_color, label) {
  counts <- counts %>%
    column_to_rownames(var = "sample_id")
  
  # uncomment to run age analysis
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
  
  ## Get feature info
  
  
  analyze_top <- top_DE %>%
    rownames_to_column(var = "gene_id") %>%
    inner_join(pbmc_rowFeature, by = "gene_id") %>%
    mutate(delabel = ifelse(gene_name %in% c("DAAM2", "MAOA", "S100A12", "RNASE2", "HLA-DQA1", "HLA-DRA", "HLA-DMB", "HLA-DMA", "FCER1A", "LRRN3", "DEFT1P", "PCSK9"), 
                            as.character(gene_name), NA))
  
  sig_DE <-  topTable(vfit, coef = "group1", sort.by = "none",
                      number = Inf, p.value = 0.05)
  
  a <- ggplot(analyze_top, aes(x = logFC, y = -log10(adj.P.Val), label = delabel)) +
    geom_point(aes(color = diffexpressed), size = 0.8) +
    scale_color_manual(name = "", values = c("grey", up_color, down_color)) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    ggtitle("pbmc") + xlab("logFC") +
    ylab("-log10(Padj)") + my.theme  + 
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    geom_label_repel(max.overlaps = Inf, min.segment.length = 0, box.padding = 0.5, seed = 21,
                     fill = alpha(c("white"),0.5)) + ggtitle(paste0(label))
  
  # Compute gene ranks
  
  gene.ranks <- -log10(top_DE$P.Value) * sign(top_DE$logFC)
  names(gene.ranks) <- rownames(top_DE)
  gene.ranks <- sort(gene.ranks, decreasing = TRUE)
  
  # Load hallmark pathways
  hallmark <- msigdbr::msigdbr(species = "Homo sapiens",
                               category = "H")
  hallmark.list <- split(x = hallmark$ensembl_gene, f = hallmark$gs_name)
  
  # fgsea
  set.seed(19) ## Reproducibility
  hallmark.gsea <- fgseaMultilevel(
    pathways = hallmark.list,
    stats = gene.ranks,
    minSize = 15,
    maxSize = 500,
    nproc = 1
  )
  
  # export
  write_csv(analyze_top, here(output_dir, paste0(label, "_age_sex_top_table.csv")))
  write_csv(sig_DE, here(output_dir, paste0(label, "_age_sex_sig_de.csv")))
  qsave(hallmark.gsea, here(output_dir, paste0(label, "_age_sex_hallmark.qs")))
  return(list(sig_de = sig_DE, volcano = a, hallmark = hallmark.gsea, annotated_top = analyze_top))
}

pbmc_only <- de_and_pathway(pbmc_counts, pbmc_only_metadata, pbmc_only_up, pbmc_only_down, "pbmc_only_age_sex")
nasal_only <- de_and_pathway(nasal_counts, nasal_only_metadata, nasal_only_up, nasal_only_down, "nasal_only_age_sex")
pbmc_complete <- de_and_pathway(pbmc_counts, pbmc_complete_metadata, pbmc_complete_up, pbmc_complete_down, "pbmc_ct_age_sex")
nasal_complete <- de_and_pathway(nasal_counts, nasal_complete_metadata, nasal_complete_up, nasal_complete_down, "nasal_ct_age_sex")

