# Here I analyze the COMET bulk data and perform a train/test split

# Are main-idea commments

## Are details/steps

# Load packages ----

library(dplyr)
library(tibble)
library(limma)
library(DESeq2)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggeffects)
library(here)
library(qs)
library(fgsea)
library(stringr)

here::i_am("COMET.Rproj")

# Read genecounts data -----
counts_25 <- read.csv(here("genecounts/Seq155_COMET_PBMCs_counts.csv"), row.names = 1)
counts_24 <- read.csv(here("genecounts/COMET_PBMC_genecounts.csv"), row.names = 1)

mapping_dict <- setNames(counts_24[[1]], rownames(counts_24))

# Read metadata -----

metadata <- read.csv(here("metadata_COMET_validation_simplified.csv"))
strandedness <- read.csv(here("strandedness_metadata.csv"))



## Clean counts_24 ------

c.24 <- colnames(counts_24)
c.24.D0 <- c.24[grep("D0", c.24)]

# Extract patient name (the number after "HS")
patient_name <- str_extract(c.24, "(?<=HS)\\d+")

# Extract time points (the number after "D")
time_point <- as.numeric(str_extract(c.24, "(?<=D)\\d+"))

# Create a dataframe
df <- data.frame(Sample = c.24, Patient = patient_name, Time = time_point)

# Select the earliest time point for each patient
earliest_samples <- df %>%
  group_by(Patient) %>%
  filter(Time == min(Time)) %>%
  dplyr::filter(Time <= 2) %>%
  ungroup() %>% 
  pull(Sample)


hs_numbers <- as.numeric(gsub(".*HS(\\d+).*", "\\1", earliest_samples)) + 1000
seq_numbers <- as.numeric(gsub(".*RSQ(\\d+).*", "\\1", earliest_samples)) 

c.24.names <- data.frame(Original_ColNames = earliest_samples, 
                          HS_Numbers = hs_numbers, 
                          RSQ_Numbers = seq_numbers)

c.24.names.filtered <- c.24.names %>%
  group_by(HS_Numbers) %>%
  slice_max(RSQ_Numbers, with_ties = FALSE) %>%
  ungroup()

unstranded <- strandedness %>%
  mutate(new_sample_id = str_replace_all(sample_id, "-", ".")) %>%
  mutate(new_sample_id = paste0("MVIR1.", new_sample_id)) %>% 
  left_join(c.24.names.filtered, ., by = c("Original_ColNames" = "new_sample_id"))


c.24.names.filtered <- unstranded
  

## Match counts to metadata

metadata.24 <- metadata %>%
  dplyr::filter(Batch.year == 2024) %>% 
  left_join(c.24.names.filtered, by = c("Patient.ID" = "HS_Numbers")) %>%
  mutate(File.name = Original_ColNames) %>%
  select(-Original_ColNames) %>%
  select(-RSQ_Numbers) %>%
  drop_na(File.name) %>%
  dplyr::filter(File.name %in% colnames(counts_24)) %>%
  dplyr::filter(COVID.Testing.Status == "Confirmed Positive") %>%
  dplyr::filter(salmon_strandedness != "S" | is.na(salmon_strandedness)) %>%
  dplyr::select(-c(sample_id, salmon_strandedness))
  
metadata.25 <- metadata %>%
  dplyr::filter(Batch.year == 2025) %>%
  dplyr::filter(File.name %in% colnames(counts_25))

counts <- cbind(counts_24, counts_25)
meta <- rbind(metadata.24, metadata.25)


non_impacc <- read.csv(here("unique_COMET_patients.csv"))

meta <- meta %>%
  dplyr::filter(File.name %in% non_impacc$File.name)

write.csv(x = counts, file = here("counts_preQC.csv"))
write.csv(x = meta, file = here("meta_preQC.csv"))

## Count QC -------

counts <- counts[, colSums(counts>0) >= 10000]

## Keep genes with >=10 counts in >=20% of samples
keep <- rowSums(counts >= 10) >= 0.2*ncol(counts)
counts <- counts[keep, ]

meta <- meta %>%
  dplyr::filter(File.name %in% colnames(counts))

counts <- counts %>%
  dplyr::select(c(meta$File.name))


## Convert categorial variables into factors

for (i in c("Deceased","Sex.at.Birth", "Batch.year")) {
  meta[,i] <- as.factor(meta[,i])
}

write.csv(meta, here("comet_meta.csv"))
write.csv(counts, here("comet_counts.csv"))

## Set up ggplot2 theme

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )


# DE analysis -----

## Make design matrix, controlling for age, sex
design <- model.matrix(~Deceased + Age.at.Admission + Sex.at.Birth + Batch.year,
                       data = meta)

## limma-voom normalization
vwts <- voom(counts, 
             design = design,
             normalize.method = "quantile",
             plot = T)

# DE calculation -----
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)


## Extract genes
top_DE <- topTable(vfit, coef = "DeceasedYes", sort.by = "none",
                   number = Inf, p.value = 1)

top_DE$Significance <- ifelse(top_DE$adj.P.Val < 0.05, "Significant (FDR < 0.05)", "Not Significant")

top_DE$diffexpressed <- "Not significant"
top_DE$diffexpressed[top_DE$logFC > 0 & top_DE$adj.P.Val < 0.05] <- "Significantly upregulated"
top_DE$diffexpressed[top_DE$logFC < 0 & top_DE$adj.P.Val < 0.05] <- "Significantly downregulated"


top_DE$Gene.name <- mapping_dict[rownames(top_DE)]

table(top_DE$diffexpressed)
write.csv(top_DE, "bulk_COMET_topTable.csv")

pretty_plot <- ggplot(top_DE, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = diffexpressed), size = 0.8, show.legend = FALSE) +
  scale_color_manual(name = "", values = c("grey", "#6b85cd", "#d2485a"),
                     labels = c("Not significant", "Downregulated", "Upregulated")) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  theme(legend.position = NA) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression(-log[10](P[adj]))) + 
  xlim(c(-4, 4)) +
  my.theme

ggsave("comet_DE.svg", pretty_plot, width = 4.5, height = 4.5)
 

# Pathway analysis ------


reactome <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = "C2", subcategory = "CP:REACTOME")
reactome.list <- split(x = reactome$human_ensembl_gene, f = reactome$gs_description)

gene.ranks <- -log10(top_DE$P.Value) * sign(top_DE$logFC)
names(gene.ranks) <- rownames(top_DE)


set.seed(19) 


reactome.gsea <- fgseaMultilevel(
  pathways = reactome.list,
  stats = gene.ranks,
  minSize = 10,
  maxSize = 500,
  nproc = 1)

top_positive <- reactome.gsea %>%
  filter(padj < 0.05) %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(8)

top_negative <- reactome.gsea %>%
  filter(padj < 0.05) %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(8)


pathway_ready <- rbind(top_positive, top_negative)

fwrite(pathway_ready, here("comet_pathway.csv"))

g <- ggplot(data = pathway_ready, 
            aes(x = reorder(pathway, NES), y = NES)) +
  geom_segment(aes(xend = reorder(pathway, NES), yend = 0)) +
  geom_point(aes(fill = ifelse(NES > 0, "Upregulated", "Downregulated")),
             size = 4, shape = 21, stroke = 1) +  
  scale_fill_manual(
    name = "",
    values = c("Upregulated" = "#d2485a", "Downregulated" = "#6b85cd")) +
  my.theme + 
  xlab("") +
  ylim(-3.5, 3.5) +
  ylab("Normalized enrichment score") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))  # Complete box around plot


qsave(pathway_ready, here("pathway_data.qs"))
ggsave(here("pathways.svg"), g, width = 6.8, height = 6)


## genex boxplot
count <- read.csv(here("comet_counts.csv"))
count <- count %>%
  column_to_rownames("X")

count_cpm <- edgeR::cpm(count, log=TRUE)

genes <- c("ENSG00000112149", "ENSG00000146122", "ENSG00000129244", "ENSG00000152463")

count_cpm <- count_cpm %>%
  t(.) %>%
  as.data.frame(.) %>%
  dplyr::select(genes) %>%
  rownames_to_column(var = "File.name") %>%
  inner_join(., meta, by = "File.name")

df <- count_cpm[, c("Deceased", genes)]


plotting_df <- reshape2::melt(df, id.vars = "Deceased", variable.name = "gene", value.name = "expression")
plotting_df$gene <- mapping_dict[match(plotting_df$gene, names(mapping_dict))]



pretty_plot <- ggplot(plotting_df, aes(x = factor(gene, levels = c("CD83", "ATP1B2", "DAAM2", "OLAH")), y = expression, fill = factor(Deceased))) +
  geom_boxplot(position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("#6b85cd", "#d2485a"),
                    labels = c("No" = "No mortality",
                               "Yes" = "Mortality")) +
  labs(x = "Gene", y = "Normalized expression", fill = "Trajectory") +
  my.theme +
  theme(legend.position = "none")

ggsave(here("box_plot_genex_all_data.svg"), pretty_plot, height = 4.5, width = 4.5)
