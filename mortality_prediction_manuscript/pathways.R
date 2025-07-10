## Here I make figures 2b and 2d

## Load packages for data wrangling, manipulation, plotting ======

library("tidyverse")
library("here")
library("janitor")
library("readr")
library("qs")
library("ggrepel")
library(ggplot2)
library("reshape2")

## Load pathway analysis packages
library(fgsea)


## Set up plotting specs ===========

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

## Load in data ================
top_DE <- read_csv("03_bucket_DEG/pbmc_ct_age_sex_age_sex_top_table.csv")

## Load in pathways =============
reactome <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = "C2", 
                             subcategory = "CP:REACTOME")
reactome.list <- split(x = reactome$ensembl_gene, f = reactome$gs_description)

## Pathway analysis (2b) =============


gene.ranks <- -log10(top_DE$P.Value) * sign(top_DE$logFC)
names(gene.ranks) <- top_DE$gene_id
gene.ranks <- sort(gene.ranks, decreasing = TRUE)


set.seed(19) ## Reproducibility
reactome.gsea <- fgseaMultilevel(
  pathways = reactome.list,
  stats = gene.ranks,
  minSize = 10,
  maxSize = 500,
  nproc = 1
)

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

g <- ggplot(data = pathway_ready, 
            aes(x = reorder(pathway, NES), y = NES)) +
  geom_segment(aes(xend = reorder(pathway, NES), yend = 0)) +
  geom_point(aes(fill = ifelse(NES > 0, "Upregulated", "Downregulated")),
             size = 4, shape = 21, stroke = 1) +  
  scale_fill_manual(
    name = "",
    values = c("Upregulated" = "#d2485a", "Downregulated" = "#6b85cd")) +
  my.theme + 
  ylim(c(-3.0, 3.0)) + 
  xlab("") +
  ylab("Normalized enrichment score") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  theme(legened.position = "none")# Complete box around plot

ggsave("p_gsea_pathway_2b.svg", g, width = 8.5, height = 5)
fwrite(pathway_ready, "reactome_results_pbmc.csv")

## Pathway analysis (2d) =============

top_DE.n <- read_csv("03_bucket_DEG/nasal_ct_age_sex_age_sex_top_table.csv")

gene.ranks.n <- -log10(top_DE.n$P.Value) * sign(top_DE.n$logFC)
names(gene.ranks.n) <- top_DE.n$gene_id
gene.ranks.n <- sort(gene.ranks.n, decreasing = TRUE)


set.seed(19) ## Reproducibility
reactome.gsea.n <- fgseaMultilevel(
  pathways = reactome.list,
  stats = gene.ranks.n,
  minSize = 10,
  maxSize = 500,
  nproc = 1
)

top_positive.n <- reactome.gsea.n %>%
  filter(padj < 0.05) %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(8)

top_negative.n <- reactome.gsea.n %>%
  filter(padj < 0.05) %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(8)

pathway_ready.n <- rbind(top_positive.n, top_negative.n)

n <- ggplot(data = pathway_ready.n, 
            aes(x = reorder(pathway, NES), y = NES)) +
  geom_segment(aes(xend = reorder(pathway, NES), yend = 0)) +
  geom_point(aes(fill = ifelse(NES > 0, "Upregulated", "Downregulated")),
             size = 4, shape = 21, stroke = 1) +  
  scale_fill_manual(
    name = "",
    values = c("Upregulated" = "#d2485a", "Downregulated" = "#6b85cd")) +
  my.theme + 
  ylim(c(-3.0, 3.0)) + 
  xlab("") +
  ylab("Normalized enrichment score") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  theme(legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))  # Complete box around plot

fwrite(pathway_ready.n, "reactome_results_NS.csv")

ggsave("n_gsea_pathway_2d.svg", n, width = 8.5, height = 5)
