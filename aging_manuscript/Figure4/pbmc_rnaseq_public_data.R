# Here I will analyze the DE genes from https://doi.org/10.1038/ncomms9570, Supplementary Data 1
# And I will compare their enriched pathways to IMPACC's

# Source: 20231017_public_data.R

library(tidyverse)
library(limma)
library(fgsea)

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0,1,1,0), "cm")
  )

# Load data -----
## Import IMPACC results -----
# DE result (from Figure4/pbmc_rnaseq.Rmd)
res.impacc <- read.csv(
  "bulk_RNAseq_pbmc_visit1_DE_no-viralload.csv",
  row.names=1
)

# GSEA result (from Figure4/pbmc_rnaseq.Rmd)
reactome.gsea.res.impacc <- data.table::fread(
  "bulk_RNAseq_pbmc_visit1_admit-age_sex-ordinalscore_reactome-gsea_fdr1.0.tsv",
  sep="\t", sep2=c("", " ", "")
)

## Import healthy control meta-analysis result -----
# Download from https://doi.org/10.1038/ncomms9570, Supplementary Data 1
res.hc <- read.csv(
  "Supp Data 1_main_table.csv"
)
# Adjusted P
res.hc$Padj <- p.adjust(res.hc$P, method = "BH")

# Remove duplicated gene symbol, keep the one that has the lowest adjusted P value
res.hc <- res.hc[order(res.hc$Padj, decreasing=FALSE),]
res.hc <- res.hc[!duplicated(res.hc$NEW.Gene.ID),]

## Reactome pathways -----
reactome <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = "C2", subcategory = "CP:REACTOME")
reactome.list <- split(x = reactome$gene_symbol, f = reactome$gs_name)
reactome.list.ensembl <- split(x = reactome$ensembl_gene, f = reactome$gs_name)

# Compare IMPACC with healthy controls -----
## fGSEA Reactome on z-score -----
gene.ranks <- res.hc$Zscore
names(gene.ranks) <- res.hc$NEW.Gene.ID
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
reactome.gsea.res.hc <- fgseaMultilevel(
  pathways = reactome.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

# Plot
ggplot(data = reactome.gsea.res.hc, aes(x=reorder(pathway,NES), y=NES)) +
  geom_segment(aes(xend=reorder(pathway,NES)), yend=0) +
  geom_point(aes(fill=as.factor(NES>0), color=as.factor(padj<0.05)),
             pch=21, size=2, stroke=1) +
  scale_color_manual(values=c("white","black")) +
  labs(x="reactome pathways", y="NES", 
       title="Pathways upregulated with age") +
  ylim(-3,3) +
  coord_flip() + theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0,1,1,0), "cm"),
    legend.key = element_rect(fill = "gray90")
  )

## Check overlapping pathways -----
print(sprintf("Pathways up with age in HC: %d",
              sum(reactome.gsea.res.hc$padj<0.05 & reactome.gsea.res.hc$NES>0)))
print(sprintf("Pathways up with age in in IMPACC: %d",
              sum(reactome.gsea.res.impacc$padj<0.05 & reactome.gsea.res.impacc$NES>0)))
print(sprintf("Pathways down with age in HC: %d",
              sum(reactome.gsea.res.hc$padj<0.05 & reactome.gsea.res.hc$NES<0)))
print(sprintf("Pathways down with age in in IMPACC: %d",
              sum(reactome.gsea.res.impacc$padj<0.05 & reactome.gsea.res.impacc$NES<0)))

# List of pathways up in HC & IMPACC
up.both <- intersect(
  reactome.gsea.res.hc[reactome.gsea.res.hc$padj<0.05 & reactome.gsea.res.hc$NES>0,"pathway"],
  reactome.gsea.res.impacc[reactome.gsea.res.impacc$padj<0.05 & reactome.gsea.res.impacc$NES>0,"pathway"]
)
# List of pathways up in IMPACC only
up.impacc <- setdiff(
  reactome.gsea.res.impacc[reactome.gsea.res.impacc$padj<0.05 & reactome.gsea.res.impacc$NES>0,"pathway"],
  reactome.gsea.res.hc[reactome.gsea.res.hc$padj<0.05 & reactome.gsea.res.hc$NES>0,"pathway"]
)
# List of pathways up in HC only
up.hc <- setdiff(
  reactome.gsea.res.hc[reactome.gsea.res.hc$padj<0.05 & reactome.gsea.res.hc$NES>0,"pathway"],
  reactome.gsea.res.impacc[reactome.gsea.res.impacc$padj<0.05 & reactome.gsea.res.impacc$NES>0,"pathway"]
)

# List of pathways down in HC & IMPACC
down.both <- intersect(
  reactome.gsea.res.hc[reactome.gsea.res.hc$padj<0.05 & reactome.gsea.res.hc$NES<0,"pathway"],
  reactome.gsea.res.impacc[reactome.gsea.res.impacc$padj<0.05 & reactome.gsea.res.impacc$NES<0,"pathway"]
)
# List of pathways down in IMPACC only
down.impacc <- setdiff(
  reactome.gsea.res.impacc[reactome.gsea.res.impacc$padj<0.05 & reactome.gsea.res.impacc$NES<0,"pathway"],
  reactome.gsea.res.hc[reactome.gsea.res.hc$padj<0.05 & reactome.gsea.res.hc$NES<0,"pathway"]
)
# List of pathways down in HC only
down.hc <- setdiff(
  reactome.gsea.res.hc[reactome.gsea.res.hc$padj<0.05 & reactome.gsea.res.hc$NES<0,"pathway"],
  reactome.gsea.res.impacc[reactome.gsea.res.impacc$padj<0.05 & reactome.gsea.res.impacc$NES<0,"pathway"]
)

# Merge gsea results
merged.gsea <- rbind(
  reactome.gsea.res.hc %>% mutate(type="HC"),
  reactome.gsea.res.impacc %>% mutate(type="IMPACC")
)

# Plot select pathways
up.pathways <- c(
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING","REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_INTERLEUKIN_2_FAMILY_SIGNALING","REACTOME_INTERLEUKIN_10_SIGNALING",
  "REACTOME_NEUTROPHIL_DEGRANULATION","REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
  "REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE",
  "REACTOME_CASPASE_ACTIVATION_VIA_DEATH_RECEPTORS_IN_THE_PRESENCE_OF_LIGAND",
  "REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
  "REACTOME_DEATH_RECEPTOR_SIGNALLING","REACTOME_THE_NLRP3_INFLAMMASOME",
  "REACTOME_TRAF6_MEDIATED_IRF7_ACTIVATION"
)
merged.gsea.toplot <- merged.gsea %>%
  subset(pathway %in% up.pathways)

# Reformat pathway names
merged.gsea.toplot$pathway <- gsub("^REACTOME_","", merged.gsea.toplot$pathway)
merged.gsea.toplot$pathway <- gsub("_"," ", merged.gsea.toplot$pathway)
merged.gsea.toplot$pathway <- stringr::str_to_sentence(merged.gsea.toplot$pathway)

# Fix capitalization of pathway names
my.dict <- list("tlr"="TLR",
                "Toll like receptor "="",
                "Caspase activation via death receptors in the presence of ligand" = "Caspase activation via depth receptors",
                "Caspase activation via extrinsic apoptotic signalling pathway" = "Caspase activation via apoptotic signalling")
for (i in names(my.dict)) {
  merged.gsea.toplot$pathway <- gsub(i, my.dict[[i]],
                                       merged.gsea.toplot$pathway)
}

# Dotplot
p <- ggplot(data = merged.gsea.toplot,
            aes(x = factor(type, levels=c("IMPACC","HC")), 
                y = reorder(pathway, NES),
                color = NES)) +
  geom_point(color = "white", size = 3) +
  geom_point(shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = . %>% subset(padj < 0.05), size = 4) +
  scale_y_discrete(limits=rev) +
  scale_color_gradientn(
    colors = c("dodgerblue4","#f7f7f7","darkred"),
    limits = c(-4,4),
    breaks = c(-4,0,4),
    name = "") +
  labs(
    x = "", y = "",
    title = "Pathways upregulated\nwith age"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black", angle=45, hjust=1),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0,1,1,0), "cm"),
    panel.spacing = unit(0.5, "cm"),
    legend.key = element_rect(fill = "gray90"),
    panel.grid.minor=element_blank()
  )

## Check overlapping DE genes -----
# Filter IMPACC DE results to remove duplicated gene symbols, and keep the one
# with the lowest P value
res.impacc2 <- res.impacc[order(res.impacc$adj.P.Val, decreasing=FALSE),]
res.impacc2 <- res.impacc2[!duplicated(res.impacc2$gene_name),]

# Merge IMPACC DE result with the healthy control result
DE.impacc.hc <- merge(
  res.impacc2[,c("gene_id","gene_name","logFC","adj.P.Val")],
  res.hc[,c("NEW.Gene.ID","Zscore","Padj")],
  by.x="gene_name", by.y="NEW.Gene.ID",
  all=FALSE, sort=FALSE
)

# Figure S5A
# Plot slope
p <- ggplot(data=DE.impacc.hc,
            aes(x=logFC, y=Zscore)) +
  geom_point(color="gray", size=0.4) +
  geom_point(data= . %>% subset((adj.P.Val<0.05) & (Padj<0.05)), size=0.4) +
  geom_text(data= . %>%
              subset(gene_name %in% c("IFI30","IRF5","IRF9","TLR2","TLR4","TLR6","TLR8")),
            aes(label=gene_name)) + # genes up in both cohorts
  geom_text(data= . %>%
              subset(gene_name %in% c("CXCR3","CXCR5","HLA-DOA","HLA-DOB")),
            aes(label=gene_name)) + # genes down in both cohorts
  geom_hline(yintercept=0, color="black", 
             linetype="dashed", linewidth=0.6) +
  geom_vline(xintercept=0, color="black",
             linetype="dashed", linewidth=0.6) +
  labs(x="COVID-19 patients: Slope",
       y="Healthy control: Meta-analysis Z-score") +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  "DE_impacc-vs-hc.png",
  plot=p,
  dpi=600,
  width=3.5, height=3.5, units="in"
)

# Genes that go up with age in both studies
up.both.DE <- DE.impacc.hc %>%
  subset((adj.P.Val<0.05) & (Padj<0.05) &
           (logFC>0) & (Zscore>0))
# Genes that go down with age in both studies
down.both.DE <- DE.impacc.hc %>%
  subset((adj.P.Val<0.05) & (Padj<0.05) &
           (logFC<0) & (Zscore<0))
# Genes that go up with age in IMPACC, go down with age in rotterdam
up.impacc.DE <- DE.impacc.hc %>%
  subset((adj.P.Val<0.05) & (Padj<0.05) &
           (logFC>0) & (Zscore<0))
# Genes that go up with age in rotterdam, go down with age in IMPACC
down.impacc.DE <- DE.impacc.hc %>%
  subset((adj.P.Val<0.05) & (Padj<0.05) &
           (logFC<0) & (Zscore>0))
