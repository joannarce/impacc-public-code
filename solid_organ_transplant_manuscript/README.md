# README #

## Introduction

This repository contains the code used in the IMPACC solid organ transplant paper. The code is organized by omics data. The code for each figure is indiciated below.

The paths and patient IDs used in the code were intended for use with the data on IMPACC's own internal server. Therefore, they are different from the deposited data.

## Nasal swab SARS-CoV-2 viral abundance

Figure 2a: ./nasal_viral_load/nasal_viral_load.Rmd/viralload_visit1.pdf
<br>
Figure 2b: ./nasal_viral_load/nasal_viral_load.Rmd/viralload_organ_visit1.pdf
<br>
Figure 2c: ./nasal_viral_load/nasal_viral_load.Rmd/viralload_longitudinal_gamm4.pdf

Supplementary Figure 1: ./nasal_viral_load/nasal_viral_load.Rmd/viralload_mngs-vs-qpcr.pdf

## Serum antibodies

Figure 3C: ./serum_antibodies/serum_antibodies.Rmd #line-189

Figure 3D: ./serum_antibodies/serum_antibodies.Rmd #line-237

Supplementary Figure 2: ./serum_antibodies/serum_antibodies.Rmd #line-301

## Blood CyTOF

Figure 3A: ./blood_cytof/blood_cytof.Rmd #line-230

Figure 3B: ./blood_cytof/blood_cytof.Rmd #line-268

Figure 7A: ./blood_cytof/blood_cytof.Rmd #line-391

Supplementary Figure 3: ./blood_cytof/blood_cytof.Rmd #line-306

## Serum soluble inflammation proteins (Olink)

Figure 4a: ./serum_protein_olink/serum_protein_olink.Rmd/olink_visit1_significant_no-viralload.pdf
<br>
Figure 4b: ./serum_protein_olink/serum_protein_olink.Rmd/olink_visit1_significant_examples.pdf
<br>
Figure 4c: ./serum_protein_olink/serum_protein_olink.Rmd/olink_longitudinal_CXCL11.pdf

Figure 7b: ./nasal_viral_load/nasal_viral_load.Rmd/olink_visit1_TG.pdf

Supplementary Figure 5a: ./serum_protein_olink/serum_protein_olink.Rmd/olink_visit1_significant_yes-vs-no-viralload.pdf
<br>
Supplementary Figure 5b: ./serum_protein_olink/serum_protein_olink.Rmd/olink_visit1_all_vs_viralload.pdf
<br>
Supplementary Figure 5c: ./serum_protein_olink/serum_protein_olink.Rmd/olink_visit1_CXCL8_vs_viralload.pdf
<br>
Supplementary Figure 5d: ./serum_protein_olink/serum_protein_olink.Rmd/olink_longitudinal_examples_viralload.pdf

## PBMC transcriptomics

Figure 5A: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-311

Figure 5B: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-273

Figure 5C: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-323

Figure 5D: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-677

Figure 5E: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-738

Figure 7C: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-423

Supplementary Figure 5A: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-273

Supplementary Figure 5B: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-677

Supplementary Figure 7A: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-506

Supplementary Figure 7C: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-515

## Nasal transcriptomics

Figure 6A: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-273

Figure 6B: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-634

Figure 6C: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-695

Figure 5D: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line

Figure 5E: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-

Figure 7D: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-379

Supplementary Figure 5C: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-273

Supplementary Figure 5D: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-634

Supplementary Figure 7B: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-463

Supplementary Figure 7D: ./pbmc_transcriptomics/pbmc_transcriptomics.Rmd #line-472

## Required software

The codes were run in R v4.0.3. Below is the sessionInfo():

```
R version 4.0.3 (2020-10-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8   
 [6] LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C        
[11] LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] gamm4_0.2-6     mgcv_1.8-33     nlme_3.1-149    lme4_1.1-34     Matrix_1.2-18   patchwork_1.1.2 ggplot2_3.4.0     dplyr_1.0.10      pheatmap_1.0-12     clusterProfiler_4.6-2     ReactomePA_1.42-0     RColorBrewer_1.1-3      biomaRt_2.54-1      edgeR_3.40-2      limma_3.54-2

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.9       nloptr_2.0.3     pillar_1.8.1     compiler_4.0.3   vipor_0.4.5      tools_4.0.3     
 [7] boot_1.3-25      lifecycle_1.0.3  tibble_3.1.8     gtable_0.3.1     lattice_0.20-41  pkgconfig_2.0.3 
[13] rlang_1.0.6      DBI_1.1.3        cli_3.6.0        rstudioapi_0.14  beeswarm_0.4.0   xfun_0.39       
[19] withr_2.5.0      knitr_1.43       generics_0.1.3   vctrs_0.5.1      ggeffects_1.1.4  grid_4.0.3      
[25] tidyselect_1.2.0 glue_1.6.2       R6_2.5.1         fansi_1.0.4      ggbeeswarm_0.7.2 minqa_1.2.5     
[31] magrittr_2.0.3   scales_1.2.1     splines_4.0.3    MASS_7.3-53      assertthat_0.2.1 colorspace_2.0-3
[37] utf8_1.2.3       munsell_0.5.0   
```

Installing the required packages takes approximately 30 minutes.