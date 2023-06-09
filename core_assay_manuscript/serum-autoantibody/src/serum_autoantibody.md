Serum auto-antibody against IFN’s data analysis for core-assay
manuscript
================
07 June, 2023

### Load libraries

``` r
suppressPackageStartupMessages(library(package = "knitr"))
suppressPackageStartupMessages(library(package = "impute"))
suppressPackageStartupMessages(library(package = "igraph"))
suppressPackageStartupMessages(library(package = "ggnetwork"))
suppressPackageStartupMessages(library(package = "ordinal"))
suppressPackageStartupMessages(library(package = "qvalue"))
suppressPackageStartupMessages(library(package = "ggbeeswarm"))
suppressPackageStartupMessages(library(package = "ggeffects"))
suppressPackageStartupMessages(library(package = "pvca"))
suppressPackageStartupMessages(library(package = "GSA"))
suppressPackageStartupMessages(library(package = "lme4"))
suppressPackageStartupMessages(library(package = "nlme"))
suppressPackageStartupMessages(library(package = "tidyverse"))
suppressPackageStartupMessages(library(package = "gtsummary"))

# set session options
options(show_col_types = FALSE)
VarCorr <- lme4::VarCorr
```

### Load codebase.R, load data matrices and set path to your local output directory

``` r
# load codebase.R, serum-rbd-abtiters data and clinical information
# source("~/bitbucket/impacc/Analysis/data_analysis_template_codebase.R")
source("../../Codebase/codebase_v2.R")
data_env <- load_IMPACC_datasets(DATA_VERSION                = "2022-01-01", 
                                 KEEP_COVID19_POS            = TRUE, 
                                 FILTER_BY_CORE_ASSAY_COHORT = TRUE,
                                 ASSAY_NAME = "serum_autoantibody")
for (n in grep(pattern = "serum_autoantibody|clinical_data", 
               names(data_env), 
               value   = TRUE)) {
  assign(n, value = data_env[[n]])
}


# Setup the assay specific output directory
out_dir = "../output"
```

 

 

#### Combine clinical data and serum auto-antobody data

 

 

#### Get numbers and stats for Table 2

    ##         
    ##           no yes
    ##   Female 177   3
    ##   Male   291  18

    ##          
    ##            no yes
    ##   [19,35]  44   1
    ##   [36,50] 108   2
    ##   [51,65] 172   7
    ##   [66,80] 118   9
    ##   [81,95]  26   2

    ##    
    ##      no yes
    ##   1 101   0
    ##   2 129   2
    ##   3  97   5
    ##   4  94  10
    ##   5  47   4

### Session info

``` r
sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 22.04.2 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## attached base packages:
    ## [1] parallel  grid      stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] doParallel_1.0.16    doRNG_1.8.2          rngtools_1.5        
    ##  [4] missForest_1.4       itertools_0.1-3      iterators_1.0.13    
    ##  [7] foreach_1.5.1        randomForest_4.6-14  ComplexHeatmap_2.6.2
    ## [10] corrr_0.4.3          gridExtra_2.3        rlang_1.1.1         
    ## [13] RColorBrewer_1.1-2   cowplot_1.1.1        ggpubr_0.4.0        
    ## [16] pals_1.7             gtsummary_1.6.0      forcats_0.5.1       
    ## [19] stringr_1.4.0        dplyr_1.0.9          purrr_0.3.4         
    ## [22] readr_2.1.2          tidyr_1.2.0          tibble_3.1.7        
    ## [25] tidyverse_1.3.1      nlme_3.1-148         lme4_1.1-27.1       
    ## [28] Matrix_1.2-18        GSA_1.03.2           pvca_0.1.0          
    ## [31] ggeffects_1.1.3      ggbeeswarm_0.6.0     qvalue_2.22.0       
    ## [34] ordinal_2019.12-10   ggnetwork_0.5.10     ggplot2_3.4.0       
    ## [37] igraph_1.4.2         impute_1.64.0        knitr_1.39          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] minqa_1.2.4         colorspace_2.0-2    rjson_0.2.20       
    ##  [4] ggsignif_0.6.2      ellipsis_0.3.2      rio_0.5.27         
    ##  [7] circlize_0.4.13     GlobalOptions_0.1.2 fs_1.5.2           
    ## [10] dichromat_2.0-0     clue_0.3-59         rstudioapi_0.13    
    ## [13] fansi_0.4.1         lubridate_1.7.10    xml2_1.3.3         
    ## [16] codetools_0.2-16    splines_4.0.2       jsonlite_1.7.2     
    ## [19] nloptr_1.2.2.2      Cairo_1.5-12.2      gt_0.6.0           
    ## [22] broom_0.8.0         cluster_2.1.0       dbplyr_2.1.1       
    ## [25] png_0.1-7           mapproj_1.2.7       compiler_4.0.2     
    ## [28] httr_1.4.4          backports_1.2.0     assertthat_0.2.1   
    ## [31] fastmap_1.1.0       cli_3.6.1           htmltools_0.5.2    
    ## [34] tools_4.0.2         gtable_0.3.0        glue_1.6.2         
    ## [37] reshape2_1.4.4      maps_3.3.0          Rcpp_1.0.8         
    ## [40] carData_3.0-4       cellranger_1.1.0    vctrs_0.6.2        
    ## [43] broom.helpers_1.7.0 xfun_0.31           openxlsx_4.2.4     
    ## [46] rvest_1.0.0         lifecycle_1.0.3     rstatix_0.7.0      
    ## [49] MASS_7.3-51.6       scales_1.2.1        hms_1.1.0          
    ## [52] yaml_2.2.1          curl_4.3            stringi_1.5.3      
    ## [55] ucminf_1.1-4        S4Vectors_0.28.1    BiocGenerics_0.36.1
    ## [58] boot_1.3-25         zip_2.2.0           shape_1.4.6        
    ## [61] pkgconfig_2.0.3     matrixStats_0.59.0  evaluate_0.15      
    ## [64] lattice_0.20-41     tidyselect_1.1.1    plyr_1.8.6         
    ## [67] magrittr_2.0.3      R6_2.5.0            IRanges_2.24.1     
    ## [70] generics_0.1.2      DBI_1.1.1           pillar_1.7.0       
    ## [73] haven_2.4.1         foreign_0.8-80      withr_2.5.0        
    ## [76] abind_1.4-5         modelr_0.1.8        crayon_1.4.1       
    ## [79] car_3.0-11          utf8_1.1.4          tzdb_0.4.0         
    ## [82] rmarkdown_2.9       GetoptLong_1.0.5    readxl_1.3.1       
    ## [85] data.table_1.14.0   reprex_2.0.0        digest_0.6.27      
    ## [88] numDeriv_2016.8-1.1 stats4_4.0.2        munsell_0.5.0      
    ## [91] beeswarm_0.4.0      vipor_0.4.5

<!--
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8    
 [5] LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C             
 [9] LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

attached base packages:
[1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gtsummary_1.5.2       doParallel_1.0.17     doRNG_1.8.2           rngtools_1.5.2       
 [5] missForest_1.4        itertools_0.1-3       iterators_1.0.14      foreach_1.5.2        
 [9] randomForest_4.6-10   ComplexHeatmap_2.11.1 corrr_0.4.3           gridExtra_2.3        
[13] rlang_1.0.1           RColorBrewer_1.1-2    cowplot_1.1.1         ggpubr_0.4.0         
[17] pals_1.7              forcats_0.5.1         stringr_1.4.0         dplyr_1.0.8          
[21] purrr_0.3.4           readr_2.1.2           tidyr_1.2.0           tibble_3.1.6         
[25] tidyverse_1.3.1       nlme_3.1-147          lme4_1.1-28           Matrix_1.2-18        
[29] GSA_1.03.2            pvca_0.1.0            ggeffects_1.1.1       ggbeeswarm_0.6.0     
[33] qvalue_2.22.0         ordinal_2019.12-10    ggnetwork_0.5.10      ggplot2_3.3.5        
[37] igraph_1.3.0          impute_1.64.0         knitr_1.37           

loaded via a namespace (and not attached):
 [1] minqa_1.2.4         colorspace_2.0-3    ggsignif_0.6.3      rjson_0.2.21       
 [5] ellipsis_0.3.2      circlize_0.4.14     GlobalOptions_0.1.2 fs_1.5.2           
 [9] dichromat_2.0-0     clue_0.3-60         rstudioapi_0.13     DT_0.15            
[13] fansi_1.0.2         lubridate_1.8.0     xml2_1.3.2          codetools_0.2-16   
[17] splines_4.0.2       jsonlite_1.7.3      nloptr_2.0.0        gt_0.4.0           
[21] broom_0.7.12        cluster_2.1.0       dbplyr_2.1.1        png_0.1-7          
[25] mapproj_1.2.8       compiler_4.0.2      httr_1.4.2          backports_1.4.1    
[29] assertthat_0.2.1    fastmap_1.1.0       cli_3.2.0           htmltools_0.5.2    
[33] tools_4.0.2         gtable_0.3.0        glue_1.6.1          reshape2_1.4.4     
[37] maps_3.4.0          Rcpp_1.0.8          carData_3.0-5       cellranger_1.1.0   
[41] vctrs_0.3.8         broom.helpers_1.6.0 xfun_0.29           rvest_1.0.2        
[45] lifecycle_1.0.1     rstatix_0.7.0       MASS_7.3-51.6       scales_1.1.1       
[49] hms_1.1.1           yaml_2.3.5          sass_0.4.0          stringi_1.7.6      
[53] ucminf_1.1-4        S4Vectors_0.28.1    checkmate_2.0.0     BiocGenerics_0.36.1
[57] boot_1.3-25         shape_1.4.6         commonmark_1.7      pkgconfig_2.0.3    
[61] matrixStats_0.61.0  evaluate_0.14       lattice_0.20-41     htmlwidgets_1.5.4  
[65] tidyselect_1.1.2    plyr_1.8.6          magrittr_2.0.2      R6_2.5.1           
[69] IRanges_2.24.1      generics_0.1.2      DBI_1.1.2           pillar_1.7.0       
[73] haven_2.4.3         withr_2.4.3         abind_1.4-5         modelr_0.1.8       
[77] crayon_1.5.0        car_3.0-12          utf8_1.2.2          tzdb_0.2.0         
[81] rmarkdown_2.11      GetoptLong_1.0.5    readxl_1.3.1        reprex_2.0.1       
[85] digest_0.6.29       numDeriv_2016.8-1.1 stats4_4.0.2        munsell_0.5.0      
[89] beeswarm_0.4.0      vipor_0.4.5
-->
