
load_blood_cell_markers <- function(setting = "protein", transcript_rowfeatures, sheet = 1){
  require(readxl)
  expression_of_cells <- read_excel("/scratch/metagenomics/age_analysis/databases/41590_2011_BFni2067_MOESM29_ESM.xlsx", sheet = sheet)
  colnames(expression_of_cells) <- expression_of_cells[4,]
  expression_of_cells <- expression_of_cells[-1:-4,]
  if(sheet == 1){
    expression_of_cells[,5:10] <- apply(expression_of_cells[,5:10], 2, as.numeric)
    expression_of_cells$`Highly Expressed in PBMC type` <- gsub("T-cells", "T_cells", gsub("B-cells", "B_cells", expression_of_cells$`Highly Expressed in PBMC type`))
    GSEA_ready_mapping <- expression_of_cells %>%
      separate_rows(`Highly Expressed in PBMC type`, sep = "-")
    GSEA_ready_mapping <- GSEA_ready_mapping[GSEA_ready_mapping$`Highly Expressed in PBMC type` != "",c(11,3)] %>%
      distinct()
  } else if (sheet ==2){
    expression_of_cells[,5:8] <- apply(expression_of_cells[,5:8], 2, as.numeric)
    GSEA_ready_mapping <- expression_of_cells %>% separate_rows(`Highly Expressed in B cells (PBMC) and highly expressed in B cell subsets`, sep = "-")
    GSEA_ready_mapping <- GSEA_ready_mapping[GSEA_ready_mapping$`Highly Expressed in B cells (PBMC) and highly expressed in B cell subsets` != "",c(9,3)] %>%
      distinct()
  }
  
  if(setting == "gene"){
    transcript_rowfeatures <- transcript_rowfeatures %>%
      rownames_to_column() %>%
      as.data.frame()
    GSEA_ready_mapping <- merge(GSEA_ready_mapping, transcript_rowfeatures, by.x = "Gene Symbol", by.y = "gene_name", all.y = T)
    GSEA_ready_mapping <- GSEA_ready_mapping[,2:3]
  }
  colnames(GSEA_ready_mapping) <- c("cluster", "feature")
  return(GSEA_ready_mapping)
}


load_nasal_cell_markers <- function(broad = T, granulocytes = T, setting = "gene", transcript_rowfeatures){
  ziegler_full <- read_excel("/scratch/metagenomics/age_analysis/databases/Nasal_upper_airway_gene_markers_scRNAseq.xlsx", sheet = 1)
  ordovas <- read_excel("/scratch/metagenomics/age_analysis/databases/Nasal_upper_airway_gene_markers_scRNAseq.xlsx", sheet = 2)
  ziegler <- read_excel("/scratch/metagenomics/age_analysis/databases/Nasal_upper_airway_gene_markers_scRNAseq.xlsx", sheet = 3)
  colnames(ordovas) <- gsub("gene_name", "gene", colnames(ordovas))
  colnames(ziegler) <- gsub("gene_name", "gene", colnames(ziegler))
  unique(ziegler$cluster)
  unique(ziegler_full$cluster)
  unique(ordovas$cluster)
  
  if(broad){
    if(granulocytes){
      active_dataset <- as.data.frame(rbind(ziegler[,7:8], ordovas[,7:8]))
    } else {
      active_dataset <- as.data.frame(ziegler[,7:8])
    }
  } else {
    if(granulocytes){
      active_dataset <- as.data.frame(rbind(ziegler_full[,6:7], ordovas[,7:8]))
    } else {
      active_dataset <- as.data.frame(ziegler_full[,6:7])
    }
  }
  
  if(setting == "gene"){
    transcript_rowfeatures <- transcript_rowfeatures %>%
      rownames_to_column() %>%
      as.data.frame()
    active_dataset <- merge(active_dataset, transcript_rowfeatures, by.x = "gene", by.y = "gene_name", all.y = T)
    GSEA_ready_mapping <- active_dataset[,2:3]
  } else {
    GSEA_ready_mapping <- active_dataset
  }
  colnames(GSEA_ready_mapping) <- c("term", "feature")
  return(GSEA_ready_mapping)
}

load_omnipath_ligand <- function(lr_network = T, sig_network = T,
                                 gr_network = F, secreted_only = F,
                                 steps_in_gene_set = 1, report_step=F,
                                 transcript_rowfeatures){
  require(OmnipathR)
  
  lr_Network_Omnipath_Unique <- import_ligrecextra_interactions() %>%
    dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>%
    dplyr::filter(from != to) %>% ## No self signal
    separate_rows(from, sep = "_") %>%
    separate_rows(to, sep = "_") %>%
    mutate(dataframe = "lr_network")
  InterCell_Annotations <- import_omnipath_intercell() ## Annotations for ligands/receptors
  if(secreted_only){
    lr_Network_Omnipath_Unique <- lr_Network_Omnipath_Unique %>%
      filter(from %in% InterCell_Annotations$genesymbol[InterCell_Annotations$secreted])
  }
  
  sig_Network_Omnipath <-     
    import_post_translational_interactions(exclude = "ligrecextra") %>% 
    dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>% 
    dplyr::filter(consensus_direction == "1") %>% 
    dplyr::distinct(from, to, .keep_all = TRUE) %>%
    separate_rows(from, sep = "_") %>%
    separate_rows(to, sep = "_") %>%
    dplyr::distinct() %>%
    mutate(dataframe = "sig_network")
  
  
  gr_Interactions_Omnipath <- 
    import_dorothea_interactions(dorothea_levels = c('A','B','C')) %>%  
    # dplyr::select(source_genesymbol, target_genesymbol) %>% 
    dplyr::rename(from=source_genesymbol, to=target_genesymbol)%>%
    mutate(dataframe = "gr_network")
  
  if(gr_network){
    connection_df <- rbind.data.frame(lr_Network_Omnipath_Unique[,c(3:4, 6, 7, 16)],
                                      sig_Network_Omnipath[,c(3:4, 6, 7, 16)],
                                      gr_Interactions_Omnipath[,c(3:4, 6, 7, 17)])
  } else {
    connection_df <- rbind.data.frame(lr_Network_Omnipath_Unique[,c(3:4, 6, 7, 16)],
                                      sig_Network_Omnipath[,c(3:4, 6, 7, 16)])
  }
  
  ### Run provided number of steps
  gene_set <- data.frame()
  for(j in connection_df$from[connection_df$dataframe == "lr_network"]){
    active_network <- j
    for(i in 1:steps_in_gene_set){
      active_network <- c(active_network, connection_df$to[connection_df$from %in% active_network])
    }
    if(report_step){
      gene_set <- rbind.data.frame(gene_set, data.frame(term = j, feature = active_network[-1], step = i)) 
    } else {
      gene_set <- rbind.data.frame(gene_set, data.frame(term = j, feature = active_network[-1])) 
    }
    
  }
  
  if(!missing(transcript_rowfeatures)){
    gene_set <- merge(gene_set, transcript_rowfeatures %>% rownames_to_column("ensembl_gene_id"),
                      by.x = "feature", by.y = "gene_name", all.x = T)
    gene_set <- gene_set[,c("term", "ensembl_gene_id")]
    colnames(gene_set) <- c("term", "feature")
  }
  
  filter <- c("", "2")
  gene_set <- gene_set[!gene_set$term %in% filter & !gene_set$feature %in% filter,]
  return(gene_set)
  
}


circle_plot <- function(results, title){
  fun_color_range <- colorRampPalette(c("blue", "white", "red"))
  my_colors <- fun_color_range(100)
  coord_polar_free <- coord_polar()
  coord_polar_free$is_free <- function() TRUE
  
  tmp.data <- results %>%
    mutate(ID = gsub("_", " ", ID)) %>%
    mutate(signed_log10_padj = -log10(p.adjust) * (NES/abs(NES))) 
  
  ggplot(tmp.data, aes(x = factor(ID), y= NES, alpha = p.adjust <= 0.05)) +
    geom_bar(stat="identity", width = 0.98, aes(fill = signed_log10_padj), color = "black") +
    coord_polar() +
    theme_bw() +
    scale_fill_distiller(palette = "RdBu", 
                         limits = c( -1*max(abs(tmp.data$signed_log10_padj), na.rm = T),
                                     max(abs(tmp.data$signed_log10_padj),na.rm = T) )
    ) + 
    #facet_wrap() + xlab("") +
    scale_alpha_manual(values = c("FALSE" = 0.1,"TRUE"= 1)) +
    theme(axis.text.x = element_text(size=10))
}




omnipath_connection_pull <- function(features, gr_network = F, double_step = F){
  lr_Network_Omnipath_Unique <- import_ligrecextra_interactions() %>%
    dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>%
    dplyr::filter(from != to) %>% ## No self signal
    separate_rows(from, sep = "_") %>%
    separate_rows(to, sep = "_") %>%
    mutate(dataframe = "lr_network")
  
  #InterCell_Annotations <- import_omnipath_intercell() ## Annotations for ligands/receptors
  
  sig_Network_Omnipath <-     
    import_post_translational_interactions(exclude = "ligrecextra") %>% 
    dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>% 
    dplyr::filter(consensus_direction == "1") %>% 
    dplyr::distinct(from, to, .keep_all = TRUE) %>%
    separate_rows(from, sep = "_") %>%
    separate_rows(to, sep = "_") %>%
    dplyr::distinct() %>%
    mutate(dataframe = "sig_network")
  
  
  gr_Interactions_Omnipath <- 
    import_dorothea_interactions(dorothea_levels = c('A','B','C')) %>%  
    # dplyr::select(source_genesymbol, target_genesymbol) %>% 
    dplyr::rename(from=source_genesymbol, to=target_genesymbol)%>%
    mutate(dataframe = "gr_network")
  
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  genes_of_interest <- lr_Network_Omnipath_Unique[lr_Network_Omnipath_Unique$from %in% features, ]
  genes_of_interest$order <- "first-order"
  sig_network_selected <- sig_Network_Omnipath[sig_Network_Omnipath$from %in% c(genes_of_interest$to, features), ]
  sig_network_selected$order <- "second-order"
  if(double_step){
    sig_network_selected2 <- sig_Network_Omnipath[sig_Network_Omnipath$from %in% c(sig_network_selected$to), ]
    sig_network_selected2$order <- "third-order"
    sig_network_selected <- rbind.data.frame(sig_network_selected, sig_network_selected2)
  }
  gr_network_selected <- gr_Interactions_Omnipath[gr_Interactions_Omnipath$from %in% c(sig_network_selected$to, genes_of_interest$to, features),]
  
  if(gr_network){
    final_selected <- rbind.data.frame(genes_of_interest[,c(3:4, 6, 7, 16)],
                                       sig_network_selected[,c(3:4, 6, 7, 16)],
                                       gr_network_selected[,c(3:4, 6, 7, 17)])
  } else {
    final_selected <- rbind.data.frame(genes_of_interest[,c(3:4, 6, 7, 16)],
                                       sig_network_selected[,c(3:4, 6, 7, 16)])
  }
  
  final_selected <- final_selected %>%
    distinct(from, to, .keep_all = TRUE)
  
  final_selected <- final_selected[!grepl("_", final_selected$from),]
  final_selected <- final_selected[!grepl("_", final_selected$to),]
  return(final_selected)
}




