##function does not run successfully from codebase


## Goal: translate row/column names of a data frame to readable names from row files
## used for pbmc/nt/metabolites
name_translator = function(data_env, DF, 
                                 omic_short_name="pbmc", 
                                 omic_name=NULL,
                                 omic_name_ref = NULL,
                                 mod = "row"){
  if(is.null(omic_name) & is.null(omic_name_ref)){
    stop("no omic name or omic_name_ref provided.")
  }else{
    if(is.null(omic_name)){
      omic_name = omic_name_ref[omic_short_name]
    }
  }
  row_feature_file = paste0(omic_name,"_rowfeature")
  rowfeature_DF = data_env0[[row_feature_file]]
  readable_columns = c("pbmc_transcriptomics" = "gene_name",
                       "nasal_transcriptomics" = "gene_name",
                       "plasma_metabolomics_global" = "CHEMICAL_NAME")
  #DF = impacc_analysis_record$projection_evaluation_train$projection_coef[[4]][,1,drop = F]
  if(omic_name%in%c("pbmc_transcriptomics","nasal_transcriptomics", "plasma_metabolomics_global")){
    if(mod == "row"){
      rownames(DF) = rowfeature_DF[rownames(DF),readable_columns[omic_name]]
    }else{
      colnames(DF) = rowfeature_DF[colnames(DF),readable_columns[omic_name]]
    }
    
  }else{
    stop("unsupported omic_name.")
  }
  return(DF)
}



#migrated/adjusted from kegg_pipeline/05162022/annotation_pipeline.R/convert_mciaObj_to_annotation_formatV2
## Goal: organize by mod
convert_mciaObj_to_annotation_format <- function(projection_coef, omic_names){
  
  remap_projection_coef_master <- list()
  if("data.frame" %in% class(projection_coef[[1]])){
    for(d in 1:length(projection_coef)){
      projection_coef[[d]] = data.frame(projection_coef[[d]])
      colnames(projection_coef[[d]]) = as.character(1:ncol(projection_coef[[d]]))
    }
  }
  for(index in 1: ncol(projection_coef[[1]])){
    mod_name = paste0("mod", index)
    remap_mcia_object <- list()
    for(d in 1:length(projection_coef)){
      remap_mcia_object[[d]] = projection_coef[[d]][,index,drop = F]
    }
    names(remap_mcia_object) = omic_names
    remap_projection_coef_master[[mod_name]] <- remap_mcia_object
  }
  return( remap_projection_coef_master)
}



#MGHtest
mgh_test_wrapper = function(geneList,
                            universe,
                            TERM2GENE,
                            minGSSize = 2,
                            maxGSSize = 500,
                            pvalueCutoff = 2,
                            p.adjust.methods="BH"){
  TERM2GENE=TERM2GENE[TERM2GENE[,2]%in%universe,,drop = F]
  genes = names(geneList)
  genes = intersect(genes, universe)
  ##keep only terms satisfying the size restriction
  size_counts = table(TERM2GENE[,1])
  terms_keep = names(size_counts[size_counts>=minGSSize & size_counts <=maxGSSize])
  TERM2GENE = TERM2GENE[TERM2GENE[,1]%in%terms_keep,,drop = F]
  pathway_names = unique(TERM2GENE[,1])
  genes_all = c(genes, setdiff(universe,genes))
  result_table = data.frame(ID = pathway_names)
  result_table[["Description"]] = pathway_names
  result_table[["setSize"]] = NA
  result_table[["pvalue"]] = NA
  result_table[["p.adj"]] = NA
  result_table[["rank"]] = NA
  result_table[["coreSetSiz"]] = NA
  result_table[["core_enrichment"]] = NA
  m = length(genes)
  for(pathway_id in 1:length(pathway_names)){
    pathway_name0 = pathway_names[pathway_id]
    term_genes = TERM2GENE[TERM2GENE[,1] ==pathway_name0,2]
    result_table[["setSize"]][pathway_id] = length(term_genes)
    mgh_vec = rep(0,length(genes_all))
    names(mgh_vec) = genes_all
    mgh_vec[term_genes]=1
    mHGtest_res = mHG::mHG.test(lambdas = mgh_vec, n_max = m)
    result_table[["pvalue"]][pathway_id] = mHGtest_res$p.value
    result_table[["rank"]][pathway_id] = mHGtest_res$n
    result_table[["coreSetSiz"]][pathway_id] = mHGtest_res$b
    result_table[["core_enrichment"]][pathway_id] = paste0(names(mgh_vec)[mgh_vec==1][1:mHGtest_res$b],collapse = "#")
  }
  result_table=result_table[order(result_table[["pvalue"]]),]
  result_table[["p.adj"]] = stats::p.adjust(p =result_table[["pvalue"]] , method = p.adjust.methods)
  ##turn to GSEA class object?
  return(result_table)
}



#Composite gene creation
comoposite_features_creation_unit = function(omic_assay_DF, 
                                        projection_coef_DF,
                                        leading_edges,
                                        sep = "/",
                                        is_vector = FALSE){
  #leading_edges can be either "/" (sep) seperated or as a vector of features
  if(!is_vector){
    leading_edges = strsplit(leading_edges, sep)[[1]]
    leading_edges = intersect(rownames(projection_coef_DF),leading_edges)
  }
  if(length(leading_edges)>0){
    w =  abs(projection_coef_DF[leading_edges,1])
    x = as.matrix(omic_assay_DF[,leading_edges])
    z = x %*% matrix(w, nrow = length(w), ncol = 1)
    z = as.data.frame(z) %>% `rownames<-`(rownames(omic_assay_DF))
  }else{
    z = NULL
  }
  
  return(z)
}


automatic_kegg_hallmark_sub_composite_run = function(test_type,base_obj,  pval=0.05){
  assay_short_names = names(base_obj$omic_names)
  if(test_type=="kegg_pos"){
    tmp_obj = base_obj$MODmghRes[["kegg_pos_impacc"]]
  }else if(test_type=="kegg_neg"){
    tmp_obj = base_obj$MODmghRes[["kegg_neg_impacc"]]
  }else if(test_type == "hallmark_pos"){
    tmp_obj = base_obj$MODmghRes[["hallmark_pos_impacc"]]
  }else if(test_type == "hallmark_neg"){
    tmp_obj = base_obj$MODmghRes[["hallmark_neg_impacc"]]
  }else if(test_type == "sub_pos"){
    tmp_obj = base_obj$MODmghRes[["subpathway_pos_impacc"]]
  }else if(test_type == "sub_neg"){
    tmp_obj = base_obj$MODmghRes[["subpathway_neg_impacc"]]
  }else{
    stop("do not support test type")
  }
  collected_terms = c()
  if(length(tmp_obj)>0){
    for(j in 1:length(tmp_obj)){
      idx = which(tmp_obj[[j]]$pvalue<=pval)
      collected_terms=c(collected_terms, tmp_obj[[j]]$ID[idx])
    }
  }

  collected_terms = unique(collected_terms)
  term_members = list()
  term_scores = list()
  term_scores[["train"]] = list()
  term_scores[["test"]] = list()
  if(length(collected_terms)>0){
    for(l in 1:length(collected_terms)){
      ava_omic=c()
      term0 = collected_terms[l]
      for(omic_name0 in names(tmp_obj)){
        tmp1 = tmp_obj[[omic_name0]]
        idx = which(tmp1$ID==term0)
        if(length(idx)>0){
          if(tmp1[idx,"pvalue"]<0.25){
            tmp2 = c(tmp1[idx,"core_enrichment"])
            names(tmp2) = term0
            term_members[[omic_name0]]=c(term_members[[omic_name0]],tmp2)
            ava_omic=c( ava_omic, omic_name0)
          }
        }
      }
      tmp3 = base_obj$composite_features_creation(
        result_table_list =tmp_obj, term = term0, core_enrichment_name = "core_enrichment",phase = "train", sep = "#")
      tmp4 = base_obj$composite_features_creation(result_table_list =tmp_obj,
                                                  term = term0,core_enrichment_name = "core_enrichment",phase = "test", sep = "#")
      for(omic_name0 in ava_omic){
        tmp2 =  tmp3[[omic_name0]]
        colnames(tmp2) = term0
        tmp22 =  tmp4[[omic_name0]]
        colnames(tmp22) = term0
        if(is.null(term_scores[["train"]][[omic_name0]])){
          term_scores[["train"]][[omic_name0]]= tmp2
          term_scores[["test"]][[omic_name0]] = tmp22
        }else{
          term_scores[["train"]][[omic_name0]]=cbind(term_scores[["train"]][[omic_name0]], tmp2)
          term_scores[["test"]][[omic_name0]]=cbind(term_scores[["test"]][[omic_name0]], tmp22)
        }
      }
    }
  }
 
  return(list(term_members=term_members, term_scores=term_scores))
  
}



plot_enrichment_dot = function(self, select_pathways = NULL, 
                               plot_title="mod", 
                               thr=0.1, 
                               datasets=c("ppt","ppgd","so","pmg","nt","pbmc"), 
                               test_type = "mhg", 
                               databases = c("kegg","hallmark"), 
                               top_nPathw = NULL,
                               auto_ordering = TRUE
                               ){
  if(is.null(select_pathways) & is.null(top_nPathw)) {
    stop("Specify either select_pathways or top_nPathw")
  }
  if(!is.null(select_pathways) & !is.null(top_nPathw)) {
    warning("Both select_pathways and top_nPathw are supplied; top_nPathw is ignored")
  }
  # Make gsea and mhg result formats identical
  if(test_type == "gsea") {
    test_res = self$MODgseaRes
  }
  if(test_type == "mhg") {
    test_res = self$MODmghRes
  }
  
  # Gather top_nPathw pathways for plotting, if requested.
  if(! is.null(top_nPathw)) {
    top_pathways = c()
    for(db in databases){
      db = paste0(db,"_mod")
      for(d in datasets) {
        tab = test_res[[db]][[d]]$result
        if(!is.null(tab)) {
          tab = tab[ tab$qvalues < thr, ]
          tmp_top_pathways = slice_min(tab, pvalue, n = top_nPathw, with_ties = FALSE) %>% pull(Description)
          #tmp_top_pathways = paste0(str_extract(db, "^."), "_", tmp_top_pathways)
          top_pathways = c(top_pathways, tmp_top_pathways)
        }
      }
    }
    pathways = unique(top_pathways)
  }
  
  if(!is.null(select_pathways)) {
    pathways = select_pathways
  }
  
  # Gather enrichment results
  res_df = data.frame()
  for(db in databases){
    db = paste0(db,"_mod")
    for(d in datasets) {
      tab = test_res[[db]][[d]]$result
      pathw_idx = which(tab$Description%in%pathways)
      if( length(pathw_idx) > 0 )
        res_df = rbind(res_df, 
                       data.frame(tab[pathw_idx,
                                      c("Description","NES","qvalues")], 
                                  db=db,
                                  d=d) %>%
                         setNames(c("pathway","score","sig.val","database","dataset")) 
        )
    }
  }
  
  # Calculate joint p-value
  calc_join_qval <- function(z){
    driving_assays = rep(1,length(z)) ## Currently, the order of datasets or which datasets are driving doesn't matter, so using 1 for all assays for which the pvalue exist. But if specific datasets need to be driving, make sure to match the order of dataset in z with the driving_assays vector.
    ii = !is.na(z)
    t = qnorm(z[ii])
    #t = sum(t*(1:length(t)))/sqrt(sum((1:length(t))^2))
    t = sum(t*driving_assays[ii])/sqrt(sum(driving_assays[ii]^2))
    #t = sum(t*rep(1, length(datasets)))/sqrt(sum(rep(1, length(datasets))^2))
    pnorm(t, lower.tail = T)
    #t = -2*sum(log(z[ii][driving_assays[ii]==1]))
    #pchisq(t, df=2*length(ii),lower.tail = F)
  }
  joint_sig_vals = res_df %>% group_by(pathway, database) %>% 
    #summarise(sig.val = calc_join_qval(sig.val)) %>%
    summarise(sig.val = survcomp::combine.test(sig.val)) %>%
    mutate(dataset = "joint", score=NA) %>%
    relocate(score, .after=pathway) %>%
    relocate(sig.val, .after=score)
  res_df = rbind(res_df, joint_sig_vals)
  

  # Prepare df for plotting
  res_df$log10.sig.val = -log10(res_df$sig.val)
  #res_df$pathway = paste0( str_extract(res_df$database, "^."), res_df$pathway)
  # Use the order of pathways based on the joint sig.val
  res_df$pathway = factor(
    res_df$pathway, 
    levels= rev(res_df %>% 
                  dplyr::filter(dataset=="joint") %>% 
                  arrange(desc(database), sig.val) %>% pull(pathway) 
    )
  )
  if(auto_ordering == FALSE & !is.null(select_pathways)) # If pathways are selected for plotting and auto_ordering is off, use the user provided order of pathways
    res_df$pathway = factor(res_df$pathway, levels=rev(pathways))
  res_df$dataset = factor(res_df$dataset, levels=c(datasets, "joint") )
  
  # Trim the pathway names to 50 characters
  levels(res_df$pathway) = substr(levels(res_df$pathway), 1, 50)
  
  # Make a plot
  p = ggplot(data =res_df, aes(x = dataset, y = pathway, size=log10.sig.val ,col =score)) +
    geom_point(color="white") + 
    geom_point(shape=1) + 
    geom_point(data=. %>% dplyr::filter(log10.sig.val > -log10(thr) ) ) + 
    ggtitle(plot_title) + 
    scale_color_distiller(palette = "RdBu", 
                          limits = c( -1*max(abs(res_df$score), na.rm = T),
                                      max(abs(res_df$score),na.rm = T) )
    ) + 
    theme_classic() +
    theme(panel.grid.major = element_line(color="grey", size=0.1), 
          axis.text.x = element_text(angle=90, hjust=1)) 
  #print(p)
  return(p)
}



plot_enrichment_dot_v2 = function(self, select_pathways = NULL, 
                               plot_title="mod", 
                               thr=0.1, 
                               datasets=c("ppt","ppgd","so","pmg","nt","pbmc"), 
                               test_type = "mhg", 
                               databases = c("kegg","hallmark"), 
                               #top_nPathw = NULL,
                               auto_ordering = TRUE,
                               nested_ordering_with_databases = TRUE, # If set to TRUE and auto_ordering is TRUE, the rows are ordered by databases first and then by joint pval. If FALSE and auto_ordering is TRUE, the rows are ordered by joint pval.
                               use_cauchy_jointp = TRUE,
                               pathway_groups = NULL # This option takes priority over nested_ordering_with_databases. Named vector where names represent pathways in select_pathways. If auto_ordering is TRUE, select_pathways and pathway_groups are provided, the ordering is performed first by the groups and then by joint pval.
){
  #if(is.null(select_pathways) & is.null(top_nPathw)) {
  #  stop("Specify either select_pathways or top_nPathw")
  #}
  #if(!is.null(select_pathways) & !is.null(top_nPathw)) {
  #  warning("Both select_pathways and top_nPathw are supplied; top_nPathw is ignored")
  #}
  if(use_cauchy_jointp == TRUE & test_type != "mhg") {
    stop("use_cauchy_jointp works only for mgh test")
  }
  # Make gsea and mhg result formats identical
  #if(test_type == "gsea") {
  #  test_res = self$MODgseaRes
  #}
  if(test_type == "mhg") {
    test_res = self$MODmghRes
  }
  
  # Gather all enrichment results: for all databases and datasets, and also for the joint pvalues. Get both the pvalues and adjusted pvalues.
  all_res = data.frame()
  for(db in databases) {
    if(db %in% c("metaSub")) {
      tmp_df = test_res$metabSub_mod$pmg$result %>% 
        setNames(c("pathway","p","setSize","rank","p.adj","core_enrichment","score","qvalues")) %>%
        dplyr::select("pathway","score","p","p.adj","core_enrichment") %>% 
        mutate(database = "metaSub", 
               dataset = "pmg", 
               log10.p = -log10(p), 
               log10.p.adj = -log10(p.adj)
               )
      tmp_joint_df = tmp_df %>% mutate(score=NA, core_enrichment = NA, dataset = "joint")
      all_res = rbind(all_res, 
                      tmp_df,
                      tmp_joint_df
      )
    } else {
      pos_list = test_res[[paste0(db, "_pos_impacc")]]
      neg_list = test_res[[paste0(db, "_neg_impacc")]]
      agg_mHG = enrichment_aggregation(pos_list=pos_list, neg_list = neg_list)
      for(d in names(agg_mHG$max_list)) {
        tmp_df = agg_mHG$max_list[[d]] %>% dplyr::select(ID,sign,pvalue,p.adj, core_enrichment) %>% 
          setNames(c("pathway","score","p","p.adj","core_enrichment")) %>%
          mutate(score = ifelse(score == "+", 1, -1),
                 database = db,
                 dataset = d,
                 log10.p = -log10(p),
                 log10.p.adj = -log10(p.adj)
                 )
        all_res = rbind(all_res, 
                        tmp_df
                        )
      }
      tmp_joint_df = agg_mHG$pval_mat %>% 
        rownames_to_column("pathway") %>% 
        dplyr::select(pathway, joint, joint_adj) %>% setNames(c("pathway","p","p.adj")) %>%
        mutate(score = NA, 
               database = db, 
               dataset = "joint", 
               core_enrichment = NA,
               log10.p = -log10(p),
               log10.p.adj = -log10(p.adj)
               )
      all_res = rbind(all_res, 
                      tmp_joint_df
      )
    }
  }
  
  # Select significant pathways based on the adjusted joint p value.
  pathways = all_res %>% dplyr::filter(dataset == "joint" & p.adj < thr) %>% pull(pathway) %>% unique()
  
  # Use the user provided pathways
  if(!is.null(select_pathways)) {
    pathways = select_pathways
  }
  
  # Subset the results for "pathways"
  all_res = all_res %>% dplyr::filter(pathway %in% pathways)
  
  all_res$pathway = factor(
    all_res$pathway, 
    levels= rev(all_res %>% 
                  dplyr::filter(dataset=="joint") %>% 
                  arrange(p.adj) %>% pull(pathway) 
    )
  )
  
  if(nested_ordering_with_databases) {
    all_res$pathway = factor(
      all_res$pathway, 
      levels= rev(all_res %>% 
                    dplyr::filter(dataset=="joint") %>% 
                    arrange(desc(database), p.adj) %>% pull(pathway) 
      )
    )
  }
  if(!is.null(pathway_groups) & !is.null(select_pathways)) {
    all_res$pathway = factor(
      all_res$pathway, 
      levels= rev(all_res %>% 
                    mutate(my_groups = pathway_groups[as.character(pathway)]) %>%
                    dplyr::filter(dataset=="joint") %>% 
                    arrange(my_groups, p.adj) %>% pull(pathway) 
      )
    )
  }
  if(auto_ordering == FALSE & !is.null(select_pathways)) # If pathways are selected for plotting and auto_ordering is off, use the user provided order of pathways
    all_res$pathway = factor(all_res$pathway, levels=rev(pathways))
  
  all_res$dataset = factor(all_res$dataset, levels=c(datasets, "joint") )
  
  # Trim the pathway names to 50 characters
  levels(all_res$pathway) = substr(levels(all_res$pathway), 1, 50)
  

  # Make a plot
  p = ggplot(data =all_res, aes(x = dataset, y = pathway, size=log10.p.adj, col=score)) +
    geom_point(color="white") + 
    geom_point(shape=1) + 
    geom_point(data=. %>% dplyr::filter(log10.p.adj > -log10(thr) ) ) + 
    ggtitle(plot_title) + 
    scale_color_distiller(palette = "RdBu", 
                          limits = c( -1*max(abs(all_res$score), na.rm = T),
                                      max(abs(all_res$score),na.rm = T) )
    ) + 
    theme_classic() +
    theme(panel.grid.major = element_line(color="grey", size=0.1), 
          axis.text.x = element_text(angle=90, hjust=1))
  #print(p)
  return(p)
}



plot_enrichment_dot_lg = function(self, select_pathways = NULL, 
                               plot_title="mod", 
                               thr=0.1, 
                               datasets=c("ppt","ppgd","so","pmg","nt","pbmc"), 
                               test_type = "mhg", 
                               databases = c("kegg","hallmark"), 
                               top_nPathw = NULL,
                               auto_ordering = TRUE
){
  if(is.null(select_pathways) & is.null(top_nPathw)) {
    stop("Specify either select_pathways or top_nPathw")
  }
  if(!is.null(select_pathways) & !is.null(top_nPathw)) {
    warning("Both select_pathways and top_nPathw are supplied; top_nPathw is ignored")
  }
  # Make gsea and mhg result formats identical
  if(test_type == "gsea") {
    test_res = self$MODgseaRes
  }
  if(test_type == "mhg") {
    test_res = self$MODmghRes
  }
  
  # Gather top_nPathw pathways for plotting, if requested.
  if(! is.null(top_nPathw)) {
    top_pathways = c()
    for(db in databases){
      db = paste0(db,"_mod")
      for(d in datasets) {
        tab = test_res[[db]][[d]]$result
        if(!is.null(tab)) {
          tab = tab[ tab$qvalues < thr, ]
          tmp_top_pathways = slice_min(tab, pvalue, n = top_nPathw, with_ties = FALSE) %>% pull(Description)
          #tmp_top_pathways = paste0(str_extract(db, "^."), "_", tmp_top_pathways)
          top_pathways = c(top_pathways, tmp_top_pathways)
        }
      }
    }
    pathways = unique(top_pathways)
  }
  
  if(!is.null(select_pathways)) {
    pathways = select_pathways
  }
  
  # Gather enrichment results
  res_df = data.frame()
  for(db in databases){
    db = paste0(db,"_mod")
    for(d in datasets) {
      tab = test_res[[db]][[d]]$result
      pathw_idx = which(tab$Description%in%pathways)
      if( length(pathw_idx) > 0 )
        res_df = rbind(res_df, 
                       data.frame(tab[pathw_idx,
                                      c("Description","NES","qvalues")], 
                                  db=db,
                                  d=d) %>%
                         setNames(c("pathway","score","sig.val","database","dataset")) 
        )
    }
  }
  
  # Calculate joint p-value
  calc_join_qval <- function(z){
    driving_assays = rep(1,length(z)) ## Currently, the order of datasets or which datasets are driving doesn't matter, so using 1 for all assays for which the pvalue exist. But if specific datasets need to be driving, make sure to match the order of dataset in z with the driving_assays vector.
    ii = !is.na(z)
    t = qnorm(z[ii])
    #t = sum(t*(1:length(t)))/sqrt(sum((1:length(t))^2))
    t = sum(t*driving_assays[ii])/sqrt(sum(driving_assays[ii]^2))
    #t = sum(t*rep(1, length(datasets)))/sqrt(sum(rep(1, length(datasets))^2))
    pnorm(t, lower.tail = T)
    #t = -2*sum(log(z[ii][driving_assays[ii]==1]))
    #pchisq(t, df=2*length(ii),lower.tail = F)
  }
  joint_sig_vals = res_df %>% group_by(pathway, database) %>% 
    #summarise(sig.val = calc_join_qval(sig.val)) %>%
    summarise(sig.val = survcomp::combine.test(sig.val)) %>%
    mutate(dataset = "joint", score=NA) %>%
    relocate(score, .after=pathway) %>%
    relocate(sig.val, .after=score)
  res_df = rbind(res_df, joint_sig_vals)
  
  
  # Prepare df for plotting
  res_df$log10.sig.val = -log10(res_df$sig.val)
  #res_df$pathway = paste0( str_extract(res_df$database, "^."), res_df$pathway)
  # Use the order of pathways based on the joint sig.val
  res_df$pathway = factor(
    res_df$pathway, 
    levels= rev(res_df %>% 
                  dplyr::filter(dataset=="joint") %>% 
                  arrange(desc(database), sig.val) %>% pull(pathway) 
    )
  )
  if(auto_ordering == FALSE & !is.null(select_pathways)) # If pathways are selected for plotting and auto_ordering is off, use the user provided order of pathways
    res_df$pathway = factor(res_df$pathway, levels=rev(pathways))
  res_df$dataset = factor(res_df$dataset, levels=c(datasets, "joint") )
  
  # Trim the pathway names to 50 characters
  levels(res_df$pathway) = substr(levels(res_df$pathway), 1, 50)
  
  # Make a plot
  p = ggplot(data =res_df, aes(x = dataset, y = pathway, size=log10.sig.val ,col =score)) +
    geom_point(color="white") + 
    geom_point(shape=1) + 
    geom_point(data=. %>% dplyr::filter(log10.sig.val > -log10(thr) ) ) + 
    ggtitle(plot_title) + 
    scale_color_distiller(palette = "RdBu", 
                          limits = c( -1*max(abs(res_df$score), na.rm = T),
                                      max(abs(res_df$score),na.rm = T) )
    ) + 
    theme_classic() +
    theme(panel.grid.major = element_line(color="grey", size=0.1), 
          axis.text.x = element_text(angle=90, hjust=1)) 
  #print(p)
  return(p)
}


# Modify the GSEA and MHG results to same format.
# For MHG, the pos and neg tables are combined and the duplicate pathways between them are removed by selecting the pathway with smallest p.adj.
streamline_enrichment_results <- function(self, 
                                       datasets=c("ppt","ppgd","so","pmg","nt","pbmc"), 
                                       databases = c("kegg","hallmark")
                                       ){
  # Make gsea and mhg result formats identical
  # GSEA
  for(db in databases) {
    self$MODgseaRes[[paste0(db,"_mod")]] = list()
    for(d in datasets) {
      if(!is.null(self$MODgseaRes[[db]][[d]]))
        self$MODgseaRes[[paste0(db,"_mod")]][[d]] = 
          list(result = self$MODgseaRes[[db]][[d]]@result)
    }
  }

  # MHG
  for(db in databases) {
    self$MODmghRes[[paste0(db,"_mod")]] = list()
    for(d in datasets) {
      if(!is.null(self$MODmghRes[[paste0(db,"_pos_impacc")]][[d]])) {
        pos_res = self$MODmghRes[[paste0(db,"_pos_impacc")]][[d]] %>% 
          dplyr::select("Description","pvalue","setSize","rank","p.adj","core_enrichment") %>% 
          mutate(NES = 1, qvalues = p.adj)
        neg_res = self$MODmghRes[[paste0(db,"_neg_impacc")]][[d]] %>% 
          dplyr::select("Description","pvalue","setSize","rank","p.adj","core_enrichment") %>% 
          mutate(NES = -1, qvalues = p.adj)
        res = rbind(pos_res, neg_res) %>% 
          group_by(Description) %>% 
          slice_min(pvalue, n=1, with_ties = F) %>%
          #dplyr::filter( pvalue == min(pvalue) ) %>%
          ungroup()
        self$MODmghRes[[paste0(db,"_mod")]][[d]] = list(result = res)
      }
    }
  }
}

# Combine the (modified) results from all databases
combine_enrichment_results = function(obj, databases = c("kegg","hallmark","smpdb","immglob"), datasets = c("ppt","ppgd","so","pmg","nt","pbmc"), enrichMethod = "MODmghRes") {
  obj[[enrichMethod]]$all_mod = list()
  for(d in datasets) {
    obj[[enrichMethod]]$all_mod[[d]] = list(result = data.frame())
    for(db in databases) {
      if(! is.null(obj[[enrichMethod]][[paste0(db,"_mod")]][[d]]$result) ) {
        obj[[enrichMethod]]$all_mod[[d]]$result = rbind(
          obj[[enrichMethod]]$all_mod[[d]]$result,
          obj[[enrichMethod]][[paste0(db,"_mod")]][[d]]$result %>% mutate(database = db)
        )
      }
    }
  }
}



## include_all_features: include all features of pathway that are present in the data instead of using only the leading edge features. Only kegg and hallmark databases are available; use projection_pval_cutoff to set a pval cutoff for selecting only significantly associated features.
## Don't use top_nfeat_by_projection with include_all_features.
## If a number is specified for top_nfeat_by_projection, only the top top_nfeat_by_projection leading edge features are retained.
## version2 == TRUE: use new dataset names and do not scale the triangles.
term_feature_network <- function(myobj, padj_thr = 0.1, 
                                 REDUCED_GRAPH = TRUE, 
                                 min_path_count = 1, 
                                 all_enrich_res_slot = "MODmghRes", 
                                 pathways, title="", 
                                 return_graph=F, 
                                 include_all_features = FALSE, 
                                 projection_fwer_cutoff = 0.001,
                                 top_nfeat_by_projection = NULL,
                                 useDejavuFonts = FALSE,
                                 excludeText = FALSE,
                                 version2 = FALSE,
                                 select_features = NULL) {
  if(is.numeric(top_nfeat_by_projection) & is.character(select_features)) {
    warning("WARNING:Are you sure you want to use both top_nfeat_by_projection and select_features parameters?")
  }
  enrich_res_list = myobj[[all_enrich_res_slot]]$all_mod
  
  edges = data.frame()
  pathway_scores = data.frame()
  for(d in names(enrich_res_list)) {
    res_sub = enrich_res_list[[d]]$result %>% dplyr::filter(Description %in% pathways & qvalues < padj_thr)
    if(nrow(res_sub) > 0) {
      for(i in 1:nrow(res_sub)) {
        path = res_sub[i,]$Description
        feat = str_split(res_sub[i,]$core_enrichment, "#")[[1]]
        is_leading_edge = rep(TRUE, length(feat))
        if(is.numeric(top_nfeat_by_projection)) {
          tmp_nfeat = min(length(feat), top_nfeat_by_projection)
          tmp_feat_coeff = myobj$projection_coef[[d]][feat,] %>% setNames(feat)
          feat = names(tail(sort(abs(tmp_feat_coeff)) , tmp_nfeat))
          is_leading_edge = rep(TRUE, length(feat))
        }
        if(is.character(select_features)) {
          feat = feat[ feat %in% select_features ]
          if(length(feat) == 0)
            next
          is_leading_edge = rep(TRUE, length(feat))
        }
        if(include_all_features) {
          feat = myobj$MODgseaRes$kegg[[d]]@geneSets[[path]]
          if(!is.null(myobj$MODgseaRes$hallmark[[d]])) {
            feat = c(feat, myobj$MODgseaRes$hallmark[[d]]@geneSets[[path]])
          }
          feat = unique(feat)
          idx = which(rownames(myobj$projection_coef[[d]]) %in% feat)
          names = rownames(myobj$projection_coef[[d]])[idx]
          FWER = myobj$projection_pvalue[[d]][idx,] * nrow(myobj$projection_pvalue[[d]])
          feat = names[ 
            FWER < projection_fwer_cutoff
          ]
          is_leading_edge = feat %in% str_split(res_sub[i,]$core_enrichment, "#")[[1]]
        }
        feat_coeff = myobj$projection_coef[[d]][feat,]
        feat_coeff_dir = sign(myobj$projection_coef[[d]][feat,])
        feat_pval = myobj$projection_pvalue[[d]][ match(feat, rownames(myobj$projection_coef[[d]])) ,]
        FWER = feat_pval * nrow(myobj$projection_pvalue[[d]])
        feat_pval_is_sig = FWER < projection_fwer_cutoff
        path_idx = which(path == pathways)
        mod_d = paste0(path_idx, "_", d)
        score = res_sub[i,]$NES
        qval = res_sub[i,]$qvalues
        pathway_scores = rbind(
          pathway_scores,
          data.frame(name = mod_d, score = score, score_dir = sign(score), qval = qval, log10qval = -log10(qval))
        )
        edges = rbind(
          edges,
          data.frame(x = path, y = mod_d, feat_coeff = 0, feat_coeff_dir = 0, feat_pval = 1, is_leading_edge = F, feat_pval_is_sig = F),
          #data.frame(x = mod_d, y = paste0(d, "_", feat))
          data.frame(x = mod_d, y = feat, feat_coeff = feat_coeff, feat_coeff_dir = feat_coeff_dir, feat_pval = feat_pval, is_leading_edge = is_leading_edge, feat_pval_is_sig = feat_pval_is_sig)
        )
      }
    }
  }
  
  
  # count the number of times a feature appears (counting all occurances even if the feature has same names in multiple assays).
  feature_count = data.frame(ids = c(edges$x, edges$y)) %>% mutate(ids = gsub("^\\d+_","",ids)) %>% dplyr::filter(! (ids %in% pathways | ids %in% omic_short_names) ) %>% pull(ids) %>% table() %>% as.data.frame() %>% setNames(c("feature","n")) %>% column_to_rownames("feature")
  
  if(REDUCED_GRAPH) {
    edges = edges %>% dplyr::filter(! x %in% pathways)
    edges$x = gsub("^\\d+_","",edges$x)
    
    # Keep only those features that were leading edge in min_path_count pathways
    edges = group_by(edges, x, y) %>% mutate(redun = n()) %>% dplyr::filter(redun >= min_path_count)
    
  }
  
  nodes = data.frame(name = unique(c(edges$x, edges$y))) %>%
    mutate(id = 1:nrow(.)) %>%
    dplyr::select(id, everything()) %>%
    mutate(score_dir = pathway_scores$score_dir[ match(name, pathway_scores$name) ],
           log10qval = pathway_scores$log10qval[ match(name, pathway_scores$name) ])
  
  edges$ori.x = edges$x
  edges$ori.y = edges$y
  edges$x = nodes$id[ match(edges$x, nodes$name) ]
  edges$y = nodes$id[ match(edges$y, nodes$name) ]
  
  
  
  g = graph_from_data_frame(edges, directed = TRUE, nodes)
  
  V(g)$degree <- degree(g, mode = "in")
  
  set.seed(1234)
  lay = create_layout(g, layout = "fr")
  #if(REDUCED_GRAPH) {
  lay$name2 = gsub("^\\d+_","",lay$name)
  if(version2) {
    lay$name2 = lay %>% 
      mutate(
        name2 = case_when(
          name2 == "ppt" ~ "PPT", 
          name2=="ppgd" ~ "PPG", 
          name2=="so" ~ "SPT", 
          name2=="pmg"~"PMG", 
          name2=="nt"~"NGX", 
          name2=="pbmc"~"PGX", 
          TRUE~name2)
        ) %>% 
      pull(name2)
  }
  #} else {
  #  lay$name2 = lay$name
  #}
  ds_names = omic_short_names
  if(version2) {
    ds_names = c("PPT","PPG","SPT","PMG","NGX","PGX")
  }
  lay$node_class = case_when(
    lay$name2 %in% pathways ~ "path",
    lay$name2 %in% ds_names ~ "dataset",
    TRUE ~ "feature"
  )
  
  text_col = pals::kelly(11)[c(10, 3,4,5,6,8,11)]
  lay$node_text_group = case_when(
    lay$name2 %in% pathways ~ "path",
    lay$name2 %in% omic_short_names ~ lay$name2,
    TRUE ~ "feature"
  )
  
  
  lay$node_size1 = ifelse(lay$node_class == "path", 4,
                          ifelse(lay$node_class == "dataset", 3, 2))
  lay$node_size2 = ifelse( is.na(lay$log10qval), lay$node_size1, lay$log10qval+1 )

  idx = lay$name2 %in% rownames(feature_count)
  lay$node_size2[idx] = feature_count[ lay$name2[idx], "n"]
  if(version2) {
    lay$node_size2 = lay$node_size1
  }
  
  set_graph_style()
  
  # unset graph theme 
  # unset_graph_style()
  
  # add node names
  p = ggraph(lay) + 
    geom_edge_link(
      aes(color = as.factor(feat_coeff_dir)), 
      arrow = arrow(type = "open", length = unit(2, 'mm')),
      end_cap = circle(2, 'mm'),
      alpha = 0.8
    ) + 
    #scale_edge_color_distiller(palette = "RdBu", 
    #                      limits = c( -1*max(abs(edges$feat_coeff_dir), na.rm = T),
    #                                  max(abs(edges$feat_coeff_dir),na.rm = T) )
    #) +
    scale_edge_color_manual(values = c(pals::brewer.rdbu(11)[c(2)], "grey60", pals::brewer.rdbu(11)[c(10)]) %>% setNames(c(1,0,-1)) ) +
    ggtitle(title)
  if(REDUCED_GRAPH) {
    p = p + geom_node_point(aes(shape=node_class, size=node_size2, fill = node_class )) +
      scale_shape_manual(values = c(21,22,24) %>% setNames(c("dataset","path","feature")) ) +
      scale_fill_manual(values = brewer.set2(8)[c(5,6,8)] %>% setNames(c("feature","dataset","path")), na.value = "grey")
  } else {
    p = p + geom_node_point(aes(shape=node_class, size=node_size2, fill = as.character(score_dir) )) +
      scale_shape_manual(values = c(21,22,24) %>% setNames(c("dataset","path","feature")) ) +
      scale_fill_manual(values = c("1"="red", "-1"="blue"), na.value = "grey")
  }
  p = p + 
    #geom_node_text(aes(label = name2, color = node_text_group), repel=TRUE) +
    #scale_color_manual(values = c(text_col, "black") %>% setNames(c("path",omic_short_names, "feature")) )
    scale_color_manual(values = c(text_col[c(1,3)], "black") %>% setNames(c("path","feature","dataset")) ) +
    scale_size_continuous(range = c(1, 6))

  if(!excludeText) {
    p = p + geom_node_text(aes(label = name2, color = node_class), repel=TRUE)
  }
  if(useDejavuFonts) {
    p = p + theme_graph(base_family = "DejaVu Sans")
  }
  
  if(return_graph)
    return(list(g=g, lay=lay, nodes = nodes, edges = edges, plot = p))
  return(p)
}



trans_omics_PATHscores = function(base_obs, selected_kegg_pathways_pos,
                                  selected_kegg_pathways_neg,
                                  selected_sub_pathways_pos,
                                  selected_sub_pathways_neg,
                                  selected_hallmark_pathways_pos,
                                  selected_hallmark_pathways_neg, 
                                  others = NULL,
                                  phase = "train", adjust = TRUE){
  Xlist_kegg = list()
  Xlist_hallmark = list()
  Xlist_sub = list()
  Xlist_others = list()
  assay_short_names = names(base_obj$omic_names)
  if(phase == "train"){
    composite_scores =  base_obj$composite_features_train
    olinks =  base_obj$assays.train$serum_olink
    clin_mat = base_obj$clin.train
    rownams = rownames(base_obj$assays.train$plasma_proteomics_targeted)
  }else{
    composite_scores =  base_obj$composite_features_test
    olinks =  base_obj$assays.test$serum_olink
    clin_mat = base_obj$clin.test
    rownams = rownames(base_obj$assays.test$plasma_proteomics_targeted)
  }
  W =model.matrix(~as.character(clin_mat$trajectory_group)+clin_mat$event_type+clin_mat$discretized_admit_age_quantile+clin_mat$sex-1)
  
  for(omic_name0 in assay_short_names){
    if(omic_name0!="so"){
      tmp = composite_scores$KEGG_POS[[omic_name0]]
      if(!is.null(tmp)){
        tmp1 = tmp[,colnames(tmp)%in%intersect(colnames(tmp),selected_kegg_pathways_pos),drop = F]
        if(ncol(tmp1)>0){
          colnames(tmp1)=paste0(colnames(tmp1),"+")
        }
        tmp = composite_scores$KEGG_NEG[[omic_name0]]
        tmp2 = tmp[,colnames(tmp)%in%intersect(colnames(tmp),selected_kegg_pathways_neg),drop = F]
        if(ncol(tmp2)>0){
          colnames(tmp2)=paste0(colnames(tmp2),"-")
        }
        Xlist_kegg[[omic_name0]] = cbind(tmp1, tmp2)
      }

    }
  }
  Xlist_kegg[["so"]] =olinks[,rownames(base_obj$short_list_features$unsigned_mod$so)]
  for(omic_name0 in assay_short_names){
    if(omic_name0!="so"){
      tmp =  composite_scores$HALLMARK_POS[[omic_name0]]
      tmp1 = NULL
      tmp2 = NULL
      tmp1 = tmp[,colnames(tmp)%in%intersect(colnames(tmp),selected_hallmark_pathways_pos),drop = F]
      if(!is.null(tmp1)){
        if(ncol(tmp1)>0){
          colnames(tmp1)=paste0(colnames(tmp1),"+")
        }
      }
      tmp =  composite_scores$HALLMARK_NEG[[omic_name0]]
      tmp2 = tmp[,colnames(tmp)%in%intersect(colnames(tmp),selected_hallmark_pathways_neg),drop = F]
      if(!is.null(tmp2)){
        if(ncol(tmp2)>0){
          colnames(tmp2)=paste0(colnames(tmp2),"-")
        }
      }
      if(!is.null(tmp1) & !is.null(tmp2)){
        Xlist_hallmark[[omic_name0]] = cbind(tmp1, tmp2)
      }else if(is.null(tmp1)){
        Xlist_hallmark[[omic_name0]]=tmp2
      }else{
        Xlist_hallmark[[omic_name0]]=tmp1
      }
    }
  }
  Xlist_others = list()
  for(omic_name0 in assay_short_names){
    if(omic_name0!="so"){
    tmp =  composite_scores$SUB_POS[[omic_name0]]
    tmp1 = NULL
    tmp2 = NULL
    tmp1 = tmp[,colnames(tmp)%in%intersect(colnames(tmp),selected_sub_pathways_pos),drop = F]
    if(!is.null(tmp1)){
      if(ncol(tmp1)>0){
        colnames(tmp1)=paste0(colnames(tmp1),"+")
      }
    }
    tmp =  composite_scores$SUB_NEG[[omic_name0]]
    tmp2 = tmp[,colnames(tmp)%in%intersect(colnames(tmp),selected_sub_pathways_neg),drop = F]
    if(!is.null(tmp2)){
      if(ncol(tmp2)>0){
        colnames(tmp2)=paste0(colnames(tmp2),"-")
      }
    }
    if(!is.null(tmp1) &!is.null(tmp2)){
      Xlist_sub[[omic_name0]] = cbind(tmp1, tmp2)
    }else if(is.null(tmp1)){
      Xlist_sub[[omic_name0]]=tmp1
    }else{
      Xlist_sub[[omic_name0]]=tmp2
    }
    }
  }
  if(!is.null(others)){
    for(pathway_name0 in others){
      
      for(omic_name0 in assay_short_names){
        if(omic_name0!="so"){
        tmp =  composite_scores[[pathway_name0]][[omic_name0]]
        if(!is.null(tmp)){
          if(is.null(Xlist_others[[omic_name0]])){
            Xlist_others[[omic_name0]] = tmp
          }else{
            Xlist_others[[omic_name0]] = cbind(Xlist_others[[omic_name0]],tmp)
          }
          
        }
      }
      }
    }
  }
  id_full = apply(is.na(W),1,sum)==0
  for(d in 1:length(Xlist_kegg)){
    id_full = (id_full & apply(is.na(Xlist_kegg[[d]]),1,sum)==0)
  }
  for(d in 1:length(Xlist_sub)){
    id_full = (id_full & apply(is.na(Xlist_sub[[d]]),1,sum)==0)
  }
  for(d in 1:length(Xlist_hallmark)){
    id_full = (id_full & apply(is.na(Xlist_hallmark[[d]]),1,sum)==0)
  }
  if(length(Xlist_others)>0){
    for(d in 1:length(Xlist_others)){
      id_full = (id_full & apply(is.na(Xlist_others[[d]]),1,sum)==0)
    }
  }

  rownams= rownams[id_full]
  W1 = apply(W[id_full,],2,function(z) qqnorm(z, plot.it = F)$x)
  Xlist_kegg1 = list()
  for(d in 1:length(Xlist_kegg)){
    Xlist_kegg1[[d]]=Xlist_kegg[[d]][id_full,,drop = F]
  }
  Xlist_sub1 = list()
  for(d in 1:length(Xlist_sub)){
    Xlist_sub1[[d]]=Xlist_sub[[d]][id_full,,drop = F]
  }
  Xlist_hallmark1 = list()
  for(d in 1:length(Xlist_hallmark)){
    Xlist_hallmark1[[d]]=Xlist_hallmark[[d]][id_full,,drop = F]
  }
  Xlist_others1 = list()
  if(length(Xlist_others)>0){
    for(d in 1:length(Xlist_others)){
      Xlist_others1[[d]]=Xlist_others[[d]][id_full,,drop = F]
    }
  }

  names(Xlist_kegg1)=names(Xlist_kegg)
  names(Xlist_hallmark1)=names(Xlist_hallmark)
  names(Xlist_sub1)=names(Xlist_sub)
  names(Xlist_others1)=names(Xlist_others)
  for(d in 1:length(Xlist_kegg)){
    Xlist_kegg1[[d]]=apply(Xlist_kegg1[[d]],2,function(z) qqnorm(z, plot.it = F)$x)
    if(length(Xlist_kegg1[[d]])>0){
      if(is.null(dim(Xlist_kegg1[[d]]))){
        Xlist_kegg1[[d]]=data.frame(Xlist_kegg1[[d]])
        colnames(Xlist_kegg1[[d]]) = colnames(Xlist_kegg[[d]])
      }
    }
  }
  for(d in 1:length(Xlist_hallmark)){
    Xlist_hallmark1[[d]]=apply(Xlist_hallmark1[[d]],2,function(z) qqnorm(z, plot.it = F)$x)
  }
  for(d in 1:length(Xlist_sub)){
    Xlist_sub1[[d]]=apply(Xlist_sub1[[d]],2,function(z) qqnorm(z, plot.it = F)$x)
  }
  if(length(Xlist_others)>0){
  for(d in 1:length(Xlist_others)){
    Xlist_others1[[d]]=apply(Xlist_others1[[d]],2,function(z) qqnorm(z, plot.it = F)$x)
  }
  }
  if(adjust){
    for(d in 1:length(Xlist_kegg)){
      if(!is.null(dim(Xlist_kegg1[[d]] ))){
        Xlist_kegg1[[d]] = lm(as.matrix(Xlist_kegg1[[d]])~W1)$residuals
      }
      
    }
    for(d in 1:length(Xlist_hallmark)){
      if(!is.null(dim(Xlist_hallmark1[[d]] ))){
      Xlist_hallmark1[[d]] = lm(as.matrix(Xlist_hallmark1[[d]])~W1)$residuals
      }
    }
    for(d in 1:length(Xlist_sub)){
      if(!is.null(dim(Xlist_sub1[[d]] ))){
      Xlist_sub1[[d]] = lm(as.matrix(Xlist_sub1[[d]])~W1)$residuals
      }
    }
    if(length(Xlist_others)>0){
      for(d in 1:length(Xlist_others)){
        if(!is.null(dim(Xlist_others1[[d]] ))){
          Xlist_others1[[d]] = lm(as.matrix(Xlist_others1[[d]])~W1)$residuals
        }
      }
    }

  }
  Xlist=Xlist_kegg1
  Xlist_full = Xlist_kegg
  for(omic_name0 in names(Xlist_kegg)){
    if(omic_name0!=0){
      if(omic_name0%in%names(Xlist_hallmark1)){
        tmp = Xlist_hallmark1[[omic_name0]]
        tmp1 = Xlist_hallmark[[omic_name0]]
        if(length(tmp)>0){
          if(is.null(dim(tmp))){
            tmp = data.frame(tmp)
            colnames(tmp) = colnames(tmp1)
          }
        Xlist[[omic_name0]]=cbind(Xlist[[omic_name0]],tmp)
        Xlist_full[[omic_name0]] =cbind(Xlist_full[[omic_name0]],tmp1)
        }
      }
      if(omic_name0%in%names(Xlist_sub1)){
        tmp = Xlist_sub1[[omic_name0]]
        tmp1 = Xlist_sub[[omic_name0]]
        if(length(tmp)>0){
          if(is.null(dim(tmp))){
            tmp = data.frame(tmp)
            colnames(tmp) = colnames(tmp1)
          }
        idx = which(colnames(tmp)%in%colnames(Xlist_sub1[[omic_name0]]))
        if(length(idx)>0){
          colnames(tmp)[idx] = paste0(colnames(tmp)[idx],"(sub)")
          colnames(tmp1)[idx] = paste0(colnames(tmp1)[idx],"(sub)")
        }
        Xlist[[omic_name0]]=cbind(Xlist[[omic_name0]],tmp)
        Xlist_full[[omic_name0]] =cbind(Xlist_full[[omic_name0]],tmp1)
        }
      }
      if(omic_name0%in%names(Xlist_others1)){
        tmp = Xlist_others1[[omic_name0]]
        tmp1 =Xlist_others[[omic_name0]]
        if(length(tmp)>0){
        if(is.null(dim(tmp))){
          tmp = data.frame(tmp)
          colnames(tmp) = colnames(tmp1)
        }
        idx = which(colnames(tmp)%in%colnames(Xlist_others1[[omic_name0]]))
        if(length(idx)>0){
          colnames(tmp)[idx] = paste0(colnames(tmp)[idx],"(sub)")
          colnames(tmp1)[idx] = paste0(colnames(tmp1)[idx],"(sub)")
        }
        Xlist[[omic_name0]]=cbind(Xlist[[omic_name0]],tmp)
        Xlist_full[[omic_name0]] =cbind(Xlist_full[[omic_name0]],tmp1)
        }
      }
    }
  }
  X1 = Xlist[["pmg"]]
  X1_full = Xlist_full[["pmg"]]
  X2 = NULL
  X2_full = NULL
  assay_anno = c()
  for(d in 1:length(Xlist)){
    if(names(Xlist)[d]!="pmg"){
      if(is.null(X2)){
        X2=Xlist[[d]]
        X2_full=as.matrix(Xlist_full[[d]])
      }else{
        X2=cbind(X2, Xlist[[d]])
        X2_full=cbind(X2_full, as.matrix(Xlist_full[[d]]))
      }

      assay_anno=c(assay_anno,rep(names(Xlist)[d],ncol(Xlist[[d]])))
    }
  }
  X21 = X2[,assay_anno!="so"]
  X21_full = X2_full[,assay_anno!="so"]
  assay_anno21 = assay_anno[assay_anno!="so"]
  X22 = X2[,assay_anno == "so"]
  X22_full = X2_full[,assay_anno=="so"]
  colnames(X21) = paste0(colnames(X21),"_",assay_anno21 )
  colnames(X21_full) = paste0(colnames(X21_full),"_",assay_anno21 )
  
  correlationA = cor(X1,X21)
  correlationApval = correlationA
  correlationB = cor(X1,X22)
  correlationBpval = correlationB
  for(i in 1:nrow(correlationA)){
    for(j in 1:ncol(correlationA)){
      tmp = cor.test(X1[,i],X21[,j])
      correlationApval[i,j]=tmp$p.value
    }
  }
  for(i in 1:nrow(correlationB)){
    for(j in 1:ncol(correlationB)){
      tmp = cor.test(X1[,i],X22[,j])
      correlationBpval[i,j]=tmp$p.value
    }
  }
  rownames(X1) = rownams
  rownames(X22) = rownams
  rownames(X21) = rownams
  return(list(Xpmg = X1, Xso = X22, Xothers = X21,
              Xpmg_full = X1_full, Xso_full = X22_full, Xothers_full = X21_full,
              correlation_pmg_so = correlationB,
              correlation_pmg_so_pval = correlationBpval,
              correlation_pmg_others = correlationA,
              correlation_pmg_others_pval = correlationApval,
              Xlist =Xlist, Xlist_full = Xlist_full
              ))
}


# linear_effects_simple = function(feature_mat, clin_mat, adjust = FALSE){
#   feature_mat0=feature_mat
#   idx1 = which(clin_mat$event_type=="Visit 1")
#   feature_mat1=feature_mat[idx1,]
#   logP_signed = matrix(NA, ncol = 3, nrow = ncol(feature_mat))
#   colnames(logP_signed) = c("ordinal","mortality","slope5|4")
#   rownames(logP_signed) = colnames(feature_mat)
#   z = clin_mat$event_date- clin_mat$symptom_date
#   options(na.action='na.pass')
#   W = model.matrix(~sqrt(z)+clin_mat$discretized_admit_age_quantile+clin_mat$sex-1)
#   if(adjust){
#     feature_mat0 = residuals(lm(as.matrix(feature_mat0)~W, na.action=na.exclude))
#     W1 = W[idx1,]
#     feature_mat1= residuals(lm(as.matrix(feature_mat1)~W1, na.action=na.exclude))
#   }
#   for(j in 1:ncol(feature_mat)){
#     idx1 = which(clin_mat$event_type=="Visit 1")
#     y = clin_mat$trajectory_group[idx1]
#     tmp = cor.test(y,feature_mat1[,j], method = "pearson", use = "pairwise.complete")
#     logP_signed[j,1] = -log(tmp$p.value,base = 10) * sign(tmp$estimate)
#     
#     y = clin_mat$trajectory_group
#     idx1 = which(clin_mat$event_type=="Visit 1"&y%in%c(4,5))
#     y = clin_mat$trajectory_group[idx1]
#     tmp = cor.test(y, feature_mat0[idx1,j], method = "pearson", use = "pairwise.complete")
#     logP_signed[j,2] = -log(tmp$p.value,base = 10) * sign(tmp$estimate)
#     
#     y = clin_mat$trajectory_group
#     idx1 = which(clin_mat$event_type%in%paste0("Visit ",c(1:6))&clin_mat$event_date<=30& y%in%c(4,5))
#     t =sqrt(clin_mat$event_date[idx1])
#     y = y[idx1]
#     participant = clin_mat$participant_id[idx1]
#     x = qqnorm(feature_mat[idx1,j],plot.it=F)$x
#     delta = t
#     delta[y==4] = -delta[y==4]
#     dat_j = data.frame(endpoints=y, event_date = t, delta = delta, value = x, participant=participant)
#     dat_j=dat_j[apply(is.na(dat_j),1,sum)==0,]
#     formula_use = formula("value~ s(event_date, bs = 'cr')+delta+endpoints")
#     fit <- try(gamm4::gamm4(formula_use, data = dat_j, random = ~(1|participant)))
#     tmp = anova(fit$gam)
#     logP_signed[j,3]=sign(tmp$p.coeff[2])*(-log( tmp$pTerms.pv[1],base=10))
#   }
#   return(logP_signed)
# }




replace_rownames_with_event_id <- function(assay, clinical_data, event_ids_selected){
  tmp <- merge(assay, clinical_data[,c("sample_id", "event_id")], by.x = 0, by.y = "sample_id", all.x = TRUE)
  rownames(tmp) <- tmp$event_id
  tmp <- tmp[,colnames(assay)]
  if(!missing(event_ids_selected)){
    tmp <- tmp[rownames(tmp) %in% event_ids_selected, ]
  }
  return(tmp)
}

preprocess_bld_cytof_child_granulo <- function(counts) {
  bld_cytof_child_nongran <- counts %>%
    rownames_to_column(var = "event_id") %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Basophil`,
                  #-`Neutrophil (CD16low)`,
                  #-`Neutrophil (CD16hi)`,
                  #-`Eosinophil`,
                  -`Tier1_Undefined`,
                  -`Tier2_undefined_Undefined`
    ) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_child_nongran <- bld_cytof_child_nongran/rowSums(bld_cytof_child_nongran)
  bld_cytof_child_nongran <- log1p(bld_cytof_child_nongran)
  
  bld_cytof_child_nongran <- bld_cytof_child_nongran %>%
    apply(2, scale) %>%
    data.frame(check.names = F) %>%
    `rownames<-`(row.names(bld_cytof_child_nongran))
  
  return(value = bld_cytof_child_nongran)
}

preprocess_bld_cytof_child_wo_granulo <- function(counts) {
  bld_cytof_child_nongran <- counts %>%
    rownames_to_column(var = "event_id") %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Tier1_Undefined`,
                  -`Tier2_undefined_Undefined`,
                  -`Neutrophil (CD16low)`,
                  -`Neutrophil (CD16hi)`,
                  -`Eosinophil`,
                  -`Basophil`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_child_nongran <- bld_cytof_child_nongran/rowSums(bld_cytof_child_nongran)
  bld_cytof_child_nongran <- log1p(bld_cytof_child_nongran)
  
  bld_cytof_child_nongran <- bld_cytof_child_nongran %>%
    apply(2, scale) %>%
    data.frame(check.names = F) %>%
    `rownames<-`(row.names(bld_cytof_child_nongran))
  
  return(value = bld_cytof_child_nongran)
}

preprocess_bld_cytof_parent_wo_granulo <- function(counts, rowfeature) { 
  ## Log1p normalize and scale the major parent frequnencies
  bld_cytof_parent_counts_processed <- counts %>%
    rownames_to_column(var = "event_id") %>%
    pivot_longer(cols = -event_id, names_to = "rowname") %>%
    merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
    group_by(event_id, `Broad Parent Population`) %>%
    mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
    dplyr::summarize(value = sum(value)) %>%
    pivot_wider(names_from = `Broad Parent Population`) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Neutrophil`,
                  -`Eosinophil`,
                  -`Basophil`,
                  -`Undefined`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_parent_counts_processed <-   bld_cytof_parent_counts_processed/rowSums(bld_cytof_parent_counts_processed)
  
  # Log-transform 
  bld_cytof_parent_counts_processed <- log1p(bld_cytof_parent_counts_processed)
  
  bld_cytof_parent_counts_processed <- bld_cytof_parent_counts_processed %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(bld_cytof_parent_counts_processed))
  
  return(value = bld_cytof_parent_counts_processed)
}

preprocess_bld_cytof_parent_granulo <- function(counts, rowfeature) { 
  ## Log1p normalize and scale the major parent frequnencies
  bld_cytof_parent_counts_processed <- counts %>%
    rownames_to_column(var = "event_id") %>%
    pivot_longer(cols = -event_id, names_to = "rowname") %>%
    merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
    group_by(event_id, `Broad Parent Population`) %>%
    mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
    dplyr::summarize(value = sum(value)) %>%
    pivot_wider(names_from = `Broad Parent Population`) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  #-`Neutrophil`,
                  #-`Eosinophil`,
                  #-`Basophil`,
                  -`Undefined`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_parent_counts_processed <-   bld_cytof_parent_counts_processed/rowSums(bld_cytof_parent_counts_processed)
  
  # Log-transform 
  bld_cytof_parent_counts_processed <- log1p(bld_cytof_parent_counts_processed)
  
  bld_cytof_parent_counts_processed <- bld_cytof_parent_counts_processed %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(bld_cytof_parent_counts_processed))
  
  return(value = bld_cytof_parent_counts_processed)
}

preprocess_bld_cytof_parent_granulo_separate <- function(counts, rowfeature) { 
  ## Log1p normalize and scale the major parent frequnencies
  bld_cytof_parent_counts_processed <- counts %>%
    rownames_to_column(var = "event_id") %>%
    pivot_longer(cols = -event_id, names_to = "rowname") %>%
    merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
    group_by(event_id, `Broad Parent Population`) %>%
    mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
    dplyr::summarize(value = sum(value)) %>%
    pivot_wider(names_from = `Broad Parent Population`) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Neutrophil`,
                  -`Eosinophil`,
                  -`Basophil`,
                  -`Undefined`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_parent_counts_processed <-   bld_cytof_parent_counts_processed/rowSums(bld_cytof_parent_counts_processed)
  
  # Log-transform 
  bld_cytof_parent_counts_processed <- log1p(bld_cytof_parent_counts_processed)
  
  bld_cytof_parent_counts_processed <- bld_cytof_parent_counts_processed %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(bld_cytof_parent_counts_processed))
  
  ### Now perform normalization with granulocytes included
  ### (but will only take the granulocyte values)
  
  bld_cytof_parent_counts_processed_granulo <- counts %>%
    rownames_to_column(var = "event_id") %>%
    pivot_longer(cols = -event_id, names_to = "rowname") %>%
    merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
    group_by(event_id, `Broad Parent Population`) %>%
    mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
    dplyr::summarize(value = sum(value)) %>%
    pivot_wider(names_from = `Broad Parent Population`) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  #-`Neutrophil`,
                  #-`Eosinophil`,
                  #-`Basophil`,
                  -`Undefined`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_parent_counts_processed_granulo <-   bld_cytof_parent_counts_processed_granulo/rowSums(bld_cytof_parent_counts_processed_granulo)
  
  # Log-transform 
  bld_cytof_parent_counts_processed_granulo <- log1p(bld_cytof_parent_counts_processed_granulo)
  
  bld_cytof_parent_counts_processed_granulo <- bld_cytof_parent_counts_processed_granulo %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(bld_cytof_parent_counts_processed_granulo))
  
  bld_cytof_parent_counts_processed_final <- cbind.data.frame(bld_cytof_parent_counts_processed,
                                                              bld_cytof_parent_counts_processed_granulo[,
                                                                                                        !colnames(bld_cytof_parent_counts_processed_granulo) %in%
                                                                                                          colnames(bld_cytof_parent_counts_processed)])
  return(value = bld_cytof_parent_counts_processed_final)
}

enrichment_aggregation = function(pos_list, neg_list = NULL, thr = 0.2){
  pathway_names = c()
  max_list = list()
  for(j in 1:length(pos_list)){
    rownames(pos_list[[j]]) = pos_list[[j]]$ID
    pathway_names=c( pathway_names,pos_list[[j]]$ID, neg_list[[j]]$ID)
    if(!is.null(neg_list)){
      rownames(neg_list[[j]]) = neg_list[[j]]$ID
    }
    pathway_names=c( pathway_names,neg_list[[j]]$ID)
    
  }
  pathway_names = unique(  pathway_names)
  for(j in 1:length(pos_list)){
    mat0 =  pos_list[[j]]
    if(!is.null(neg_list)){
      mat0$sign = "+"
      mat1 = neg_list[[j]]
      for(l in 1:nrow( mat0)){
        p0 = mat0$ID[l]
        l1 = which(mat1$ID==p0)
        if(mat1$pvalue[l1]<(mat0$pvalue[l])){
          mat0[l,1:ncol(mat1)] = mat1[l1,]
          mat0$sign[l] = "-"
        }
      }
      #mat0$pvalue =mat0$pvalue*2
      mat0$pvalue=ifelse(mat0$pvalue>1.0,1.0, mat0$pvalue)
    }
    
    max_list[[j]] = mat0
    max_list[[j]] =  max_list[[j]][order( max_list[[j]]$pvalue),]
    max_list[[j]]$p.adj = p.adjust(max_list[[j]]$pvalue, method = "BH")
  }
  names( max_list) = names(pos_list)
  pval_mat = matrix(NA, ncol = length(pos_list), nrow = length(pathway_names))
  for(l in 1:length(pathway_names)){
    p0 = pathway_names[l]
    for(j in 1:length(max_list)){
      mat0 = max_list[[j]]
      j0 = which(rownames(mat0) == p0)
      if(length(j0)>0){
        pval0 = mat0$pvalue[j0]
        pval_mat[l,j] = pval0
      }
    }
  }
  pval_mat = data.frame(pval_mat)
  rownames(pval_mat) = pathway_names
  colnames(pval_mat) = names(pos_list)
  pval_mat_rd = pval_mat
  pval_mat_rd[is.na(pval_mat_rd)] = -1
  pval_mat_rd[pval_mat_rd>thr] = (thr+1.0)/2
  pval_mat_rd[pval_mat_rd == -1] = NA
  pval_tan = tan((0.5- pval_mat_rd) * pi)
  pval_combine = apply(pval_tan,1,mean, na.rm = T)
  pval_combine = pcauchy(pval_combine, location = 0, scale = 1, lower.tail = F, log.p = F)
  pval_combine[is.na(pval_combine)] = 1.0
  pval_mat$joint = pval_combine
  pval_mat=pval_mat[order(pval_mat$joint),]
  pval_mat$joint_adj = p.adjust(pval_mat$joint,method = "BH")
  return(list(max_list = max_list, pval_mat=pval_mat))
  
}

enrichment_aggregation1 = function(pos_list, neg_list = NULL, thr = 0.2, type = "mhg"){
  pathway_names = c()
  max_list = list()
  if( type == "gsea"){
    for(j in 1:length(pos_list)){
    pos_list[[j]] =pos_list[[j]]@result
    }
  }
  for(j in 1:length(pos_list)){
    rownames(pos_list[[j]]) = pos_list[[j]]$ID
    pathway_names=c( pathway_names,pos_list[[j]]$ID, neg_list[[j]]$ID)
    if(!is.null(neg_list)){
      rownames(neg_list[[j]]) = neg_list[[j]]$ID
    }
    pathway_names=c( pathway_names,neg_list[[j]]$ID)
    
  }
  pathway_names = unique(  pathway_names)
  for(j in 1:length(pos_list)){
    mat0 =  pos_list[[j]]
    if(!is.null(neg_list)){
      mat0$sign = "+"
      mat1 = neg_list[[j]]
      for(l in 1:nrow( mat0)){
        p0 = mat0$ID[l]
        l1 = which(mat1$ID==p0)
        if(mat1$pvalue[l1]<(mat0$pvalue[l])){
          mat0[l,1:ncol(mat1)] = mat1[l1,]
          mat0$sign[l] = "-"
        }
      }
      #mat0$pvalue =mat0$pvalue*2
      mat0$pvalue=ifelse(mat0$pvalue>1.0,1.0, mat0$pvalue)
    }else if(type == "gsea"){
      mat0$sign = mat0$NES
    }
    
    max_list[[j]] = mat0
    max_list[[j]] =  max_list[[j]][order( max_list[[j]]$pvalue),]
    max_list[[j]]$p.adj = p.adjust(max_list[[j]]$pvalue, method = "BH")
  }
  names( max_list) = names(pos_list)
  pval_mat = matrix(NA, ncol = length(pos_list), nrow = length(pathway_names))
  direction_mat = matrix(NA, ncol = length(pos_list), nrow = length(pathway_names))
  for(l in 1:length(pathway_names)){
    p0 = pathway_names[l]
    for(j in 1:length(max_list)){
      mat0 = max_list[[j]]
      j0 = which(rownames(mat0) == p0)
      if(length(j0)>0){
        pval0 = mat0$pvalue[j0]
        pval_mat[l,j] = pval0
        if(!is.null(mat0$sign)){
          direction_mat[l,j]=mat0$sign[j0]
        }
      }
    }
  }
  pval_mat = data.frame(pval_mat)
  rownames(pval_mat) = pathway_names
  colnames(pval_mat) = names(pos_list)
  direction_mat = data.frame(direction_mat)
  rownames(direction_mat) = pathway_names
  colnames(direction_mat) = names(pos_list)
  pval_mat_rd = pval_mat
  pval_mat_rd[is.na(pval_mat_rd)] = -1
  pval_mat_rd[pval_mat_rd>thr] = (thr+1.0)/2
  pval_mat_rd[pval_mat_rd == -1] = NA
  pval_tan = tan((0.5- pval_mat_rd) * pi)
  pval_combine = apply(pval_tan,1,mean, na.rm = T)
  pval_combine = pcauchy(pval_combine, location = 0, scale = 1, lower.tail = F, log.p = F)
  pval_combine[is.na(pval_combine)] = 1.0
  pval_mat$joint = pval_combine
  pval_mat$joint_adj = p.adjust(pval_mat$joint,method = "BH")
  return(list(max_list = max_list, pval_mat=pval_mat, direction_mat = direction_mat))
  
}

condensed_trajectory_plot = function(expr_list,
                                     model_DF,
                                     base_obj,
                                     max_date = 28,
                                     phase = "train",
                                     event_date = "admission", cex_all = 15){
  if(phase == "train"){
    clinDF = base_obj$clin.train
  }else{
    clinDF = base_obj$clin.test
  }
  plotDF = data.frame(visit = clinDF$event_type,
                      trajectory_group = clinDF$trajectory_group,
                      baseline = clinDF$resp_status_v1,
                      ADdate =clinDF$event_date,
                      participant_id = clinDF$participant_id,
                      enrollment_site =clinDF$enrollment_site,
                      sex =clinDF$sex,
                      discretized_admit_age_quantile = clinDF$discretized_admit_age_quantile,
                      DFSO = clinDF$event_date - clinDF$symptom_date)
  
  if(event_date=="admission"){
    plotDF$event_date =  plotDF$ADdate
    xlabel = "Days from Admission"
  }else{
    plotDF$event_date =  plotDF$DFSO
    xlabel = "DFSO"
  }
  for(j in 1:length(expr_list)){
    name_j = names(expr_list)[j]
    plotDF[[name_j]] = expr_list[[name_j]][,1]
  }
  plotDF=plotDF[!is.na(plotDF$trajectory_group) & !is.na(plotDF$event_date) & clinDF$event_type %in%paste0("Visit ",1:6),]
  plotDF = plotDF[plotDF$event_date <= max_date,]
  plotDF[names(expr_list)] = apply(plotDF[names(expr_list)],2,function(z){
    z1 = z[!is.na(z)]
    z1 = qqnorm(z1,plot.it=F)$x
    z[!is.na(z)] = z1
    return(z)
  })
  plotDF$sex = factor(plotDF$sex)
  plotDF$discretized_admit_age_quantile = factor(plotDF$discretized_admit_age_quantile)
  plotDF$participant_id = factor(plotDF$participant_id)
  plotDF$enrollment_site = factor(plotDF$enrollment_site)
  plotDF$trajectory_group = ordered(factor(plotDF$trajectory_group))
  plotDF = tidyr::gather(plotDF, key = "name",
                         value = "value",
                         -visit, -trajectory_group, -baseline, -ADdate, -DFSO,
                         -participant_id, -enrollment_site,-sex,
                         -discretized_admit_age_quantile, -event_date)
  colors = c("#639A21", "#39828C",  "#6371AD","#BD7D31", "#9C3418")
  plotDF=plotDF[!is.na(plotDF$value),]
  digit=0
  s=min(cex_all/length(expr_list),3.5)
  plot_list = list()
  tmp_list = list()
  for(j in 1:length(expr_list)){
    name_j = names(expr_list)[j]
    plotExample <- plotDF[plotDF$name == name_j,] ## plotExample would be a plotDF object example in 
    tmp_list[[j]]= plot_model(plotExample, model_loop = model_DF, modelType = "smoothSpline",
                              endpoint = "trajectory_group", 
                              CI = T, facet = F, individual_points = F, knot_lines = F, individual_trendlines = F, individual_paths = F, xlabel = xlabel) 
    p_shape = formatC(model_DF[name_j,"p.slope"], format = "e", digits = digit)
    p_average =  formatC(model_DF[name_j,"p.intercept"], format = "e", digits = digit)
    text1 = paste0("all.shp/ave:", p_shape,"/", p_average)
    p_shape = formatC(model_DF[name_j,"p.slope_4v5"], format = "e", digits = digit)
    p_average =  formatC(model_DF[name_j,"p.intercept_4v5"], format = "e", digits = digit)
    text2 = paste0("4|5.shp/ave:", p_shape,"/", p_average)
  }
  plot_join = ggpubr::ggarrange(plotlist = tmp_list, nrow=1, common.legend = TRUE, legend="bottom")
  
  plot_join=annotate_figure(plot_join, top = text_grob(term, 
                                                       color = "black", face = "bold", size = 12))
  
  plot_join
  
  return(plot_join)
}

condensed_trajectory_plotI = function(expr_list, TG, base_obj,
                                      model_DF = NULL,
                                      plot_name = NULL,
                                      max_date = 28,
                                      phase = "train",
                                      event_date = "admission", cex_all = 15,
                                      colors =  c("#639A21", "#39828C",  "#6371AD","#BD7D31", "#9C3418")){
  if(phase == "train"){
    clinDF = base_obj$clin.train
  }else{
    clinDF = base_obj$clin.test
  }
  plotDF = data.frame(visit = clinDF$event_type,
                      trajectory_group = TG,
                      baseline = clinDF$resp_status_v1,
                      ADdate =clinDF$event_date,
                      participant_id = clinDF$participant_id,
                      enrollment_site =clinDF$enrollment_site,
                      sex =clinDF$sex,
                      discretized_admit_age_quantile = clinDF$discretized_admit_age_quantile,
                      DFSO = clinDF$event_date - clinDF$symptom_date)
  
  if(event_date=="admission"){
    plotDF$event_date =  plotDF$ADdate
    xlabel = "Days from Admission"
  }else{
    plotDF$event_date =  plotDF$DFSO
    xlabel = "DFSO"
  }
  for(j in 1:length(expr_list)){
    name_j = names(expr_list)[j]
    plotDF[[name_j]] = expr_list[[name_j]][,1]
  }
  plotDF=plotDF[!is.na(plotDF$trajectory_group) & !is.na(plotDF$event_date) & clinDF$event_type %in%paste0("Visit ",1:6),]
  plotDF = plotDF[plotDF$event_date <= max_date,]
  plotDF[names(expr_list)] = apply(plotDF[names(expr_list)],2,function(z){
    z1 = z[!is.na(z)]
    z1 = qqnorm(z1,plot.it=F)$x
    z[!is.na(z)] = z1
    return(z)
  })
  plotDF$sex = factor(plotDF$sex)
  plotDF$discretized_admit_age_quantile = factor(plotDF$discretized_admit_age_quantile)
  plotDF$participant_id = factor(plotDF$participant_id)
  plotDF$enrollment_site = factor(plotDF$enrollment_site)
  plotDF$trajectory_group = ordered(factor(plotDF$trajectory_group))
  plotDF = tidyr::gather(plotDF, key = "name",
                         value = "value",
                         -visit, -trajectory_group, -baseline, -ADdate, -DFSO,
                         -participant_id, -enrollment_site,-sex,
                         -discretized_admit_age_quantile, -event_date)
  
  plotDF=plotDF[!is.na(plotDF$value),]
  digit=0
  s=min(cex_all/length(expr_list),3.5)
  if(is.null(model_DF)){
    model_DF = model_loop(plotDF, modelType = "smoothSpline")
  }
  plot_list = list()
  tmp_list = list()
  
  
  for(j in 1:length(expr_list)){
    name_j = names(expr_list)[j]
    plotExample <- plotDF[plotDF$name == name_j,] ## plotExample would be a plotDF object example in 
    tmp_list[[j]]= plot_model(plotExample, model_loop = model_DF, modelType = "smoothSpline",
                              endpoint = "trajectory_group", 
                              CI = T, facet = F, individual_points = F, knot_lines = F,
                              individual_trendlines = F, individual_paths = F, xlabel = xlabel, colors = colors) 
  }
  plot_join = ggpubr::ggarrange(plotlist = tmp_list, nrow=1, common.legend = TRUE, legend="bottom")
  
  plot_join=annotate_figure(plot_join, top = text_grob(plot_name, 
                                                       color = "black", face = "bold", size = 12))
  
  
  return(list(plot = plot_join, modelDF =model_DF))
}


scale_tidy <- function(x){
  qqnorm(x,plot.it = F)$x
}
model_cell_pop <- function(tmp.long, active_factor = 'Factor1'){
  out.p <- data.frame()
  for(active_name in unique(tmp.long$name)){
    
    tmp <- tmp.long %>%
      filter(name == active_name) 
    
    lme.formula.test <- formula(paste0(active_factor, " ~", "value + sex+  + discretized_admit_age_quantile"))
    lme.formula.base <- formula(paste0(active_factor, " ~", " sex+ discretized_admit_age_quantile"))
    fit <- lme(fixed = lme.formula.test,
               random = ~1|enrollment_site, data = tmp %>% as.data.frame() )
    
    out.p<- rbind.data.frame(out.p, data.frame("Feature" = active_name,
                                               "Factor" = active_factor,
                                               coefficient = fit$coefficients$fixed[2],
                                               "p.model" =  anova(fit, type = "marginal")$`p-value`[2]
    ))
    
  }
  return(out.p)
}
cor_loop <- function(mat1, mat2, method = "spearman"){
  out <- data.frame()
  for(column1 in colnames(mat1)){
    for(column2 in colnames(mat2)){
      tmp.test <- cor.test(mat1[,column1], mat2[,column2], method = "spearman")
      out <- rbind(out, cbind(column1, column2, tmp.test$estimate, tmp.test$p.value))
    }
  }
  rownames(out) <- 1:nrow(out); colnames(out) <- c("Factor", "Feature", "Rho", "p-value")
  out$Rho <- as.numeric(out$Rho); out$`p-value` <- as.numeric(out$`p-value`)
  out$p.adj <- p.adjust(out$`p-value`, method = "fdr")
  return(as.data.frame(out))
}
preprocess_bld_cytof_child_seperate_granulo <- function(counts){
  bld_cytof_child_counts_processed <- preprocess_bld_cytof_child_wo_granulo(counts)
  
  bld_cytof_child_gran <- counts %>%
    rownames_to_column(var = "event_id") %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Tier1_Undefined`,
                  -`Tier2_undefined_Undefined`#,
                  #-`Neutrophil (CD16low)`,
                  #-`Neutrophil (CD16hi)`,
                  #-`Eosinophil`,
                  #-`Basophil`
    ) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_child_gran <- bld_cytof_child_gran/rowSums(bld_cytof_child_gran)
  bld_cytof_child_gran <- log1p(bld_cytof_child_gran)
  
  bld_cytof_child_gran <- bld_cytof_child_gran %>%
    apply(2, scale) %>%
    data.frame(check.names = F) %>%
    `rownames<-`(row.names(bld_cytof_child_gran))
  
  final <- cbind.data.frame(bld_cytof_child_counts_processed,
                            bld_cytof_child_gran[,!colnames(bld_cytof_child_gran) %in% colnames(bld_cytof_child_counts_processed)])
  return(final)
}

preprocess_bld_cytof_child_granulo <- function(counts) {
  bld_cytof_child_nongran <- counts %>%
    rownames_to_column(var = "event_id") %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Basophil`,
                  #-`Neutrophil (CD16low)`,
                  #-`Neutrophil (CD16hi)`,
                  #-`Eosinophil`,
                  -`Tier1_Undefined`,
                  -`Tier2_undefined_Undefined`
    ) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_child_nongran <- bld_cytof_child_nongran/rowSums(bld_cytof_child_nongran)
  bld_cytof_child_nongran <- log1p(bld_cytof_child_nongran)
  
  bld_cytof_child_nongran <- bld_cytof_child_nongran %>%
    apply(2, scale) %>%
    data.frame(check.names = F) %>%
    `rownames<-`(row.names(bld_cytof_child_nongran))
  
  return(value = bld_cytof_child_nongran)
}

preprocess_bld_cytof_child_wo_granulo <- function(counts) {
  bld_cytof_child_nongran <- counts %>%
    rownames_to_column(var = "event_id") %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Tier1_Undefined`,
                  -`Tier2_undefined_Undefined`,
                  -`Neutrophil (CD16low)`,
                  -`Neutrophil (CD16hi)`,
                  -`Eosinophil`,
                  -`Basophil`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_child_nongran <- bld_cytof_child_nongran/rowSums(bld_cytof_child_nongran)
  bld_cytof_child_nongran <- log1p(bld_cytof_child_nongran)
  
  bld_cytof_child_nongran <- bld_cytof_child_nongran %>%
    apply(2, scale) %>%
    data.frame(check.names = F) %>%
    `rownames<-`(row.names(bld_cytof_child_nongran))
  
  return(value = bld_cytof_child_nongran)
}

preprocess_bld_cytof_parent_wo_granulo <- function(counts, rowfeature) { 
  ## Log1p normalize and scale the major parent frequnencies
  bld_cytof_parent_counts_processed <- counts %>%
    rownames_to_column(var = "event_id") %>%
    pivot_longer(cols = -event_id, names_to = "rowname") %>%
    merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
    group_by(event_id, `Broad Parent Population`) %>%
    mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
    dplyr::summarize(value = sum(value)) %>%
    pivot_wider(names_from = `Broad Parent Population`) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Neutrophil`,
                  -`Eosinophil`,
                  -`Basophil`,
                  -`Undefined`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_parent_counts_processed <-   bld_cytof_parent_counts_processed/rowSums(bld_cytof_parent_counts_processed)
  
  # Log-transform 
  bld_cytof_parent_counts_processed <- log1p(bld_cytof_parent_counts_processed)
  
  bld_cytof_parent_counts_processed <- bld_cytof_parent_counts_processed %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(bld_cytof_parent_counts_processed))
  
  return(value = bld_cytof_parent_counts_processed)
}

preprocess_bld_cytof_parent_granulo <- function(counts, rowfeature) { 
  ## Log1p normalize and scale the major parent frequnencies
  bld_cytof_parent_counts_processed <- counts %>%
    rownames_to_column(var = "event_id") %>%
    pivot_longer(cols = -event_id, names_to = "rowname") %>%
    merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
    group_by(event_id, `Broad Parent Population`) %>%
    mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
    dplyr::summarize(value = sum(value)) %>%
    pivot_wider(names_from = `Broad Parent Population`) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  #-`Neutrophil`,
                  #-`Eosinophil`,
                  #-`Basophil`,
                  -`Undefined`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_parent_counts_processed <-   bld_cytof_parent_counts_processed/rowSums(bld_cytof_parent_counts_processed)
  
  # Log-transform 
  bld_cytof_parent_counts_processed <- log1p(bld_cytof_parent_counts_processed)
  
  bld_cytof_parent_counts_processed <- bld_cytof_parent_counts_processed %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(bld_cytof_parent_counts_processed))
  
  return(value = bld_cytof_parent_counts_processed)
}

preprocess_bld_cytof_parent_granulo_separate <- function(counts, rowfeature) { 
  ## Log1p normalize and scale the major parent frequnencies
  bld_cytof_parent_counts_processed <- counts %>%
    rownames_to_column(var = "event_id") %>%
    pivot_longer(cols = -event_id, names_to = "rowname") %>%
    merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
    group_by(event_id, `Broad Parent Population`) %>%
    mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
    dplyr::summarize(value = sum(value)) %>%
    pivot_wider(names_from = `Broad Parent Population`) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  -`Neutrophil`,
                  -`Eosinophil`,
                  -`Basophil`,
                  -`Undefined`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_parent_counts_processed <-   bld_cytof_parent_counts_processed/rowSums(bld_cytof_parent_counts_processed)
  
  # Log-transform 
  bld_cytof_parent_counts_processed <- log1p(bld_cytof_parent_counts_processed)
  
  bld_cytof_parent_counts_processed <- bld_cytof_parent_counts_processed %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(bld_cytof_parent_counts_processed))
  
  ### Now perform normalization with granulocytes included
  ### (but will only take the granulocyte values)
  
  bld_cytof_parent_counts_processed_granulo <- counts %>%
    rownames_to_column(var = "event_id") %>%
    pivot_longer(cols = -event_id, names_to = "rowname") %>%
    merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
    group_by(event_id, `Broad Parent Population`) %>%
    mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
    dplyr::summarize(value = sum(value)) %>%
    pivot_wider(names_from = `Broad Parent Population`) %>%
    dplyr::select(-contains(match = "Debris"), 
                  -contains(match = "Multiplet"),
                  -`Platelets`, 
                  -`RBC`, 
                  #-`Neutrophil`,
                  #-`Eosinophil`,
                  #-`Basophil`,
                  -`Undefined`) %>%
    column_to_rownames(var = "event_id")
  bld_cytof_parent_counts_processed_granulo <-   bld_cytof_parent_counts_processed_granulo/rowSums(bld_cytof_parent_counts_processed_granulo)
  
  # Log-transform 
  bld_cytof_parent_counts_processed_granulo <- log1p(bld_cytof_parent_counts_processed_granulo)
  
  bld_cytof_parent_counts_processed_granulo <- bld_cytof_parent_counts_processed_granulo %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(bld_cytof_parent_counts_processed_granulo))
  
  bld_cytof_parent_counts_processed_final <- cbind.data.frame(bld_cytof_parent_counts_processed,
                                                              bld_cytof_parent_counts_processed_granulo[,
                                                                                                        !colnames(bld_cytof_parent_counts_processed_granulo) %in%
                                                                                                          colnames(bld_cytof_parent_counts_processed)])
  return(value = bld_cytof_parent_counts_processed_final)
}

load_blood_cell_markers <- function(setting = "protein", transcript_rowfeatures, sheet = 1){
  require(readxl)
  expression_of_cells <- read_excel("/scratch/data-integration/pipelines/Interpretation/Reference_databases/41590_2011_BFni2067_MOESM29_ESM.xlsx", sheet = sheet)
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
  ziegler_full <- read_excel("/scratch/data-integration/pipelines/Interpretation/Reference_databases/Nasal_upper_airway_gene_markers_scRNAseq.xlsx", sheet = 1)
  ordovas <- read_excel("/scratch/data-integration/pipelines/Interpretation/Reference_databases/Nasal_upper_airway_gene_markers_scRNAseq.xlsx", sheet = 2)
  ziegler <- read_excel("/scratch/data-integration/pipelines/Interpretation/Reference_databases/Nasal_upper_airway_gene_markers_scRNAseq.xlsx", sheet = 3)
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


CyTOF_normalize <- function(counts, population = "childen", granulocytes = FALSE, rowfeature, log1p = TRUE){
  exclude <- c("Debris (MDIPA)", "Debris (CD45low)", "Debris","Multiplet", "Multiplet (Mon_NK)", "Multiplet (Neu_B)","Multiplet (Neu_Mon)","Multiplet (Neu_NK)", 
               "Multiplet (Neu_T)", "Multiplet (B_T)", "Multiplet (Mon_B)", "RBC", "Platelets", "Tier1_Undefined", "Tier2_undefined_Undefined", "Undefined")
  exclude_all <- c(exclude, "Neutrophil (CD16low)", "Neutrophil (CD16hi)", "Neutrophil", "Eosinophil", "Basophil")
  if(population %in% c("children", "child")){
    tmp.data <- counts %>%
      rownames_to_column(var = "event_id") %>%
      dplyr::select(!any_of(exclude_all)) %>%
      column_to_rownames(var = "event_id")
  } else if (population == "parent"){
    tmp.data <- counts %>%
      rownames_to_column(var = "event_id") %>%
      pivot_longer(cols = -event_id, names_to = "rowname") %>%
      merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
      group_by(event_id, `Broad Parent Population`) %>%
      mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
      dplyr::summarize(value = sum(value)) %>%
      pivot_wider(names_from = `Broad Parent Population`) %>%
      dplyr::select(!any_of(exclude_all)) %>%
      column_to_rownames(var = "event_id")
  }
  tmp.data <-   tmp.data/rowSums(tmp.data)
  
  if(log1p){
    tmp.data <- log1p(tmp.data)
    
    tmp.data <- tmp.data %>%
      apply(2, scale) %>%
      data.frame(check.names = F) %>%
      `rownames<-`(row.names(tmp.data))
  }
  
  
  if(granulocytes){
    if(population %in% c("children", "child")){
      tmp.data.granulo <- counts %>%
        rownames_to_column(var = "event_id") %>%
        dplyr::select(!any_of(exclude)) %>%
        column_to_rownames(var = "event_id")
    } else if (population == "parent"){
      tmp.data.granulo <- counts %>%
        rownames_to_column(var = "event_id") %>%
        pivot_longer(cols = -event_id, names_to = "rowname") %>%
        merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
        group_by(event_id, `Broad Parent Population`) %>%
        mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
        dplyr::summarize(value = sum(value)) %>%
        pivot_wider(names_from = `Broad Parent Population`) %>%
        dplyr::select(!any_of(exclude)) %>%
        column_to_rownames(var = "event_id")
    }
    tmp.data.granulo <-   tmp.data.granulo/rowSums(tmp.data.granulo)
    
    if(log1p){
      tmp.data.granulo <- log1p(tmp.data.granulo)
      
      tmp.data.granulo <- tmp.data.granulo %>%
        apply(2, scale) %>%
        data.frame(check.names = F) %>%
        `rownames<-`(row.names(tmp.data.granulo))
    }
    
    ## Add only the columns not found in the non-granulo data 
    tmp.data <- cbind.data.frame(tmp.data,tmp.data.granulo[, !colnames(tmp.data.granulo) %in% colnames(tmp.data)])
    
  }
  
  return(tmp.data)
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
