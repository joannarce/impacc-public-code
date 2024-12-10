#Authors: Ravi, Jeremy, Anna, Cole, Pramod, Leying
pkgs = c("tidyverse",
         "ggplot2",
         "pals",
         "stringr",
         "ggpubr",
         "cowplot",
         "RColorBrewer",
         "rlang",
         "gridExtra",
         "corrr",
         "ComplexHeatmap")
for(i in 1:length(pkgs)){
  require(pkgs[i], character.only = TRUE)
}

predict = stats::predict # otherwise masked by MOFA2::predict
library(qvalue)
library(lme4)
library(ordinal)
library(nlme)
library(gamm4)
library(mgcv)
library(ggsignif)
###Small functions
### Assign elements of a list to variables
assign_env_vars <- function(env) {
  lnames <- names(list)
  if( ! is.null( lnames ) ) {
    for( n in lnames ) {
      assign( n, list[[n]] )
    }
  }
}


get.event_id.from.sample_id <- function(s_id, clinical_data){
  dict <- distinct(clinical_data, sample_id, event_id) %>% filter(sample_id==s_id)
  if(nrow(dict) > 0){
    return(dict$event_id[1])
  } else {
    return(NA)
  }
}

get.training_test.from.sample_id <- function(s_id, clinical_data){
  dict <- distinct(clinical_data, sample_id, training_test) %>% filter(sample_id==s_id)
  if(nrow(dict) > 0){
    return(dict$training_test[1])
  } else {
    return(NA)
  }
}




coefficient_name_translator = function(data_env, projection_coef_factor, omic_names){
  for(j in 1:length(omic_names)){
    omic_name = omic_names[j]
    if(omic_name %in% c("pbmc_transcriptomics","nasal_transcriptomics", "plasma_metabolomics_global")){
      projection_coef_factor[[j]] = name_translator(data_env = data_env,
                                                  DF = projection_coef_factor[[j]], 
                                                  omic_name=omic_name,
                                                  mod = "row")
    }
  }
  return(projection_coef_factor)
}

cutlist = function(projection_coef, projection_pval, mag_thr = 0.2, FWER_thr = 1e-2, direction = NULL){
  short_list_coef = projection_coef
  for(j in 1:length(projection_coef)){
    vec1 = projection_pval[[j]][,1]
    vec2 = projection_coef[[j]][,1]
    if(is.null(direction)){
      vec2 = abs(vec2)
    }else if(direction == -1){
      vec2 = -vec2
    }else if(direction == 1){
      vec2 = vec2
    }else{
      vec2 = abs(vec2)
    }
    FWER = vec1 * length(vec1)
    short_list_coef[[j]] = projection_coef[[j]][FWER<FWER_thr&vec2>mag_thr,,drop = F]
  }
  return(short_list_coef)
}

runGSEA = function(projection_coef, Term2GeneList){
  gsea_result = list()
  for(j in 1:length(projection_coef)){
    if(names(projection_coef)[j] %in% names(Term2GeneList)){
      geneList = projection_coef[[j]][,1]
      names(geneList) = rownames(projection_coef[[j]])
      geneList = sort(geneList, decreasing = T)
      Term2Gene = Term2GeneList[[names(projection_coef)[j]]]
      Term2Gene=Term2Gene[Term2Gene[,2]%in%names(geneList),]
      gsea_result[[names(projection_coef)[j]]] = clusterProfiler::GSEA(geneList = geneList,
                                                                       minGSSize = 3,
                                                                       maxGSSize = 500,
                                                                       pvalueCutoff = 2,
                                                                       TERM2GENE = Term2Gene)
    }
  }
  return(gsea_result)
}

runMGH = function(short_list_coef, Term2GeneList, projection_coef, backround = "impacc"){
  enrichr_result = list()
  for(j in 1:length(short_list_coef)){
    if(names(short_list_coef)[j] %in% names(Term2GeneList)){
      print(names(short_list_coef)[j])
      geneList = short_list_coef[[j]][,1]
      names(geneList) = rownames(short_list_coef[[j]])
      geneList = abs(geneList)
      geneList = sort(geneList,decreasing = T)
      Term2Gene = Term2GeneList[[names(short_list_coef)[j]]]
      universe = NULL
      if(backround == "impacc"){
        j0 = which(names(projection_coef)==names(short_list_coef)[j])
        universe = rownames(projection_coef[[j0]])
      }else{
        universe = unique(Term2Gene[,2])
      }
      #minGSSize = term overlapping with universe>=2
      if(length( geneList)>0){
        enrichr_result[[names(short_list_coef)[j]]] = mgh_test_wrapper(geneList = geneList,
                                                                       minGSSize = 2,
                                                                       maxGSSize = 500,
                                                                       pvalueCutoff = 2,
                                                                       universe = universe,
                                                                       TERM2GENE = Term2Gene,
                                                                       p.adjust.methods = "BH")
      }
      
    }
  }
  return(enrichr_result)
}
###########################Data alignment and MOFA#######

data_prepare_DR = function(omics_used, data_env, clinical_data, rna_keep = 5000){
  # assay data
  datasets.mofa = list()
  #### Identify outlying samples and filter out them.
  for( i in 1:length(omics_used) ) {
    o_name = omics_used[i]
    count_df <- data_env[[omics_used[i]]]
    datasets.mofa[[o_name]] = count_df
  }
  datasets <- datasets.mofa
  train.datasets <- list()
  test.datasets <- list()
  datasets.names <- names(datasets.mofa)
  for(i in 1:length(datasets)){
    dataset.name <- datasets.names[i]
    dataset <- datasets[[i]]
    print(dataset.name)
    ##filter and transformation
    ###filtering transcriptomics
    sds = apply(dataset, 2, sd)
    dataset = dataset[,!is.na(sds)]
    sds = apply(dataset, 2, sd)
    nas = apply(dataset!=0, 2, sum)
    if(grepl("transcriptomics",dataset.name) & !is.null(rna_keep)){
      thr_sd =sort(sds,decreasing = T)[rna_keep]
    }else{
      thr_sd = -Inf
    }
    dataset = dataset[,sds>=thr_sd & nas >= 10]
    
    # Change sample name to event_id
    train_test = sapply(rownames(dataset), function(z) get.training_test.from.sample_id(z, clinical_data=clinical_data))
    train.dataset = dataset[which(train_test == "training"),]
    test.dataset = dataset[which(train_test == "test"),]
    
    tmp =  sapply(rownames(train.dataset), function(z) get.event_id.from.sample_id(z, clinical_data=clinical_data))
    train.dataset <- train.dataset[!is.na(tmp),]
    rownames(train.dataset) <- tmp
    train.datasets[[i]] <- train.dataset
    
    tmp =  sapply(rownames(test.dataset), function(z) get.event_id.from.sample_id(z, clinical_data=clinical_data))
    test.dataset <- test.dataset[!is.na(tmp),]
    rownames(test.dataset) <- tmp
    test.datasets[[i]] <- test.dataset
  }
  
  # Name the training/testing
  names(train.datasets) <- names(test.datasets) <- names(datasets)
  
  return(list(train = train.datasets, test = test.datasets))
}





#######################Marginal test of Modules
mixed_ordinal = function(my.formula0, my.formula1, data_use){
  fit = ordinal::clmm(my.formula1, data = data_use)
  av =  ordinal:::anova.clm(ordinal::clmm(my.formula0, data = data_use),
             ordinal::clmm(my.formula1, data = data_use))
  res = rep(0, 3)
  names(res) = c("AIC", "pval")
  res[1] = av$AIC[2]
  res[2] = (av$`Pr(>Chisq)`)[2]
  res[3] = ifelse(fit$coefficients[5] > 0, "Severe", "Mild")
  return(res)
}

mixed_pairwise_coef = function(my.formula0, my.formula1, data_use){
  endpoints0 = sort(unique(data_use$endpoints))
  pair_names = c()
  res_table = c()
  for(i in 1:(length(endpoints0)-1)){
    for(j in (i+1):(length(endpoints0))){
      pair_names = c(pair_names, paste0(endpoints0[i],"|", endpoints0[j]))
      data_tmp = data_use[data_use$endpoints == endpoints0[i] | data_use$endpoints == endpoints0[j],]
      fit = lme4::lmer(my.formula1,data = data_tmp)
      res_table = c(res_table, ifelse(fit@beta[2] > 0, "Severe", "Mild"))
    }
  }
  names(res_table) = pair_names
  return(res_table)
}

mixed_pairwise = function(my.formula0, my.formula1, data_use){
  endpoints0 = sort(unique(data_use$endpoints))
  pair_names = c()
  res_table = c()
  for(i in 1:(length(endpoints0)-1)){
    for(j in (i+1):(length(endpoints0))){
      pair_names = c(pair_names, paste0(endpoints0[i],"|", endpoints0[j]))
      data_tmp = data_use[data_use$endpoints == endpoints0[i] | data_use$endpoints == endpoints0[j],]
      av = anova(lme4::lmer(my.formula0,data = data_tmp),
                 lme4::lmer(my.formula1,data = data_tmp))
      res_table = c( res_table, av$`Pr(>Chisq)`[2])
    }
  }
  names(res_table) = pair_names
  return(res_table)
}


visit1_analysis_module_func = function(data_use, module_names, endpoints = "endpoints"){
  # Association test
  # a model with only intercept and random effect across sites
  res_table_ordinal = data.frame(matrix(0, ncol = 4, nrow = length(module_names)))
  colnames(res_table_ordinal) = c("AIC", "pval", "qval", "direction")
  rownames(res_table_ordinal) = c(module_names)
  for(j in 1:length(module_names)){
    print(j)
    my.formula0 = as.formula(paste0(endpoints,'~ 1+(1|sites)+age+sex'))
    my.formula1 = as.formula(paste0(endpoints,'~',module_names[j],'+(1|sites)+age+sex'))
    res_table_ordinal[j,c(1,2,4)]=mixed_ordinal(my.formula0,  my.formula1, data_use)
  }
  res_table_ordinal$AIC = as.numeric(res_table_ordinal$AIC)
  res_table_ordinal$pval = as.numeric(res_table_ordinal$pval)
  res_table_ordinal$qval =  p.adjust(res_table_ordinal$pval, method = "BH")
  for(j in 1:length(module_names)){
    print(j)
    my.formula0 = paste0(module_names[j],"~ (1|sites)+age+sex")
    my.formula1 = paste0(module_names[j],"~  endpoints+(1|sites)+age+sex")
    if(j == 1){
      tmp = mixed_pairwise(my.formula0, my.formula1, data_use)
      tmp_coef = mixed_pairwise_coef(my.formula0, my.formula1, data_use)
      res_table_pairwise_pvalue = data.frame(matrix(0, ncol = length(tmp), nrow = length(module_names)))
      res_table_pairwise_coef= data.frame(matrix(0, ncol = length(tmp), nrow = length(module_names)))
      res_table_pairwise_pvalue[j,] = tmp
      res_table_pairwise_coef[j,] = tmp_coef
      colnames(res_table_pairwise_pvalue) = names(tmp)
      colnames(res_table_pairwise_coef) = names(tmp_coef)
    }else{
      res_table_pairwise_pvalue[j,] = mixed_pairwise(my.formula0, my.formula1, data_use)
      res_table_pairwise_coef[j,] <- mixed_pairwise_coef(my.formula0, my.formula1, data_use)
    }
  }
  
  rownames(res_table_pairwise_pvalue) = module_names
  rownames(res_table_pairwise_coef) = module_names
  res_table_pairwise_qvalue = apply(res_table_pairwise_pvalue, 2, function(x){
    stats::p.adjust(x,method = "BH")})
  if(is.null(dim(res_table_pairwise_qvalue))){
    res_table_pairwise_qvalue = matrix(res_table_pairwise_qvalue,nrow = 1)
  }
  
  res_table_pairwise_qvalue = data.frame(res_table_pairwise_qvalue)
  colnames(res_table_pairwise_qvalue) = colnames(res_table_pairwise_pvalue)
  rownames(res_table_pairwise_qvalue) = rownames(res_table_pairwise_pvalue)
  return(list(res_table_ordinal =res_table_ordinal,
              res_table_pairwise_coef=res_table_pairwise_coef,
              res_table_pairwise_pvalue=res_table_pairwise_pvalue,
              res_table_pairwise_qvalue=res_table_pairwise_qvalue))
}


pairwise_format <- function(pairwise_table){
  pairwise_table12 <- data.frame(module = rownames(pairwise_table),group1="1",group2="2",q.val = pairwise_table$`1|2`)
  pairwise_table23 <- data.frame(module = rownames(pairwise_table),group1="2",group2="3",q.val = pairwise_table$`2|3`)
  pairwise_table34 <- data.frame(module = rownames(pairwise_table),group1="3",group2="4",q.val = pairwise_table$`3|4`)
  pairwise_table45 <- data.frame(module = rownames(pairwise_table),group1="4",group2="5",q.val = pairwise_table$`4|5`)
  pairwise_table13 <- data.frame(module = rownames(pairwise_table),group1="1",group2="3",q.val = pairwise_table$`1|3`)
  pairwise_table24 <- data.frame(module = rownames(pairwise_table),group1="2",group2="4",q.val = pairwise_table$`2|4`)
  pairwise_table35 <- data.frame(module = rownames(pairwise_table),group1="3",group2="5",q.val = pairwise_table$`3|5`)
  pairwise_table14 <- data.frame(module = rownames(pairwise_table),group1="1",group2="4",q.val = pairwise_table$`1|4`)
  pairwise_table25 <- data.frame(module = rownames(pairwise_table),group1="2",group2="5",q.val = pairwise_table$`2|5`)
  pairwise_table15 <- data.frame(module = rownames(pairwise_table),group1="1",group2="5",q.val = pairwise_table$`1|5`)
  pairwise_table_all <- rbind(pairwise_table12,pairwise_table23,pairwise_table34,pairwise_table45,pairwise_table13,pairwise_table24,pairwise_table35,pairwise_table14,pairwise_table25,pairwise_table15)
  pairwise_table_all$q.val <- signif(pairwise_table_all$q.val,3)
  pairwise_table_all$q.signif <- cut(pairwise_table_all$q.val, 
                                     breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                     labels = c("****", "***", "**", "*", "ns"))
  pairwise_table_all$group1 = paste0("TG",pairwise_table_all$group1)
  pairwise_table_all$group2 = paste0("TG",pairwise_table_all$group2)
  return(pairwise_table_all)
}

ordinal_format <- function(ordinal_table){
  ordinal_table$module <- rownames(ordinal_table)
  ordinal_table$qval <- signif(ordinal_table$qval,3)
  ordinal_table$q.signif <- cut(ordinal_table$qval, 
                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                labels = c("****", "***", "**", "*", "ns"))
  ordinal_table$modp <- paste0(ordinal_table$module, ": ", ordinal_table$qval)
  ordinal_table$modsig <- paste0(ordinal_table$module, ": ", ordinal_table$q.signif) 
  
  return(ordinal_table)
}



visit1_box_plot_func = function(data_use,  res_table_pairwise, res_table,
                                modules = NULL, ncol = 4, show_legend = FALSE,
                                label = "q.signif", show_significance = TRUE,
                                colors = c("#639A21", "#39828C",  "#6371AD","#BD7D31", "#9C3418")){
  
  pairwise_res_table_formatted= pairwise_format(pairwise_table = res_table_pairwise)
  ordinal_res_table_formatted= ordinal_format(ordinal_table = res_table)
  ordinal_res_table_formatted = ordinal_res_table_formatted[,c(4,6,7)]
  data_plot <- data_use[,colnames(data_use)%in%rownames(res_table)]
  if(is.null(dim(data_plot))){
    data_plot = data.frame(data_plot)
    colnames(data_plot)=rownames(res_table)
  }
  data_plot$endpoints = data_use$endpoints
  data_plot<- data_plot%>%gather(key = 'module', value = 'expression', -endpoints)
  data_plot$expression = as.numeric(data_plot$expression)
  # Adding p-values:
  data_plot$p.val <- sapply(data_plot$module, function(mod){
    return(res_table[which(rownames(res_table) == mod), 2])
  })
  # Adding q-values:
  data_plot$q.val <- sapply(data_plot$module, function(mod){
    return(res_table[which(rownames(res_table) == mod), 3])
  })
  
  module_pvalue <- ordinal_res_table_formatted[,2]
  names(module_pvalue) <- rownames(ordinal_res_table_formatted)
  module_psig <- ordinal_res_table_formatted[,3]
  names(module_psig) <- rownames(ordinal_res_table_formatted)
  if(!is.null(modules)){
    data_plot_01 <- data_plot[which(data_plot$module %in% modules),]
    res_table_pairwise_all01 <- pairwise_res_table_formatted[which(pairwise_res_table_formatted$module%in%modules),]
    module_pvalue = module_pvalue[modules]
    module_psig = module_psig[modules]
  }else{
    data_plot_01 <- data_plot
    res_table_pairwise_all01 <- pairwise_res_table_formatted
  }
  data_plot_01$endpoints = factor(paste0("TG",data_plot_01$endpoints),levels = c("TG1","TG2","TG3","TG4","TG5"))
  
  mod01_box0 <- ggboxplot(data_plot_01, x = "endpoints", y = "expression",
                          color = "endpoints", palette =colors,
                          add = "jitter") +
    facet_wrap(~ module,labeller = labeller(module_pvalue), ncol = min(ncol, length(unique(data_plot_01$module))))+
    theme_bw()
  if(!show_legend){
    mod01_box0<- mod01_box0 + theme(legend.position="none")
  }
  if(show_significance){
    ll1 = (res_table_pairwise_all01$group1=="TG4")&(res_table_pairwise_all01$group2=="TG5")
    ll2 = (res_table_pairwise_all01$group1=="TG5")&(res_table_pairwise_all01$group2=="TG4")
    res_table_pairwise_all02 = res_table_pairwise_all01[ll1|ll2,]
    res_table_pairwise_all02 =  res_table_pairwise_all02[res_table_pairwise_all02$q.signif!="ns",]
    if(nrow(res_table_pairwise_all02)>0){
      mod01_box0<- mod01_box0 +  stat_pvalue_manual(data = res_table_pairwise_all02, 
                                                    y.position = quantile(data_plot_01$expression,0.99), step.increase = 0.05, step.group.by = "module", bracket.size = 0.2,  remove.bracket  = F,
                                                    label = "q.signif") 
    }

  }
  return(mod01_box0)
}


run_visit1_module_analysis = function(impacc_analysis_record, phase = "train",
                                      show_legend = FALSE, 
                                      show_significance = TRUE,
                                      flipping = NULL){
  irank = impacc_analysis_record$rank_selected
  module_names = paste0("factor",1:irank)
  if(is.null(flipping)){
    flipping = rep(F, irank)
  }
  if(phase == "train"){
    factors_used = apply(impacc_analysis_record$mcia.factors.train,2,function(z) qqnorm(z,plot.it = F)$x)
    idx = which(impacc_analysis_record$clin.train$event_type %in% c("Visit 1"))
    clinical_mat.temp =impacc_analysis_record$clin.train[idx,]
    y = impacc_analysis_record$clin.train$trajectory_group[idx]
  }else{
    factors_used = apply(impacc_analysis_record$mcia.factors.test,2,function(z) qqnorm(z,plot.it = F)$x)
    idx = which(impacc_analysis_record$clin.test$event_type %in% c("Visit 1"))
    clinical_mat.temp =impacc_analysis_record$clin.test[idx,]
    y = impacc_analysis_record$clin.test$trajectory_group[idx]
  }
  for(j in 1:irank){
    if(flipping[j]){
      factors_used[,j] = factors_used[,j] * (-1)
    }
  }

  factors_used = factors_used[,1:irank]
  x = factors_used[idx,]
  colnames(x) = module_names
  data_use = data.frame(x)
  colnames(data_use) = module_names
  data_use$endpoints = y
  data_use$sites = clinical_mat.temp$enrollment_site
  data_use$age = clinical_mat.temp$discretized_admit_age_quantile
  data_use$sex =  clinical_mat.temp$sex
  data_use$control = factor(1:nrow(data_use))
  # Only use samples that have valid endpoints (not NA/missing)
  data_use = data_use[!is.na(data_use$endpoints),]
  data_use$endpoints = factor(data_use$endpoints) 
  visit1_mod_train_result = visit1_analysis_module_func(data_use = data_use,module_names =  module_names, endpoints = "endpoints")
  
  ###plot individual and pooled figures
  figures_visit1_module = list()
  for(j in 1:irank){
    figures_visit1_module[[j]] = visit1_box_plot_func(data_use,
                                                      res_table_pairwise=visit1_mod_train_result$res_table_pairwise_pvalue,
                                                      res_table=visit1_mod_train_result$res_table_ordinal,
                                                      modules = module_names[j], ncol = 4, show_legend = show_legend,
                                                      label = "q.signif", show_significance =   show_significance,
                                                      colors = c("#639A21", "#39828C",  "#6371AD","#BD7D31", "#9C3418"))
  }
  names(figures_visit1_module) = module_names
  figures_visit1_module[["pooled"]] = visit1_box_plot_func(data_use,
                                                           res_table_pairwise=visit1_mod_train_result$res_table_pairwise_pvalue,
                                                           res_table=visit1_mod_train_result$res_table_ordinal,
                                                           modules = module_names, ncol = 4, show_legend = show_legend,
                                                           label = "q.signif", show_significance =  show_significance ,
                                                           colors = c("#639A21", "#39828C",  "#6371AD","#BD7D31", "#9C3418"))
  
  return(list(visit1_mod_train_result=visit1_mod_train_result,
              figures_visit1_module=figures_visit1_module))
  
}

run_visit1_analysis_inputresponse = function(x,  impacc_analysis_record, fname = NULL,phase = "train",
                                             show_legend = FALSE, 
                                             show_significance = TRUE){
  x =qqnorm(x,plot.it = F)$x[,1]
  idx0 = which(!is.na(x))
  if(phase == "train"){
    idx = intersect(which(impacc_analysis_record$clin.train$event_type %in% c("Visit 1")),idx0)
    clinical_mat.temp =impacc_analysis_record$clin.train[idx,]
    y = impacc_analysis_record$clin.train$trajectory_group[idx]
  }else{
    idx = intersect(which(impacc_analysis_record$clin.test$event_type %in% c("Visit 1")),idx0)
    clinical_mat.temp =impacc_analysis_record$clin.test[idx,]
    y = impacc_analysis_record$clin.test$trajectory_group[idx]
  }
  x=x[idx]
  x = data.frame(x)
  if(is.null(fname)){
    fname = "feature"
  }
  colnames(x) = fname
  data_use = x
  data_use$endpoints = y
  data_use$sites = clinical_mat.temp$enrollment_site
  data_use$age = clinical_mat.temp$discretized_admit_age_quantile
  data_use$sex =  clinical_mat.temp$sex
  data_use$control = factor(1:nrow(data_use))
  data_use = data_use[!is.na(data_use$endpoints),]
  data_use$endpoints = factor(data_use$endpoints) 
  visit1_mod_train_result = visit1_analysis_module_func(data_use = data_use,module_names =  fname, endpoints = "endpoints")
  
  ###plot individual and pooled figures
  figures_visit1_module = visit1_box_plot_func(data_use,
                                               res_table_pairwise=visit1_mod_train_result$res_table_pairwise_pvalue,
                                               res_table=visit1_mod_train_result$res_table_ordinal,
                                               modules = fname, ncol = 4, show_legend = show_legend,
                                               label = "q.signif", show_significance =   show_significance,
                                               colors = c("#639A21", "#39828C",  "#6371AD","#BD7D31", "#9C3418"))

  return(list(FIG = figures_visit1_module, test_result = visit1_mod_train_result))
  
}

data_prepare_trajectory_inputfeature = function(x, feature_names, impacc_analysis_record, phase = "train",event_date_thr=42){
  factors_used = x
  if(phase == "train"){
    clinical_subsets_train =impacc_analysis_record$clin.train
  }else{
    clinical_subsets_train =impacc_analysis_record$clin.test
  }
  modules_scores_train = data.frame(factors_used)
  rownames(modules_scores_train) = clinical_subsets_train$event_id
  colnames(modules_scores_train) = feature_names
  
  #filtering
  used_ids_train = rownames(clinical_subsets_train[clinical_subsets_train$event_type %in% paste0("Visit ", 1:6) &
                                                     clinical_subsets_train$event_date <= event_date_thr,])
  modules_filtered_train =modules_scores_train[used_ids_train,]
  clinical_subsets_train = clinical_subsets_train[used_ids_train,]
  
  clinical_subsets_train$visits = factor(clinical_subsets_train$event_type, labels = 1:6)
  inputDF_train <- cbind(modules_filtered_train, clinical_subsets_train[c("trajectory_group", "respiratory_status_day14", "visits","event_date","enrollment_site", "participant_id","sex","admit_age", "discretized_admit_age_quantile")])
  inputDF_train <- inputDF_train %>% mutate(outcomeD14  = cut(respiratory_status_day14, 
                                                              breaks = c(1, 2, 4, 6, 7), 
                                                              include.lowest = TRUE)) %>%
    pivot_longer(cols = -c("event_date", "participant_id", "visits",
                           "enrollment_site", "trajectory_group",
                           "respiratory_status_day14", "outcomeD14",
                           "sex", "admit_age", "discretized_admit_age_quantile")) %>%
    filter(!is.na(value) & !is.na(trajectory_group))
  
  inputDF_train$trajectory_group <- as.factor(inputDF_train$trajectory_group)
  inputDF_train$name <- as.factor(inputDF_train$name)
  inputDF_train$sex <- factor(inputDF_train$sex, levels = c("Female", "Male"))
  inputDF_train$discretized_admit_age_quantile <- as.factor(inputDF_train$discretized_admit_age_quantile)
  inputDF_train$participant_id <- as.factor(inputDF_train$participant_id)
  inputDF_train$enrollment_site <- as.factor(inputDF_train$enrollment_site)
  return(inputDF_train)
}

data_prepare_trajectory = function(impacc_analysis_record, phase = "train", flipping = NULL){
  # Data preparation
  irank = impacc_analysis_record$rank_selected
  module_names = paste0("factor",1:irank)
  if(phase == "train"){
    factors_used = apply(impacc_analysis_record$mcia.factors.train,2,function(z) qqnorm(z,plot.it = F)$x)
    clinical_subsets_train =impacc_analysis_record$clin.train
  }else{
    factors_used = apply(impacc_analysis_record$mcia.factors.test,2,function(z) qqnorm(z,plot.it = F)$x)
    clinical_subsets_train =impacc_analysis_record$clin.test
  }
  factors_used = factors_used[,1:irank]
  if(is.null(flipping)){
    flipping = rep(F, irank)
  }
  for(j in 1:irank){
    if(flipping[j]){
      factors_used[,j] = factors_used[,j] * (-1)
    }
  }
  modules_scores_train = data.frame(factors_used)
  rownames(modules_scores_train) = clinical_subsets_train$event_id
  colnames(modules_scores_train) = module_names
  
  #filtering
  used_ids_train = rownames(clinical_subsets_train[clinical_subsets_train$event_type %in% paste0("Visit ", 1:6) &
                                                     clinical_subsets_train$event_date <= 42,])
  modules_filtered_train =modules_scores_train[used_ids_train,]
  clinical_subsets_train = clinical_subsets_train[used_ids_train,]
  
  clinical_subsets_train$visits = factor(clinical_subsets_train$event_type, labels = 1:6)
  inputDF_train <- cbind(modules_filtered_train, clinical_subsets_train[c("trajectory_group", "respiratory_status_day14", "visits","event_date","enrollment_site", "participant_id","sex","admit_age", "discretized_admit_age_quantile")])
  inputDF_train <- inputDF_train %>% mutate(outcomeD14  = cut(respiratory_status_day14, 
                                                              breaks = c(1, 2, 4, 6, 7), 
                                                              include.lowest = TRUE)) %>%
    pivot_longer(cols = -c("event_date", "participant_id", "visits",
                           "enrollment_site", "trajectory_group",
                           "respiratory_status_day14", "outcomeD14",
                           "sex", "admit_age", "discretized_admit_age_quantile")) %>%
    filter(!is.na(value) & !is.na(trajectory_group))
  
  inputDF_train$trajectory_group <- as.factor(inputDF_train$trajectory_group)
  inputDF_train$name <- as.factor(inputDF_train$name)
  inputDF_train$sex <- factor(inputDF_train$sex, levels = c("Female", "Male"))
  inputDF_train$discretized_admit_age_quantile <- as.factor(inputDF_train$discretized_admit_age_quantile)
  inputDF_train$participant_id <- as.factor(inputDF_train$participant_id)
  inputDF_train$enrollment_site <- as.factor(inputDF_train$enrollment_site)
  return(inputDF_train)
}


linear_effects_simple=function(feature_mat, clin_mat, adjust = FALSE){
  feature_mat0=feature_mat
  idx1 = which(clin_mat$event_type=="Visit 1")
  feature_mat1=feature_mat[idx1,,drop = F]
  logP_signed = matrix(NA, ncol = 3, nrow = ncol(feature_mat))
  colnames(logP_signed) = c("ordinal","mortality","slope5|4")
  rownames(logP_signed) = colnames(feature_mat)
  z = clin_mat$event_date- clin_mat$symptom_date
  options(na.action='na.pass')
  #W = model.matrix(~sqrt(z)+clin_mat$discretized_admit_age_quantile+clin_mat$sex+as.character(clin_mat$resp_status_v1)-1, )
  W = model.matrix(~sqrt(z)+clin_mat$discretized_admit_age_quantile+clin_mat$sex-1)
  W1 = W[idx1,]
  ww = apply(is.na(W),1,sum)==0
  if(adjust){
    ll_v1= (apply(is.na(W1),1,sum)==0)
    ll = (apply(is.na(W),1,sum)==0)
    for(j in 1:ncol(feature_mat0)){
      ll1 = which((!is.na(feature_mat0[,j])) & ll) 
      feature_mat0[-ll1,j]=NA
      feature_mat0[ll1,j] = lm(as.matrix(feature_mat0[ll1,j])~W[ll1,], na.action = 'na.pass')$residuals
      
      ll1 = which((!is.na(feature_mat0[idx1,j])) & ll_v1) 
      feature_mat1[-ll1,]=NA
      feature_mat1[ll1,j] = lm(as.matrix(feature_mat0[idx1[ll1],j])~W1[ll1,], na.action = 'na.pass')$residuals
    }
  }
  for(j in 1:ncol(feature_mat)){
    idx1 = which(clin_mat$event_type=="Visit 1")
    y = clin_mat$trajectory_group[idx1]
    tmp = cor.test(y,feature_mat1[,j], method = "spearman", use = "pairwise.complete")
    logP_signed[j,1] = -log(tmp$p.value,base = 10) * sign(tmp$estimate)
    
    y = clin_mat$trajectory_group
    idx1 = which(clin_mat$event_type=="Visit 1"&y%in%c(4,5) & (!is.na(feature_mat0[,j])))
    y = clin_mat$trajectory_group[idx1]
    tmp = cor.test(y, feature_mat0[idx1,j], method = "spearman", use = "pairwise.complete")
    tmpWilcox = wilcox.test(x=feature_mat0[idx1,j][y==4], y=feature_mat0[idx1,j][y==5])
    logP_signed[j,2] = -log(tmpWilcox$p.value,base = 10) * sign(tmp$estimate)
    
    y = clin_mat$trajectory_group
    idx1 = which(clin_mat$event_type%in%paste0("Visit ",c(1:6))&clin_mat$event_date<=30& y%in%c(4,5))
    t =sqrt(clin_mat$event_date[idx1])
    y = y[idx1]
    participant = clin_mat$participant_id[idx1]
    x = qqnorm(feature_mat[idx1,j],plot.it=F)$x
    delta = t
    delta[y==4] = -delta[y==4]
    dat_j = data.frame(endpoints=y, event_date = t, delta = delta, value = x, participant=participant)
    dat_j = dat_j[apply(is.na(dat_j),1,sum)==0,]
    formula_use = formula("value~ s(event_date, bs = 'cr')+delta+endpoints")
    fit <- try(gamm4::gamm4(formula_use, data = dat_j, random = ~(1|participant)))
    tmp = anova(fit$gam)
    logP_signed[j,3]=sign(tmp$p.coeff[2])*(-log( tmp$pTerms.pv[1],base=10))
  }
  return(logP_signed)
}


simple_pathwaytest_wrapper = function(composite_mat, base_obj, phase = "others", adjust = F){
  if(phase == "train"){
    clin_dat = base_obj$clin.train
  }else{
    clin_dat =base_obj$clin.test
  }
  W =model.matrix(~as.character(clin_dat$trajectory_group)+clin_dat$event_type+clin_dat$discretized_admit_age_quantile+clin_dat$sex-1)
  #id1 = !is.na(composite_mat[,1])
  #clin_mat = clin_dat[id1,,drop = F]
  #composite_mat=composite_mat[id1,]
  clin_mat = clin_dat
  fitted = linear_effects_simple(feature_mat=composite_mat, clin_mat=clin_mat,adjust = adjust)
  return(fitted)
}

mgcv_global_coef = function(inputDF, age_sex = T, endpoints = "trajectory_group"){
  features_names = sort(unique(inputDF$name))
  
  res_table = data.frame(matrix(NA, ncol = 1, nrow = length(features_names)))
  colnames(res_table) = c("p.intercept.direction")
  
  for(i in 1:length(features_names)){
    print(i)
    feature_name = features_names[i]
    tmp_input = inputDF[inputDF$name==feature_name,]
    tmp_input[[endpoints]] =ordered(as.factor(tmp_input[[endpoints]]))
    tmp_input$participant_id = as.factor(tmp_input$participant_id)
    if(age_sex){
      formula_use = formula(paste0("value~ s(event_date, bs = 'cr')+
              s(event_date, bs = 'cr', by =", endpoints, ")+", endpoints,
                                   "+ sex + discretized_admit_age_quantile"))
    }else{
      formula_use = formula(paste0("value~ s(event_date, bs = 'cr')+
              s(event_date, bs = 'cr', by =", endpoints, ")+",endpoints))
    }
    fit <- try(gamm4::gamm4(formula_use, data = tmp_input, random = ~(1|enrollment_site/participant_id)))
    if(class(fit)!="try-error"){
      res_table[i,1] = ifelse(fit$gam$coefficients["trajectory_group.L"] > 0, "Severe", "Mild")
    }
  }
  rownames(res_table) = features_names
  return(res_table)
}



model_loop_coef <- function(inputDF, age_sex =TRUE, knots = c(1, 4, 7, 14, 21),
                            endpoint = "trajectory_group", old_p_corrections = FALSE){
  inputDF$participant_id = as.factor(inputDF$participant_id)
  ## Because of the pairwise function requires there to be a column for "name"
  ## the rownames are copied to a new column here, ultimately this column is removed
  inputDF_s = inputDF
  inputDF_s$event_date_transformed = inputDF_s$event_date
  s_tmp <- mgcv_global_coef(inputDF_s, endpoint, age_sex = age_sex)
  
  lmDF <- inputDF %>%
    group_by(name) %>% 
    do(p.slope.direction = {
      fit <- try(lme(fixed = as.formula(paste0("value ~ event_date * ", endpoint,
                                               "+sex+discretized_admit_age_quantile")), data = .,
                     random =  ~1|enrollment_site/participant_id))
      if(class(fit)[1]!="try-error"){
        result = ifelse(fit$coefficients$fixed["event_date:trajectory_group5"] > 0, "Severe", "Mild")
      } else {
        result = NA
      }
    }) %>% ungroup() %>%
    mutate(p.slope.direction = unlist(p.slope.direction),
           p.intercept.direction  = s_tmp$p.intercept.direction) %>% as.data.frame()
  
  message("\n Initial calculation of p.slope and p.intercept complete. Beginning pairwise comparisons.")
  #### Pairwise comparisons #### 
  ## Loop the pairwise comparisons ##
  for(i in seq_along(rownames(lmDF))){ 
    tmp.out = matrix(NA,ncol = ncol(combn(unique(inputDF[[endpoint]]), 2)), nrow =2)
    tmp.out <- pairwise_coef(inputDF[inputDF$name == lmDF$name[i],],
                             variable = endpoint, fixedKnots = F)
    if(i==1){
      values.out <- c(tmp.out[1,], tmp.out[2,])
      names(values.out) <- c(paste0(rownames(tmp.out)[1],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])),
                             paste0(rownames(tmp.out)[2],"_", colnames(tmp.out[,1:ncol(tmp.out), drop=FALSE])))
    } else {
      values.out <- rbind(values.out, c(tmp.out[1,], tmp.out[2,]))
    }
  }
  
  lmDF_pairwise <- as.data.frame(cbind(lmDF, values.out))
  rownames(lmDF_pairwise) <- lmDF_pairwise$name
  lmDF_pairwise <- lmDF_pairwise[,colnames(lmDF_pairwise) != "name"]
  return(lmDF_pairwise)
}



pairwise_coef <- function(data, variable = "trajectory_group", fixedKnots = F, sex_age = T,
                          knots = c(1,4,7,14,21)) {
  for(i in seq_along(unique(data$name))){
    ## retrieve comparisons
    active_data <- data[data$name == unique(data$name)[i], ]
    pairwise_comparisons1 <- combn(levels(active_data[[variable]]), 2)
    
    levels_encode = sort(unique(active_data[[variable]]))
    pairwise_comparisons2 <- combn(levels_encode, 2)
    
    compout <- matrix(NA, ncol(pairwise_comparisons1), nrow = 2)
    rownames(compout) = c("p.slope", "p.intercept")
    colnames(compout) <- paste0(pairwise_comparisons1[1, ], "v", pairwise_comparisons1[2, ])
    
    for(j in 1:ncol(pairwise_comparisons1)){
      data_tmp = active_data[active_data[[variable]] == pairwise_comparisons1[1, j] |
                               active_data[[variable]] == pairwise_comparisons1[2, j],]
      formula_use= formula(paste0("value ~ event_date * ",variable, "+sex +discretized_admit_age_quantile"))
      lmfit= lme(fixed = formula_use, random = ~1|enrollment_site/participant_id, data = data_tmp)
      compout[1,j] = ifelse(lmfit$coefficients$fixed[9] > 0, "Severe", "Mild")
      
      tmp_data = active_data[active_data[[variable]]==pairwise_comparisons2[1,j] | 
                               active_data[[variable]]==pairwise_comparisons2[2,j],]
      tmp_data$participant_id = as.factor(tmp_data$participant_id )
      tmp_data$enrollment_site = as.factor(tmp_data$enrollment_site )
      tmp_data[[variable]] = ordered(as.factor(as.character(tmp_data[[variable]])), 
                                     levels = c(as.character(pairwise_comparisons2[1,j]), 
                                                as.character(pairwise_comparisons2[2,j])))
      formula_use_s = formula(paste0("value~ s(event_date, bs = 'cr') + 
                                     s(event_date, bs = 'cr', by =", variable, ")+", variable,
                                     "+ sex + discretized_admit_age_quantile"))
      sfit <- try(gamm4::gamm4(formula_use_s, data = tmp_data, random = ~(1|enrollment_site/participant_id)))
      compout[2,j] = ifelse(sfit$gam$coefficients["trajectory_group.L"] > 0, "Severe", "Mild")
    }
    
    if(length(unique(data$name)) > 1){
      if(i ==1){
        compout_final <- list(compout)
        names(compout_final)[i] <- unique(data$name)[i]
      } else {
        compout_final <- c(compout_final, list(compout))
        names(compout_final)[i] <- unique(data$name)[i]
      }
    } else {
      compout_final <- compout
    } 
  }
  return(compout_final)
}


#A1-A4: table of interest, B1-B4: table for filtering
#test="one-sided": if the direction of A - B does not match: p-val = 1. 
module_test_table_creation_visit1_helper = function(A1, A2, A3, A4, B1, B2, B3, B4, 
                                                    global_qval = 0.05, pairwise_qval = 0.05,
                                                    test = "two-sided",
                                                    columns_pairwise_interest = columns_pairwise_interest, columns_pairwise_redefine = columns_pairwise_redefine){
  #ordinal
  idx_ordinal = which(B1$qval <=global_qval )
  res_table = data.frame("Test" = rep("ordinal",length(idx_ordinal)), 
                         "Factor" = paste0("Factor",idx_ordinal),
                         "pval" = A1$pval[idx_ordinal], 
                         "adj_pval" = A1$qval[idx_ordinal],
                         "direction" = A1$direction[idx_ordinal])
  for(j in 1:length(columns_pairwise_interest)){
    col0 = columns_pairwise_interest[j]
    idx = which(B3[[col0]]<=pairwise_qval)
    if(length(idx)>0){
      tmp = data.frame("Test" = rep(columns_pairwise_redefine[j],length(idx)), 
                       "Factor" = paste0("Factor", idx),
                       "pval" = A2[[col0]][idx], 
                       "adj_pval" = A3[[col0]][idx],
                       "direction" = A4[[col0]][idx])
      if(test == "one-sided"){
        B4direction = B4[[col0]][idx]
        for(l in 1:nrow(tmp)){
          if(B4direction[l]==tmp$direction[l]){
            tmp$pval[l] = tmp$pval[l]/2
            tmp$adj_pval[l] = tmp$adj_pval[l]/2
          }else{
            tmp$pval[l] = 1.0
            tmp$qval[l] = 1.0
          }
        }
      }
      res_table=rbind(res_table,tmp)
    }
    
  }
  return(res_table)
}

module_test_table_creation_lonitudinaltest_helper = function(smooth_pvalue,  direction,
                                                             smooth_pvalue_filter,  direction_filter,
                                                             global_qval = 0.05,
                                                             pairwise_qval = 0.05,
                                                             test = "two-sided"){
  global_smooth_pvalue =smooth_pvalue[,c("p.slope","p.intercept","adjp.slope","adjp.intercept")]
  global_smooth_pvalue[,3] = p.adjust(global_smooth_pvalue[,1])
  global_smooth_pvalue[,4] = p.adjust(global_smooth_pvalue[,2])
  
  global_smooth_pvalue_filter =smooth_pvalue_filter[,c("p.slope","p.intercept","adjp.slope","adjp.intercept")]
  global_smooth_pvalue_filter[,3] = p.adjust(global_smooth_pvalue_filter[,1])
  global_smooth_pvalue_filter[,4] = p.adjust(global_smooth_pvalue_filter[,2])
  #ordinal
  idx_ordinal_slope = which(global_smooth_pvalue_filter$adjp.slope <=global_qval )
  idx_ordinal_intercept = which(global_smooth_pvalue_filter$adjp.intercept <=global_qval )
  res_table = data.frame("Test" = c(rep("slope_global",length(idx_ordinal_slope)),
                                    rep("intercept_global",length(idx_ordinal_intercept))),
                         "Factor" = c(paste0("Factor",idx_ordinal_slope),
                                      paste0("Factor", idx_ordinal_intercept)),
                         "pval" = c(global_smooth_pvalue$p.slope[idx_ordinal_slope], global_smooth_pvalue$p.intercept[idx_ordinal_intercept]),
                         "adj_pval" = c(global_smooth_pvalue$adjp.slope[idx_ordinal_slope], global_smooth_pvalue$adjp.intercept[idx_ordinal_intercept])
  )
  
  res_table$direction = NA
  tmp1 = intersect(colnames(smooth_pvalue), colnames(direction))
  direction = direction[tmp1]
  tmp2 = sapply(strsplit(tmp1, "_"), function(z) paste0(z[1], ".adj","_", z[2]))
  smooth_pvalue_selected = smooth_pvalue[tmp1]
  smooth_qvalue_selected = smooth_pvalue[tmp2]
  direction_filter = direction_filter[tmp1]
  smooth_pvalue_filter_selected = smooth_pvalue_filter[tmp1]
  smooth_qvalue_filter_selected = smooth_pvalue_filter[tmp2]  
  
  for(j in 1:ncol(smooth_pvalue_filter_selected)){
    col0 = smooth_qvalue_filter_selected[,j]
    col0name = colnames(smooth_qvalue_filter_selected)[j]
    idx = which(col0<=pairwise_qval)
    if(length(idx)>0){
      tmp = data.frame("Test" = rep(col0name,length(idx)), 
                       "Factor" = paste0("Factor", idx),
                       "pval" = smooth_pvalue_selected[,j][idx], 
                       "adj_pval" = smooth_qvalue_selected[,j][idx],
                       "direction" = direction[,j][idx])
      if(test == "one-sided"){
        direction1 =direction_filter[idx,j]
        for(l in 1:nrow(tmp)){
          if(direction1[l]==tmp$direction[l]){
            tmp$pval[l] = tmp$pval[l]/2
            tmp$adj_pval[l] = tmp$adj_pval[l]/2
          }else{
            tmp$pval[l] = 1.0
            tmp$adj_pval[l] = 1.0
          }
        }
      }
      res_table = rbind(res_table,tmp)
    }
  }
  return(res_table)
}


composite_mat_subsetting = function(composite_score_mat_list,selected_pathways){
  pathway_vec = c()
  omic_vec = c()
  composiste_mat = NULL
  for(d in 1:length(composite_score_mat_list)){
    tmp_X = composite_score_mat_list[[d]]
    intersect_pathways = intersect(selected_pathways,colnames(tmp_X ))
    if(length(intersect_pathways)>0){
      tmp_X = tmp_X[,intersect_pathways]
      pathway_vec=c(pathway_vec,intersect_pathways)
      omic_vec = c(omic_vec,rep(names(composite_score_mat_list)[d], length(intersect_pathways)))
      if(is.null(composiste_mat)){
        composiste_mat  = as.matrix(tmp_X)
      }else{
        composiste_mat  = cbind( composiste_mat,as.matrix(tmp_X))
      }
    }
  }
  return(list(composiste_mat =composiste_mat , omic_vec=omic_vec,
              pathway_vec = pathway_vec))
}

pathway_selection = function(signed_log10_pathway, thr = 0.05){
  selected_pathways = NULL
  thr = -log(thr, base = 10)
  idx1 = which(apply(signed_log10_pathway>=thr,1,sum)>=2)
  idx2 = which(apply(signed_log10_pathway< -thr,1,sum)>=2)
  idx = c(idx1, idx2)
  if(length(idx)>0){
    idx = sort(unique(idx))
    selected_pathways=signed_log10_pathway[idx,,drop = F]
  }
  return(selected_pathways)
}

TG45_selection = function(composite_mat,base_obj, sign, thr = 0.05, phase = "train"){
  tmp = simple_pathwaytest_wrapper(composite_mat = composite_mat,
                                   base_obj = base_obj,
                                   phase = phase, adjust = F)
  tmp1 = simple_pathwaytest_wrapper(composite_mat =composite_mat,
                                    base_obj = base_obj,
                                    phase = phase, adjust = T)
  ###stability selection
  signed_log10_pathway = cbind(tmp[,c(2,3)], tmp1[,2])
  selected1 = (sign*signed_log10_pathway[,1])>abs(log(thr, base = 10))
  selected2 = rep(F, length(selected1))
  tmp2 = pathway_selection(signed_log10_pathway, thr = thr)
  if(length(tmp2)>0){
    selected2[rownames(tmp1)%in%rownames(tmp2)] = T
  }
  w = data.frame(cbind(tmp[,c(1, 2,3)], tmp1[,2]))
  w$selected = (selected1 & selected2)
  colnames(w) = c("severity", "mortality", "slope_adj",   "mortality_adj","selected")
  return(w)
  
  
}

time_robust_pval_construction = function(pval_mat, sign_mat){
  z = rep(NA, nrow(pval_mat))
  for(i in 1:nrow(pval_mat)){
    z[i] = pval_mat[i,1]
    tmp_z = pval_mat[i,-1]
    tmp_z[sign_mat[i,-1]!=sign_mat[i,1]] = 1.0
    tmp_z = min(tmp_z)
    tmp_z = 1.0 - (1.0-tmp_z)^2
    z[i] = sqrt(z[i]*max(z[i], tmp_z))
  }
  return(z)
}