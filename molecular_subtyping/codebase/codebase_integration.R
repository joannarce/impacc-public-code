packages <-
  c("pals",
    "stringr",
    "ggpubr",
    "cowplot",
    "RColorBrewer",
    "rlang",
    "biomaRt",
    "limma",
    "ComplexHeatmap",
    "GSEABase",
    "edgeR",
    "pbmcapply",
    "org.Hs.eg.db",
    "ReactomePA",
    "igraph",
    "tidyverse",
    "metaboliteIDmapping")
for (n in 1:length(packages)) {
  suppressMessages(library(packages[n], character.only = TRUE))
}
#################################################
# data process II
##################################################
#Maintainer: Jeremy
#
#
get.event_type.from.part_event <- function(g_id){
  dict <- distinct(clinical_data, event_type, part.event) %>% filter(part.event==g_id)
  if(nrow(dict) > 0){
    return(dict$event_type)
  } else {
    return(NA)
  }
}


get.part_event.from.sample_id <- function(g_id, clinical_data){
  dict <- dplyr::distinct(clinical_data, sample_id, part.event) %>% filter(sample_id == g_id)
  if(nrow(dict) > 0){
    return(dict$part.event)
  } else {
    return(NA)
  }
}

get.event_type.from.sample_id <- function(g_id, clinical_data){
  dict <- distinct(clinical_data, sample_id, event_type) %>% filter(sample_id == g_id)
  if(nrow(dict) > 0){
    return(dict$event_type)
  } else {
    return(NA)
  }
}

get.phase.from.sample_id <- function(g_id, clinical_data){
  dict <- distinct(clinical_data, sample_id, phase) %>% filter(sample_id == g_id)
  if(nrow(dict) > 0){
    return(dict$phase)
  } else {
    return(NA)
  }
}

get.dates.average <- function(clinical_data, clinical_used){
  n = nrow(clinical_used)
  event_dates = rep(NA, n)
  sym_dates = rep(NA, n)
  for(i in 1:n){
    v = rownames(clinical_used)[i]
    ll= which(clinical_data$part.event == v)
    v1 = clinical_data$event_date[ll]
    event_dates[i] = mean(v1, na.rm = T)
    v2 = clinical_data$symptom_date[ll]
    if(sum(!is.na(v2))>0){
      sym_dates[i] = event_dates[i]-mean(v2, na.rm = T)
    }
  }
  clinical_used$event_date = event_dates
  clinical_used$symptom_dates = sym_dates
  return(clinical_used)
}


data_prepare_DR = function(omics_used, data_env, clinical_data, rna_keep = 5000){
  #clinical data preparation
  clinical_df <- clinical_data %>% distinct(part.event, sex, ethnicity, admit_age,  discretized_admit_age_quantile,
                                            respiratory_status,respiratory_status_day14, respiratory_status_day28,
                                            trajectory_group, race, enrollment_site,
                                            discharge_dates, death_date)
  # make sex one hot encoded:
  clinical_df$sex_binomial <- ifelse(clinical_df$sex == "Male", 1, 0)
  # make hispanic_or_latino one hot encoded:
  clinical_df$hispanic_or_latino <- ifelse(clinical_df$ethnicity == "Hispanic or Latino", 1, 0)
  # make race one hot encoded
  clinical_df <- clinical_df %>% mutate(value = 1)  %>% spread(race, value,  fill = 0)
  # Sort:
  clinical_df <- clinical_df %>% dplyr::select(starts_with("respiratory"), trajectory_group,sex_binomial, everything())
  clinical_df <- clinical_df %>% mutate(respiratory_status_day14_endpoint = sapply(clinical_df$respiratory_status_day14, function(r){
    if(is.na(r)){
      return(NA)
    } else if(r < 3){
      return(0)
    } else if (r < 5){
      return(1)
    } else if (r < 7){
      return(2)
    } else{
      return(3)
    }
  })
  ) %>% mutate(respiratory_status_day28_endpoint = sapply(clinical_df$respiratory_status_day28, function(r){
    if(is.na(r)){
      return(NA)
    } else if(r < 3){
      return(0)
    } else if (r < 5){
      return(1)
    } else if (r < 7){
      return(2)
    } else{
      return(3)
    }
  })
  ) %>%  dplyr::select(part.event,respiratory_status, respiratory_status_day14, respiratory_status_day14_endpoint,
                       respiratory_status_day28, respiratory_status_day28_endpoint, trajectory_group,
                       everything())
  # turn into matrix
  clinical_matrix <- as.matrix(dplyr::select(clinical_df, -part.event, -sex, -ethnicity))
  rownames(clinical_matrix) <- clinical_df$part.event
  
  # assay data
  datasets.mofa = list()
  #### Identify outlying samples and filter out them.
  for( i in 1:length(omics_used) ) {
    o_name = omics_used[i]
    count_df <- data_env[[omics_used[i]]]
    datasets.mofa[[o_name]] = count_df
  }
  datasets <- datasets.mofa
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
    # a = qnorm((c(1:nrow(dataset))-0.5)/nrow(dataset))
    # for(j in 1:ncol(dataset)){
    #   b = rank(dataset[,j])
    #   dataset[,j] = a[b]
    # }
    # Change sample name to part_event
    tmp =  sapply(rownames(dataset), function(z) get.part_event.from.sample_id(z, clinical_data=clinical_data))
    dataset = dataset[!is.na(tmp),]
    rownames(dataset) <- tmp[!is.na(tmp)]
    datasets[[i]] <- dataset
  }
  # Combine:
  part.event.tab = c()
  for(i in 1:length(datasets)){
    part.event.tab = c(part.event.tab, rownames(datasets[[i]]))
  }
  part.event.tab <- table( part.event.tab)
  
  
  #participants.to.use <- names(part.event.tab)[which(part.event.tab == 4)]
  participants.to.use = names(part.event.tab)
  
  #unique part.event
  clinical_used = clinical_df
  rownames(clinical_used) =  clinical_df$part.event
  clinical_used = clinical_used[participants.to.use,]
  datasets_combined = list()
  for(i in 1:length(datasets)){
    datasets_combined[[i]] = matrix(NA, ncol = ncol(datasets[[i]]), nrow= length(participants.to.use))
    rownames(datasets_combined[[i]]) = participants.to.use
    datasets_combined[[i]] = data.frame(datasets_combined[[i]])
    rownames(datasets_combined[[i]]) = participants.to.use
    colnames(datasets_combined[[i]]) = colnames(datasets[[i]])
    datasets_combined[[i]][rownames(datasets[[i]]),] = datasets[[i]]
  }
  #unique column names
  for(i in 1:length(datasets_combined)){
    dataset <- datasets_combined[[i]]
    colnames(dataset) <- paste0(omics_used[i], "_", colnames(dataset))
    datasets_combined[[i]] <- dataset
  }
  missed = matrix(0, ncol = length(datasets_combined), nrow = nrow(datasets_combined[[1]]))
  for(d in 1:length(datasets_combined)){
    a = apply(is.na(datasets_combined[[d]]),1,sum)
    missed[a>0,d]=1
  }
  tmp = as.integer(sapply(strsplit(clinical_used$part.event,"[.]"),function(z) z[3]))
  missing_assays = apply(missed[tmp==1,],1,sum)
  clinical_used=get.dates.average(clinical_data, clinical_used)
  print(table(  missing_assays))
  print(sapply(datasets, function(z) dim(z)))
  sample_metadata <- data.frame(
    traj = as.factor(clinical_used$trajectory_group),
    day14 = clinical_used$respiratory_status_day14_endpoint,
    event_date = clinical_used$event_date,
    symptom_date = clinical_used$symptom_date,
    vists = as.factor(as.integer(sapply(strsplit(clinical_used$part.event,"[.]"),function(z) z[[3]]))),
    escalation = sapply(strsplit(clinical_used$part.event,"[.]"),function(z) z[[2]]),
    participant =(as.character(sapply(strsplit(clinical_used$part.event,"[.]"),function(z) z[[1]]))),
    enrollment_site = clinical_used$enrollment_site,
    death_date = clinical_used$death_date,
    discharge_dates = clinical_used$discharge_dates,
    sex = clinical_used$sex,
    age = clinical_used$admit_age,
    age_discrete = clinical_used$discretized_admit_age_quantile
  )
  names(datasets_combined) =  omics_used
  return(list(datasets_combined = datasets_combined, 
              sample_metadata = sample_metadata, omics_used = omics_used))
}

mofa_run = function(datasets_combined, num_factors = 50,
                    seed = 42, maxiter = 1000, drop_factor_threshold = -1, 
                    convergence_mode = "fast", startELBO = 2, stochastic=F,
                    mofa_file = NULL){
  # Run MOFA+:
  x.mofa <- list()
  for(d in 1:length(datasets_combined)){
    x.mofa[[d]] = t(datasets_combined[[d]])
  }
  names(x.mofa) = names(datasets_combined)
  MOFAobject <- create_mofa(x.mofa)
  # MOFA+ specific parameters:
  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  train_opts <- get_default_training_options(MOFAobject)
  model_opts$num_factors =num_factors
  train_opts$convergence_mode = convergence_mode
  train_opts$startELBO = startELBO
  train_opts$maxiter = maxiter
  train_opts$stochastic = stochastic
  train_opts$seed = seed
  train_opts$drop_factor_threshold =drop_factor_threshold
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  MOFAobject.trained <- run_mofa(MOFAobject,use_basilisk = TRUE, save_data = F)
  return(MOFAobject.trained)
}

##########################################
# Dimension reduction (only if extra process needed)
##########################################


#################################################
# Association test
# TODO: ordering the figures based on mod number in the plot
##################################################
#Maintainer: Leqi Xu
#' modules: NULL=plot all; 
#' label: q.signif or  q.val
association_plot = function(data_use,  res_table_pairwise, res_table,
                            modules = NULL, ncol = 4, 
                            label = "q.signif",
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
  #data_plot <- data_use%>%gather(key = 'module', value = 'expression', -endpoints, -control, -sites)
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
  names(module_pvalue) <- ordinal_res_table_formatted$module
  module_psig <- ordinal_res_table_formatted[,3]
  names(module_psig) <-ordinal_res_table_formatted$module
  if(!is.null(modules)){
    data_plot_01 <- data_plot[which(data_plot$module %in% modules),]
    res_table_pairwise_all01 <- pairwise_res_table_formatted[which(pairwise_res_table_formatted$module%in%modules),]
    module_pvalue = module_pvalue[modules]
    module_psig = module_psig[modules]
  }else{
    data_plot_01 <- data_plot
    res_table_pairwise_all01 <- pairwise_res_table_formatted
  }
  data_plot_01$endpoints = factor(data_plot_01$endpoints, labels = c("TG1", "TG2", "TG3", "TG4", "TG5"))
  if(length(row.names(res_table_pairwise_all01[res_table_pairwise_all01$q.signif!="ns",])) == 0){
    mod01_box0 <- ggboxplot(data_plot_01, x = "endpoints", y = "expression",
                            color = "endpoints", palette =colors,
                            add = "jitter") + 
      facet_wrap(~ module,labeller = labeller(module = module_pvalue), ncol = min(ncol, length(unique(data_plot_01$module))))+
      theme_bw()
  } else{
    mod01_box0 <- ggboxplot(data_plot_01, x = "endpoints", y = "expression",
                            color = "endpoints", palette =colors,
                            add = "jitter") +
    stat_pvalue_manual(data = res_table_pairwise_all01[res_table_pairwise_all01$q.signif!="ns",], 
                       y.position = quantile(data_plot_01$expression,0.99), step.increase = 0.05, step.group.by = "module", bracket.size = 0.2,  remove.bracket  = F,
                       label = "q.signif") +
    facet_wrap(~ module,labeller = labeller(module = module_pvalue), ncol = min(ncol, length(unique(data_plot_01$module))))+
    theme_bw()
  }
  return(mod01_box0)
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
                                     breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
                                     labels = c("*****","****", "***", "**", "*", "ns"))
  pairwise_table_all$group1 = paste0("TG",pairwise_table_all$group1)
  pairwise_table_all$group2 = paste0("TG",pairwise_table_all$group2)
  return(pairwise_table_all)
}

ordinal_format <- function(ordinal_table){
  ordinal_table$module <- rownames(ordinal_table)
  ordinal_table$qval <- signif(ordinal_table$qval,3)
  ordinal_table$q.signif <- cut(ordinal_table$qval, 
                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
                                labels = c("*****","****", "***", "**", "*", "ns"))
  ordinal_table$modp <- paste0(ordinal_table$module, ": ", ordinal_table$qval)
  ordinal_table$modsig <- paste0(ordinal_table$module, ": ", ordinal_table$q.signif) 
  
  return(ordinal_table)
}
#################################################
# Prediction
##################################################
#Maintainer: Anna, Casey,Jeremy Gygi

ordinalreg_l1selection = function(dat.tr, module_names, response_name = "traj", ordering = NULL, nfolds = NULL){
  tmp <- ordinalNet::ordinalNet(as.matrix(dat.tr[module_names]), y =dat.tr[[response_name]] , family = "cumulative",link = "logit", parallelTerms = TRUE, nonparallelTerms = FALSE, standardize = F)
  idxs = coef(tmp, criteria = "aic")[-(1:(length(unique(dat.tr$traj))-1))]
  aa = paste(names(idxs),collapse="+")
  tmp <- MASS::polr(formula(paste0(response_name,"~",aa)), data = dat.tr, Hess=TRUE)
  return(tmp)
}



prediction_hierachical = function(dat.te, module_names,layer1model,  layer2model, layer2_level = list(c(1),c(2,3), c(4,5)),
                                  ordinal = T){
  if(length(layer2_level)>2){
    if(ordinal){
      pred.layer1 = predict(layer1model, dat.te, type = "prob")
    }else{
      pred.layer1 = predict(layer1model, as.matrix(dat.te[module_names]), type = "response", s = "lambda.min")
    }
  }else{
    pred.layer1 = predict(layer1model, as.matrix(dat.te[module_names]), type = "response", s = "lambda.min")
    pred.layer1 = cbind(1-pred.layer1, pred.layer1)
  }
 
  ii = sapply(layer2_level, length)
  ii1 = which(ii>1)
  K2 = length((ii1))
  if(K2>0){
    pred.layer2s = matrix(0, ncol = K2, nrow = nrow(dat.te))
    for(k in 1:K2){
      pred.layer2s[,k] = predict(layer2model[[k]], as.matrix(dat.te[module_names]), type = "response", s = "lambda.min")
    }
    pred.final =  matrix(0, nrow = nrow(pred.layer1), ncol = length(unlist(layer2_level)))
    for(k in 1:length(layer2_level)){
      pred.final[,layer2_level[[k]]] = pred.layer1[,k]
    }
    for(k in 1:K2){
      pred.final[,layer2_level[[ii1[k]]][1]] = pred.layer1[,ii1[k]]*(1-pred.layer2s[,k])
      pred.final[,layer2_level[[ii1[k]]][2]] = pred.layer1[,ii1[k]]*pred.layer2s[,k]
    }
  }else{
    pred.final = pred.layer1
  }


  return(pred.final)
}

hireg = function(dat.tr, module_names, response_name = "traj", nfolds = 20,
                 layer2_level = list(c(1),c(2,3), c(4,5)), seed = 2022, 
                 update.response = TRUE, ordinal= T, alpha = 1){
  #create layer 1 ,2 response
  dat_sub.tr = list()
  ii = sapply(layer2_level, length)
  ii1 = which(ii>1)
  K2 = length((ii1))
  if(K2>0){
    for(k in 1:K2){
      dat_sub.tr[[k]] = dat.tr[dat.tr[[response_name]] %in%layer2_level[[ii1[k]]],]
      dat_sub.tr[[k]][[response_name]] = as.integer( dat_sub.tr[[k]][[response_name]])
      aa = dat_sub.tr[[k]][[response_name]]
      dat_sub.tr[[k]][[response_name]][aa == layer2_level[[ii1[k]]][1]] = 0
      dat_sub.tr[[k]][[response_name]][aa == layer2_level[[ii1[k]]][2]] = 1
      dat_sub.tr[[k]][[response_name]] = as.integer( dat_sub.tr[[k]][[response_name]])
    }
  }

  dat.tr.mutated = dat.tr
  if(update.response){
    if(K2>0){
      for(k in 1:K2){
        dat.tr.mutated[[response_name]][dat.tr[[response_name]]%in%layer2_level[[ii1[k]]]] = ii1[k]
      }
    }

  }
  dat.tr.mutated[[response_name]] = factor(as.integer(dat.tr.mutated[[response_name]]))
  # #prediction layer 1
  if(length(layer2_level)>2){
    if(ordinal){
      model_ordinal=ordinalreg_l1selection(dat.tr = dat.tr.mutated, module_names = module_names, response_name = response_name)
    }else{
      model_ordinal = cv.glmnet(x = as.matrix(dat.tr.mutated[module_names]), y =dat.tr.mutated[[response_name]], family = "multinomial", nfolds = nfolds, seed = seed, alpha = alpha, standardize = F)
    }
  }else{
    model_ordinal=cv.glmnet(x = as.matrix(dat.tr.mutated[module_names]), y =dat.tr.mutated[[response_name]], family = "binomial", nfolds = nfolds, seed = seed, alpha = alpha, standardize = F)
  }
  #perform a layer 2 task
  model_layer2 = list()
  if(K2>0){
    for(k in 1:K2){
      model_layer2[[k]] = cv.glmnet(x = as.matrix(dat_sub.tr[[k]][module_names]), y =dat_sub.tr[[k]][[response_name]], family = "binomial", nfolds = nfolds, seed = seed, alpha = alpha, standardize = F)
    }
  }

  
  return(list(layer1model =  model_ordinal, layer2model=model_layer2, layer2_level = layer2_level))
}


ordinalreg_l1selectionV2 = function(dat.tr, module_names, response_name = "traj",  nfolds = NULL){
  model <- ordinalNet::ordinalNet(as.matrix(dat.tr[module_names]), y =dat.tr[[response_name]] , family = "cumulative",link = "logit", 
                                  parallelTerms = TRUE, nonparallelTerms = FALSE, standardize = F, nLambda = 100, lambdaMinRatio=1e-4)
  
  return(model)
}

hiregV2 = function(dat.tr, module_names, response_name = "traj", nfolds = 20,
                   layer2_level = list(c(1),c(2), c(3),c(4,5)), seed = 2022, 
                   update.response = TRUE, ordinal= T, alpha = 1){
  #global ordinal model
  set.seed(seed)
  model_ordinal=ordinalreg_l1selectionV2(dat.tr = dat.tr, module_names = module_names, response_name = response_name, nfolds = nfolds)
  #create layer 1 ,2 response
  dat_sub.tr = list()
  ii = sapply(layer2_level, length)
  ii1 = which(ii>1)
  K2 = length((ii1))
  if(K2>0){
    for(k in 1:K2){
      dat_sub.tr[[k]] = dat.tr[dat.tr[[response_name]] %in%layer2_level[[ii1[k]]],]
      dat_sub.tr[[k]][[response_name]] = as.integer( dat_sub.tr[[k]][[response_name]])
      aa = dat_sub.tr[[k]][[response_name]]
      dat_sub.tr[[k]][[response_name]][aa == layer2_level[[ii1[k]]][1]] = 0
      dat_sub.tr[[k]][[response_name]][aa == layer2_level[[ii1[k]]][2]] = 1
      dat_sub.tr[[k]][[response_name]] = as.integer( dat_sub.tr[[k]][[response_name]])
    }
  }
  #perform a layer 2 task
  model_layer2 = list()
  if(K2>0){
    for(k in 1:K2){
      model_layer2[[k]] = cv.glmnet(x = as.matrix(dat_sub.tr[[k]][module_names]), y =dat_sub.tr[[k]][[response_name]], family = "binomial", nfolds = nfolds, seed = seed, alpha = alpha, standardize = F)
    }
  }
  return(list(layer1model =  model_ordinal, layer2model=model_layer2, layer2_level = layer2_level))
}

prediction_hierachicalV2 = function(dat.te, module_names,layer1model,  layer2model, layer2_level = list(c(1),c(2,3), c(4,5))){
  pred.layer1 = predict(layer1model, as.matrix(dat.te[module_names]), type = "response", criteria = "aic")
  #merge prediction and aggregate
  ii = sapply(layer2_level, length)
  ii1 = which(ii>1)
  K2 = length((ii1))
  if(K2>0){
    pred.layer2s = matrix(0, ncol = K2, nrow = nrow(dat.te))
    for(k in 1:K2){
      pred.layer2s[,k] = predict(layer2model[[k]], as.matrix(dat.te[module_names]), type = "response", s = "lambda.min")
    }
    pred.final =  matrix(0, nrow = nrow(pred.layer1), ncol = length(unlist(layer2_level)))
    for(k in 1:length(layer2_level)){
      if(length(layer2_level[[k]])==1){
        pred.final[,layer2_level[[k]]] = pred.layer1[,layer2_level[[k]]]
      }else{
        pred.final[,layer2_level[[k]]] = apply(pred.layer1[,layer2_level[[k]],drop = F],1,sum)
      }
    }
    for(k in 1:K2){
      pred.final[,layer2_level[[ii1[k]]][1]] = pred.final[,layer2_level[[ii1[k]]][1]]*(1-pred.layer2s[,k])
      pred.final[,layer2_level[[ii1[k]]][2]] = pred.final[,layer2_level[[ii1[k]]][2]]*pred.layer2s[,k]
    }
  }else{
    pred.final = pred.layer1
  }
  return(pred.final)
}

# #prediction_hierachicalV2 do not perform any grouping, grouping is done at a post-training fasion later
# prediction_hierachicalV3 = function(dat.te, module_names, model_ordinal ,   model_logit_list){
#   pred1 = predict(model_ordinal, dat.te, type = "probs")
#   pred2 = predict(object = model_logit_list, newx = as.matrix(dat.te[module_names]), type = "response")
#   pred.prob5 = pred2[,3,1]
#   pred.prob.others = pred1*(1.0-pred.prob5)
#   pred.prob.all = cbind(pred.prob.others, pred.prob5) 
#   return(list(pred.prob.all=pred.prob.all, pred.ordinal = pred1, pred.logit = pred2))
# }
# 
# #prediction_hierachicalV2 do not perform any grouping, grouping is done at a post-training fasion later
# prediction_hierachicalV2 = function(dat.te, module_names, model_ordinal ,   model_logit_list){
#   pred1 = predict(model_ordinal, dat.te, type = "probs")
#   pred2 = sapply(model_logit_list, function(z) predict(z, as.matrix(dat.te[module_names]), type = "response", s = "lambda.min"))
#   pred.prob5 = pred2*pred1/(1.0-pred1+pred2*pred1)
#   tmp1 = pred.prob5 + pred.prob5*(1.0-pred1)
#   pred.prob5 = apply(pred.prob5 * tmp1,1,sum)/apply(tmp1,1,sum)
#   pred.prob.others = pred1*(1.0-pred.prob5)
#   pred.prob.all = cbind(pred.prob.others, pred.prob5) 
#   return(list(pred.prob.all=pred.prob.all, pred.ordinal = pred1, pred.logit = pred2))
# }
# 
# prediction_prob_combined = function(pred.ordinal, pred.logit){
#   pred1 = pred.ordinal
#   pred2 =  pred.logit
#   pred.prob5 = apply(pred1 * pred2,1,sum)/apply(1.0-pred2+pred1 * pred2,1,sum)
#   pred.prob.others = apply(pred1,2,function(z) z*(1.0-pred.prob5))
#   it = 0
#   pred.prob.others = pred1*(1.0-pred.prob5)
#   while(it < 10){
#     tmp1 = pred.prob5 +apply(pred.prob.others,2,function(z) z*(1.0-pred.prob5))
#     pred.prob5 = pred2*pred1/(1.0-pred2+pred2*pred1)
#     pred.prob5 = apply(pred.prob5 * tmp1,1,sum)/apply(tmp1,1,sum)
#     pred.prob.others = apply(pred1,2,function(z) z*(1.0-pred.prob5))
#     it = it +1
#   }
#   pred.prob.all = cbind(pred.prob.others, pred.prob5) 
#   return(pred.prob.all=pred.prob.all)
# }


balance_predictions = function(probs){
  if(is.null(dim(probs))){
    probs = log(probs/(1.0 -probs)) - log(mean(probs)/(1-mean(probs)))
    probs = 1.0/(exp(-probs)+1)
  }else{
    K = ncol(probs)
    probs = t(apply(probs, 1, function(z) z/sum(z)))
    log_probs = log(probs)
    for(k in 1:K){
      log_probs[,k] = log_probs[,k]-log(mean( probs[,k]))+log(1/K)
    }
    log_probs = log_probs - apply(log_probs,1,function(z) matrixStats::logSumExp(z))
    probs = exp(log_probs)
  }
  return(probs)
}

balance_predictionsV2 = function(probs, train_ratio = NULL,balance_ratio = NULL){
  if(!is.null(train_ratio)|!is.null(balance_ratio)){
    if(is.null(dim(probs))){
      probs = log(probs/(1.0 -probs)) - log(train_ratio[2]/(1-train_ratio[2]))+log(balance_ratio[2]/(1-balance_ratio[2]))
      probs = 1.0/(exp(-probs)+1)
    }else{
      K = ncol(probs)
      probs = t(apply(probs, 1, function(z) z/sum(z)))
      log_probs = log(probs)
      for(k in 1:K){
        log_probs[,k] = log_probs[,k]-log(train_ratio[k])+log(balance_ratio[k])
      }
      log_probs = log_probs - apply(log_probs,1,function(z) matrixStats::logSumExp(z))
      probs = exp(log_probs)
    }
  }
  return(probs)
}

balanced_errors = function(z, balanced = T){
  if(balanced){
    s = 0
    for(i in 1:nrow(z)){
      s = sum(z[i,-i])/sum(z[i,])+s
    }
    s = s/nrow(z)
  }else{
    s = sum(z[row(z)!=col(z)])/sum(z)
  }
  
  return(s)
}

cvfit_evalution = function(probs, y, measure ="sens", x.measure ="fpr" ){
  cv_predictions =   probs
  if(is.null(dim( cv_predictions))){
    cv_predictions = cbind(1-cv_predictions,cv_predictions)
  }
  ypred =apply(cv_predictions, 1,which.max)
  confusion_integration = array(0, c(ncol(cv_predictions),ncol(cv_predictions)))
  for(k in 1:ncol(cv_predictions)){
    for(k1 in 1:ncol(cv_predictions)){
      confusion_integration[k,k1] = sum(ypred == k1&y==k)
    }
  }
  y_extended = matrix(0, nrow = length(y), ncol = ncol(probs))
  for(k in 1:ncol(probs)){
    y_extended[y==k,k]=1
  }
  res_return = list()
  res_return$berr = round(balanced_errors(confusion_integration, balanced = T),3)
  res_return$err = round(balanced_errors(confusion_integration, balanced = F),3)
  res_return$mse =mean((y_extended-probs)^2)
  res_return$extremal = mean(abs(ypred-y)>=2)
  #res_return$roc_lists = list()
  res_return$confusion= confusion_integration
  # for(k in 1:ncol(confusion_integration)){
  #   tmp = ROCR::prediction(cv_predictions[,k], ifelse(y==k,1,0))
  #   perf = performance(tmp ,measure,x.measure)
  #   res_return$roc_lists[[k]]=perf
  # }
  return(res_return)
}


make_pred_result_tab = function(pred.test, y.mutated, xmiss, balance = T){
  prediction_cv0 = pred.test
  train_ratio = table(y.mutated)/length(y.mutated)
  balance_ratio = rep(1,length(train_ratio))/length(train_ratio)
  if(balance){
    prediction_cv =balance_predictionsV2(prediction_cv0, train_ratio=train_ratio, balance_ratio = balance_ratio)
  }else{
    prediction_cv = prediction_cv0
  }
  K=ncol(pred.test)
  prediction_result_tableI = data.frame(matrix(NA, ncol = (K+1), nrow = (K+1)*3))
  prediction_result_tableII = data.frame(matrix(NA, ncol = (K+1), nrow = (K+1)*3))
  colnames(prediction_result_tableI) <- colnames(prediction_result_tableII) <-c("",paste0("pred",1:K))
  prediction_result_tableI[c(c(2:(K+1)),c((K+1+2):(2*K+2)),c((2*K+2+2):(3*K+3))),1] <-
    prediction_result_tableII[c(c(2:(K+1)),c((K+1+2):(2*K+2)),c((2*K+2+2):(3*K+3))),1] <-
    rep(paste0("TG",1:K),3)
  
  n = nrow(prediction_cv)
  tmp1 = cvfit_evalution(prediction_cv, y= y.mutated, measure ="sens", x.measure ="fpr" )
  prediction_result_tableI[1,1] <-prediction_result_tableII[1,1] <- paste0("no filtering:",n, ", berr:", tmp1$berr, ", mse:", round(tmp1$mse,3))
  prediction_result_tableI[c(2:(K+1)),-1] = tmp1$confusion
  prediction_result_tableII[c(2:(K+1)),-1]<-t(apply(tmp1$confusion,1,function(z) z/sum(z)))

  ll_complete = which(apply(xmiss,1,sum)<=3)
  if(balance){
    prediction_cv =balance_predictionsV2(prediction_cv0[ll_complete,], train_ratio=train_ratio, balance_ratio = balance_ratio)
  }else{
    prediction_cv = prediction_cv0[ll_complete,]
  }
  n = nrow(prediction_cv)
  prediction_result_tableI[(K+1+1),1] <-prediction_result_tableII[6,1] <- paste0("half observed:",length(ll_complete), ", berr:", tmp1$berr, ", mse:", round(tmp1$mse,3))
  prediction_result_tableI[c((K+1+2):(2*K+2)),-1] = tmp1$confusion
  prediction_result_tableII[c((K+1+2):(2*K+2)),-1]<-t(apply(tmp1$confusion,1,function(z) z/sum(z)))
  
  ll_complete = which(apply(xmiss,1,sum)<=0)
  if(balance){
    prediction_cv = balance_predictionsV2(prediction_cv0[ll_complete,], train_ratio=train_ratio, balance_ratio = balance_ratio)
  }else{
    prediction_cv = prediction_cv0[ll_complete,]
  }
  
  tmp1 = cvfit_evalution(prediction_cv, y= y.mutated[ll_complete], measure ="sens", x.measure ="fpr" )
  prediction_result_tableI[(2*K+2+1),1] <-prediction_result_tableII[11,1] <- paste0("all observed:",length(ll_complete), ", berr:", tmp1$berr, ", mse:", round(tmp1$mse,3))
  prediction_result_tableI[(2*K+2+2):(3*K+3),-1] = tmp1$confusion
  prediction_result_tableII[(2*K+2+2):(3*K+3),-1]<-t(apply(tmp1$confusion,1,function(z) z/sum(z)))
  
  return(list(prediction_result_tableI=prediction_result_tableI, prediction_result_tableII=prediction_result_tableII))
}


##################################
# Subtyping
##################################
make_data_long = function(factor_scores, sample_meta, visit_max, mean_removal = F){
  participant_ids = unique(sample_meta$participant[sample_meta$vists==1 & sample_meta$escalation=="Visit"])
  sample_meta_participant_level = unique(sample_meta[,c("participant", "traj", "sex", "discharge_dates", "death_date" ,"age", "age_discrete" ,"enrollment_site")])
  rownames(sample_meta_participant_level) = sample_meta_participant_level$participant
  sample_meta_participant_level = sample_meta_participant_level[participant_ids,]
  X = data.frame(matrix(NA, nrow = length(participant_ids), ncol =  visit_max*ncol(factor_scores)))
  rownames(X) =  participant_ids
  for(v in 1:visit_max){
    ll = which(sample_meta$vists==v&sample_meta$escalation=="Visit")
    sub_factors = factor_scores[ll,]
    sub_sample_meta = sample_meta[ll,]
    rownames(sub_factors) =sub_sample_meta$participant
    ll1 = intersect(rownames(sub_factors),rownames(X))
    ll2 = (ncol(factor_scores)*(v-1)+1):(ncol(factor_scores)*v)
    X[ll1,ll2] =  sub_factors[ll1,]
  }
  if(mean_removal){
    for(j in 1:ncol(X)){
      ll = !is.na(X[,j])
      X[ll,j] = X[ll,j] - mean(X[ll,j])
    }
  }
  return(list(X = X, sample_meta_participant_level=sample_meta_participant_level))
  ###back trimming until 
}

#######################################
#Trajectory analysis
######################################
#Maintainer: Cole Maguire



#######################################
#Intepretation
######################################
#Maintainer: Leqi Xu and Pramod
Read_KEGG <- function(pathway = "/data/resources/databases/KEGG/current/ko00001.keg"){
  kegg_raw = readLines(pathway)
  kegg_raw = kegg_raw[-c(1:3,56205:56208)]
  level = substr(kegg_raw,1,1)
  
  cate_p = which(level %in% c("A"))
  subcate_p = which(level %in% c("B"))
  path_p = which(level %in% c("C"))
  gene_p = which(level %in% c("D"))
  
  ##gene data info (D)
  gene_data <- data.frame(all = substring(kegg_raw[gene_p],7))
  gene_data <- gene_data %>%
    mutate(gene_kegg_id = substr(all,1,7)) %>% 
    mutate(gene_all = substring(all,9)) %>% 
    separate(gene_all,c("gene_name","enzyme_name_all"),sep=";",extra = "merge") %>%
    separate(enzyme_name_all,c("enzyme_name","enzyme_id_raw"),sep = "\\[EC") %>%
    separate(enzyme_id_raw,c("enzyme_id_raw"),sep = "\\]") %>%
    mutate(enzyme_id = paste0("EC",enzyme_id_raw)) %>%
    dplyr::select(gene_kegg_id,gene_name,enzyme_name, enzyme_id)
  gene_data$enzyme_id <- str_replace_all(gene_data$enzyme_id,"ECNA","NA")
  
  ##path data info (C)
  path_data <- data.frame(all = substring(kegg_raw[path_p],6))
  path_data <- path_data %>%
    mutate(path_kegg_id = substr(all,1,6)) %>%
    mutate(path_all = substring(all,7)) %>%
    separate(path_all,c("path_name","path_id_raw"),sep = "\\[PATH") %>%
    separate(path_id_raw,c("path_id_raw"),sep = "\\]") %>%
    mutate(path_id = paste0("PATH",path_id_raw)) %>%
    dplyr::select(path_kegg_id,path_name,path_id)
  path_data$path_id <- str_replace_all(path_data$path_id,"PATHNA","NA")
  
  ##path_gene info (C-D)
  path_repeat_num <- c()
  path_kegg_id_repeat <- c()
  for(i in 1:(length(path_p)-1)){
    path_repeat_num <- c(path_repeat_num, length(which(gene_p > path_p[i] & gene_p < path_p[i+1])))
    path_kegg_id_repeat <- c(path_kegg_id_repeat, rep(path_data$path_kegg_id[i],path_repeat_num[i]))
  }
  path_repeat_num <- c(path_repeat_num, length(which(gene_p > path_p[i+1])))
  path_kegg_id_repeat <- c(path_kegg_id_repeat, rep(path_data$path_kegg_id[i+1],path_repeat_num[i+1]))
  
  path_gene <- data.frame(path_kegg_id = path_kegg_id_repeat, gene_kegg_id = gene_data$gene_kegg_id)
  
  ##cate info (A)
  cate_data <- data.frame(all = substring(kegg_raw[cate_p],2))
  cate_data <- cate_data %>%
    mutate(cate_kegg_id = substr(all,1,6)) %>%
    mutate(cate_name = substring(all,7)) %>%
    dplyr::select(cate_kegg_id,cate_name)
  
  ##subcate info (B)
  subcate_data <- data.frame(all = substring(kegg_raw[subcate_p],3))
  subcate_p <- subcate_p[which(subcate_data$all != "")]
  subcate_data <- data.frame(all = subcate_data[which(subcate_data$all != ""),])
  subcate_data <- subcate_data %>%
    mutate(subcate_kegg_id = substr(all,1,6)) %>%
    mutate(subcate_name = substring(all,7)) %>%
    dplyr::select(subcate_kegg_id,subcate_name)
  
  ##cate_subcate info (A-B)
  cate_repeat_num <- c()
  cate_kegg_id_repeat <- c()
  for(i in 1:(length(cate_p)-1)){
    cate_repeat_num <- c(cate_repeat_num, length(which(subcate_p > cate_p[i] & subcate_p < cate_p[i+1])))
    cate_kegg_id_repeat <- c(cate_kegg_id_repeat, rep(cate_data$cate_kegg_id[i],cate_repeat_num[i]))
  }
  cate_repeat_num <- c(cate_repeat_num, length(which(subcate_p > cate_p[i+1])))
  cate_kegg_id_repeat <- c(cate_kegg_id_repeat, rep(cate_data$cate_kegg_id[i+1],cate_repeat_num[i+1]))
  
  cate_subcate <- data.frame(cate_kegg_id = cate_kegg_id_repeat, subcate_kegg_id = subcate_data$subcate_kegg_id)
  
  ##subcate_path info (B-C)
  subcate_repeat_num <- c()
  subcate_kegg_id_repeat <- c()
  for(i in 1:(length(subcate_p)-1)){
    subcate_repeat_num <- c(subcate_repeat_num, length(which(path_p > subcate_p[i] & path_p < subcate_p[i+1])))
    subcate_kegg_id_repeat <- c(subcate_kegg_id_repeat, rep(subcate_data$subcate_kegg_id[i],subcate_repeat_num[i]))
  }
  subcate_repeat_num <- c(subcate_repeat_num, length(which(path_p > subcate_p[i+1])))
  subcate_kegg_id_repeat <- c(subcate_kegg_id_repeat, rep(subcate_data$subcate_kegg_id[i+1],subcate_repeat_num[i+1]))
  
  subcate_path <- data.frame(subcate_kegg_id = subcate_kegg_id_repeat, path_kegg_id = path_data$path_kegg_id)
  
  ##cate_subcate_path infor (A-B-C)
  cate_subcate_path <- merge(cate_subcate,subcate_path)[,c(2,1,3)]
  
  kegg_data <- list(gene_data = gene_data, path_data = path_data, subcate_data = subcate_data, 
                    cate_data = cate_data, path_gene = path_gene, cate_subcate_path = cate_subcate_path)
  return(kegg_data)
}

nameMap_LG <- function(factor_loadings, row_feature_files, data_env0){
  mapped_name_frames = list()
  for(d in 1:length(factor_loadings)){
    ##Protein: current name, protein name, UNIPROT, NA
    ##Olink: current name, protein name, UNIPROT, NA
    ##Metabolite: current name, metabolite name, KEGG
    ##gene: current name, gene name, entrezID
    tmp1 = data_env0[[row_feature_files[d]]]
    mapped_name_frames[[d]] = data.frame(matrix(NA, ncol = 4, nrow = nrow(factor_loadings[[d]])))
    mapped_name_frames[[d]][,1] <-rownames(mapped_name_frames[[d]]) <-rownames(factor_loadings[[d]])
    if(grepl("proteomics",row_feature_files[d])){
      tmp1 = tmp1[mapped_name_frames[[d]][,1], ,drop = F]
      mapped_name_frames[[d]][,2] = rownames(tmp1)
      mapped_name_frames[[d]][,3] = tmp1[,1]
      colnames(mapped_name_frames[[d]]) = c("ori", "protein", "uniprotID", "empty")
    }else if(grepl("olink",row_feature_files[d])){
      tmp1$protein=rownames(tmp1)
      rownames(tmp1) = tmp1$Gene_name
      tmp1 = tmp1[mapped_name_frames[[d]][,1], ,drop = F]
      mapped_name_frames[[d]][,2] = mapped_name_frames[[d]][,1]
      mapped_name_frames[[d]][,3] = tmp1$Uniprot.ID
      colnames(mapped_name_frames[[d]]) = c("ori", "gene", "uniprotID", "empty")
    }else if(grepl("metabolomics",row_feature_files[d])){
      tmp1 = tmp1[mapped_name_frames[[d]][,1], ,drop = F]
      mapped_name_frames[[d]][,2] =  tmp1$CHEMICAL_NAME
      mapped_name_frames[[d]][,3] = tmp1$HMDB
      mapped_name_frames[[d]][,4] = tmp1$KEGG
      colnames(mapped_name_frames[[d]]) = c("ori", "CHEMICAL_NAME",  "HMDB", "KEGG")
    }else if(grepl("transcriptomics", row_feature_files[d])){
      tmp1 = tmp1[mapped_name_frames[[d]][,1], ,drop = F]
      mapped_name_frames[[d]][,2] = tmp1$gene_name
      ###map to entrez id
      tmp2=clusterProfiler::bitr(geneID = mapped_name_frames[[d]][,1],
                                fromType = "ENSEMBL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db,
                                drop = F)
      
      ###map to multiple NA
      tmp21 = table(tmp2[,1])
      tmp22= names(tmp21[tmp21==1])
      tmp2 = tmp2[tmp2[,1]%in%tmp22,]
      tmp21 = table(tmp2[,2])
      tmp22 = names(tmp21[tmp21==1])
      tmp2 = tmp2[tmp2[,2]%in%tmp22,]
      mapped_name_frames[[d]][ tmp2[,1],3] = tmp2[,2]
      colnames(mapped_name_frames[[d]]) = c("ori", "gene", "ENTREZID", "empty")
    }else{
      stop("unsupported data types")
    }
  }
  return(mapped_name_frames)
  
}

subpath_metabo_extraction = function(){
  tmp1 = data_env0[["plasma_metabolomics_global_rowfeature"]]
  tmp1 = tmp1[!is.na(tmp1$SUB_PATHWAY),]
  tmp1$MET_ID = rownames(tmp1)
  subpaths =  unique(tmp1$SUB_PATHWAY)
  superpaths =  unique(tmp1$SUPER_PATHWAY)
  subpath_mat = NULL
  superpath_mat = NULL
  for(i in 1:length(superpaths)){
    idxs = which(tmp1$SUPER_PATHWAY==superpaths[i])
    tmp2 = tmp1[idxs, c("SUPER_PATHWAY", "MET_ID")]
    if(is.null(superpath_mat)){
      superpath_mat= tmp2
    }else{
      superpath_mat= rbind(superpath_mat, tmp2)
    }
  }
  for(i in 1:length(subpaths)){
    idxs = which(tmp1$SUB_PATHWAY==subpaths[i])
    tmp2 = tmp1[idxs, c("SUB_PATHWAY", "MET_ID")]
    if(is.null(subpath_mat)){
      subpath_mat = tmp2
    }else{
      subpath_mat = rbind(subpath_mat, tmp2)
    }
  }
  return(list(subpath_mat = subpath_mat, superpath_mat = superpath_mat))
}

TERM2GENE_creation_LG = function(databases = c("kegg")){
  pathways <- getMultiOmicsFeatures( dbs = databases, layer ="all",
                                     returnTranscriptome = "ENSEMBL",
                                     returnProteome = "UNIPROT",
                                     returnMetabolome = "HMDB",
                                     useLocal = FALSE)
  TERM2GENE = list()
  for(i in 1:3){
    print(paste0("#########", i, "#########"))
    tmp1 = pathways[[i]]
    pathway_names = names(tmp1)
    pathway_names = sapply(strsplit( pathway_names, ") "), function(z) z[[2]])
    for(j in 1:length(tmp1)){
      tmp21 = tmp1[[j]]
      if(length(tmp21)>0){
        tmp22 = cbind(rep(pathway_names[j],length(tmp21)), tmp21)
        colnames(tmp22) = c("term", "gene")
        if(is.null(TERM2GENE[i][[1]])){
          TERM2GENE[[i]] =tmp22
        }else{
          TERM2GENE[[i]] = rbind(TERM2GENE[[i]], tmp22)
        }
      }
    }
  }
  check1 = sapply(pathways, function(z) sapply(z, function(z) length(z[!is.na(z)])))
  ll = which(apply(check1,1,sum)>0)
  for(i in 1:length(pathways)){
    pathways[[i]] = pathways[[i]][ll]
  }
  check1 = check1[ll,]
  names(TERM2GENE) = names(pathways)
  return(TERMlist = pathways,TERM2GENE= TERM2GENE, summary = check1)
}

GSEAanalysis_LG = function(factor_loadings_one, omic_short_names, mapped_names, kegg_reftab){
  kegg_test_one = list()
  for(d in 1:length(omic_short_names)){
    print(omic_short_names[d])
    gene_stats = factor_loadings_one[[d]][,1]
    gene_names = rownames(factor_loadings_one[[d]])
    names(gene_stats) = gene_names
    gene_stats = sort(gene_stats, decreasing = T)
    map0 = mapped_names[[d]][names(gene_stats),]
    if(omic_short_names[d]%in%c("ppt", "ppgd", "so")){
      gene_names_transformed = map0$uniprotID
      idxs = which(!is.na(gene_names_transformed))
      gene_stats = gene_stats[idxs]
      gene_stats1 = gene_stats
      names(gene_stats1) =gene_names_transformed[idxs] 
      tmp1 = clusterProfiler::gseKEGG(geneList =gene_stats1,organism = "hsa", keyType="uniprot",    pAdjustMethod = "fdr", minGSSize = 1, maxGSSize = 2000, pvalueCutoff = 1.2)
    }else if(omic_short_names[d] %in% c("nt","pbmc")){
      gene_names_transformed = map0$ENTREZID
      idxs = which(!is.na(gene_names_transformed))
      gene_stats = gene_stats[idxs]
      gene_stats1 = gene_stats
      names(gene_stats1) =gene_names_transformed[idxs] 
      tmp1 = clusterProfiler::gseKEGG(geneList =gene_stats1,organism = "hsa", keyType="ncbi-geneid",    pAdjustMethod = "fdr", minGSSize = 1, maxGSSize = 2000, pvalueCutoff = 1.2)
    }else if(omic_short_names[d] == "pmg"){
      gene_names_transformed = map0$HMDB
      idxs = which(!is.na(gene_names_transformed) & !grepl(",", gene_names_transformed))
      gene_stats1 = gene_stats[idxs]
      gene_names_transformed=gene_names_transformed[idxs] 
      gene_names_transformed_uni = table(gene_names_transformed)
      idxs = which(gene_names_transformed%in%names( gene_names_transformed_uni[  gene_names_transformed_uni==1]))
      gene_stats1 = gene_stats1[idxs]
      names(gene_stats1) =gene_names_transformed[idxs] 
      pos_list = gene_stats1[gene_stats1>=sort(gene_stats1, decreasing = T)[50]]
      # tmp11 = clusterProfiler::enricher(
      #                 names( pos_list), 
      #                 TERM2GENE=kegg_reftab$TERM2GENE$metabolome, 
      #                 minGSSize = 1,
      #                 pAdjustMethod = "BH",pvalueCutoff =0.5)
      # neg_list = gene_stats1[gene_stats1<=sort(gene_stats1, decreasing = F)[50]]
      # tmp12 = clusterProfiler::enricher(names(neg_list), 
      #                 TERM2GENE=kegg_reftab$TERM2GENE$metabolome, 
      #                 minGSSize = 1,pAdjustMethod = "BH",pvalueCutoff =0.5)
      tmp1 =clusterProfiler::GSEA(geneList =gene_stats1, TERM2GENE= kegg_reftab$TERM2GENE$metabolome,  minGSSize = 1,maxGSSize = 2000,by = "fgsea", pAdjustMethod = "none",pvalueCutoff =.999,verbose=TRUE)
      
    }
    kegg_test_one[[d]] = tmp1
    print(kegg_test_one[[d]])
  }
  names(kegg_test_one) =omic_short_names
  return(kegg_test_one)
}



GSEAanalysis_general_RKP = function(factor_loadings_one, omic_short_names, mapped_names, reftab){
  test_one = list()
  for(d in 1:length(omic_short_names)){
    print(omic_short_names[d])
    if(is.na(omic_short_names[d]))
      next
    gene_stats = factor_loadings_one[[d]][,1]
    gene_names = rownames(factor_loadings_one[[d]])
    names(gene_stats) = gene_names
    gene_stats = sort(gene_stats, decreasing = T)
    map0 = mapped_names[[d]][names(gene_stats),]
    if(omic_short_names[d]%in%c("ppt", "ppgd", "so")){
      gene_names_transformed = map0$uniprotID
      idxs = which(!is.na(gene_names_transformed))
      gene_stats = gene_stats[idxs]
      gene_stats1 = gene_stats
      names(gene_stats1) =gene_names_transformed[idxs] 
      TERM2GENE = reftab$TERM2GENE$proteome
    }else if(omic_short_names[d] %in% c("nt","pbmc")){
      gene_names_transformed = map0$ENTREZID
      idxs = which(!is.na(gene_names_transformed))
      gene_stats = gene_stats[idxs]
      gene_stats1 = gene_stats
      names(gene_stats1) =gene_names_transformed[idxs] 
      TERM2GENE = reftab$TERM2GENE$transcriptome
    }else if(omic_short_names[d] == "pmg"){
      gene_names_transformed = map0$HMDB
      idxs = which(!is.na(gene_names_transformed) & !grepl(",", gene_names_transformed))
      gene_stats1 = gene_stats[idxs]
      gene_names_transformed=gene_names_transformed[idxs] 
      gene_names_transformed_uni = table(gene_names_transformed)
      idxs = which(gene_names_transformed%in%names( gene_names_transformed_uni[  gene_names_transformed_uni==1]))
      gene_stats1 = gene_stats1[idxs]
      names(gene_stats1) =gene_names_transformed[idxs]
      TERM2GENE = reftab$TERM2GENE$metabolome
    }
    test_one[[d]] = clusterProfiler::GSEA(geneList =gene_stats1, TERM2GENE= TERM2GENE,  minGSSize = 1,maxGSSize = 2000,by = "fgsea", pAdjustMethod = "none",pvalueCutoff =.999,verbose=TRUE)
    print(test_one[[d]])
  }
  #names(test_one) = omic_short_names[!is.na(omic_short_names)]
  names(test_one) = omic_short_names[1:length(test_one)]
  return(test_one)
}


short_feature_list = function(factor_loadings_one, factor_loading_pval_one, factor_regression_one, omic_short_names, mapped_names,
                              qval_thr = 1e-2, max_counts = c(50, 100, 50, 100, 100, 100), driver_counts = c(20, 50, 20, 50, 50, 50)){
  gene_list_downstsream_top = list()
  gene_list_driven_top = list()
  ps = sapply(factor_loadings_one,nrow)
  pss = c(0,cumsum(ps))
  for(d in 1:length(omic_short_names)){
    map0 = mapped_names[[d]]
    gene_list_driven_top[[d]] = data.frame(factor_regression_one[[d]][order(-abs(factor_regression_one[[d]])),,drop = F])
    gene_list_driven_top[[d]] = cbind(map0[rownames(gene_list_driven_top[[d]]),c(1,2)],gene_list_driven_top[[d]])
    if(omic_short_names[d]%in%c("ppt","ppgd","so" )){
      colnames(gene_list_driven_top[[d]]) = c("original_name","protein_name","reg.coef")
    }else if(omic_short_names[d]%in%c("nt","pbmc")){
      colnames(gene_list_driven_top[[d]]) = c("original_name","gene_name","reg.coef")
    }else if(omic_short_names[d]=="pmg"){
      colnames(gene_list_driven_top[[d]]) = c("original_name","chemical_name","reg.coef")
    }
    gene_list_driven_top[[d]] = gene_list_driven_top[[d]][1:driver_counts[d],,F]
    idx1 = which(factor_loading_pval_one[[d]]<=qval_thr/length(factor_loading_pval_one[[d]]))
    if(length(idx1)>0){
      gene_list_downstsream_top[[d]] = factor_loadings_one[[d]][idx1,,drop = F]
      gene_list_downstsream_top[[d]] = gene_list_downstsream_top[[d]][order(-abs(gene_list_downstsream_top[[d]])),,drop = F]
      if(length(gene_list_downstsream_top[[d]][,1])>max_counts[d]){
        gene_list_downstsream_top[[d]] = gene_list_downstsream_top[[d]][1:max_counts[d],,drop= F]
      }
      gene_list_downstsream_top[[d]] = cbind(map0[rownames(gene_list_downstsream_top[[d]]),c(1,2)],gene_list_downstsream_top[[d]] )
      if(omic_short_names[d]%in%c("ppt","ppgd","so" )){
        colnames(gene_list_downstsream_top[[d]]) = c("original_name","protein_name","proj.coef")
      }else if(omic_short_names[d]%in%c("nt","pbmc")){
        colnames(gene_list_downstsream_top[[d]]) = c("original_name","gene_name","proj.coef")
      }else if(omic_short_names[d]=="pmg"){
        colnames(gene_list_downstsream_top[[d]]) = c("original_name","chemical_name","proj.coef")
      }
    }else{
      gene_list_downstsream_top[[d]] = NULL
    }
  }
  names(gene_list_driven_top) <- names(gene_list_downstsream_top) <- omic_short_names
  return(list(gene_list_driven_top = gene_list_driven_top, gene_list_downstsream_top = gene_list_downstsream_top))
}


data_heatmap_prepare_func_LG = function(data_subsets_imputed, top_annotations_list){
  data_prepared_plot = list()
  for(d in 1:length(data_subsets_imputed)){
    data_d = data_subsets_imputed[[d]]
    data_d = data_d[,order(top_annotations_list[[d]][,1])]
    data_d$mod = factors_one_visit1
    data_d$traj =traj_group
    data_d$type = omic_short_names[d]
    data_d = data_d[order(data_d$traj,data_d$mod),]
    data_prepared_plot[[d]] = data_d
  }
  return(data_prepared_plot)
}

data_heatmap_plot_func_LG = function(data_gather0, top_annotations_list, col_fun1, col_fun2, column_fontsize = 7){
  data_gather = NULL
  annotations_row = NULL
  annotations_column = NULL
  data_type = NULL
  for(d in 1:length(data_gather0)){
    tmp = data_gather0[[d]]
    annotations_row = tmp[,c("traj","mod")]
    if(is.null(data_gather)){
      data_type = rep(tmp[1,c("type")],ncol(tmp)-3)
      data_gather =tmp[,!(colnames( tmp)%in%c("traj","mod","type"))]
      annotations_column = data.frame(mod.coef =sort(top_annotations_list[[d]][,1]),
                                      traj.cor = top_annotations_list[[d]][order(top_annotations_list[[d]][,1]),2])
      
    }else{
      data_type =c(data_type, rep(tmp[1,c("type")],ncol(tmp)-3))
      data_gather = cbind(data_gather,tmp[,!(colnames( tmp)%in%c("traj","mod","type"))])
      annotations_column = rbind(annotations_column,data.frame(mod.coef =sort(top_annotations_list[[d]][,1]),
                                                               traj.cor = top_annotations_list[[d]][order(top_annotations_list[[d]][,1]),2]))
      
    }
  } 
  row_annotations =  rowAnnotation(traj =annotations_row$traj, mod =annotations_row$mod, col = list(traj=brewer.reds(5) %>% setNames(1:5),mod=col_fun2))
  top_annotations = columnAnnotation(mod.coef =annotations_column$mod.coef,
                                     traj.cor = annotations_column$traj.cor, 
                                     col = list(mod.coef=col_fun1, traj.cor=col_fun1))
  p <-
    Heatmap(
      data_gather, 
      column_title_rot = 90,
      column_title_gp = gpar(fontsize=10),
      show_row_names = FALSE,
      cluster_rows  = F,
      cluster_columns   = F,
      column_split = data_type,
      row_split = annotations_row$traj,
      #show_column_names = FALSE,
      left_annotation = row_annotations,
      top_annotation =  top_annotations,
      column_names_gp = grid::gpar(fontsize = column_fontsize)
    )
  return(p)
}

#' @param min_num minimum number of pathways to show
#' @param max_num maximum number of pathways to show
enrichment_gseaplot1 = function(kegg_test_one, mod_name, omic_names, thr = 0.05, min_num=3, max_num = 50,  driving_assays = NULL){
  pathway_selected_ids = c()
  pathway_names = c()
  if(is.null(  driving_assays)){
    driving_assays = rep(1,length(kegg_test_one))
  }
  for(d in 1:length(kegg_test_one)){
    A = kegg_test_one[[d]]
    if(class(A) == "gseaResult"){
      tab <- A@result
      pathways <- tab$Description[tab$core_enrichment!=""]
      qvals = tab$qvalues[tab$core_enrichment!=""]
      descriptions = tab$Description[tab$core_enrichment!=""]
      pathway_selected_ids = c(pathway_selected_ids,pathways)
      pathway_names = c(pathway_names, descriptions)
    }
  }
  pathway_selected_ids = unique(pathway_selected_ids)
  pathway_names = unique(pathway_names)
  pathway.matrix  = list()
  pathway.matrix.weighted = list()
  qvals = matrix(NA, ncol = length( kegg_test_one), nrow = length(pathway_selected_ids))
  signs = matrix(NA, ncol = length( kegg_test_one), nrow = length(pathway_selected_ids))
  for(d in 1:length( kegg_test_one)){
    A = kegg_test_one[[d]]
    if(class(A) == "gseaResult"){
      tab <- A@result
      pathways <-pathway_selected_ids
      pathway.list <- list()
      qvals_assay = c()
      nes_assay = c()
      for(pathway in pathways){
        is = which(tab$Description == pathway & tab$core_enrichment!="")
        if(length(is)>0){
          genes <- tab$core_enrichment[is]
          genes <- stringr::str_split(genes, "/")[[1]]
          pathway.list[[pathway]] <- genes
          qvals_assay = c(qvals_assay,tab$pvalue[is])
          nes_assay = c(nes_assay, tab$NES[is])
        }
      }
      pathway.matrix[[d]] <- matrix(0, nrow = length(unique(unlist(pathway.list))), ncol = length(pathway_selected_ids))
      pathway.matrix.weighted[[d]] =  pathway.matrix[[d]]
      rownames(pathway.matrix[[d]]) <-rownames(pathway.matrix.weighted[[d]])<- unique(unlist(pathway.list))
      colnames(pathway.matrix[[d]]) <- colnames(pathway.matrix.weighted[[d]]) <-pathways
      for(j in 1:length(pathway.list)){
        j0 = which(pathway_selected_ids==names(pathway.list)[j])
        pathway.matrix[[d]][pathway.list[[j]],j0] = 1
        qvals[j0,d] = qvals_assay[j]
        signs[j0,d] = sign(nes_assay[j])
        pathway.matrix.weighted[[d]][pathway.list[[j]],j0] = A@geneList[pathway.list[[j]]]
      }
    }
  }
  qvals = data.frame(qvals)
  colnames(qvals) = omic_names
  qvals$joint = apply(qvals,1,function(z){
    ii = !is.na(z)
    t = qnorm(z[ii])
    t = sum(t*driving_assays[ii])/sqrt(sum(driving_assays[ii]^2))
    pnorm(t, lower.tail = T)
    #t = -2*sum(log(z[ii][driving_assays[ii]==1]))
    #pchisq(t, df=2*length(ii),lower.tail = F)
  })
  rownames(qvals)=pathway_names
  library(qvalue)
  try1 <-try(qvalue(qvals$joint))
  if(class(try1)=="try-error"){
    qvals$joint.adjust = qvalue::qvalue(qvals$joint, pi0 = 1, fdr.level=0.05)$qvalues
  }else{
    qvals$joint.adjust =  try1$qvalues
  }
 
  o = order(qvals$joint)
  qvals = qvals[o,,drop = F]
  pathway_names = pathway_names[o]
  signs =signs[o,,drop = F]
  idx2 =which(qvals$joint.adjust<=thr)
  if(length(idx2)<min_num){
    idx2 = which(qvals$joint <= sort(qvals$joint)[min_num])
  }else if(length(idx2)>max_num){
    idx2 = which(qvals$joint.adjust <= sort(qvals$joint.adjust)[max_num])
  }
  pathway_names = pathway_names[idx2]
  qvals = qvals[idx2,,drop = F]
  signs = signs[idx2,,drop = F]
  for(d in 1:length(kegg_test_one)){
    A = kegg_test_one[[d]]
    if(class(A) == "gseaResult"){
      pathway.matrix[[d]] = pathway.matrix[[d]][,o][,idx2]
      pathway.matrix[[d]] = pathway.matrix[[d]][apply(abs(pathway.matrix[[d]]),1,sum)!=0,]
      pathway.matrix.weighted[[d]] = pathway.matrix.weighted[[d]][,o][,idx2]
      pathway.matrix.weighted[[d]] = pathway.matrix.weighted[[d]][apply(abs(pathway.matrix.weighted[[d]]),1,sum)!=0,]
    }
  }
  signed_logpvals = abs(log(qvals,base = 10))
  signed_logpvals[,1:length(kegg_test_one)]  = signed_logpvals[,1:length(kegg_test_one)]*signs
  signed_logpvals = data.frame(signed_logpvals)
  signed_logpvals = signed_logpvals[,apply(!is.na(signed_logpvals),2,sum)>0]
  signed_logpvals$pathway = rownames(signed_logpvals)
  signed_logpvals = gather(signed_logpvals, key = "type", val = "signed.log10pval", -pathway)
  signed_logpvals$ES = ifelse(signed_logpvals$signed.log10pval<0,"neg", "pos")
  signed_logpvals$log10pval = abs(signed_logpvals$signed.log10pval)
  signed_logpvals$type = factor(signed_logpvals$type, levels = unique(signed_logpvals$type))
  signed_logpvals$ES[!signed_logpvals$type%in%omic_names] = "unsigned"
  signed_logpvals$ES = factor(signed_logpvals$ES, levels = c("unsigned", "pos","neg"))
  signed_logpvals$pathway = factor(signed_logpvals$pathway, levels = rownames(qvals)[nrow(qvals):1])
  signed_logpvals = signed_logpvals[signed_logpvals$type!="joint",]
  fig1 = ggplot(data =signed_logpvals, aes(x =type, y = pathway, size=log10pval,col =ES))+geom_point()+ggtitle(mod_name)+ scale_color_manual(values=c("black","red","blue"))
  return(list(pval_tb = qvals, pval_fig = fig1, pathway.matrix = pathway.matrix,  pathway.matrix.weighted = pathway.matrix.weighted,
              pathway.names = pathway_names))
}


### GSEAobject to overlap plot
#    Jeremy Gygi
#  Input: GSEA object from clusterProfiler::GSEA(...)
#  Output: Table with rows as union of all leading edge genes for each pathway,
#                     cols as each pathway
#                     entries as 1 (if gene is in leading edge for pathway) or 0 otherwise
#
# Necessary libraries
# ggplot2
# reshape2
# gplots
# ggdendro
# dplyr
# RColorBrewer
#
# Parameters:
# GSEAobject - object returned from GSEA function
# num_groups - number of cuts to make into the tree of pathways (will be represented by different colors)
# num_min_overlap - minimum number of overlapping features to be shown (defaults to 10)
# max_num_pathways - max number of pathways to be shown (ranked by p.adjust). Defaults to 10
# dendrogram.width - how wide should the dendrogram be in proportion to the plot? 1 would be equal size. Defaults to .15
#
# EXAMPLE: (remove "#. " from lines)
#. # Attach GSEAobject here:
#. load("/scratch/data-integration/pipelines/Interpretation/gsea_test_mcia_subset.RData")
#. GSEAobject <- kegg_test
#. # Parameters:
#. val <- 4
#. module <- c("Module 1", "Module 4", "Module 6", "Module 24")[val]
#. dataset <- "pbmc_transcriptomics"
#. # Get results table
#. GSEAobject <- GSEAobject[[val]][[dataset]]
#. # Make the plot:
#. generate_GSEA_table_plot(GSEAobject, max_num_pathways = 13)
generate_GSEA_table_plot <- function(GSEAobject, num_groups = 5, num_min_overlap = 10,
                                     max_num_pathways = 10, dendrogram.width = .15){
  tab <- GSEAobject
  # shorten if too many
  if(nrow(tab) > max_num_pathways){
    tab <- tab[1:max_num_pathways,]
  }
  # Make list of leading edge genes per pathway:
  pathways <- tab$Description
  pathway.list <- list()
  for(pathway in pathways){
    genes <- tab$core_enrichment[which(tab$Description == pathway)]
    genes <- stringr::str_split(genes, "/")[[1]]
    pathway.list[[pathway]] <- genes
  }
  
  # Make matrix:
  pathway.matrix <- matrix(0, nrow = length(unique(unlist(pathway.list))), ncol = length(pathway.list))
  rownames(pathway.matrix) <- unique(unlist(pathway.list))
  colnames(pathway.matrix) <- pathways
  for(i in 1:nrow(pathway.matrix)){
    for(j in 1:ncol(pathway.matrix)){
      cur.path <- pathways[j]
      cur.gene <- rownames(pathway.matrix)[i]
      if(cur.gene %in% pathway.list[[cur.path]]){
        pathway.matrix[i,j] <- 1
      }
    }
  }
  
  
  # Parameters:
  # l - a named list of event_ids (or whatever you want to get the overlap for)
  l <- pathway.list
  groups <- dplyr::arrange(tab, p.adjust)$Description
  color_vec <- viridis::magma(n = length(l), begin = .4)
  title <- module
  
  # Get the intersection df from the gplots::venn function
  v.table <- gplots::venn(l, show.plot = FALSE, simplify = TRUE, intersections = FALSE)
  # Set the attributes to NULL, make data.frame (otherwise you can't melt it to plot the points underneath)
  attr(v.table, "class") <- NULL
  v.table <- as.data.frame(v.table)
  # Filter out the case with 000 (we didn't supply a 'universe' parameter, so '000' can be removed)
  v.table <- dplyr::arrange(dplyr::filter(v.table, num >= num_min_overlap), num)
  
  ### Clustering:
  hc <- hclust(dist(t(pathway.matrix)))
  dend <- as.dendrogram(hc)
  dend_data <- ggdendro::dendro_data(dend)
  pathway.groups <- cutree(hc, k = num_groups)[dend_data$labels$label]
  
  ### TODO: Cut tree into how_many_param (let user decide)
  ### Add dot to end of dendrogram for 
  ### color = abs(log10(p.adjusted)) (positive = red, neg = blue)
  
  segment_data <- with(ggdendro::segment(dend_data), data.frame(x = y, y = x, xend = yend, yend = xend))
  #segment_data$x <- segment_data$x + length(values$static_names) + .5
  #segment_data$xend <- segment_data$xend + length(values$static_names) + .5
  gene_pos_table <- with(dend_data$labels, data.frame(y_center = x, PATHWAY = as.character(label), height = 1))
  # Rearrange columns:
  v.table <- v.table[c("num", dend_data$labels$label)]
  groups <- dend_data$labels$label
  # add p.values:
  p.values <- sapply(groups, function(g){
    ES <- tab$NES[which(tab$Description == g)]
    p.val <- -log10(tab$p.adjust[which(tab$Description == g)])
    if(ES >= 0){
      return(p.val)
    }
    return(-p.val)
  })
  
  # Add 'comb' to v.table (rownames are "000"..."111")
  v.table <- dplyr::mutate(v.table, comb = rownames(v.table))
  # Melt so we can plot the individual points (can use pivot_longer if preferred)
  v.table.melt <- reshape2::melt(v.table)
  # Convert group names into indices, remove the 'num' column
  v.table.melt <- dplyr::filter(v.table.melt, variable != "num")
  v.table.melt$ypos <- sapply(v.table.melt$variable, function(var){
    return(-which(groups == var))
  })
  # Shift the ypos up .3 (looks better in my opinion)
  v.table.melt$ypos <- v.table.melt$ypos + .3
  # Remove all points that shouldn't be there (i.e. 3rd point for '110', middle point for '101...')
  v.table.melt <- dplyr::filter(v.table.melt, value > 0)
  
  # Scaling parameters (used to make sure the points aren't too close together)
  # Feel free to edit these if they look weird in the plot
  text.scale <- max(v.table$num)/20
  y.scale <- max(v.table$num)/10
  
  # Group colors:
  group_colors <- RColorBrewer::brewer.pal(num_groups, "Set1")
  
  # Generate plots:
  g.overlap <- ggplot2::ggplot(v.table) +
    ggplot2::geom_bar(ggplot2::aes(x = factor(comb, levels = rev(comb)), y = num), fill = "#F0F0F0", color = "black", stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(x = factor(comb, levels = rev(comb)), y = num + text.scale, label = num)) +
    ggplot2::geom_hline(yintercept = y.scale * (seq(-1, -length(l)) + .3), size = .2) +
    ggplot2::geom_line(data = v.table.melt, ggplot2::aes(x = factor(comb, levels = rev(v.table$comb)), y = ypos * y.scale), size = 2) +
    ggplot2::geom_point(data = v.table.melt, ggplot2::aes(x = factor(comb, levels = rev(v.table$comb)), y = ypos * y.scale, fill = variable), size = 4, shape = 21) +
    ggplot2::scale_y_continuous(breaks = y.scale * (seq(-1, -length(l)) + .3), labels = groups,
                                limits = c(y.scale * -length(l) + .3, max(v.table$num) + text.scale)) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(NULL) +
    # Better theme is cowplot::theme_map(), but I didn't want to give too many dependencies
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8, hjust = 1, color = group_colors[pathway.groups]),
                   title = ggplot2::element_text(size = 10, face = "plain"),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = color_vec) +
    ggplot2::guides(size = FALSE, color = FALSE, fill = FALSE)
  
  
  g.dendrogram <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = segment_data$x, 
                                       y = y.scale * (-segment_data$y + .3), 
                                       xend = segment_data$xend, 
                                       yend = y.scale * (-segment_data$yend + .3)),
                          color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(x = .5,
                                     y = y.scale * (seq(-1, -length(l)) + .3),
                                     fill = p.values), shape = ifelse(p.values >= 0, 24, 25), size = 3) +
    ggplot2::ylim(c(y.scale * -length(l) + .3, max(v.table$num) + text.scale)) +
    ggplot2::scale_fill_gradient2(low = "darkblue", high = "darkred") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(vjust = .8, hjust = .5))
  
  g.dendrogram.legend <- cowplot::get_legend(g.dendrogram)
  g.dendrogram <- g.dendrogram + ggplot2::guides(fill = "none")
  g.total <- cowplot::plot_grid(g.overlap, g.dendrogram, nrow = 1, align = "h", rel_widths = c(1, dendrogram.width))
  return(cowplot::plot_grid(g.total, g.dendrogram.legend, ncol = 1, rel_heights = c(1, .15)))
}




generate_GSEA_table_plotV3 <- function(pathway.matrices, pathway.matrices.weighted, pvals, num_groups  = 5,align = T, fig_name = NULL, 
                                       dist.method = "binary",dist.method2 = "euclidean", max_num_pathways = 10, font_size = 8, width = 160, height = 160, 
                                       cor_max = 0.5, thr = 0){
  # shorten if too many
  if( (ncol(pathway.matrices[[1]]) > max_num_pathways)){
    pvals = pvals[1:max_num_pathways]
    for(d in 1:length(pathway.matrices)){
      pathway.matrices[[d]] = pathway.matrices[[d]][,1: max_num_pathways]
      pathway.matrices[[d]] = ifelse(abs(pathway.matrices[[d]])<=thr, 0, pathway.matrices[[d]])
      pathway.matrices.weighted[[d]] = pathway.matrices.weighted[[d]][,1: max_num_pathways]
    }
  }
  #align = T: match features in nt and pbmc
  if(align){
    tmp1 = rownames(pathway.matrices[[5]])
    tmp2 = rownames(pathway.matrices[[6]])
    tmp3 = union(tmp1, tmp2)
    tmp_mat1 = data.frame(array(0, dim = c(length(tmp3), ncol(pathway.matrices[[5]]))))
    tmp_mat2 = data.frame(array(0, dim = c(length(tmp3), ncol(pathway.matrices[[5]]))))
    colnames(tmp_mat1) <- colnames(tmp_mat2) <-  colnames(pathway.matrices[[5]])
    rownames(tmp_mat1) <- rownames(tmp_mat2) <-   tmp3
    tmp_mat1[rownames(pathway.matrices[[5]]),] = pathway.matrices[[5]]
    tmp_mat2[rownames(pathway.matrices[[6]]),] =pathway.matrices[[6]]
    pathway.matrices[[5]] =tmp_mat1
    pathway.matrices[[6]] =tmp_mat2
    tmp_mat1[rownames(pathway.matrices.weighted[[5]]),] = pathway.matrices.weighted[[5]]
    tmp_mat2[rownames(pathway.matrices.weighted[[6]]),] = pathway.matrices.weighted[[6]]
    pathway.matrices.weighted[[5]] =tmp_mat1
    pathway.matrices.weighted[[6]] =tmp_mat2
  }
  # Make list of leading edge genes per pathway:
  combined_pathway.matrix = NULL
  combined_pathway.matrix.weighted = NULL
  for(d in 1:length(pathway.matrices)){
    if(is.null(combined_pathway.matrix)){
      combined_pathway.matrix = pathway.matrices[[d]]
      combined_pathway.matrix.weighted = pathway.matrices.weighted[[d]]
    }else{
      combined_pathway.matrix = rbind(combined_pathway.matrix, pathway.matrices[[d]])
      combined_pathway.matrix.weighted = rbind(combined_pathway.matrix.weighted, pathway.matrices.weighted[[d]])
    }
  }
  ### Clustering: rows
  hc <- hclust(dist(t(combined_pathway.matrix),method = dist.method))
  dend <- as.dendrogram(hc)
  dend_data <- ggdendro::dendro_data(dend)
  pathway.groups <- cutree(hc, k = num_groups)[dend_data$labels$label]
  ### Clustering: columns
  hc_cluster_results =list()
  for(d in 1:4){
    if(nrow(pathway.matrices.weighted[[d]])>3){
      hc.tmp <- hclust(dist(pathway.matrices.weighted[[d]],method = dist.method2))
      hc_cluster_results[[d]] = hc.tmp
      pathway.matrices.weighted[[d]] =pathway.matrices.weighted[[d]][hc.tmp$order, ]
      pathway.matrices[[d]] =pathway.matrices[[d]][hc.tmp$order, ]
    }
  }
  if(align){
    hc.tmp <- hclust(dist(cbind(pathway.matrices.weighted[[5]],pathway.matrices.weighted[[6]]),method = dist.method2))
    for(d in c(5,6)){
      hc_cluster_results[[d]] = hc.tmp
      pathway.matrices.weighted[[d]] =pathway.matrices.weighted[[d]][hc.tmp$order, ]
      pathway.matrices[[d]] =pathway.matrices[[d]][hc.tmp$order, ]
    }
  }else{
    for(d in 5:6){
      if(nrow(pathway.matrices.weighted[[d]])>3){
        hc.tmp <- hclust(dist(pathway.matrices.weighted[[d]],method = dist.method2))
        hc_cluster_results[[d]] = hc.tmp
        pathway.matrices.weighted[[d]] =pathway.matrices.weighted[[d]][hc.tmp$order, ]
        pathway.matrices[[d]] =pathway.matrices[[d]][hc.tmp$order, ]
      }
    }
  }
  combined_pathway.matrix.weighted.reordered = NULL
  data_types = NULL
  for(d in 1:length(pathway.matrices)){
    if(is.null(combined_pathway.matrix)){
      combined_pathway.matrix.weighted.reordered = pathway.matrices.weighted[[d]]
      data_types = rep(omic_short_names[d], nrow(pathway.matrices.weighted[[d]]))
    }else{
      combined_pathway.matrix.weighted.reordered = rbind(combined_pathway.matrix.weighted.reordered, pathway.matrices.weighted[[d]])
      data_types = c(data_types,rep(omic_short_names[d], nrow(pathway.matrices.weighted[[d]])))
    }
  }
  #heatmap
  col_fun_tmp = colorRamp2(c(0, min(max(-log(pvals, base = 10)),10)), c("white", "red"))
  col_fun_tmp0 = unique(c(as.vector(RColorBrewer::brewer.pal(8, "Dark2")),as.vector(RColorBrewer::brewer.pal(12, "Paired"))))
                    #[c(1,2,5,6,7,10)], RColorBrewer::brewer.pal(9, "BuGn")[c(9)],RColorBrewer::brewer.pal(9, "Blues")[c(9)])
  
  col_fun_tmp1 =col_fun_tmp0[1:num_groups]
  names(col_fun_tmp1) = 1:num_groups
  llorder = order(pathway.groups,hc$order)
  row_annotations =  rowAnnotation(groups = pathway.groups[llorder], 
                                   log10p = -log(pvals[llorder ], base = 10),
                                   col = list(groups=col_fun_tmp1,
                                              log10p=col_fun_tmp), width = unit(.3, "cm"),
                                   show_legend = c(FALSE, FALSE)) 
  
  ll = which(data_types%in%omic_short_names[1:6])
  p1 <-
    Heatmap(
      t(combined_pathway.matrix.weighted.reordered)[llorder,ll], 
      cluster_rows = F, cluster_columns = F,
      show_row_names = T,
      show_column_names = F,
      column_title_rot = 90,
      column_title_gp = gpar(fontsize= font_size),
      column_gap = unit(.8, "mm"),
      column_split = data_types[ll],
      show_heatmap_legend = F,
      column_names_gp = grid::gpar(fontsize =font_size),
      row_names_gp = gpar(fontsize =font_size,col = col_fun_tmp1[pathway.groups[llorder]]),
      left_annotation = row_annotations,
      col =colorRamp2(c(-cor_max , 0, cor_max ), c("blue", "gray", "red")),
      width = unit(width, "mm"),
      height = unit(height, "mm")
    )
    return(p1)
}


generate_GSEA_table_plotV4 <- function(pathway.matrices, pathway.matrices.weighted, pvals, num_groups  = 5,align = T, fig_name = NULL, 
                                       dist.method = "binary",dist.method2 = "euclidean", hclust.method = "complete",max_num_pathways = 10, font_size = 8, width = 160, height = 160, 
                                       cor_max = 0.5, thr = 0){
  # shorten if too many
  if( (ncol(pathway.matrices[[1]]) > max_num_pathways)){
    pvals = pvals[1:max_num_pathways]
    for(d in 1:length(pathway.matrices)){
      pathway.matrices[[d]] = pathway.matrices[[d]][,1: max_num_pathways]
      pathway.matrices[[d]] = ifelse(abs(pathway.matrices[[d]])<=thr, 0, pathway.matrices[[d]])
      pathway.matrices.weighted[[d]] = pathway.matrices.weighted[[d]][,1: max_num_pathways]
    }
  }
  ##keep only useful genes
  for(d in 1:length(pathway.matrices)){
    ll_keep = which(apply(abs(pathway.matrices[[d]])>0,1,sum)>0)
    if(length(ll_keep)==0){
      pathway.matrices.weighted[[d]] =NULL
      pathway.matrices[[d]] = NULL
    }else{
      pathway.matrices[[d]] = pathway.matrices[[d]][ll_keep,,drop = F]
      pathway.matrices.weighted[[d]] = pathway.matrices.weighted[[d]][ll_keep,,drop = F]
    }
  }
  #align = T: match features in nt and pbmc
  if(align){
    tmp1 = rownames(pathway.matrices[[5]])
    tmp2 = rownames(pathway.matrices[[6]])
    tmp3 = union(tmp1, tmp2)
    tmp_mat1 = data.frame(array(0, dim = c(length(tmp3), ncol(pathway.matrices[[5]]))))
    tmp_mat2 = data.frame(array(0, dim = c(length(tmp3), ncol(pathway.matrices[[5]]))))
    colnames(tmp_mat1) <- colnames(tmp_mat2) <-  colnames(pathway.matrices[[5]])
    rownames(tmp_mat1) <- rownames(tmp_mat2) <-   tmp3
    tmp_mat1[rownames(pathway.matrices[[5]]),] = pathway.matrices[[5]]
    tmp_mat2[rownames(pathway.matrices[[6]]),] =pathway.matrices[[6]]
    pathway.matrices[[5]] =tmp_mat1
    pathway.matrices[[6]] =tmp_mat2
    tmp_mat1[rownames(pathway.matrices.weighted[[5]]),] = pathway.matrices.weighted[[5]]
    tmp_mat2[rownames(pathway.matrices.weighted[[6]]),] = pathway.matrices.weighted[[6]]
    pathway.matrices.weighted[[5]] =tmp_mat1
    pathway.matrices.weighted[[6]] =tmp_mat2
  }
  # Make list of leading edge genes per pathway:
  combined_pathway.matrix = NULL
  combined_pathway.matrix.weighted = NULL
  for(d in 1:length(pathway.matrices)){
    if(is.null(combined_pathway.matrix)){
      combined_pathway.matrix = pathway.matrices[[d]]
      combined_pathway.matrix.weighted = pathway.matrices.weighted[[d]]
    }else{
      combined_pathway.matrix = rbind(combined_pathway.matrix, pathway.matrices[[d]])
      combined_pathway.matrix.weighted = rbind(combined_pathway.matrix.weighted, pathway.matrices.weighted[[d]])
    }
  }
  ### Clustering: rows
  hc <- hclust(dist(t(combined_pathway.matrix),method = dist.method),method = hclust.method)
  dend <- as.dendrogram(hc)
  dend_data <- ggdendro::dendro_data(dend)
  pathway.groups <- cutree(hc, k = num_groups)[dend_data$labels$label]
  ### Clustering: columns
  hc_cluster_results =list()
  for(d in 1:4){
    if(nrow(pathway.matrices.weighted[[d]])>3){
      hc.tmp <- hclust(dist(pathway.matrices.weighted[[d]],method = dist.method2),method = hclust.method)
      hc_cluster_results[[d]] = hc.tmp
      pathway.matrices.weighted[[d]] =pathway.matrices.weighted[[d]][hc.tmp$order, ]
      pathway.matrices[[d]] =pathway.matrices[[d]][hc.tmp$order, ]
    }
  }
  if(align){
    hc.tmp <- hclust(dist(cbind(pathway.matrices.weighted[[5]],pathway.matrices.weighted[[6]]),method = dist.method2),method = hclust.method)
    for(d in c(5,6)){
      hc_cluster_results[[d]] = hc.tmp
      pathway.matrices.weighted[[d]] =pathway.matrices.weighted[[d]][hc.tmp$order, ]
      pathway.matrices[[d]] =pathway.matrices[[d]][hc.tmp$order, ]
    }
  }else{
    for(d in 5:6){
      if(nrow(pathway.matrices.weighted[[d]])>3){
        hc.tmp <- hclust(dist(pathway.matrices.weighted[[d]],method = dist.method2),method = hclust.method)
        hc_cluster_results[[d]] = hc.tmp
        pathway.matrices.weighted[[d]] =pathway.matrices.weighted[[d]][hc.tmp$order, ]
        pathway.matrices[[d]] =pathway.matrices[[d]][hc.tmp$order, ]
      }
    }
  }
  combined_pathway.matrix.weighted.reordered = NULL
  data_types = NULL
  for(d in 1:length(pathway.matrices)){
    if(is.null(combined_pathway.matrix)){
      combined_pathway.matrix.weighted.reordered = pathway.matrices.weighted[[d]]
      data_types = rep(omic_short_names[d], nrow(pathway.matrices.weighted[[d]]))
    }else{
      combined_pathway.matrix.weighted.reordered = rbind(combined_pathway.matrix.weighted.reordered, pathway.matrices.weighted[[d]])
      data_types = c(data_types,rep(omic_short_names[d], nrow(pathway.matrices.weighted[[d]])))
    }
  }
  #heatmap
  col_fun_tmp = colorRamp2(c(0, min(max(-log(pvals, base = 10)),10)), c("white", "red"))
  col_fun_tmp0 = unique(c(as.vector(RColorBrewer::brewer.pal(8, "Dark2")),as.vector(RColorBrewer::brewer.pal(12, "Paired"))))
  
  col_fun_tmp1 =col_fun_tmp0[1:num_groups]
  names(col_fun_tmp1) = 1:num_groups
  #llorder = order(pathway.groups,hc$order)
  llorder = hc$order
  row_annotations =  rowAnnotation(groups = pathway.groups, 
                                   log10p = -log(pvals[llorder], base = 10),
                                   col = list(groups=col_fun_tmp1,
                                              log10p=col_fun_tmp), width = unit(.3, "cm"),
                                   show_legend = c(FALSE, FALSE)) 
  
  ll = which(data_types%in%omic_short_names[1:6])
  p1 <-
    Heatmap(
      t(combined_pathway.matrix.weighted.reordered)[llorder,ll], 
      cluster_rows = F, cluster_columns = F,
      show_row_names = T,
      show_column_names = F,
      column_title_rot = 90,
      column_title_gp = gpar(fontsize= font_size),
      column_gap = unit(.8, "mm"),
      column_split = data_types[ll],
      show_heatmap_legend = F,
      column_names_gp = grid::gpar(fontsize =font_size),
      row_names_gp = gpar(fontsize =font_size,col = col_fun_tmp1[pathway.groups]),
      left_annotation = row_annotations,
      col =colorRamp2(c(-cor_max , 0, cor_max ), c("blue", "gray", "red")),
      width = unit(width, "mm"),
      height = unit(height, "mm")
    )
  return(p1)
}

#############################
#Others
###################################

visit_summary_plot = function(clinical_data_use){
  clinical_data_use$visits = paste0("Visit", clinical_data_use$visits)
  plotDF <- clinical_data_use %>%
    dplyr::filter(grepl(pattern = "Visit", event_type)) %>%
    add_column(flag = 1) %>%
    dplyr::select(participant_id,visits,flag) %>%
    distinct() %>%
    pivot_wider(names_from = visits, values_from = flag) %>%
    dplyr::select(-participant_id) %>%
    group_by_all() %>%
    summarize(n = n()) %>%
    arrange(desc(n)) %>%
    rowid_to_column() %>%
    pivot_longer(cols = -c(rowid, n), names_to = "event_type") %>%
    dplyr::filter(!is.na(value))
  
  ggplot(data = plotDF,
         mapping = aes(x = event_type, y = rowid)) +
    geom_point(mapping = aes(size = n)) +
    geom_line(mapping = aes(group = rowid)) +
    scale_y_continuous(breaks = 1:max(plotDF$rowid), 
                       labels = distinct(dplyr::select(plotDF,rowid, n))$n) +
    labs(y = "Number of participants") +
    theme_bw() +
    theme(axis.ticks.y = element_blank(),
          axis.text.x  = element_text(angle = 45, hjust = 1))
}






##########################
## Data preprocessing
##########################
##### transcriptomics

# Keep only protein-coding genes
preprocess_transcriptomics = function(assay, data_env){
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  rowfeature = data_env[[ data_env$assay_files_info[[assay]]$rowfeature_var_name ]]
  genes_pc <- rowfeature %>%
    # get rownames as an explicit column
    tibble::rownames_to_column(var = "gene_id") %>%
    # remove row names
    tibble::remove_rownames() %>%
    # keep only gene ids with valid hgnc symbols
    tidyr::drop_na(gene_name) %>%
    # keep only protein coding genes
    dplyr::filter(gene_biotype == "protein_coding") %>%
    # remove redundant genes if any and make unique
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    # add gene ids as row names
    column_to_rownames(var = "gene_id")
  
  
  # sample_id %in% rownames(nasal_transcriptomics_counts)
  
  # define good samples
  # No need of this step now as I already filtered for good samples while loading the data
  
  # goodSamples <- nasal_transcriptome_QCMetrics_data %>%
  #   filter(passQC == "Pass") %>%
  #   filter(sample_id %in% rownames(nasal_transcriptomics_counts)) %>%
  #   pull(sample_id)
  
  # get counts data only for pc genes
  counts1 <- counts[,rownames(genes_pc)] %>%
    t()
  counts1 <- counts1[rownames(genes_pc), ]  # 19835  1086
  
  # compute CPM values
  dge_list <- edgeR::DGEList(counts = counts1, genes = genes_pc)
  dge_list <- edgeR::calcNormFactors(dge_list)
  
  # expressed gene (CPM >=1) should be present in at least 10% of samples
  cut.filter <- 0.1
  keepRows <- rowSums(round(edgeR::cpm(dge_list$counts)) >= 1) >= cut.filter*ncol(counts1)
  #table(keepRows)
  curDGE <- dge_list[keepRows,]
  curDGE <- edgeR::calcNormFactors(curDGE) # 15957   1086
  voomCounts_1 <- limma::voom(curDGE, design = NULL, plot = FALSE, save.plot = FALSE)
  processed_data = t(voomCounts_1$E)
  
  #Scale the processed (normalized and log-transformed) data
  processed_data <- apply(processed_data, 2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(processed_data))
  
  return(processed_data)
}


##proteomics
preprocess_proteomics <- function(assay, data_env){
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  
  counts <- 2^counts #unlog2 the counts
  counts_return <- counts #copy. This is the one we will change and return
  
  NormalisationFactor <- median(rowSums(counts, na.rm = TRUE)) #calculate row(sample) sums. Take the median of all sums.
  
  #run for loop over each over the rows (samples)
  for(i in 1:nrow(counts)) {
    #get sample i, apply normalisation factor on it. Change sample i in the normalized df.
    counts_return[i,] <- counts[i,] * (NormalisationFactor / sum(counts[i,], na.rm = TRUE))
  }
  #remove samples with more than 50% missing features? remove samples with more 50% missing samples?
  
  #Impute with half min value per Protein
  for(i in 1:ncol(counts_return)) {       # for-loop over columns
    halfmin <- min(counts_return[,i], na.rm = TRUE)/2
    counts_return[,i][is.na(counts_return[,i])] <- halfmin
  }
  
  # Log-transform and scale data
  counts_return <- log1p(counts_return) %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(counts_return))
  
  return(counts_return)
}


##proteomics
preprocess_proteomics_LG <- function(assay, data_env){
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  
  counts <- 2^counts #unlog2 the counts
  counts_return <- counts #copy. This is the one we will change and return
  
  NormalisationFactor <- median(rowSums(counts, na.rm = TRUE)) #calculate row(sample) sums. Take the median of all sums.
  
  #run for loop over each over the rows (samples)
  for(i in 1:nrow(counts)) {
    #get sample i, apply normalisation factor on it. Change sample i in the normalized df.
    counts_return[i,] <- counts[i,] * (NormalisationFactor / sum(counts[i,], na.rm = TRUE))
  }
  #remove samples with more than 95% missing features? remove samples with more 50% missing samples?
  counts_return = counts_return[,apply(!is.na(counts_return),2,sum)>(0.01*nrow(counts_return))]
  #Impute with half min value per Protein
  for(i in 1:ncol(counts_return)) {       # for-loop over columns
    halfmin <- min(counts_return[,i], na.rm = TRUE)/2
    counts_return[,i][is.na(counts_return[,i])] <- halfmin
  }
  
  # Log-transform and scale data
  counts_return <- log1p(counts_return) %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(counts_return))
  
  return(counts_return)
}

##olink
preprocess_serum_olink <- function(assay, data_env) {
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  counts_processed <- counts[rowMeans(is.na(counts)) < 1, ] %>%
    as.matrix() %>%
    impute.knn() %>%
    .$data
  # Scaling the data (log-transformation is not performed since there are negative values)
  counts_processed <- apply(counts_processed, 2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(counts_processed))
  return(counts_processed)
}


#cytof
#' @title preprocesseed blood cytof children populations
#' @description remove debris, red blood cells and undefined populations and normalize on total number of events per sample 
#' @author Slim Fourati
#' @author Brian Lee
#' @author Jingjing Qi
#' @examples \dontrun{bld_cytof_parent_counts_processed <- preprocess_bld_cytof(counts)}
preprocess_bld_cytof <- function(assay, data_env) {
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  output_counts <-  counts %>%
    rownames_to_column(var = "sample_id") %>%
    select(-contains(match = "Debris"), 
           -contains(match = "Multiplet"),
           -`Platelets`, 
           -`RBC`, 
           -`Tier1_Undefined`,
           -`Tier2_undefined_Undefined`) %>%
    column_to_rownames(var = "sample_id")
  output_counts <- output_counts/rowSums(output_counts)
  
  # Log-transform and scale data
  output_counts <- log1p(output_counts) %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(output_counts))
  
  return(value = output_counts)
}

#' @title preprocesseed blood cytof parent populations
#' @description remove debris, red blood cells and undefined populations, summarize parent populations and normalize on total number of events per sample 
#' @author Slim Fourati
#' @author Brian Lee
#' @author Jingjing Qi
#' @examples \dontrun{bld_cytof_parent_counts_processed <- preprocess_bld_cytof_parent(counts, rowfeature)}
preprocess_bld_cytof_parent <- function(assay, data_env) {
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  rowfeature = data_env[[ data_env$assay_files_info[[assay]]$rowfeature_var_name ]]
  bld_cytof_parent_counts_processed <- counts %>%
    rownames_to_column(var = "sample_id") %>%
    pivot_longer(cols = -sample_id, names_to = "rowname") %>%
    merge(y = rownames_to_column(rowfeature), by = "rowname") %>%
    group_by(sample_id, `Broad Parent Population`) %>%
    mutate(`Broad Parent Population` = gsub(pattern = "\\+", replacement = "", `Broad Parent Population`)) %>%
    summarize(value = sum(value)) %>%
    pivot_wider(names_from = `Broad Parent Population`) %>%
    select(-contains(match = "Debris"), 
           -contains(match = "Multiplet"),
           -`Platelets`, 
           -`RBC`, 
           -`Undefined`) %>%
    column_to_rownames(var = "sample_id")
  bld_cytof_parent_counts_processed <-   bld_cytof_parent_counts_processed/rowSums(bld_cytof_parent_counts_processed)
  
  # Log-transform and scale data
  bld_cytof_parent_counts_processed <- log1p(bld_cytof_parent_counts_processed) %>% 
    apply(2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(bld_cytof_parent_counts_processed))
  
  return(value = bld_cytof_parent_counts_processed)
}

preprocess_default <- function(assay, data_env) {
  counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  # Imputation using missForest
  missForest_imputation(data_env$assay_files_info[[assay]], data_env, alpha = data_env$alpha, random_seed = data_env$random_seed, cores = data_env$impute_cores)
  # Return the imputed count data; the imputed count data frame is stored in data_env by missForest_imputation function
  output_counts = data_env[[ data_env$assay_files_info[[assay]]$count_var_name ]]
  
  # Scaling the data (log-transformation is not performed since there could be negative values for some datasets. If log-transformation is needed for a specific assay, create a new function for that assay and log-transform only for that assay)
  output_counts <- apply(output_counts, 2, scale) %>% 
    data.frame(check.names = F) %>% 
    `rownames<-`(row.names(output_counts))
  
  return( output_counts )
}

preprocess_metabolomics <- function(assay, data_env){
  counts = preprocess_default(assay, data_env)
  return(counts)
}

preprocess_nasal_viralload <- function(assay, data_env){
  counts = preprocess_default(assay, data_env)
  return(counts)
}

preprocess_serum_rbd_abtiters <- function(assay, data_env){
  counts = preprocess_default(assay, data_env)
  return(counts)
}

preprocess_serum_sarscov2_abtiters <- function(assay, data_env){
  counts = preprocess_default(assay, data_env)
  return(counts)
}

preprocess_metagenomics <- function(assay, data_env){
  counts = preprocess_default(assay, data_env)
  return(counts)
}





### Main function for data preprocessing
#' @title preprocess data
#' @description preprocess the data using assay specific preprocessing functions.
#' @author Ravi Patel
#' @param data_env The data environment produced by load_IMPACC_datasets / load_IMPACC_Publication_datasets functions.
#' @param alpha Optional, features with "alpha" or more missing values are discarded
#' @param random_seed Optional, random seed used by missForest
#' @param impute_cores Optional, the number of cores to use for imputation.
#' @examples \dontrun{ preprocessed_data <- preprocess_data(data_env) }
preprocess_data <- function(data_env, alpha=0.3, random_seed=2021,  impute_cores=8) {
  # List of preprocessing functions
  preprocess_funcs = list(
    "plasma_proteomics_targeted" = preprocess_proteomics,
    "plasma_proteomics_global_dda" = preprocess_proteomics,
    "serum_olink" = preprocess_serum_olink,
    "nasal_viralload" = preprocess_nasal_viralload,
    "serum_rbd_abtiters" = preprocess_serum_rbd_abtiters,
    "serum_sarscov2_abtiters" = preprocess_serum_sarscov2_abtiters,
    "nasal_transcriptomics" = preprocess_transcriptomics,
    "plasma_metabolomics_global"= preprocess_metabolomics,
    "bld_cytof" = preprocess_bld_cytof,
    #"ea_transcriptomics" = preprocess_transcriptomics,
    "pbmc_transcriptomics" = preprocess_transcriptomics,
    #"ea_metagenomics" = preprocess_metagenomics,
    "nasal_metagenomics" = preprocess_metagenomics
  )
  
  
  # Load required packages
  library(impute)
  
  # Load the imputation-related parameters into data_env for an easy access in the downstream functions.
  data_env$alpha = alpha
  data_env$random_seed = random_seed
  data_env$impute_cores = impute_cores
  
  # Preprocess the data
  preprocessed_data = list()
  for(assay in names(data_env$assay_files_info) ) {
    # Preprocess the data for the assays that have preprocessing functions defined.
    if( ! is.null(preprocess_funcs[[assay]]) ) {
      cat("Preprocessing", assay, "\n")
      # Preprocess the data
      preprocessed_data[[ assay ]] = preprocess_funcs[[assay]](assay, data_env)
    }
  }
  return(preprocessed_data)
}

preprocess_data_LG <- function(data_env, alpha=0.3, random_seed=2021,  impute_cores=8) {
  # List of preprocessing functions
  preprocess_funcs = list(
    "plasma_proteomics_targeted" = preprocess_proteomics_LG,
    "plasma_proteomics_global_dda" = preprocess_proteomics_LG,
    "serum_olink" = preprocess_serum_olink,
    "nasal_viralload" = preprocess_nasal_viralload,
    "serum_rbd_abtiters" = preprocess_serum_rbd_abtiters,
    "serum_sarscov2_abtiters" = preprocess_serum_sarscov2_abtiters,
    "nasal_transcriptomics" = preprocess_transcriptomics,
    "plasma_metabolomics_global"= preprocess_metabolomics,
    "bld_cytof" = preprocess_bld_cytof,
    "ea_transcriptomics" = preprocess_transcriptomics,
    "pbmc_transcriptomics" = preprocess_transcriptomics,
    "ea_metagenomics" = preprocess_metagenomics,
    "nasal_metagenomics" = preprocess_metagenomics
  )
  
  
  # Load required packages
  library(impute)
  
  # Load the imputation-related parameters into data_env for an easy access in the downstream functions.
  data_env$alpha = alpha
  data_env$random_seed = random_seed
  data_env$impute_cores = impute_cores
  
  # Preprocess the data
  preprocessed_data = list()
  for(assay in names(data_env$assay_files_info) ) {
    # Preprocess the data for the assays that have preprocessing functions defined.
    if( ! is.null(preprocess_funcs[[assay]]) ) {
      cat("Preprocessing", assay, "\n")
      # Preprocess the data
      preprocessed_data[[ assay ]] = preprocess_funcs[[assay]](assay, data_env)
    }
  }
  return(preprocessed_data)
}



pairwise_overexpression=function(X, subtype_pred, FWER = 0.01){
  grp_id0 = sort(unique(subtype_pred))
  sigificant_features = list()
  for(g1 in 1:length(grp_id0)){
    sigificant_features[[g1]] = list()
    for(g2 in 1:length(grp_id0)){
      sigificant_features[[g1]][[g2]] =matrix(NA, ncol = 2, nrow = ncol(X))
      if(g1!=g2){
        #g1 > g2
        ll1 = which(subtype_pred==g1)
        ll2 = which(subtype_pred==g2)
        z = c(rep(1,length(ll1)),rep(0,length(ll2)))
        for(j in 1:ncol(X)){
          tmp1 = cor.test(X[c(ll1,ll2),j], z, method = "spearman", alternative = "greater")
          sigificant_features[[g1]][[g2]][j,1] = tmp1$estimate
          sigificant_features[[g1]][[g2]][j,2] = tmp1$p.value
        }
        colnames(sigificant_features[[g1]][[g2]])=c("rho","pval")
      }
    }
  }
  significant_all_top = list()
  for(g1 in 1:length(grp_id0)){
    tmp1 =  sigificant_features[[g1]]
    tmp2 = sapply(tmp1[-g1], function(z) z[,1])
    tmp3 = sapply(tmp1[-g1], function(z) z[,2])
    ll1 =  which(apply(tmp3,1,min)<=FWER/(ncol(X)*length(grp_id0)))
    ll2 = which(apply(tmp2,1,min)>=0)
    ll0 = intersect(ll1,ll2)
    if(length(ll0)>0){
      significant_all_top[[g1]] =ll0
    }else{
      significant_all_top[[g1]] =c(NA)
    }
    
  }
  return(list(sigificant_features=sigificant_features,significant_all_top=significant_all_top))
}
