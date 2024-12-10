library(R6)
## I've commented out to make sure I'm using the version in 0.Codebase instead - Cole
# source("~/impacc/Analysis/data_analysis_template_codebase.R")
# source("~/impacc/Analysis/data-integration/scripts/CleanCode/codebase_basic.R")
# source("~/impacc/Analysis/data-integration/scripts/CleanCode/codebase_mod.R")
# source("~/impacc/Analysis/data-integration/scripts/CleanCode/codebase_prediction.R")
# source("~/impacc/Analysis/data-integration/scripts/CleanCode/codebase_interpretation.R")

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
         "ComplexHeatmap",
         "MOFA2",
         "glmnet",
         "MASS",
         "ordinalNet",
         "ggpubr",
         "dplyr",
         "ggeffects")

for(pkg in pkgs){
  require(pkg, character.only = TRUE)
}
predict = stats::predict

copy_interpretation_class = function(old_class_obj,impacc_analysis_record){
  new_class_obj = impacc_interpretation_pipeline$new(mod_name = old_class_obj$mod_name,
                                                     impacc_analysis_record = impacc_analysis_record,
                                                     omic_names =  old_class_obj$omic_names,
                                                     projection_coef =  old_class_obj$projection_coef, 
                                                     projection_pval =  old_class_obj$projection_pval,
                                                     regression_coef =  old_class_obj$regression_coef)
  nams = intersect(names(old_class_obj),names(new_class_obj))
  for(inam in 1:length(nams)){
    nam = nams[inam]
    if(!("function" %in% class(old_class_obj[[nam]]))){
      new_class_obj[[nam]] = old_class_obj[[nam]]
    }
  }
  return(new_class_obj)
}


copy_computational_class = function(impacc_analysis_record){
  new_class_obj =  impacc_analysis$new(rank_max=impacc_analysis_record$rank_max, rank_grids=impacc_analysis_record$rank_grids, 
                                       nest_cv_folds= impacc_analysis_record$nest_cv_folds,
                                       inner_loop_folds=impacc_analysis_record$inner_loop_folds,
                                       random_seed=impacc_analysis_record$random_seed)
  nams = intersect(names(impacc_analysis_record),names(new_class_obj))
  for(inam in 1:length(nams)){
    nam = nams[inam]
    if(!("function" %in% class(impacc_analysis_record[[nam]]))){
      new_class_obj[[nam]] = impacc_analysis_record[[nam]]
    }
  }
  return(new_class_obj)
}
#Update
###LG: 20221021: add finialized marginal test results for modules
impacc_analysis = R6Class(classname = "impacc_computational_pipeline",
                          public =list(
                            rank_max = 40,
                            rank_grids = c(c(2:10),15, 20, 30, 40),
                            rank_selected = NULL,
                            deviance_loss = NULL,
                            model_mcia_list = NULL,
                            model_mcia = NULL,
                            model_clin = NULL,
                            model_ensemble = NULL,
                            random_seed = 2022,
                            nest_cv_folds = 50,
                            inner_loop_folds = 5,
                            nest_cv_foldid = NULL,
                            MofaOutPut_trainPath =NULL,
                            MofaOutPut_testPath = NULL,
                            MCIA_block_prep='lambda_all',
                            MCIA_deflat_method = "globalScore",
                            mcia_model.train = NULL,
                            x.mcia.train = NULL,
                            x.mcia.test = NULL,
                            clin.train = NULL,
                            clin.test = NULL,
                            mcia.factors.train = NULL,
                            mcia.factors.test = NULL,
                            mcia.visit1.test_plot = NULL,
                            mcia.trajectory.test_plot = NULL,
                            projection_evaluation_train = NULL,
                            layer2_level = NULL,
                            y = NULL,
                            deviance_loss_vec_cv = NULL,
                            Cor_glob_vec_cv = NULL,
                            Cor_TG45_vec_cv = NULL,
                            pred_train_ids = NULL,
                            other_model_list = list(),
                            test_prediction_list = list(),
                            
                            initialize=function(rank_max=40, rank_grids=c(c(2:10),15, 20, 30, 40), 
                                                nest_cv_folds= 50,
                                                inner_loop_folds=10,
                                                random_seed=2022){
                              self$rank_max = rank_max
                              self$rank_grids=rank_grids
                              self$nest_cv_folds = nest_cv_folds
                              self$inner_loop_folds = inner_loop_folds
                              self$random_seed = random_seed
                            },
                            mcia_model_construction = function(MofaOutPut_trainPath, MCIA_block_prep='lambda_all',
                                                         MCIA_deflat_method = "globalScore"){
                              self$MofaOutPut_trainPath =MofaOutPut_trainPath
                              self$MCIA_block_prep=MCIA_block_prep
                              self$MCIA_deflat_method = MCIA_deflat_method
                              imputed_data = readRDS(self$MofaOutPut_trainPath)
                              self$clin.train = imputed_data$clinical
                              x.mcia.train = list()
                              D = length(imputed_data$datasets_imputed)
                              for(d in 1:D){
                                x.mcia.train[[d]] = t(imputed_data$datasets_imputed[[d]])
                              }
                              x.mcia.train =nsc_prep(data_input=x.mcia.train, self$rank_max)
                              self$x.mcia.train<-processOpt(x.mcia.train,scale=FALSE,center=FALSE,self$rank_max,option=self$MCIA_block_prep)
                              self$mcia_model.train<-mbpca(self$x.mcia.train,ncomp=self$rank_max,k="all",method=self$MCIA_deflat_method,
                                                            option="uniform",center=FALSE,scale=FALSE,
                                                            unit.p=TRUE,unit.obs=TRUE,moa=FALSE)
                              gs=new_gs(self$x.mcia.train,self$mcia_model.train)
                              ###if in consistency, trucate
                              cors =  diag(cor(gs,self$mcia_model.train$t))
                              idx = which(cors<=0.99)
                              if(length(idx)>0){
                                idx = idx[1]
                                self$rank_max = idx-1
                              }
                              self$mcia.factors.train = gs[,1:self$rank_max]
                              self$projection_evaluation_train = projection_coef_func(
                                mcia_factors=self$mcia.factors.train,
                                datasets=imputed_data$datasets_imputed[1:D])
                            },
                            
                            create_nested_foldid = function(N){
                              set.seed(self$random_seed)
                              self$nest_cv_foldid= sample(rep(1:self$nest_cv_folds,ceiling(N/self$nest_cv_folds)),N)
                            },
                            run_prediction_composite = function(x, y, layer2_level = NULL){
                              if(!is.null(layer2_level)){
                                self$layer2_level = layer2_level
                              }
                              if(is.null(self$layer2_level)){
                                stop("no layer2_level provided!")
                              }
                              if(length(self$nest_cv_foldid)!=length(y)){
                                stop("fold assignment does not match the sample size!")
                              }
                              if(is.null(colnames(x))){
                                stop("x does not have column names!")
                              }
                              if(sum(colnames(x) !=colnames(data.frame(x)))){
                                stop("improper column names for x: data.frame(x) modifies the original colnames!")
                              }
                              tmp1 = sort(unlist( self$layer2_level))
                              tmp2 = sort(unique(y))
                              if(length(setdiff(tmp1, tmp2)) != 0 | length(setdiff(tmp2, tmp1)) != 0){
                                stop("layer2_level does not match the labels in y!")
                              }
                              fitted=nested_cv_prediction_composite(X=x, Y=y, module_names =colnames(x), 
                                                                       nest_cv_foldid = self$nest_cv_foldid, 
                                                                       nest_cv_folds = self$nest_cv_folds, 
                                                                       nested_cv_seed =self$random_seed, nfolds_train =self$inner_loop_folds,
                                                                       train_seed =self$random_seed,
                                                                       layer2_level = self$layer2_level)
                              return(fitted)
                            },
                            
                            run_mcia_varying_rank_composite = function(row_ids){
                              self$pred_train_ids=row_ids
                              mcia_use = self$mcia.factors.train
                              sample_meta_data = self$clin.train 
                              x_all = sapply(1:ncol(mcia_use[row_ids,]), function(i) qqnorm(mcia_use[row_ids,i],plot.it=F)$x)
                              module_names = paste0("mod",1:ncol(x_all))
                              colnames(x_all) = module_names
                              y = as.integer(sample_meta_data$trajectory_group[row_ids])
                              N = length(y)
                              self$create_nested_foldid(N)
                              y[y>=3] =y[y>=3]-1
                              self$y = y
                              layer2_level = list(c(1),c(2), c(3,4))
                              self$layer2_level = layer2_level
                              self$model_mcia_list = list()
                              for(irank in 1:length(self$rank_grids)){
                                rank_cur = self$rank_grids[irank]
                                print(rank_cur)
                                x = x_all[,1:rank_cur]
                                self$model_mcia_list[[irank]] = nested_cv_prediction_composite(X=x, Y=y, module_names =colnames(x), 
                                                                                               nest_cv_foldid = self$nest_cv_foldid, 
                                                                                               nest_cv_folds = self$nest_cv_folds, 
                                                                                               nested_cv_seed =self$random_seed, nfolds_train =self$inner_loop_folds,
                                                                                               train_seed =self$random_seed+1,
                                                                                               layer2_level = self$layer2_level)
                              }
                            },
                            
                            mcia.test_factor_construction = function(MofaOutPut_testPath){
                              self$MofaOutPut_testPath =MofaOutPut_testPath
                              imputed_data = readRDS(self$MofaOutPut_testPath)
                              self$clin.test = imputed_data$clinical
                              x.mcia.train = list()
                              D = length(imputed_data$datasets_imputed)
                              for(d in 1:D){
                                x.mcia.train[[d]] = t(imputed_data$datasets_imputed[[d]])
                              }
                              x.mcia.train =nsc_prep(data_input=x.mcia.train, self$rank_max)
                              self$x.mcia.test<-processOpt(x.mcia.train,scale=FALSE,center=FALSE,self$rank_max,option=self$MCIA_block_prep)
                              gs=new_gs(self$x.mcia.test,self$mcia_model.train)
                              self.mcia.factors.test = gs[,1:self$rank_max]
                            },
                            
                            evaluation_varying_rank = function(model_list=NULL){
                              if(is.null(model_list)){
                                model_list=self$model_mcia_list
                              }
                              self$deviance_loss_vec_cv = rep(NA,  length(model_list))
                              self$Cor_glob_vec_cv = rep(NA, length(model_list))
                              self$Cor_TG45_vec_cv = rep(NA, length(model_list))
                              for(irank in 1:length(model_list)){
                                self$deviance_loss_vec_cv[irank]=mean(deviance_loss_multinomial(y =self$y, probs =   model_list[[irank]]$nested_cv_probs))
                                y1 = self$y
                                #merge to create 1 vs (23) vs (45)
                                y1[y1==4]=y1[y1==4]-1
                                self$Cor_glob_vec_cv[irank] = cor(1-model_list[[irank]]$nested_cv_probs[,1],y1,method = "spearman")
                                # compare TG4 and TG5
                                ll1 = which(self$y%in%c(3,4))
                                y1 = self$y[ll1]
                                tmp1 = model_list[[irank]]$nested_cv_probs
                                self$Cor_TG45_vec_cv[irank] = cor(tmp1[ll1,4]/( tmp1[ll1,3]+tmp1[ll1,4]),y1,method = "spearman")
                              }
                              self$rank_selected = self$rank_grids[which.min(self$deviance_loss_vec_cv)]
                              self$model_mcia = self$model_mcia_list[[which.min(self$deviance_loss_vec_cv)]]
                             
                            },
                            
                            model_ensemble_construct = function(z1, z2, y){
                              dat = z1
                              dat$y = factor(y)
                              model_ensemble = list()
                              model_ensemble$layer1 = MASS::polr(y~. , data = dat)
                              idx = which(y%in%c(3,4))
                              y2 = y[idx]-3
                              z2 = z2[idx,]
                              model_ensemble$layer2 = glmnet(x = z2, y = y2, family="binomial", lower.limits = 0, lambda = 1e-5)
                              return(model_ensemble)
                            },
                            make_prediction=function(xnew, model){
                              dat.te = data.frame(xnew)
                              probs.test = prediction_hierachical(dat.te = dat.te, module_names = colnames(dat.te),
                                                                    layer1model=model$layer1model, 
                                                                    layer2model = model$layer2model,
                                                                    layer2_level = model$layer2_level)
                              return(probs.test)
                            },
                            make_ensemble_prediction = function(z1, z2, model){
                              layer2_level = self$layer2_level
                              ii = sapply(layer2_level, length)
                              ii1 = which(ii>1)
                              K2 = length((ii1))
                              prob1 = predict(model$layer1, newdata = z1, type = "probs")
                              prob2 = predict(model$layer2, newx = as.matrix(z2), type = "response")
                              pred.final =  pred.final =  matrix(0, nrow = nrow(z1), ncol = length(unlist(layer2_level)))
                              if(K2 > 0){
                                for(k in 1:length(layer2_level)){
                                  if(length(layer2_level[[k]])==1){
                                    pred.final[,layer2_level[[k]]] = prob1[,layer2_level[[k]]]
                                  }else{
                                    pred.final[,layer2_level[[k]]] = apply(prob1[,layer2_level[[k]],drop = F],1,sum)
                                  }
                                }
                                for(k in 1:K2){
                                  pred.final[,layer2_level[[ii1[k]]][1]] = pred.final[,layer2_level[[ii1[k]]][1]]*(1-prob2)
                                  pred.final[,layer2_level[[ii1[k]]][2]] = pred.final[,layer2_level[[ii1[k]]][2]]*prob2
                                }
                              }else{
                                pred.final = pred.layer1
                              }
                              return(pred.final)
                            },
                            coef_extract_composite = function(model, feature_names=NULL, method = "MCIA"){
                              idx = which.min(model$layer1model$aic)
                              if(is.null(feature_names)){
                                feature_names=rownames(coef(model$layer2model[[1]], s = "lambda.min"))[-1]
                              }
                              coef_glob = -model$layer1model$coefs[idx,feature_names]
                              coef_4vs5 = coef(model$layer2model[[1]], s = "lambda.min")[-1]
                              names(coef_4vs5) =feature_names
                              coef_mat = data.frame(ordinal= coef_glob, mortality.in.severe=coef_4vs5)
                              coef_mat$Feature = feature_names
                              coef_mat = gather(coef_mat, key = "model", value = "coef", -Feature)
                              coef_mat$method = method
                              return(coef_mat)
                            }

                            
                          ),
                          private = list()
                          )


run_impacc_analysis = function(MofaOutPut_trainPath,
                               MofaOutPut_testPath,
                               Clin_trainPath,
                               base_file_impacc_analysis,
                               file_impacc_analysis,
                               rank_max=40, 
                               rank_grids=c(c(2:10),15, 20, 30, 40),
                               nest_cv_folds= 50, 
                               inner_loop_folds = 5,
                               random_seed=2022){
  impacc_analysis_record_new = impacc_analysis$new(rank_max=rank_max, rank_grids=rank_grids, 
                                                   nest_cv_folds= nest_cv_folds,
                                                   inner_loop_folds=inner_loop_folds,
                                                   random_seed=random_seed)
  if(!file.exists(file_impacc_analysis)){
    
    if(file.exists(base_file_impacc_analysis)){
      impacc_analysis_record_base = readRDS(base_file_impacc_analysis)
      impacc_analysis_record_new$MofaOutPut_trainPath =impacc_analysis_record_base$MofaOutPut_trainPath
      impacc_analysis_record_new$clin.train=impacc_analysis_record_base$clin.train
      impacc_analysis_record_new$MCIA_block_prep=impacc_analysis_record_base$MCIA_block_prep
      impacc_analysis_record_new$MCIA_deflat_method =impacc_analysis_record_base$MCIA_deflat_method
      impacc_analysis_record_new$x.mcia.train =impacc_analysis_record_base$x.mcia.train
      impacc_analysis_record_new$mcia_model.train =impacc_analysis_record_base$mcia_model.train
      impacc_analysis_record_new$mcia.factors.train=impacc_analysis_record_base$mcia.factors.train
      impacc_analysis_record_new$projection_evaluation_train=impacc_analysis_record_base$projection_evaluation_train
      impacc_analysis_record_new$MofaOutPut_testPath=impacc_analysis_record_base$MofaOutPut_testPath
      impacc_analysis_record_new$clin.test=impacc_analysis_record_base$clin.test
      impacc_analysis_record_new$x.mcia.test=impacc_analysis_record_base$x.mcia.test
      impacc_analysis_record_new$mcia.factors.test=impacc_analysis_record_base$mcia.factors.test
    }else{
      impacc_analysis_record_new$mcia_model_construction(MofaOutPut_trainPath=MofaOutPut_trainPath)
      impacc_analysis_record_new$mcia.test_factor_construction(MofaOutPut_testPath=MofaOutPut_testPath)
    }
    row_visit1<- which(impacc_analysis_record_new$clin.train$event_type =="Visit 1")
    impacc_analysis_record_new$run_mcia_varying_rank_composite(row_ids = row_visit1)
    impacc_analysis_record_new$evaluation_varying_rank(impacc_analysis_record_new$model_mcia_list)
    ##clin model
    y =impacc_analysis_record_new$y
    x_clin = read.csv(file = Clin_trainPath)[,-1]
    impacc_analysis_record_new$model_clin = impacc_analysis_record_new$run_prediction_composite(x=x_clin, y)
    saveRDS(impacc_analysis_record_new, file = file_impacc_analysis)
  }else{
    impacc_analysis_record = readRDS(file_impacc_analysis)
    nams = names(impacc_analysis_record)
    for(inam in 1:length(nams)){
      nam = nams[inam]
      if(!("function" %in% class(impacc_analysis_record[[nam]]))){
        impacc_analysis_record_new[[nam]] = impacc_analysis_record[[nam]]
      }
    }
  }
  return(impacc_analysis_record_new)
  
}

#' This function helps us to perform the non-essential rearrangement
#' 1. Flip the module directions based on signs
organize_impacc_analysis = function(impacc_analysis_record, 
                                    sign_flips = NULL){
  if(!is.null(sign_flips)){
    for(j in 1:length(sign_flips)){
      if(sign_flips[j]){
        for(l in 1:length(impacc_analysis_record$mcia_model.train$pb)){
          #flip the signs in MCIA construction
          impacc_analysis_record$mcia_model.train$pb[[l]][,j] = -impacc_analysis_record$mcia_model.train$pb[[l]][,j]
          #flip the signs impacc_analysis_record$projection_evaluation_train$projection_coef[[l]][,j] 
          impacc_analysis_record$projection_evaluation_train$projection_coef[[l]][,j] = -impacc_analysis_record$projection_evaluation_train$projection_coef[[l]][,j]
          
        }
        impacc_analysis_record$mcia.factors.test[,j] = -impacc_analysis_record$mcia.factors.test[,j]
        impacc_analysis_record$mcia.factors.train[,j] = -impacc_analysis_record$mcia.factors.train[,j]
      }
    }
    ####rerun the full prediction model
    for(j in 1:length(impacc_analysis_record$model_mcia_list)){
      tmp = impacc_analysis_record$model_mcia_list[[j]]
      rank1 = impacc_analysis_record$rank_grids[j]
      mods = paste0("mod",1:rank1)
      model = tmp$model_full
      sign_flips1 =ifelse(sign_flips, -1,1)
      if(rank1 > length(sign_flips)){
        sign_flips1 = c(sign_flips1, rep(1, rank1-length(sign_flips)))
      }else{
        sign_flips1 = sign_flips1[1:rank1]
      }
      model$layer1model$coefs[,paste0("mod",1:rank1)] = t(t(model$layer1model$coefs[,paste0("mod",1:rank1)]) * sign_flips1)
      model$layer2model[[1]]$glmnet.fit$beta = model$layer2model[[1]]$glmnet.fit$beta*sign_flips1
      impacc_analysis_record$model_mcia_list[[j]]$model_full = model
    }
    impacc_analysis_record$model_mcia = impacc_analysis_record$model_mcia_list[[which.min(impacc_analysis_record$deviance_loss_vec_cv)]]
  }
  return(impacc_analysis_record)
}


impacc_analysis_object_dropdata = function(impacc_analysis_record){
  impacc_analysis_record$clin.train = NA
  impacc_analysis_record$clin.test= NA
  impacc_analysis_record$mcia.factors.train = NA
  impacc_analysis_record$mcia.factors.test = NA
  return(impacc_analysis_record)
  
}

impacc_interpretation_pipeline = R6Class(classname = "impacc_mod_interpretation",
                                         lock_objects = FALSE,
                                public = list(
                                  mod_name = NULL,
                                  clin.train = NULL,
                                  clin.test = NULL,
                                  mcia.factors.train = NULL,
                                  mcia.factors.test = NULL,
                                  MofaOutPut_trainPath = NULL,
                                  MofaOutPut_testPath = NULL,
                                  projection_coef = NULL,
                                  projection_pvalue = NULL,
                                  regression_coef = NULL,
                                  omic_names = NULL,
                                  assays.train = NULL,
                                  assays.test = NULL,
                                  short_list_features = list(),
                                  others = list(),
                                  composite_features_train = list(),
                                  original_features_train = list(),
                                  composite_features_test = list(),
                                  original_features_test = list(),
                                  composite_features_members = list(),
                                  MODtrajectory_tests = NULL,
                                  MODvisit1_tests = NULL,
                                  composite_trajectory_tests = list(),
                                  composite_visit_tests = list(),
                                  params = list(),
                                  MODgseaRes = list(),
                                  MODmghRes = list(),
                                  MODenrichrRes = list(),
                                  Metabo_Cytokine = list(),
                                  Metabo_Gene = list(),
                                  Cytokine_Gene = list(),
                                  initialize=function(mod_name,
                                                      impacc_analysis_record,
                                                      omic_names,
                                                      projection_coef, 
                                                      projection_pval,
                                                      regression_coef
                                                      ){
                                    self$mod_name = mod_name
                                    self$clin.train = impacc_analysis_record$clin.train
                                    self$clin.test= impacc_analysis_record$clin.test
                                    self$omic_names = omic_names
                                    self$projection_coef = projection_coef
                                    self$projection_pvalue = projection_pval
                                    self$regression_coef = regression_coef
                                    self$mcia.factors.train = impacc_analysis_record$mcia.factors.train
                                    self$mcia.factors.test = impacc_analysis_record$mcia.factors.test
                                    imputed_data = readRDS(impacc_analysis_record$MofaOutPut_trainPath)
                                    self$assays.train = imputed_data$datasets
                                    imputed_data = readRDS(impacc_analysis_record$MofaOutPut_testPath)
                                    self$assays.test= imputed_data$datasets
                                    self$MofaOutPut_trainPath = impacc_analysis_record$MofaOutPut_trainPath
                                    self$MofaOutPut_testPath = impacc_analysis_record$MofaOutPut_testPath
                                    
                                  },
                                  coefficient_name_translator = function(data_env){
                                    for(j in 1:length(self$omic_names)){
                                      omic_name = self$omic_names[j]
                                      if(omic_name %in% c("pbmc_transcriptomics","nasal_transcriptomics", "plasma_metabolomics_global")){
                                        self$projection_coef[[j]] = name_translator(data_env = data_env,
                                                                                    DF = self$projection_coef[[j]], 
                                                                                    omic_name=omic_name,
                                                                                    mod = "row")
                                        self$projection_pvalue[[j]] = name_translator(data_env = data_env,
                                                                                    DF = self$projection_pvalue[[j]], 
                                                                                    omic_name=omic_name,
                                                                                    mod = "row")
                                      }
                                    }
                                  },
                                  assay_name_translator = function(data_env){
                                    for(j in 1:length(self$omic_names)){
                                      omic_name = self$omic_names[j]
                                      if(omic_name %in% c("pbmc_transcriptomics","nasal_transcriptomics", "plasma_metabolomics_global")){
                                        self$assays.train[[j]] = name_translator(data_env = data_env,
                                                                                    DF = self$assays.train[[j]], 
                                                                                    omic_name=omic_name,
                                                                                    mod = "col")
                                        self$assays.test[[j]] = name_translator(data_env = data_env,
                                                                                 DF = self$assays.test[[j]], 
                                                                                 omic_name=omic_name,
                                                                                 mod = "col")
                                      }
                                    }
                                  },
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
                                  },
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
                                  },
                                  runEnrichr = function(short_list_coef, Term2GeneList, backround = "impacc"){
                                    #backround = "impacc": use our full list as background
                                    #backround = "term": all features from term as the background
                                    enrichr_result = list()
                                    for(j in 1:length(short_list_coef)){
                                      if(names(short_list_coef)[j] %in% names(Term2GeneList)){
                                        geneList = rownames(short_list_coef[[j]])
                                        Term2Gene = Term2GeneList[[names(short_list_coef)[j]]]
                                        universe = NULL
                                        if(backround == "impacc"){
                                          universe = rownames(self$projection_coef[[j]])
                                        }else{
                                          universe = unique(Term2Gene[,2])
                                        }
                                        enrichr_result[[names(short_list_coef)[j]]] = clusterProfiler::enricher(gene = geneList,
                                                                                                         minGSSize = 2,
                                                                                                         maxGSSize = 500,
                                                                                                         pvalueCutoff = 2,
                                                                                                         universe = universe,
                                                                                                         TERM2GENE = Term2Gene)
                                      }
                                    }
                                    return(enrichr_result)
                                  },
                                  #Term-TermSize-LeadingCount-LeadingEdges -mimic GSEA object
                                  runMGH = function(short_list_coef, Term2GeneList, backround = "impacc"){
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
                                          j0 = which(names(self$projection_coef)==names(short_list_coef)[j])
                                          universe = rownames(self$projection_coef[[j0]])
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
                                  },
                                  composite_features_creation = function(result_table_list,
                                                                          term,
                                                                          core_enrichment_name = "core_enrichment",
                                                                          phase = "train",
                                                                          sep = "#"
                                                                          ){
                                    composite_expressions = list()
                                    for(j in 1:length(self$projection_coef)){
                                      name_j = names(self$projection_coef)[j]
                                      if(name_j %in% names(result_table_list)){
                                        result_table_DF = result_table_list[[name_j]]
                                        if(nrow(result_table_DF)>0){
                                          idx = (result_table_DF$Description==term)
                                          if(sum(idx) > 0){
                                            leading_edges = result_table_DF[[core_enrichment_name]][idx]
                                            if(leading_edges != ""){
                                              if(phase == "train"){
                                                omic_assay_DF = self$assays.train[[self$omic_names[name_j]]]
                                              }else{
                                                omic_assay_DF = self$assays.test[[self$omic_names[name_j]]]
                                              }
                                              projection_coef_DF=self$projection_coef[[name_j]]
                                              composite_expressions[[name_j]] = comoposite_features_creation_unit(omic_assay_DF, 
                                                                                                                  projection_coef_DF,
                                                                                                                  leading_edges,
                                                                                                                  sep = sep,
                                                                                                                  is_vector = FALSE)
                                            }

                                          }
                                        }
                                      }
                                    }
                                    return(composite_expressions)
                                    
                                  },
                                  box_plot = function(
                                    plot_name,
                                    expr_list,
                                    phase = "train",
                                    baseline_control = FALSE
                                  ){
                                    colors = c("#639A21", "#39828C",  "#6371AD","#BD7D31", "#9C3418")
                                    if(phase == "train"){
                                      clinDF = self$clin.train
                                    }else{
                                      clinDF = self$clin.test
                                    }
                                    visit_ids = c(1)
                                    plotDF = data.frame(visit = clinDF$event_type,
                                                        TG = clinDF$trajectory_group,
                                                        baseline = clinDF$resp_status_v1)
                                    for(j in 1:length(expr_list)){
                                      name_j = names(expr_list)[j]
                                      plotDF[[name_j]] = expr_list[[name_j]][,1]
                                    }
                                    plotDF = plotDF[plotDF$visit%in%c(paste0("Visit ", visit_ids)),]
                                    if(baseline_control){
                                      plotDF=plotDF[!is.na(plotDF$baseline),]
                                    }
                                    plotDF[names(expr_list)] = apply(plotDF[names(expr_list)],2,function(z){
                                      z1 = z[!is.na(z)]
                                      z1 = qqnorm(z1,plot.it=F)$x
                                      if(baseline_control){
                                        x1 = model.matrix(~plotDF$baseline-1)
                                        x1 = x1[!is.na(z),]
                                        z1  = lm(z1~x1)$residuals
                                      }
                                      z[!is.na(z)] = z1
                                      return(z)
                                    })
                                    plotDF$TG = factor(plotDF$TG)
                                    plotDF$baseline = factor(plotDF$baseline)
                                    plotDF = tidyr::gather(plotDF, key = "omics",
                                                           value = "expr",
                                                           -visit, -TG, -baseline)
                                    plot_name1 = plot_name
                                    if(baseline_control){
                                      plot_name1 = paste0(plot_name1,"|baseline")
                                    }
                                    plot_name1 = paste0(plot_name1, ",", phase)
                                    p1 = ggboxplot(plotDF,
                                                   x="TG", y="expr", 
                                                   color = "TG",
                                                   palette =colors,
                                                   add = "jitter") +
                                      facet_wrap(~ omics)+ggtitle(label = plot_name1)
                                    if(visit_ids==1){
                                      my_comparisons1 = list(c("4", "5"))
                                      p1 = p1 + stat_compare_means()+stat_compare_means(comparisons =my_comparisons1,
                                                                                        label.y = c(2))
                                      
                                    }
                                    return(p1)
                                    
                                    
                                  },
                                  heatmap_plot = function(){
                                    
                                  },
                                  trajectory_test = function(expr_list,
                                                             phase = "train",
                                                             event_date = "admission",
                                                             max_date = 28){
                                    if(phase == "train"){
                                      clinDF = self$clin.train
                                    }else{
                                      clinDF = self$clin.test
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
                                    }else{
                                      plotDF$event_date =  plotDF$DFSO
                                    }
                                    for(j in 1:length(expr_list)){
                                      name_j = names(expr_list)[j]
                                      plotDF[[name_j]] = expr_list[[name_j]][,1]
                                    }
                                    plotDF=plotDF[!is.na(plotDF$trajectory_group) & !is.na(plotDF$event_date),]
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
                                    plotDF = tidyr::gather(plotDF, key = "name",
                                                                 value = "value",
                                                                 -visit, -trajectory_group, -baseline, -ADdate, -DFSO,
                                                                 -participant_id, -enrollment_site,-sex,
                                                                 -discretized_admit_age_quantile, -event_date)
                                    
                                    smooth_model_loop <- model_loop(plotDF, age_sex=T, modelType = "smoothSpline",
                                                                    endpoint = "trajectory_group")
                                    
                                    return(smooth_model_loop)
                                  },
                                  trajectory_plot = function(expr_list,
                                                             model_DF,
                                                             max_date = 28,
                                                             phase = "train",
                                                             event_date = "admission",
                                                             p_adjust = F,
                                                             signif_markers = T, knot_lines = F,  remove_NS = T,
                                                             remove_NS_threshold = 0.05,knots =  c(1, 4, 7, 14, 21),
                                                             individual_points = T, age_sex = F,group_trendline = T, individual_trendlines = T,
                                                             term = NULL
                                                             ){
                                    if(phase == "train"){
                                      clinDF = self$clin.train
                                    }else{
                                      clinDF = self$clin.test
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
                                    plot_list = list()
                                    for(j in 1:length(expr_list)){
                                      name_j = names(expr_list)[j]
                                      plotExample <- plotDF[plotDF$name == name_j,] ## plotExample would be a plotDF object example in 
                                      if(is.null(term)){
                                        title = paste0(name_j,"(",phase,")",":shape/average")
                                      }else{
                                        title = paste0(term,"-",name_j,"(",phase,")",":shape/average")
                                      }
                                      
                                      plot_list[[name_j]] = plot_model(plotExample, model_loop = model_DF, modelType = "smoothSpline",
                                                                  endpoint = "trajectory_group", 
                                                                  p_adjust = p_adjust,
                                                                  knot_lines = knot_lines,
                                                                  signif_markers = signif_markers, 
                                                                  remove_NS = remove_NS,
                                                                  remove_NS_threshold = remove_NS_threshold,
                                                                  knots = knots,
                                                                  individual_points = individual_points, 
                                                                  age_sex = age_sex,
                                                                  group_trendline = group_trendline, 
                                                                  individual_trendlines = individual_trendlines,
                                                                  xlabel = xlabel,
                                                                  title = title,
                                                                  CI = T) 
                                    }
                                    
                                    return(plot_list)
                                  },
                                  condensed_trajectory_plot = function(expr_list,
                                                                       model_DF,
                                                                       functional_pvalue, 
                                                                       functional_sign,
                                                                       
                                                                       max_date = 28,
                                                                       phase = "train",
                                                                       event_date = "admission", cex_all = 15){
                                    if(phase == "train"){
                                      clinDF = self$clin.train
                                    }else{
                                      clinDF = self$clin.test
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
                                      functional_sign0 = functional_sign[name_j]
                                      if(functional_sign0==1){
                                        functional_sign0="+"
                                      }else if(functional_sign0 == -1){
                                        functional_sign0="-"
                                      }else{
                                        functional_sign0='unsigned'
                                      }
                                      functional_pval0 =formatC(functional_pvalue[name_j], format = "e", digits = digit)
                                      text0 = paste0( functional_sign0," ", functional_pval0)
                                      tmp_list[[j]] = tmp_list[[j]]+annotate("text", x = 2, y = -2, label=text0, size= s)
                                      tmp_list[[j]] = tmp_list[[j]]+annotate("text", x=max_date/2,y=2, label= text1,  size=s) + annotate("text", x=max_date/2, y=1.6, label =text2,  size=s)
                                    }
                                    plot_join = ggpubr::ggarrange(plotlist = tmp_list, nrow=1, common.legend = TRUE, legend="bottom")
                                    
                                    plot_join=annotate_figure(plot_join, top = text_grob(term, 
                                                                                         color = "black", face = "bold", size = 12))
                                    
                                    plot_join
                                    
                                    return(plot_join)
                                  }
                                )
)


replace_rownames_with_event_id <- function(assay, clinical_data, event_ids_selected){
  tmp <- merge(assay, clinical_data[,c("sample_id", "event_id")], by.x = 0, by.y = "sample_id", all.x = TRUE)
  rownames(tmp) <- tmp$event_id
  tmp <- tmp[,colnames(assay)]
  if(!missing(event_ids_selected)){
    tmp <- tmp[rownames(tmp) %in% event_ids_selected, ]
  }
  return(tmp)
}
