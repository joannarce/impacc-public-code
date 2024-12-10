##########

ordinalreg_l1selection = function(dat.tr, module_names, response_name = "traj",  nfolds = NULL){
  model <- ordinalNet::ordinalNet(as.matrix(dat.tr[module_names]), y =dat.tr[[response_name]] , family = "cumulative",link = "logit", 
                                  parallelTerms = TRUE, nonparallelTerms = FALSE, standardize = F, nLambda = 100, lambdaMinRatio=1e-4)
  
  return(model)
}

hireg = function(dat.tr, module_names, response_name = "traj", nfolds = 20,
                   layer2_level = list(c(1),c(2), c(3),c(4,5)), seed = 2022, 
                   update.response = TRUE, ordinal= T, alpha = 1){
  #global ordinal model
  set.seed(seed)
  model_ordinal=ordinalreg_l1selection(dat.tr = dat.tr, module_names = module_names, response_name = response_name, nfolds = nfolds)
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

prediction_hierachical = function(dat.te, module_names,layer1model,  layer2model, layer2_level = list(c(1),c(2,3), c(4,5))){
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


deviance_loss_multinomial = function(y, probs, eps = 1e-10){
  K = length(unique(y))
  y.extended = matrix(0, ncol  = K, nrow = length(y))
  for(k in 1:K){
    y.extended[y==k,k] = 1
  }
  apply(-2*log(probs+eps)*y.extended,1,sum)
  
}

nested_cv_prediction_composite = function(X, Y, module_names,  layer2_level = list(c(1),c(2), c(3), c(4,5)),
                                                nest_cv_foldid = NULL, nest_cv_folds = 10, nested_cv_seed = 2022, nfolds_train = 10, train_seed = 2023){
  n = length(Y)
  dat = data.frame(X)
  dat$traj = factor(Y)
  if(is.null(nest_cv_foldid)){
    nest_cv_foldid= sample(rep(1:nfolds,ceiling(n/nest_cv_folds)),n)
  }else{
    nest_cv_folds = length(unique(nest_cv_foldid))
  }
  pred.test = matrix(0, nrow = n, ncol = length(unique(dat$traj)))
  for(fold_id in 1:nest_cv_folds){
    print(fold_id)
    dat.tr =dat[nest_cv_foldid!=fold_id,]
    dat.te =dat[nest_cv_foldid==fold_id,]
    hireg_models = hireg(dat.tr = dat.tr, module_names = module_names, response_name = "traj", nfolds =nfolds_train, layer2_level = layer2_level, seed = train_seed)
    pred.hi = prediction_hierachical(dat.te = dat.te, module_names = module_names,
                                       layer1model=hireg_models$layer1model, 
                                       layer2model = hireg_models$layer2model,
                                       layer2_level = hireg_models$layer2_level)
    pred.test[which(nest_cv_foldid==fold_id),] =  pred.hi
  }
  ##full model
  model_full = hireg(dat.tr = dat, module_names = module_names, response_name = "traj", nfolds =nfolds_train, layer2_level = layer2_level, seed = train_seed)
  return(list(model_full = model_full, nested_cv_probs = pred.test))
}



nested_cv_prediction_multinomial = function(X, Y, module_names, nest_cv_foldid = NULL, nest_cv_folds = 10, nested_cv_seed = 2022, nfolds_train = 10, train_seed = 2023){
  n = length(Y)
  Y = as.integer(Y)
  dat = data.frame(X)
  dat$traj = factor(Y)
  if(is.null(nest_cv_foldid)){
    nest_cv_foldid= sample(rep(1:nfolds,ceiling(n/nest_cv_folds)),n)
  }else{
    nest_cv_folds = length(unique(nest_cv_foldid))
  }
  pred.test = matrix(0, nrow = n, ncol = length(unique(dat$traj)))
  for(fold_id in 1:nest_cv_folds){
    print(fold_id)
    X.tr = X[nest_cv_foldid!=fold_id,]
    Y.tr = Y[nest_cv_foldid!=fold_id]
    X.te = X[nest_cv_foldid==fold_id,]
    dat.tr =dat[nest_cv_foldid!=fold_id,]
    dat.te =dat[nest_cv_foldid==fold_id,]
    set.seed(train_seed)
    lasso_model = cv.glmnet(x = X.tr, y = Y.tr, family = "multinomial", nfolds = nfolds_train, standardize = F)
    pred = predict(lasso_model,newx = X.te, s = "lambda.min", type = "response")
    pred.test[which(nest_cv_foldid==fold_id),] =   pred[,,1]
  }
  ##full model
  set.seed(train_seed)
  model_full =cv.glmnet(x = X, y = Y, family = "multinomial", nfolds = nfolds_train, standardize = F)
  return(list(model_full = model_full, nested_cv_probs = pred.test))
}
