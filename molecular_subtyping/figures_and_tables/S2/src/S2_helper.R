### Helpers for S2 ###

# Ordinal regression adjusted for demographics (one vs other comparison)
test_ordinal_demo <- function(subtype, response, demo_data) {
  id = (!is.na(subtype)) & (!is.na(response)) & (!rowSums(is.na(demo_data)))
  subtype = factor(subtype[id])
  response = factor(response[id])
  demo_data = demo_data[id,]
  nclass = length(unique(subtype))
  nlevel_resp = nlevels(response)
  
  x_ordinal = model.matrix(~subtype-1)
  x_ordinal = data.frame(x_ordinal[,1:(nclass-1)]-x_ordinal[,nclass])
  colnames(x_ordinal) = paste0("Subtype", levels(subtype)[1:(nclass-1)])
  df_ordinal = cbind(data.frame(response = response), x_ordinal, demo_data)
  model1 <- clm(response ~ ., data = df_ordinal)
  model_AIC <- AIC(model1)
  model_loglik <- logLik(model1)
  coef_ordinal <- head(summary(model1)$coefficients, 
                       ((nlevel_resp-1)+(nclass-1)))
  coef_ordinal <- tail(coef_ordinal, nclass-1)
  est_ordinal <- coef_ordinal[,1]
  p_ordinal <- coef_ordinal[,4]
  
  x_ordinal = model.matrix(~subtype-1)
  x_ordinal = data.frame(x_ordinal[,2:(nclass)]-x_ordinal[,1])
  colnames(x_ordinal) = paste0("Subtype", levels(subtype)[2:nclass])
  df_ordinal = cbind(data.frame(response = response), x_ordinal, demo_data)
  model1 <- clm(response ~ ., data = df_ordinal)
  coef_ordinal1 <- head(summary(model1)$coefficients, 
                        ((nlevel_resp-1)+(nclass-1)))
  coef_ordinal1 <- coef_ordinal1[(nlevel_resp-1)+(nclass-1),]
  est_ordinal1 <- coef_ordinal1[1]
  p_ordinal1 <- coef_ordinal1[4]
  est_ordinal = append(est_ordinal, est_ordinal1)
  names(est_ordinal) = paste0("Subtype", levels(subtype)[1:nclass])
  p_ordinal = append(p_ordinal, p_ordinal1)
  names(p_ordinal) = paste0("Subtype", levels(subtype)[1:nclass])
  return(list(est_ordinal = est_ordinal, p_ordinal = p_ordinal,
              model_AIC = model_AIC, model_loglik = model_loglik))
}