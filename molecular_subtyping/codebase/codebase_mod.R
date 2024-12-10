
####MCIA related functions

# FUNCTIONS FOR MCIA_mbpca.Rmd
# Anna Konstorum (anna.ionstorum@yale.edu)

# From: https://github.com/mengchen18/mogsa/blob/master/R/mbpca.R

library(RSpectra)
library(corpcor)

# Source: mogsa/svd.nipals.R
svd.solver <- function(
  x, nf, opts.svd.nipals=list(), opts.svds=list(), opts.fast.svd = list(), opts.svd=list()
) {
  if (miss.nf <- missing(nf))
    nf <- min(dim(x))
  
  smallnf <- nf < 3
  
  rt <- nrow(x)/ncol(x)
  if (any(is.na(x))) {
    r <- do.call(svd.nipals, c_mat(list(x=x, nf = nf), opts.svd.nipals))
    solver <- "svd.nipals"
    # cat("svd.nipals used. \n" )
  } else if (inherits(x, "Matrix") || smallnf) {
    r <- do.call(svds, c_mat(list(A=x, k= nf), opts.svds))
    solver <- "svds"
    # cat("svds used. \n" )
  } else if (length(x) >= 1e6 & abs(log10(rt)) >= 1) {
    if (!miss.nf)
      message("function 'fast.svd' used, all singular vectors will be computed")
    r <- do.call(fast.svd, c_mat(list(m=x), opts.fast.svd))
    solver <- "fast.svd"
    # cat("fast.svd used. \n" )
  } else {
    if (!miss.nf)
      message("Function 'svd' used, all singular vectors will be computed")
    r <- do.call(svd, c_mat(list(x=x), opts.svd))
    solver <- "svd"
    # cat("svd used. \n" )
  }
  if (inherits(x, "Matrix")) {
    r$u <- Matrix(r$u)
    r$v <- Matrix(r$v)
  }  
  attr(r, "solver") <- solver
  r
}

c_mat <- function(x, ...) {
  if (inherits(x, "Matrix"))
    as.numeric(x) else
      base::c(x, ...)
}

# Modified from function in mogsa library
# Added 'lambda_all' block-level processing option
processOpt <-
  function(x, center=TRUE, scale=FALSE, num_comps, option = c("lambda1", "inertia", "uniform","lambda_all")) {
    
    opt <- match.arg(option)  
    
    if (is.null(names(x)))
      names(x) <- paste("data", 1:length(x), sep = "_")
    
    x <- lapply(x, function(x) {
      r <- scale(x, center = center, scale = scale)
      if (inherits(x, "Matrix"))
        r <- Matrix(r)
      r
    })
    if (opt == "lambda1") {
      w <- sapply(x, function(xx) 1/svd.solver(xx, nf = 1)$d)
    } else if (opt == "inertia") {
      w <- sapply(x, function(xx) 1/sqrt(sum(xx^2)))
    } else if (opt == "uniform") {
      w <- rep(1, length(x))
    } else if (opt == "lambda_all"){
      w <- sapply(x, function(xx) 1/sum(svd.solver(xx, nf = num_comps)$d))
    }
    mapply(SIMPLIFY = FALSE, function(xx, ww) xx*ww, xx=x, ww=w)
  }


# Source: ade4 library

"array2ade4" <-
  function(dataset, pos=FALSE,  trans=FALSE){
    
    # Allows matrix, data.frame, ExpressionSet, marrayRaw to be read as data.frame
    if (!is.data.frame(dataset)) dataset<-getdata(dataset)
    
    if (any(is.na(dataset)))
      stop("Arraydata must not contain NA values. Use impute.knn in library(impute), KNNimpute from Troyanskaya et al., 2001 or LSimpute from Bo et al., 2004 to impute missing values\n")
    
    
    # COA needs table of positive data, will add real no to make +ve
    if(pos){
      if (any(dataset < 0)) {
        num<-round(min(dataset)-1)
        dataset<-dataset+abs(num)
      }
    }
    
    if(trans) {
      # Transpose matrix  (as BGA, CIA expects the samples to be in the rows)
      # dudi.nsc should not be transposed, use t.dudi instead to ensure row weight are equal
      # There is a horrible bug is dudi.pca/coa etc, if a dataset with vars>>cases is given
      # It can end abruptly crashing the session. This is a bug in sweep
      # There will now use t.dudi rather than transpose the data
      
      # using t convert data.frame to matrix and messes up affymetrix probe ID names
      # It changes all of the "-" to "." in probeids like AFFX-CreX-5_at
      # So save names and change the names of the transposed matrix
      
      colnam= colnames(dataset)
      rownam = rownames(dataset)               
      dataset<-t(dataset)		
      dimnames(dataset) = list(colnam, rownam) 
      
      # convert matrix to data.frame for ade4
      dataset <- as.data.frame(dataset)
      
      if (!is.data.frame(dataset)) stop("Problems checking dataset")
    }
    
    data.out<-dataset        
    return(data.out)
  }

# A.K. functions
# center and scale all feature x sample arrays in list 
center_scale<-function(x){
  data_out<-lapply(x, function(x) {
    r<-t(scale(t(x), center=TRUE, scale=TRUE))
  })
  return(data_out)}

# Preprocess arrays using nsc
nsc_prep<-function(data_input, num_comps){
  df.list<-data_input
  df.list <- lapply(df.list, as.data.frame)
  df.list <- lapply(df.list, array2ade4, pos = TRUE)
  coa.list <- lapply(df.list, dudi.nsc, scannf = FALSE, nf = num_comps)
  coa.list.t <- lapply(coa.list, ade4:::t.dudi)
  dfl <- lapply(coa.list, function(x) x$tab)
  return(dfl)
}

new_gs<-function(data_input,mcia_out){
  nblocks=length(data_input)
  num_comps=dim(mcia_out$t)[2]
  #data_prep<-lapply(data_input,function(x) as.matrix(data.frame(x)))
  blocks_out<-lapply(1:nblocks,matrix,data=NA,nrow=dim(data_input[[1]])[2],ncol=num_comps)
  names(blocks_out)<-names(data_input)
  for (i in 1:nblocks){
    name = names(data_input)[i]
    block = t(as.matrix(data.frame(data_input[name]))) %*% as.matrix(data.frame(mcia_out$pb[name]))
    blocks_out[[i]]<-block
    rownames(blocks_out[[i]])<-colnames(data_input[name])
  }
  weights = mcia_out$w
  
  gs = matrix(0,dim(data_input[[1]])[2],ncol(mcia_out$t))
  for (i in 1:num_comps){
    out_sum=0
    for (j in 1:nblocks){
      out = weights[j,i]*as.matrix(data.frame(blocks_out[[j]][,i]))
      out_sum=out_sum+out
    }
    gs[,i] = out_sum
  }
  return(gs)
}

# run mcia
# data_input is list of feature x samples arrays; samples must match

#Output:
# mcia_obj$tb - the block scores
# mcia_obj$pb - the block loadings
# mcia_obj$t - the global scores
# mcia_obj$w - the weights of block scores to construct the global score
mcia_mbpca<-function(data_input,num_comps,preprocess='nsc',block_prep='inertia',
                     deflat_method='blockLoading'){
  if (preprocess=='nsc'){
    table_out<-nsc_prep(data_input, num_comps)
  }
  else{
    table_out<-center_scale(data_input)
  }
  # block-level preprocess 
  final_out<-processOpt(table_out,scale=FALSE,center=FALSE,num_comps,option=block_prep)
  
  mcia_out<-mbpca(final_out,ncomp=num_comps,k="all",method=deflat_method,
                  option="uniform",center=FALSE,scale=FALSE,
                  unit.p=TRUE,unit.obs=TRUE,moa=FALSE)
  rownames(mcia_out$t)<-colnames(data_input[[1]])
  mcia_out<-list(data_prep=final_out,mcia_result=mcia_out)
  return(mcia_out)
}

#input: mcia model and raw datasets before MOFA imputation
projection_coef_func = function(mcia_factors,datasets){
  # datasest 
  xlist= datasets
  xlist1 = datasets
  factors_mcia =mcia_factors
  projection_coef = list()
  projection_pval = list()
  ave_variance_explained = matrix(0, nrow = ncol(factors_mcia), ncol = length(xlist))
  total_variance_explained = matrix(0, nrow = ncol(factors_mcia), ncol = length(xlist))
  colnames(ave_variance_explained) <- colnames(total_variance_explained) <- names(xlist)
  # standardize data with quantile normalization (for raw features with missing data)
  for(d in 1:length(xlist1)){
    llkeep = which(!is.na(xlist1[[d]][,1]))
    xlist[[d]] = apply(xlist[[d]],2, function(z){
      idx = which(!is.na(z))
      zobs = qqnorm(z[idx], plot.it = F)$x
      zout = z
      zout[idx] = zobs
      return(zout)
    })
    xlist1[[d]][,] = NA
    xlist1[[d]][llkeep,] = xlist[[d]][llkeep,]
  }
  # standardize factor
  factors_normalized = apply(as.matrix(factors_mcia),2,function(z) qqnorm(z,plot.it = F)$x)
  for(d in 1:length(xlist)){
    print(d)
    # use observed features
    idx_keep = which(!is.na(xlist1[[d]][,1]))
    tmp = lm(as.matrix(xlist1[[d]][idx_keep,])~ scale(as.matrix(factors_normalized[idx_keep,])))
    tmp1 = summary(tmp)
    projection_coef[[d]] =data.frame(t(tmp$coefficients[-1,]))
    projection_coef[[d]][is.na(projection_coef[[d]])] = 0
    projection_pval[[d]] = array(NA, dim = dim(projection_coef[[d]]))
    for(j in 1:ncol(xlist[[d]])){
      projection_pval[[d]][j,]=tmp1[[j]]$coefficients[-1,4]
    }
    projection_pval[[d]] = data.frame(projection_pval[[d]])
    #projection_coef[[d]]: assay d, p by K (p features, K factors)
    ave_variance_explained[,d] = apply(projection_coef[[d]]^2,2,mean)
    total_variance_explained[,d] = ave_variance_explained[,d]*nrow(projection_pval[[d]] )
    tmp_name = colnames(xlist1[[d]])
    tmp_name = sapply(strsplit(tmp_name, "_"),function(z) z[length(z)])
    rownames(projection_coef[[d]]) <- rownames(projection_pval[[d]]) <- tmp_name
  }
  projection_evaluation = list()
  projection_evaluation$projection_coef = projection_coef
  projection_evaluation$projection_pval = projection_pval
  #ave_variance_explained: (K factor by d assay)
  projection_evaluation$ave_variance_explained = ave_variance_explained
  #ave_variance_explained: (K factor by d assay)
  projection_evaluation$total_variance_explained =total_variance_explained
  return(  projection_evaluation)
}