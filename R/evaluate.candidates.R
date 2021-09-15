
#' evaluate.candidates
#'
#' This function performs for each candidate matrix, a Generalized Fused Lasso
#' (sparse fused lasso 2d or 3d) and computes various statistics and
#' information criteria related to the constructed model.
#'
#' @param candmatrices List. Output of the "candidates" function. The spicefp
#' dimension is the first element. The second contains many lists of one
#' candidate matrix and related vector with index and numbers of class
#' intervals used per predictor. The other elements of the lists are the
#' inputs of "candidates" function. If the user does not need the "candidates"
#' function for the creation of candmatrices, it is possible to build a list
#' provided that it respects the same structure as well as the names
#' of the outputs of the "candidates" function. In this case only the first two
#' elements of the list are essential: spicefp.dimension and candidates.
#' The remaining elements can be NULL.
#' @param y numerical vector. Contains the dependent variable. This vector will
#' be used as response variable in the construction of models involving each
#' candidate matrix.
#' @param penratios numeric vector with values greater than or equal to 0. It
#' represents the ratio between the regularization parameters of parsimony and
#' fusion. When penratios=0, it corresponds to the pure fusion. The higher
#' its value, the more parsimonious the model is.
#' @param nknots integer. For one value in penratios vector, it represents the
#' number of models that will be constructed for each candidate matrix. It is
#' the argument "nlam" of  \code{\link[genlasso:coef.genlasso]{coef.genlasso}}
#' function. This argument can also be NULL. In this case, the argument
#' appropriate.df must be provided.
#' @param appropriate.df (appropriate degree of freedom) NULL by default.
#' Numerical vector with values greater than or equal to 1. The degree of
#' freedom of generalized fused problem is equal the number of connected
#' components. A connected component gives information on a group of non-zero
#' coefficients sharing the same value and connected by a contiguity matrix.
#' More simply, it can be interpreted as a group of coefficients that have a
#' unique influence. When the user has a prior idea of the number of zones of
#' influence that the desired solution could contain, it is advisable to provide
#'  appropriate.df, a vector of appropriate degrees of freedom. In this case,
#'  nknots must be NULL.
#' @param ncores numbers of cores that will be used for parallel computation.
#' By default, it is equal to detectCores()-1.
#' @param penfun function with 2 arguments (dim1, dim2) when dealing with 2
#' dimensional spiceFP or 3 arguments (dim1, dim2, dim3) when dealing with 3
#' dimensional spiceFP. The argument order in the penalty function is
#' associated with the order of numbers of class intervals used per predictor
#' in the second element of candmatrices argument. NULL by default.
#' When penfun=NULL, getD2dSparse of genlasso or getD3dSparse is used according
#' to the dimension of spiceFP.
#' @param file_name character. It is the name of the file in which the
#' evaluation summary of all the candidate matrices is stored.
#' This file is saved in your working directory.
#' @param write.external.file logical. Indicates whether the result table should
#' be written as a file (txt) in your working directory. It is recommended to use
#' write.external.file=TRUE when evaluating a large number of candidate matrices
#' (more than 100) in order to keep memory available.
#'
#'
#' @details This function mainly returns statistics on the models built based
#' on the candidate matrices. For each candidate matrix,
#' length(penratios) x nknots or length(penratios) x length(appropriate.df)
#' models are constructed in order to estimate the regularization parameters
#' and to perform a variable selection. The computed statistics provide
#' information on the quality of the models. For obvious reasons of memory
#' management, the coefficients related to each of these models are not stored.
#' The statistics are stored in a file named via the argument file_name and
#' can be consulted to get an idea of the state of progress of the program.
#' The genlasso package is used for the implementation of the
#' Generalized Fused Lasso.
#'
#' @return The output is a list with :
#' \describe{
#' \item{evaluation.result}{Same as file_name. The file contains a matrix with in columns : the candidate index (Candidate_id),
#' the value of penratios used for this model (Pen_ratio), the parameter that penalizes the difference in related coefficients
#' (PenPar_fusion), the degree of freedom of the model (Df_), the residual sum of squares (RSS_), the Akaike information criterion
#' (AIC_), the Bayesian information criterion (BIC_), the Mallows' Cp
#' (Cp_), the Generalized Cross Validation (GCV_), the slope of the regression lm(\eqn{y} ~ \eqn{X \beta}) (Slope_), the ratio \eqn{var(y-X \beta)/var(y)}
#' (Var_ratio).}
#' \item{response.variable, penalty.ratios, nknots, appropriate.df, penalty.function}{Exactly the inputs y, penratios, nknots,
#' appropriate.df, penfun}
#' }
#'
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores
#' @importFrom purrr pmap map
#' @importFrom genlasso fusedlasso fusedlasso2d coef.genlasso
#' @importFrom stats sd lm var
#' @importFrom utils write.table read.table
#'
#' @examples
#' \donttest{
#' # Constructing 2 candidates for spiceFP data (temperature and Irradiance)
#' linbreaks<-function(x,n){
#'      sort(round(seq(trunc(min(x)),
#'                 ceiling(max(x)+0.001),
#'                 length.out =unlist(n)+1),
#'             1)
#'          )
#' }
#' # In this example, we will evaluate 2 candidates (each having 10
#' # temperature classes and respectively 10 and 20 irradiance classes).
#' # Only one value is used for alpha (logbreaks argument)
#' tpr.nclass=10
#' irdc.nclass=c(10,20)
#' irdc.alpha=0.005
#' p2<-expand.grid(tpr.nclass, irdc.alpha, irdc.nclass)
#' parlist.tpr<-split(p2[,1], seq(nrow(p2)))
#' parlist.irdc<-split(p2[,2:3], seq(nrow(p2)))
#' parlist.irdc<-lapply(
#'    parlist.irdc,function(x){
#'    list(x[[1]],x[[2]])}
#' )
#' m.irdc <- as.matrix(Irradiance[,-c(1)])
#' m.tpr <- as.matrix(Temperature[,-c(1)])
#' test2<-candidates(fp1=m.irdc,
#'                  fp2=m.tpr,
#'                  fun1=logbreaks,
#'                  fun2=linbreaks,
#'                  parlists=list(parlist.irdc,
#'                                parlist.tpr),
#'                  xcentering = TRUE,
#'                  xscaling = FALSE,
#'                  ncores=2)
#' # Evaluating candidates
#' # For the constructed models, only one regularization parameter ratio
#' # penratios=c(1) is used. In a real case, we will have to evaluate
#' # more candidates and regularization parameters ratio.
#' start_time_ev <- Sys.time()
#' evcand<-evaluate.candidates(candmatrices = test2,
#'                            y=FerariIndex_Difference$fi_dif,
#'                            penratios=c(1),
#'                            appropriate.df=NULL,
#'                            nknots = 100,
#'                            ncores=2,
#'                            write.external.file = FALSE)
#' duration_ev <- Sys.time() - start_time_ev
#' tab_res<-evcand$evaluation.result
#' dim(tab_res)
#' tab_res[which.min(tab_res$AIC_),]
#'
#' closeAllConnections()
#'
#'}
#' @export

evaluate.candidates <- function(candmatrices, y,
                                penratios, nknots,
                                appropriate.df=NULL,
                                ncores=parallel::detectCores()-1,
                                penfun=NULL, file_name="parametertable",
                                write.external.file=TRUE){
  try(if( !is.logical(write.external.file)) stop("write.external.file must be logical "))
  try(if(is.null(appropriate.df)+is.null(nknots)>1)
    stop("Only one of appropriate.df or nknots must be provided. The other must be NULL"))
  try(if(!is.list(candmatrices)|sum(names(candmatrices)==c("spicefp.dimension","candidates","fp1","fp2","fp3","fun1","fun2",
                                                       "fun3","parlists","xcentering","xscaling"))!=11)
    stop("candmatrices x must be a list containing 11 items named spicefp.dimension,candidates,fp1,fp2,fp3,fun1,fun2,
         fun3,parlists,xcentering,xscaling. It corresponds to an object of class xxxxxx"))
  try(if(isTRUE(write.external.file) & !is.character(file_name))
    stop("If write.external.file=TRUE then a filename (file_name) must be provided"))


fun.D=ifelse(is.null(penfun),
             ifelse(candmatrices$spicefp.dimension == 2,
                    genlasso::getD2dSparse,
                    getD3dSparse),
             penfun)
is2dim<-candmatrices$spicefp.dimension == 2
isdf<-!is.null(appropriate.df)
crit_n <- length(y)
varY <- stats::var(y)
file_name=paste0(Sys.time()," ",file_name,".txt")

#------------------------------------------------------------------- for1candidate
for1candidate <- function(candidate, dim.candidate){
  suppressWarnings({
    do.call(rbind,
            map(penratios, function(GAMMA){
              ## Find beta
              res4GAMMA<- if (is2dim){fusedlasso(y,candidate,fun.D(dim.candidate[2], dim.candidate[3]), gamma = GAMMA)} else
              {fusedlasso(y,candidate,fun.D(dim.candidate[2], dim.candidate[3], dim.candidate[4]), gamma = GAMMA)}
              Coef_gen <-if(isdf){coef.genlasso(res4GAMMA,df=appropriate.df)} else
              {coef.genlasso(res4GAMMA,nlam=min(nknots,dim(res4GAMMA$beta)[2]))}

              ## compute criteria related to each model
              # lists of lambda, coefficients and df
              lambdas_nlam <- as.list(Coef_gen$lambda)
              beta_nlam <- as.list(as.data.frame(Coef_gen$beta))
              df_est_nlam <- as.list(Coef_gen$df)
              # computing
              allpar<-do.call(rbind,
                              pmap( list(beta_gflx=beta_nlam, crit_k=df_est_nlam, LAMBDA=lambdas_nlam) ,
                                    function(beta_gflx,crit_k,LAMBDA) {
                                      ## Criteria (Efficient algorithm to select tuning parameters in
                                      ## sparse regression modeling with regularization- Kei Hirose, Shohei Tateishi and Sadanori Konishi)
                                      SCR <- sum((y-(candidate%*%beta_gflx))^2)
                                      aic_<- (crit_n*log(2*pi*varY)) + (SCR/varY) + (2*crit_k)
                                      bic_<- (crit_n*log(2*pi*varY)) + (SCR/varY) + (log(crit_n)*crit_k)
                                      aicc_<-(crit_n*log(2*pi*(SCR/crit_n)))+ crit_n - ((2*crit_n*crit_k)/(crit_n-crit_k-1))
                                      cp_<-  SCR + (2*varY*crit_k)
                                      gcv_<- (1/crit_n)*(SCR/((1-crit_k/crit_n)^2))
                                      slope_<- lm(candidate%*%beta_gflx~y)$coefficients[2] ; names(slope_)=c()
                                      ratio_var <- var(y-(candidate%*%beta_gflx))/var(y)
                                      return(c(Candidate_id=dim.candidate[1], Pen_ratio=GAMMA, PenPar_fusion=LAMBDA, Df_=crit_k,
                                               RSS_=SCR, AIC_=aic_, BIC_=bic_, AICc_=aicc_, Cp_=cp_, GCV_=gcv_, Slope_=slope_, Var_ratio=ratio_var))
                                    }))
              rownames(allpar)<-NULL ; return(allpar)
            }))
  })
}

#------------------------- apply for1candidate function for each candidate
registerDoParallel(cores=ncores)

if(isTRUE(write.external.file)){
  write.table(x=t(c("Candidate_id", "Pen_ratio", "PenPar_fusion", "Df_",
                    "RSS_", "AIC_", "BIC_", "AICc_", "Cp_", "GCV_", "Slope_", "Var_ratio")),
              file = file_name, append = FALSE, row.names=FALSE, col.names =FALSE)
  j=c()
  res_allc<-foreach(j = 1:length(candmatrices$candidates), .combine=rbind) %dopar% {
    write.table(x=for1candidate(candmatrices$candidates[[j]][[1]], candmatrices$candidates[[j]][[2]]),
                file = file_name, append = TRUE, row.names=FALSE, col.names =FALSE)
  }
  reseval<-list(evaluation.result=file_name,
                response.variable=y,
                penalty.ratios=penratios,
                nknots=nknots,
                appropriate.df=appropriate.df,
                penalty.function=fun.D)
} else {
  j=c()
  res_allc<-foreach(j = 1:length(candmatrices$candidates), .combine=rbind) %dopar% {
    for1candidate(candmatrices$candidates[[j]][[1]], candmatrices$candidates[[j]][[2]])
  }

  reseval<-list(evaluation.result= data.frame(res_allc),
                response.variable=y,
                penalty.ratios=penratios,
                nknots=nknots,
                appropriate.df=appropriate.df,
                penalty.function=fun.D)
}

return(reseval)
}



