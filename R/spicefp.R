#' spicefp
#'
#' This function is used to implement the spiceFP approach. This approach
#' transforms 2 (by default) or 3 functional predictors into candidate explonatory matrices in order to
#' identify joint classes of influence. It can take functional predictors and
#' partitioning functions as inputs in order to create candidate matrices to be
#' evaluated. The user can choose among the existing partitioning functions (as logbreaks) or provide
#' his own partitioning functions specific to the functional predictors under consideration. The user
#' can also directly provide candidate matrices already constructed as desired.
#'
#' @param y a numerical vector. Contains the dependent variable. This vector will
#' be used as response variable in the construction of models involving each
#' candidate matrix.
#' @param fp1 a numerical matrix with in columns observations of one statistical
#' individual to partition. Each column corresponds to the functional
#' predictor observation for one statistical individual. The order of the
#' statistical individuals is the same as in fp2. It is assumed that no
#' data are missing and that all functional predictors are observed on an equidistant (time) scale.
#' @param fp2 a numerical matrix with the same number of columns and rows as fp1.
#' Columns are also observations.
#' The order of the statistical individuals is the same as in fp1.
#' @param fp3 NULL by default. A numerical matrix with the same number of columns and
#' rows as fp1 and fp2.
#' The order of the statistical individuals is the same as in fp1 and fp2.
#' @param fun1 a function object with 2 arguments. First argument is fp1 and
#' the second is a list of parameters that will help to partition fp1, such as
#' the number of class intervals, etc. For example using the logbreaks function,
#' the list of parameters is equivalent to list(alpha, J). All the arguments
#' to be varied for the creation of different candidate matrices must
#' be stored in the parameter list. The other arguments must be set by default.
#' @param fun2 a function object with 2 arguments. First argument is fp2 and the
#'  second is a list of parameters.
#' @param fun3 NULL by default. Same as fun1 and fun2, a function with 2
#' arguments fp3 and a list of parameters.
#' @param parlists a list of 2 elements when fp3 and fun3 are equal to NULL or of 3
#' elements when fp3 and fun3 are provided.
#' All the elements of parlists are lists that have the same length. Each list
#' contains all the lists of parameters that have to be used to create different
#'  candidates. The first element of parlists concerns the first functional
#'  predictor fp1, the second element is relative to fp2 and the third to fp3.
#' @param xcentering TRUE by default. Defined whether or not the variables in
#' the new candidate matrices should be centered.
#' @param xscaling FALSE by default. Defined whether or not the variables in
#' the candidate matrices should be scaled.
#' @param candmatrices NULL by default. List. Output of the "candidates"
#' function. The spiceFP dimension is its first element. The second contains
#' many lists of one candidate matrix and related vector with index and numbers
#' of class intervals used per predictor. The other elements of the lists are
#' the inputs of "candidates" function. If the user does not need the
#' "candidates" function for the creation of candmatrices, it is possible to
#' build a list while making sure that it respects the same structure as well as
#' the names of the outputs of the "candidates" function. In this case, only the
#' first two elements of the list are essential: spicefp.dimension and
#' candidates. The remaining elements can be NULL.
#' @param K number of iterations of the spiceFP approach. Equal to 2 by default.
#' @param criterion character. One of "AIC_", "BIC_", "Cp_". The
#' criterion to be used in each iteration in order to identify the best
#' candidate matrix and to estimate the regulation parameters. This criterion
#' is used to perform model selection as well as variable selection.
#' @param penratios a numeric vector with values greater than or equal to 0. It
#' represents the ratio between the regularization parameters of parsimony and
#' fusion. When penratios=0, it corresponds to the pure fusion. The higher its
#' value, the more parsimonious the model is.
#' @param nknots integer. For one value in penratios vector, it represents the
#' number of models that will be constructed for each candidate matrix. It is
#' the argument "nlam" of \code{\link[genlasso:coef.genlasso]{coef.genlasso}}
#' function.
#' This argument can be also NULL. In this case, the argument appropriate.df
#' must be provided.
#' @param appropriate.df (appropriate degree of freedom) NULL by default. When used, 
#' nknots must be NULL. It is
#' the argument "df" of \code{\link[genlasso:coef.genlasso]{coef.genlasso}}
#' function. When the user has a prior idea of the number of zones of influence that the
#' solution could contain, it is advisable to provide appropriate.df,
#' a vector of appropriate degrees of freedom. 
#' appropriate.df is a numerical vector with values greater than or equal to 1.
#' The degree of freedom of generalized fused Lasso models is equal to the number of
#' connected components. A connected component gives information on a group of
#' non-zero coefficients sharing the same value and connected by a contiguity
#' matrix. More simply, it can be interpreted as a group of coefficients that
#' have a unique influence.
#' 
#' @param penfun function with 2 arguments (dim1, dim2) when dealing with 2
#' dimensional spiceFP, or with 3 arguments (dim1, dim2, dim3) when dealing with
#' 3 dimensional spiceFP. The argument order in the penalty function is
#' associated with the order of numbers of class intervals used per predictor
#' in the second element of candmatrices argument. NULL by default.
#' When penfun=NULL, getD2dSparse of genlasso or getD3dSparse is used according
#' to the dimension of spiceFP.
#' @param dim.finemesh numeric vector of length 2 or 3. This vector informs
#' about the dimension of the fine-mesh arrays (or matrices) that will be used
#' for the visualization of the sum of the coefficients selected at different
#' iterations.
#' @param file_name character vector. Of length K, it contains the list of
#' names that will be used to name the files containing informations on the
#' candidate matrix models
#' @param ncores numbers of cores that will be used for parallel computation. By default, it is equal
#' to detectCores()-1.
#' @param write.external.file logical. indicates whether the result table related
#' to each iteration should be written as a file (txt) in your working directory.
#' It is recommended to use write.external.file=TRUE when evaluating a large
#' number of candidate matrices (more than 100) in order to keep memory available.
#'
#'
#' @details  Three main steps are involved to implement spiceFP: transformation
#' of functional predictors, creation of a graph of contiguity constraints and
#' identification of the best class intervals and related regression
#' coefficients. 
#'
#' @return Returns a list with:
#'
#' \describe{
#' \item{Candidate.Matrices}{a list with candidate matrices and their
#' characteristics. same as candmatrices if it has been provided.}
#' \item{Evaluations}{List of length less than or equal to K. Each element of
#' the list contains information about an iteration. Contains the results
#' related to the evaluation of the candidate matrices. These include the name
#' of the file where the model information is stored, the best candidate matrix
#' and related coefficients, the partition vector that indexes it, the \eqn{X \beta}
#' estimation, the residuals, etc.}
#' \item{coef.NA}{List of length less than or equal to K. For each iteration,
#' it contains the coefficient vector where the coefficient value of
#' never-observed joint modalities is NA}
#' \item{coef.NA.finemeshed}{List of length less than or equal to K. For each
#' iteration, the coefficient vector is transformed into fine-mesh array or
#' matrix allowing arithmetic operations to be performed between coefficients
#' coming from different partitions}
#' \item{spicefp.coef}{fine-mesh array or matrix. Sum of the coefficients
#' selected at all iterations}
#' }
#'
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom stringr str_c
#' @importFrom parallel detectCores
#' @importFrom tidyr spread_
#' @importFrom purrr pmap map
#' @importFrom genlasso fusedlasso fusedlasso2d coef.genlasso
#' @importFrom stats sd lm var
#' @importFrom utils write.table read.table
#'
#' @examples
#'\donttest{
#' ##linbreaks: a function allowing to obtain breaks linearly
#' linbreaks<-function(x,n){
#'     sort(round(seq(trunc(min(x)),
#'                 ceiling(max(x)+0.001),
#'                 length.out =unlist(n)+1),
#'             1)
#'         )
#' }
#'
#' # In this example, we will evaluate 2 candidates with 14 temperature
#' # classes and 15 irradiance classes. The irradiance breaks are obtained
#' # according to a log scale (logbreaks function) with different alpha
#' # parameters for each candidate (0.005, 0.01).
#' ## Data and inputs
#' tpr.nclass=14
#' irdc.nclass=15
#' irdc.alpha=c(0.005, 0.01)
#' p2<-expand.grid(tpr.nclass, irdc.alpha, irdc.nclass)
#' parlist.tpr<-split(p2[,1], seq(nrow(p2)))
#' parlist.irdc<-split(p2[,2:3], seq(nrow(p2)))
#' parlist.irdc<-lapply(
#'    parlist.irdc,function(x){
#'    list(x[[1]],x[[2]])}
#' )
#' m.irdc <- as.matrix(Irradiance[,-c(1)])
#' m.tpr <- as.matrix(Temperature[,-c(1)])
#'
#' # For the constructed models, only two regularization parameter ratios
#' # penratios=c(1/25,5) are used. In a real case, we will have to evaluate
#' # more candidates and regularization parameters ratio.
#' start_time_sp <- Sys.time()
#' ex_sp<-spicefp(y=FerariIndex_Difference$fi_dif,
#'               fp1=m.irdc,
#'               fp2=m.tpr,
#'               fun1=logbreaks,
#'               fun2=linbreaks,
#'               parlists=list(parlist.irdc,
#'                             parlist.tpr),
#'               penratios=c(1/25,5),
#'               appropriate.df=NULL,
#'               nknots = 100,
#'               ncores =2,
#'               write.external.file=FALSE)
#'
#' duration_sp <- Sys.time() - start_time_sp
#' # View(ex_sp$Evaluations[[1]]$Evaluation.results$evaluation.result)
#' # View(ex_sp$Evaluations[[2]]$Evaluation.results$evaluation.result)
#' # Visualization of the coefficients
#' g<-ex_sp$spicefp.coef
#' g.x<-as.numeric(rownames(g))
#' g.y<-as.numeric(colnames(g))
#'
#' library(fields)
#' plot(c(10,2000),c(15,45),type= "n", axes = FALSE,
#'      xlab = "Irradiance (mmol/m²/s - Logarithmic scale)",
#'      ylab = "Temperature (°C)",log = "x")
#' rect(min(g.x),min(g.y),max(g.x),max(g.y), col="black", border=NA)
#' image.plot(g.x,g.y,g, horizontal = FALSE,
#'            col=designer.colors(256, c("blue","white","red")),
#'            add = TRUE)
#' axis(1) ; axis(2)
#'
#' closeAllConnections()
#'
#'}
#'
#' @export

spicefp<-function(y,
                  fp1, fp2, fp3=NULL,
                  fun1, fun2, fun3=NULL,
                  parlists,
                  xcentering=TRUE,
                  xscaling=FALSE,
                  candmatrices=NULL, K=2,
                  criterion="AIC_",
                  penratios=c(1/10,1/5,1/2,1,2,5,10),
                  nknots=50,
                  appropriate.df=NULL,
                  penfun=NULL,
                  dim.finemesh=c(1000,1000),
                  file_name=paste0("parametertable",1:2),
                  ncores=parallel::detectCores()-1,
                  write.external.file=TRUE){

try(if( !is.logical(write.external.file)) stop("write.external.file must be logical "))

try(if(sum(criterion %in% c("AIC_","BIC_","AICc_","Cp_"))!=1)
  stop("Only one criterion has to be chosen among AIC_, BIC_, AICc_, Cp_"))
try(if(length(unique(file_name))!=K)
  stop("Different filenames must be provided for iterations of spicefp.
         length(unique(file_name)) should be equal to K"))
if (is.null(candmatrices)) {
  # stop candidates
  try(if(is.null(fp3)+is.null(fun3)==1)
    stop("fp3 and fun3 must be provided both or equal to NULL if only 2 functional predictors are used"))
  try(if(is.null(fp3) & length(parlists)!=2)
    stop("The length of parlists should be equivalent to the number of functional predictors provided"))
  try(if(ncol(fp1)!=ncol(fp2))
    stop("fp1 and fp2 must have the same number of columns "))
}
if(!is.null(candmatrices)) {
  try(if(!is.list(candmatrices)|sum(names(candmatrices)==c("spicefp.dimension","candidates","fp1","fp2","fp3","fun1","fun2",
                                                           "fun3","parlists","xcentering","xscaling"))!=11)
    stop("candmatrices x must be a list containing 11 items named spicefp.dimension,candidates,fp1,fp2,fp3,fun1,fun2,
           fun3,parlists,xcentering,xscaling. It corresponds to an object of class xxxxxx"))
}
try(if(is.null(appropriate.df)+is.null(nknots)>1)
  stop("Only one of appropriate.df or nknots must be provided. The other must be NULL"))

if(is.null(candmatrices)){try(if((is.null(fp3) & length(dim.finemesh)!=2) | (!is.null(fp3) & length(dim.finemesh)!=3))
  stop("The length of the dim.finemesh vector must correspond to the spicefp dimension used"))} else  {
    try(if(length(dim.finemesh)!=candmatrices$spicefp.dimension) stop("The length of the dim.finemesh vector
                                                                      must correspond to the spicefp dimension used"))}
try(if(isTRUE(write.external.file) & !is.character(file_name))
    stop("If write.external.file=TRUE then filenames (one per iteration) must be provided"))

# Candidate construction
candmatrices <- if(is.null(candmatrices)){
                      candidates(fp1,fp2,fp3,
                                  fun1,fun2,fun3,
                                  parlists,
                                  ncores,
                                  xcentering,
                                  xscaling)
                } else {
                     candmatrices
                }

wef<-write.external.file ; rm(write.external.file)
response_var<-y
tokeep <-list()
iter_=1
check_coef_zeros <- FALSE
while(iter_ <= K & check_coef_zeros==FALSE){
  rescd <- evaluate.candidates(candmatrices,
                               y=response_var,
                               penratios,
                               nknots,
                               appropriate.df,
                               ncores,
                               penfun,
                               file_name[iter_],
                               write.external.file=wef)
  rescd_file<-if(wef){ read.table(rescd$evaluation.result, header = TRUE)} else {
    rescd$evaluation.result }

  try(if(sum(rescd_file[,"Slope_"]>0.05 & rescd_file[,"Slope_"]<0.95)==0)
    stop("At iteration ", iter_ ,", none of the candidate matrices produces a model with a slope between 0.05 and 0.95"))
  thecand.par<-as.list(rescd_file[rescd_file[,"Slope_"]>0.05 &
                                    rescd_file[,"Slope_"]<0.95,
  ][which.min(rescd_file[rescd_file[,"Slope_"]>0.05 &
                           rescd_file[,"Slope_"]<0.95, ][,criterion]),])
  rm(rescd_file)


  ## Coef for thecand.par
  thecand = candmatrices$candidates[[thecand.par$Candidate_id]][[1]]
  dim.thecand = candmatrices$candidates[[thecand.par$Candidate_id]][[2]]
  suppressWarnings({
  thecand.mdl<- if (candmatrices$spicefp.dimension==2){
                    fusedlasso(y=response_var,
                               X=thecand,
                               D=rescd$penalty.function(dim.thecand[2],
                                                        dim.thecand[3]),
                               gamma = thecand.par$Pen_ratio)
                } else {
                    fusedlasso(y=response_var,
                               X=thecand,
                               D=rescd$penalty.function(dim.thecand[2],
                                                        dim.thecand[3],
                                                        dim.thecand[4]),
                               gamma = thecand.par$Pen_ratio)
                }
  })

  thecand.coef <-coef.genlasso(thecand.mdl,lambda=thecand.par$PenPar_fusion)
  rownames(thecand.coef$beta) <- colnames(thecand)
  xbeta<- thecand%*%thecand.coef$beta
  thecand.residuals<-response_var-xbeta
  tokeep[[iter_]]<-list(Evaluation.results=rescd,
                        thecandidate.parameters=thecand.par,
                        thecandidate=thecand,
                        dim.thecandidate=dim.thecand,
                        thecandidate.model=thecand.mdl,
                        thecandidate.coefficients=thecand.coef,
                        XBeta=xbeta,
                        thecandidate.residuals=thecand.residuals)

  ## Check and next
  if(isTRUE(all.equal(c(thecand.coef$beta), rep(0,nrow( thecand.coef$beta )))))
    print(paste("Best coefficients chosen by", criterion,"at iteration",iter_,"are all equal to 0.",
                "Stop after iteration",iter_))
  if(isTRUE(all.equal(c(thecand.coef$beta), rep(0,nrow( thecand.coef$beta ))))){
    check_coef_zeros=T
  }

  iter_=iter_+1
  response_var<-thecand.residuals
}

## Treat the coef for sum
coef.NA=list()
coef.NA.finemeshed=list()
for (g in 1:length(tokeep)) {
  amop<-tokeep[[g]]
  nrowcand<-nrow(amop$thecandidate)
  tobeNA<- unlist(map(1:ncol(amop$thecandidate), function(x){
    if(sum(amop$thecandidate[,x]==0)==nrowcand){x}}))
  coefcand<-amop$thecandidate.coefficients$beta
  coefcand[tobeNA,1]<-NA
  coefcand.finemeshed <- if(candmatrices$spicefp.dimension==2){
      finemeshed2d(coefcand,
                   n.breaks1=dim.finemesh[1],
                   n.breaks2=dim.finemesh[2])$finemeshed.matrix} else {
      finemeshed3d(coefcand,
                   n.breaks1=dim.finemesh[1],
                   n.breaks2=dim.finemesh[2],
                   n.breaks3=dim.finemesh[3])$finemeshed.array}
  coef.NA[[g]]<-coefcand
  coef.NA.finemeshed[[g]]<- coefcand.finemeshed
}
spicefp_coefsum <- function(x){
  ifelse(setequal(x, rep(NA,length(x))),
         NA,
         sum(x,na.rm = T))
}
if(candmatrices$spicefp.dimension==2){
  acoef<-array(NA,c(length(coef.NA.finemeshed),dim.finemesh[1],dim.finemesh[2]))
  for (i in 1:length(coef.NA.finemeshed)) {
    acoef[i,,]<-coef.NA.finemeshed[[i]]
    }
  spicefp.coef= apply(acoef, c(2,3),FUN = function(x){spicefp_coefsum(x)})
  dimnames(spicefp.coef)=dimnames(coef.NA.finemeshed[[1]])} else {
    acoef<-array(NA,c(length(coef.NA.finemeshed),
                      dim.finemesh[1],
                      dim.finemesh[2],
                      dim.finemesh[3]))
    for (i in 1:length(coef.NA.finemeshed)) {
      acoef[i,,,]<-coef.NA.finemeshed[[i]]
      }
    spicefp.coef= apply(acoef, c(2,3,4),FUN = function(x){spicefp_coefsum(x)})
    dimnames(spicefp.coef)=dimnames(coef.NA.finemeshed[[1]])}

## Global result
return(list(Candidate.Matrices=candmatrices,
            Evaluations=tokeep,
            coef.NA=coef.NA,
            coef.NA.finemeshed=coef.NA.finemeshed,
            spicefp.coef=spicefp.coef))
}


