

#' coef_spicefp
#'
#' This function allows to obtain the coefficients of a model (involving a
#' candidate matrix and 2 regularization parameters). There are two possible options to
#' use this function: 1/ by minimizing an information criterion and selecting a number
#' of model (option by default), or 2/ directly by providing the parameters of the
#' model(s) that the user wishes to reconstruct.
#'
#' @param spicefp.result List. Outputs of the spicefp function.
#' @param iter_ integer. number of the iteration of interest.
#' @param criterion character. One of "AIC_", "BIC_", "Cp_". Can be NULL,
#' "AIC_" by default. If specified, nmodels must also be provided.
#' @param nmodels integer. Equal to 1 by default. Represents the number of best
#' models, according to the information criterion used.
#' Should be NULL if criterion = NULL.
#' @param model.parameters data.frame. NULL by default. One or more rows
#' contained in the file where the model statistics were stored. Be careful to
#' use the file related to the selected iteration. Names used in
#' model.parameters shoud be the same in the file.
#' @param dim.finemesh numeric vector of length 2 or 3. This vector informs
#' about the dimension of the fine-mesh arrays (or matrices).
#' @param ncores numbers of cores that will be used for parallel computation. By default, it is equal
#' to detectCores()-1.
#' @param write.external.file logical. indicates whether the result table related
#' to each iteration has been written as a file (txt) in your working directory.
#' This argument must be equal to the argument with the same name in the spicefp
#' function.
#'
#'
#' @details By providing criterion and nmodels, the function returns the
#' coefficients of the nmodels best models chosen by the selected information
#' criterion. When model.parameters is instead provided, it returns the
#' coefficients of the models described on each row of the data.frame.
#'
#' @return Returns a list of 2 elements:
#'
#' \describe{
#' \item{Model.parameters}{data.frame where each row contains statistics related to the
#' models of interest. Same as input if model.parameters is provided.}
#' \item{coef.list}{List of length nmodels or the number of rows in Model.parameters.
#' Each element of this list contains the  model results as provided by the genlasso
#' package, its coefficients without and with NA, a fine-mesh array with the coefficients,
#' and the estimation of \eqn{X \beta}. Coefficients with NA are coefficient vector
#' where the coefficient value of never-observed joint modalities is NA.  }
#' }
#'
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores
#' @importFrom purrr map
#' @importFrom genlasso fusedlasso coef.genlasso
#' @importFrom utils read.table
#'
#' @examples
#'
#' \donttest{
#' ##linbreaks: a function allowing to obtain equidistant breaks
#' linbreaks<-function(x,n){
#'        sort(round(seq(trunc(min(x)),
#'                 ceiling(max(x)+0.001),
#'                 length.out =unlist(n)+1),
#'             1)
#'            )
#' }
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
#'      parlist.irdc,function(x){
#'      list(x[[1]],x[[2]])}
#' )
#' m.irdc <- as.matrix(Irradiance[,-c(1)])
#' m.tpr <- as.matrix(Temperature[,-c(1)])
#'
#' # For the constructed models, only two regularization parameter ratios
#' # penratios=c(1/25,5) is used. In a real case, we will have to evaluate
#' # more candidates and regularization parameters ratio.
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
#'               write.external.file = FALSE)
#'
#' # coef_spicefp
#' ## coefficients based on the parameters of the model
#' ## focus on model selected by Mallows's Cp at iteration 1
#'
#' start_time_spc <- Sys.time()
#' results.eval.iter1<-ex_sp$Evaluations[[1]]$Evaluation.results$evaluation.result
#' c.mdl <- coef_spicefp(ex_sp, iter_=1,
#'                       criterion =NULL,
#'                       nmodels=NULL,
#'   model.parameters=results.eval.iter1[which.min(results.eval.iter1$Cp_),],
#'                       ncores = 1,
#'                       write.external.file =FALSE)
#'
#' g1<-c.mdl$coef.list$'231'$Candidate.coef.NA.finemeshed
#' g1.x<-as.numeric(rownames(g1))
#' g1.y<-as.numeric(colnames(g1))
#' duration_spc <- Sys.time() - start_time_spc
#'
#' #library(fields)
#' #plot(c(10,2000),c(15,45),type= "n", axes = FALSE,
#' #     xlab = "Irradiance (mmol/m2/s - Logarithmic scale)",
#' #     ylab = "Temperature (deg C)",log = "x")
#' #rect(min(g1.x),min(g1.y),max(g1.x),max(g1.y), col="black", border=NA)
#' #image.plot(g1.x,g1.y,g1, horizontal = FALSE,
#' #           col=designer.colors(64, c("blue","white")),
#' #           add = TRUE)
#' #axis(1) ; axis(2)
#'
#' ## Let's visualize the same model from other arguments of coef_spicefp
#' c.crit <- coef_spicefp(ex_sp, iter_=1,
#'                        criterion ="Cp_",nmodels=1,
#'                        ncores = 1,
#'                        write.external.file =FALSE)
#' g2<-c.crit$coef.list$'231'$Candidate.coef.NA.finemeshed
#' g2.x<-as.numeric(rownames(g2))
#' g2.y<-as.numeric(colnames(g2))
#' #plot(c(10,2000),c(15,45),type= "n", axes = FALSE,
#' #     xlab = "Irradiance (mmol/m2/s - Logarithmic scale)",
#' #     ylab = "Temperature (deg C)",log = "x")
#' #rect(min(g2.x),min(g2.y),max(g2.x),max(g2.y), col="black", border=NA)
#' #image.plot(g2.x,g2.y,g2, horizontal = FALSE,
#' #           col=designer.colors(64, c("blue","white")),
#' #           add = TRUE)
#' #axis(1) ; axis(2)
#' closeAllConnections()
#'
#'}
#' @export


coef_spicefp<-function(spicefp.result,
                       iter_,criterion="AIC_",
                       nmodels=1,
                       model.parameters=NULL,
                       dim.finemesh=NULL,
                       ncores=parallel::detectCores()-1,
                       write.external.file=TRUE){
  ## Warnings
  if(is.null(dim.finemesh)) {
    dim.finemesh=dim(spicefp.result$coef.NA.finemeshed[[1]])
  } else {
    dim.finemesh=dim.finemesh
  }

  try(if(length(dim.finemesh)!=length(dim(spicefp.result$coef.NA.finemeshed[[1]])))
    stop("The length of the dim.finemesh vector must correspond to the spicefp dimension used"))
  if(is.null(model.parameters)) {try(if(sum(criterion %in% c("AIC_","BIC_","AICc_","Cp_"))!=1)
    stop("Only one criterion has to be chosen among AIC_, BIC_, AICc_, Cp_"))}
  try(if(is.null(model.parameters) & (is.null(nmodels) | is.null(criterion)))
    stop("One of model.parameters or criterion and nmodels must be provided"))
  try(if((is.null(nmodels) & is.null(criterion)) & is.null(model.parameters))
    stop("One of model.parameters or criterion and nmodels must be provided"))
  if(!is.null(model.parameters)) {try(if(is.null(rownames(model.parameters)))
    stop("row names of model.parameters must be provided"))}



  rescd<-spicefp.result$Evaluations[[iter_]]$Evaluation.results
  candmatrices<-spicefp.result$Candidate.Matrices
  if(is.null(model.parameters)){
    ev_res <-if(write.external.file){read.table(rescd$evaluation.result, header = T)} else {
      rescd$evaluation.result}
    ev_res <-ev_res[order(ev_res[,criterion], decreasing = F)[1:nmodels],]
  } else {
    ev_res<-model.parameters
  }

  registerDoParallel(cores=ncores)
  ## Model
  # thecand.par = list from modelparameters
  j=c()
  res_csp<-foreach(j = 1:nrow(ev_res)) %dopar% {
    thecand.par=as.list(ev_res[j,])
  thecand = candmatrices$candidates[[thecand.par$Candidate_id]][[1]]
  dim.thecand = candmatrices$candidates[[thecand.par$Candidate_id]][[2]]
  suppressWarnings({
  thecand.mdl<-if(candmatrices$spicefp.dimension==2){
                  fusedlasso(y=rescd$response.variable,
                          X=thecand,
                          D=rescd$penalty.function(dim.thecand[2],dim.thecand[3]),
                          gamma = thecand.par$Pen_ratio)
                } else {
                  fusedlasso(y=rescd$response.variable,
                             X=thecand,
                             D=rescd$penalty.function(dim.thecand[2],dim.thecand[3],dim.thecand[4]),
                             gamma = thecand.par$Pen_ratio)}
  })
  thecand.coef <-coef.genlasso(thecand.mdl,lambda=thecand.par$PenPar_fusion)
  rownames(thecand.coef$beta) <- colnames(thecand)
  xbeta<- thecand%*%thecand.coef$beta
  # NA and finemeshed coefficients
  tobeNA<- unlist(map(1:ncol(thecand), function(x){
    if(sum(thecand[,x]==0)==nrow(thecand)){x}}))
  coefcand<-thecand.coef$beta
  coefcand[tobeNA,1]<-NA
  coefcand.finemeshed <- if(candmatrices$spicefp.dimension==2){
    finemeshed2d(coefcand, n.breaks1=dim.finemesh[1], n.breaks2=dim.finemesh[2])$finemeshed.matrix} else {
      finemeshed3d(coefcand, n.breaks1=dim.finemesh[1], n.breaks2=dim.finemesh[2], n.breaks3=dim.finemesh[3])$finemeshed.array}
  list(Candidate.model=thecand.mdl, Candidate.coef=thecand.coef, Candidate.coef.NA=coefcand,
       Candidate.coef.NA.finemeshed=coefcand.finemeshed, y.estimated=xbeta)
  }

  names(res_csp)=rownames(ev_res)
  return(list(Model.parameters=ev_res,coef.list=res_csp))
}


