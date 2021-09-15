

#' candidates
#'
#' The "candidates" function essentially provides the candidate matrices and
#' their characteristics. These candidate matrices can be constructed from 2 or
#' 3 functional predictors. 
#'
#' @param fp1 numerical matrix with in columns observations of one statistical
#'individual to partition. Each column corresponds to the functional predictor
#'observation for one statistical individual. The order of statistical
#'individuals is the same as in fp2. It is assumed that no data are
#'missing and that all functional predictors are observed on an equidistant (time) scale.
#' @param fp2 numerical matrix with the same number of columns and rows as fp1.
#' Columns are also observations. The order of statistical individuals is
#' the same as in fp1.
#' @param fp3 NULL by default. numerical matrix with the same number of columns and
#'  rows as fp1 and fp2. The order of statistical individuals is the same
#'  as in fp1 and fp2.
#' @param fun1 a function object with 2 arguments. First argument is fp1 and
#' the second is a list of parameters that will help to partition fp1, such as
#' the number of class intervals, etc. For example, the list of parameters for using 
#' the logbreaks function is equivalent to list(alpha, J). All arguments
#' to be varied for the creation of different candidate matrices must
#' be stored in the parameter list. The other arguments must be set by default.
#' @param fun2 a function object with 2 arguments. First argument is fp2 and
#' the second is a list of parameters.
#' @param fun3 NULL by default. Same as fun1 and fun2, a function with 2
#' arguments fp3 and a list of parameters.
#' @param parlists list of 2 elements when fp3 and fun3 are equal to NULL or of 3
#' elements when fp3 and fun3 are provided. All elements of parlists are
#' lists that have the same length. Each list contains all the lists of
#' parameters required to create different candidates. The first
#' element of parlists concerns the list of parameters required for fun1, the second
#' element is relative to fun2 and the third to fun3. See Example 2 below.
#' @param ncores numbers of cores that will be used for parallel computation. By default, it is equal
#' to detectCores()-1.
#' @param xcentering TRUE by default. Defined whether or not the variables in
#' the new candidate matrices should be centered.
#' @param xscaling FALSE by default. Defined whether or not the variables in
#' the candidate matrices should be scaled.
#'
#' @details  The function begins by partitioning each of the functional
#' predictors using the function and associated parameter lists. Once the class
#' intervals are obtained for each predictor, a contingency table is created
#' for each statistical individual. This table counts the components of the
#' observation variable (time for time series). The contingency table is then
#' transformed into a row vector that corresponds to a row of the candidate
#' matrix created.  The number of candidate matrices is equal to the length of
#' each element contained in parlists. For a fixed index, the functional
#' predictors (fp1, fp2, fp3), the functions (fun1, fun2, fun3) and the lists
#' of parameters associated to the index in each element of parlists
#' allow to create a single candidate matrix. In addition to constructing
#' the candidate matrices, the function associates with each matrix a vector
#' containing the index and the numbers of class intervals used per predictor.
#'
#' @return The function returns a list with:
#' \describe{
#' \item{spicefp.dimension}{the dimension of the approach. Equal to 2 if fp3=NULL and 3 if not}
#' \item{candidates}{a list that has the same length as the elements of parlists. Each element of this list contains
#' a candidate matrix and a vector with index and  the numbers of class intervals used per predictor}
#' \item{fp1, fp2, fp3, fun1, fun2, fun3, parlists, xcentering, xscaling}{same as inputs}
#' }
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom stringr str_c
#' @importFrom parallel detectCores
#' @importFrom tidyr spread_
#'
#' @examples 
#' ##linbreaks: a function allowing to obtain equidistant breaks
#' linbreaks<-function(x,n){
#'     sort(round(seq(trunc(min(x)),
#'          ceiling(max(x)+0.001),
#'          length.out =unlist(n)+1),
#'          1)
#'     )
#' }
#'
#' p<-expand.grid(c(12,15),c(15,20))
#' pl<-list(split(p[,1], seq(nrow(p))),
#'          split(p[,2], seq(nrow(p))))
#'
#' # Setting ncores=2 for this example check purpose
#' test<-candidates(fp1=matrix(rnorm(1000,52,15),ncol=10),
#'                  fp2=matrix(rpois(1000,50),ncol=10),
#'                  fun1=linbreaks,
#'                  fun2=linbreaks,
#'                  parlists=pl,
#'                  xcentering = FALSE,
#'                  xscaling = FALSE,
#'                  ncores=2)
#' str(test)
#' names(test)
#'
#' # Example 2 from the spiceFP data
#' tpr.nclass=seq(10,16,2)
#' irdc.nclass=seq(20,24,2)
#' irdc.alpha=c(0.01,0.02,0.03)
#' p2<-expand.grid(tpr.nclass, irdc.alpha, irdc.nclass)
#' parlist.tpr<-split(p2[,1], seq(nrow(p2)))
#' parlist.irdc<-split(p2[,2:3], seq(nrow(p2)))
#' parlist.irdc<-lapply(
#'   parlist.irdc,function(x){
#'     list(x[[1]],x[[2]])}
#' )
#' m.irdc <- as.matrix(Irradiance[,-c(1)])
#' m.tpr <- as.matrix(Temperature[,-c(1)])
#' test2<-candidates(fp1=m.irdc,
#'                   fp2=m.tpr,
#'                   fun1=logbreaks,
#'                   fun2=linbreaks,
#'                   parlists=list(parlist.irdc,
#'                                 parlist.tpr),
#'                   xcentering = TRUE,
#'                   xscaling = FALSE,
#'                   ncores=2)
#' length(test2$candidates)
#' class(test2$candidates)
#' #View(test2$candidates[[1]][[1]])
#' dim(test2$candidates[[1]][[1]])
#' test2$candidates[[1]][[2]]
#'
#' # Closing the connections for the example check purpose
#' closeAllConnections()
#' @export
candidates <- function(fp1,fp2,fp3=NULL,
                       fun1,fun2,fun3=NULL,
                       parlists,
                       ncores=parallel::detectCores()-1,
                       xcentering=TRUE,
                       xscaling=FALSE){

try(if(is.null(fp3)+is.null(fun3)==1)
  stop("fp3 and fun3 must be provided both or equal to NULL if only 2 functional predictors are used"))
try(if(is.null(fp3) & length(parlists)!=2)
  stop("The length of parlists should be equivalent to the number of functional predictors provided"))
try(if(ncol(fp1)!=ncol(fp2))
  stop("fp1 and fp2 must have the same number of columns "))

# Initialization
registerDoParallel(cores=ncores)
nmat <- length(parlists[[1]])
createbins <- function(b){
    paste0("[",b[1:(length(b)-1)],",",b[2:length(b)],"[")
}

pa1<-parlists[[1]]
pa2<-parlists[[2]]
if(!is.null(fp3)){
  pa3<-parlists[[3]]
}

# foreach loop
m=c()
res<-foreach(m=1:nmat) %dopar% {
  breaks1 <- fun1(c(fp1), pa1[[m]]) ; bins1 <- createbins(breaks1)
  breaks2 <- fun2(c(fp2), pa2[[m]]) ; bins2 <- createbins(breaks2)

  if(is.null(fp3)){
    # 1 histogram per individual
    ht_ind<-array(NA,
                  dim = c(ncol(fp1), length(bins1), length(bins2)),
                  dimnames=list(colnames(fp1), bins1, bins2))
    for(i in 1:ncol(fp1)){
      ht_ind[i,,]<-hist_2d(x=fp1[,i],
                           y=fp2[,i],
                           breaks_x=breaks1,
                           breaks_y=breaks2)$Hist.Values
    }

    # Transform the array into a matrix (ncol(fp1) x (length(bins1) x length(bins2)))
    n_nvvar<-spread_(as.data.frame.table(ht_ind), "Var1", "Freq")
    rownames(n_nvvar)<-str_c(n_nvvar$Var2,"_",n_nvvar$Var3)
    list(
      assign(paste0("CM",m), base::scale(t(n_nvvar[,c(3:ncol(n_nvvar))]),
                                         center = xcentering,
                                         scale = xscaling)),
      c(m, length(bins1),length(bins2))
    )
  } else {
    breaks3 <- fun3(c(fp3), pa3[[m]]) ; bins3 <- createbins(breaks3)
    # 1 histogram per individual
    ht_ind<-array(NA,
                  dim = c(ncol(fp1), length(bins1), length(bins2), length(bins3)),
                  dimnames=list(colnames(fp1), bins1, bins2, bins3))
    for(i in 1:ncol(fp1)){
      ht_ind[i,,,]<-hist_3d(x=fp1[,i],
                            y=fp2[,i],
                            z=fp3[,i],
                            breaks_x=breaks1,
                            breaks_y=breaks2,
                            breaks_z=breaks3)$Hist.Values
    }
    # Transform the array into a matrix (ncol(fp1) x (length(bins1) x length(bins2) x length(bins3)))
    n_nvvar<-spread_(as.data.frame.table(ht_ind), "Var1", "Freq")
    rownames(n_nvvar)<-str_c(n_nvvar$Var2,"_",n_nvvar$Var3,"_",n_nvvar$Var4)
    list(
      assign(paste0("CM",m), scale(t(n_nvvar[,c(4:ncol(n_nvvar))]),
                                   center = xcentering,
                                   scale = xscaling)),
      c(m, length(bins1),length(bins2),length(bins3))
    )
  }
}

toret <-list(spicefp.dimension=ifelse(is.null(fp3),2,3),
             candidates=res,
             fp1=fp1,
             fp2=fp2,
             fp3=fp3,
             fun1=fun1,
             fun2=fun2,
             fun3=fun3,
             parlists=parlists,
             xcentering=xcentering,
             xscaling=xscaling)
# class(toret)<-"CandidateMatrices"
return(toret)
}


