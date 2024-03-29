
#' meancoef
#'
#' This function can be used to compute the mean of coefficients from different
#' partitions in the context of the spicefp approach.
#'
#' @param coef.list list. The second element of the coef_spicefp function
#' outputs. It has the same name as the argument.
#' @param weight a numerical vector of weights with the same length as
#' coef.list.
#'
#' @details  Here, the fine-mesh coefficients are weighted and a weighted mean is deduced. 
#' If the user wishes, he can use as weights the slopes associated with the qualities
#' of the models concerned.
#'
#' @return Returns a list of :
#'
#' \describe{
#' \item{weighted_mean}{fine-mesh matrix or array with the weighted mean of the
#'  coefficients}
#' \item{y.estimated}{weighted estimation of \eqn{X \beta}}
#' \item{coefficients.array}{An array with all the fine-mesh coefficients that
#' will be used to compute the weighted mean}
#' \item{weight}{same as inputs}
#' }
#'
#' @importFrom stats weighted.mean
#' @importFrom purrr map_dfc
#'
#' @examples
#' \donttest{
#' ##linbreaks: a function allowing to obtain breaks linearly
#' linbreaks<-function(x,n){
#'      sort(round(seq(trunc(min(x)),
#'                 ceiling(max(x)+0.001),
#'                 length.out =unlist(n)+1),
#'             1)
#'           )
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
#'    parlist.irdc,function(x){
#'    list(x[[1]],x[[2]])}
#' )
#' m.irdc <- as.matrix(Irradiance[,-c(1)])
#' m.tpr <- as.matrix(Temperature[,-c(1)])
#'
#' # For the constructed models, only two regularization parameter ratios
#' # penratios=c(1/25,5) is used. In a real case, more candidates 
#' # and regularization parameter ratios should be evaluated.

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
#' ## Focus on the 2 best models retained by the AIC criterion at iteration 1
#' c.mdls <- coef_spicefp(ex_sp, iter_=1, criterion ="AIC_",
#'                       nmodels=2, ncores = 2,
#'                       dim.finemesh=c(1000,1000),
#'                       write.external.file = FALSE)
#'
#' # meancoef
#' # Compute the mean of the coefficients of these models
#' mean.c.mdls<-meancoef(c.mdls$coef.list,
#'                      weight = c.mdls$Model.parameters$Slope_)
#' g3<-mean.c.mdls$weighted_mean
#' g3.x<-as.numeric(rownames(g3))
#' g3.y<-as.numeric(colnames(g3))
#'
#'
#' #library(fields)
#' #plot(c(10,2000),c(15,45),type= "n", axes = FALSE,
#' #     xlab = "Irradiance (mmol/m2/s - Logarithmic scale)",
#' #     ylab = "Temperature (deg C)",log = "x")
#' #rect(min(g3.x),min(g3.y),max(g3.x),max(g3.y), col="black", border=NA)
#' #image.plot(g3.x,g3.y,g3, horizontal = FALSE,
#' #           col=designer.colors(256, c("blue","white","red")),
#' #           add = TRUE)
#' #axis(1) ; axis(2)
#'
#' closeAllConnections()
#'
#'}
#'
#' @export



meancoef<-function(coef.list, weight){
  try(if(length(weight)!=length(coef.list))
    stop("The length of weight vector must be equal to the length of  coef.list"))
  try(if(sum(is.na(weight))!=0)
    stop("weight vector must not contain missing data"))

  dim.finemesh=dim(coef.list[[1]]$Candidate.coef.NA.finemeshed)
  coef.NA.finemeshed=list()
  for (g in 1:length(coef.list)){
    coef.NA.finemeshed[[g]]<-coef.list[[g]]$Candidate.coef.NA.finemeshed
  }
  spicefp_coefweimean<-function(X,W){
    if(setequal(X, rep(NA,length(X)))) {NA} else {
      X[which(is.na(X))]<-0
      weighted.mean(X,W)
    }
  }

  if(length(dim(coef.list[[1]]$Candidate.coef.NA.finemeshed))==2){
    acoef<-array(NA,c(length(coef.NA.finemeshed),
                      dim.finemesh[1],
                      dim.finemesh[2]))
    for (i in 1:length(coef.NA.finemeshed)) {
      acoef[i,,]<-coef.NA.finemeshed[[i]]
    }
    spicefp.coef= apply(acoef,
                        c(2,3),
                        FUN = function(x){spicefp_coefweimean(X=x, W=weight)})
    dimnames(spicefp.coef)=dimnames(coef.NA.finemeshed[[1]])
  } else {
    acoef<-array(NA,c(length(coef.NA.finemeshed),
                      dim.finemesh[1],
                      dim.finemesh[2],
                      dim.finemesh[3]))
    for (i in 1:length(coef.NA.finemeshed)) {
      acoef[i,,,]<-coef.NA.finemeshed[[i]]
    }
    spicefp.coef= apply(acoef,
                        c(2,3,4),
                        FUN = function(x){spicefp_coefweimean(X=x, W=weight)})
    dimnames(spicefp.coef)=dimnames(coef.NA.finemeshed[[1]])
  }

  suppressMessages(
  y_est <-apply(map_dfc(1:length(coef.list),
                        function(x){coef.list[[x]]$y.estimated*weight[x]}),1,
                sum)/sum(weight)
  )

  return(list(weighted_mean=spicefp.coef,
              y.estimated=y_est,
              coefficients.array=acoef,
              weight=weight))
}
