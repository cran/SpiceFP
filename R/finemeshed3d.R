

#' finemeshed3d
#'
#' Function that helps to transform  a vector into a 3 dimensional array (with a fine
#' mesh). In the implementation of the spiceFP approach, it allows to transform
#' matrices of coefficients having different dimensions into matrices of the
#' same dimension in order to perform arithmetic operations. In practice, the
#' 3d array to be transformed is associated with a contingency table, which
#' implies numerical variables for which classes have been created.
#'
#'
#' @param x vector or one column matrix to scale. This vector comes from the
#' vectorization of the 3d array to be transformed. x is named using the
#' concatenation of the names of the dimension of the array to be transformed,
#' as shown in the example below.
#' @param n.breaks1 integer. Number of breaks needed for the first variable
#' The variable for which classes are in first position when constructing x's
#' names is the first variable.
#' @param n.breaks2 integer. Number of breaks needed for the second variable.
#' The variable for which classes are in second position when constructing x's
#' names is the second variable.
#' @param n.breaks3 integer. Number of breaks needed for the third variable.
#' The variable for which classes are in third position when constructing x's
#' names is the third variable.
#' @param round.breaks1 integer. Number of decimals for breaks of the first
#' variable.
#' @param round.breaks2 integer. Number of decimals for breaks of the second
#' variable.
#' @param round.breaks3 integer. Number of decimals for breaks of the third
#' variable.
#'
#' @details  This function is designed to return a 3d fine meshed array and
#' breaks associated.
#' In order to obtain a fine mesh, a high number of breaks must be fixed.
#'
#' @return Returns:
#' \describe{
#' \item{finemeshed.array}{Array of dimension n.breaks1 x n.breaks2 x n.breaks3.
#' The dimension names of finemeshed.array are the breaks created from each
#' variable and the associated n.breaks. Each value of finemeshed.array is
#' equal to the value of x indexed by the classes containing the row and column
#' names of finemeshed.array}
#' \item{finemeshed.values1}{First variable breaks}
#' \item{finemeshed.values2}{Second variable breaks}
#' \item{finemeshed.values3}{Third variable breaks}
#' }
#'
#'
#' @importFrom stringr str_split
#'
#' @examples
#' set.seed(4)
#' count_table<-hist_3d(x = rnorm(1000),
#'                      y = rnorm( 1000,5,0.1),
#'                      z = rnorm( 1000,2,1),
#'                      breaks_x = seq(-4, 4, by =1),
#'                      breaks_y = seq(2, 8, by =1),
#'                      breaks_z = seq(-3, 6, by =1))$Hist.Values
#'
#' df.x<-as.data.frame.table(count_table)
#' x<-df.x$Freq
#' names(x)<-paste0(df.x$Var1,"_",df.x$Var2,"_",df.x$Var3)
#'
#' res.fm3d<- finemeshed3d(x,10,50,100)
#' dim(res.fm3d$finemeshed.array)
#' @export


finemeshed3d <- function(x,
                         n.breaks1=10,
                         n.breaks2=1000,
                         n.breaks3=500,
                         round.breaks1=9,
                         round.breaks2=9,
                         round.breaks3=9){

  retain_c_name <- if(is.matrix(x)) {rownames(x)} else {names(x)}

  # Find the 4 boundaries of each break 2d via the name of the breaks (name
  # having a certain format)
  var1_min <- as.numeric(chartr(old = "[",
                                new = " ",
                                str_split(str_split(retain_c_name,"_",simplify = TRUE)[,1] ,
                                                    "," , simplify = TRUE)[,1]))
  var1_max <- as.numeric(chartr(old = "[",
                                new = " ",
                                str_split(str_split(retain_c_name,"_",simplify = TRUE)[,1] ,
                                                    "," , simplify = TRUE)[,2]))
  var2_min <- as.numeric(chartr(old = "[",
                                new = " ",
                                str_split(str_split(retain_c_name,"_",simplify = TRUE)[,2] ,
                                                    "," , simplify = TRUE)[,1]))
  var2_max <- as.numeric(chartr(old = "[",
                                new = " ",
                                str_split(str_split(retain_c_name,"_",simplify = TRUE)[,2] ,
                                                    "," , simplify = TRUE)[,2]))
  var3_min <- as.numeric(chartr(old = "[",
                                new = " ",
                                str_split(str_split(retain_c_name,"_",simplify = TRUE)[,3] ,
                                                    "," , simplify = TRUE)[,1]))
  var3_max <- as.numeric(chartr(old = "[",
                                new = " ",
                                str_split(str_split(retain_c_name,"_",simplify = TRUE)[,3] ,
                                                    "," , simplify = TRUE)[,2]))

  ## limits
  var1min_smr = min(var1_min) ; var1max_smr = max(var1_max)
  var2min_smr = min(var2_min) ; var2max_smr = max(var2_max)
  var3min_smr = min(var3_min) ; var3max_smr = max(var3_max)

  valeurs_var1_smr = unique(round(seq(var1min_smr,
                                      var1max_smr,
                                      length.out = n.breaks1),
                                  round.breaks1))
  valeurs_var2_smr = unique(round(seq(var2min_smr,
                                      var2max_smr,
                                      length.out = n.breaks2 ),
                                  round.breaks2))
  valeurs_var3_smr = unique(round(seq(var3min_smr,
                                      var3max_smr,
                                      length.out = n.breaks3 ),
                                  round.breaks3))

  smr_0 = array(NA, dim=c(length(valeurs_var1_smr),
                          length(valeurs_var2_smr),
                          length(valeurs_var3_smr)),
                dimnames = list(valeurs_var1_smr,
                                valeurs_var2_smr,
                                valeurs_var3_smr))

  for(i in 1:length(x)){
    smr_0[which(as.numeric(dimnames(smr_0)[[1]]) >=  var1_min[i]  &
                as.numeric(dimnames(smr_0)[[1]]) <= var1_max[i]) ,
           which(as.numeric(dimnames(smr_0)[[2]]) >=  var2_min[i]  &
                 as.numeric(dimnames(smr_0)[[2]]) <= var2_max[i]),
           which(as.numeric(dimnames(smr_0)[[3]]) >=  var3_min[i] &
                 as.numeric(dimnames(smr_0)[[3]]) <= var3_max[i])] <- x[i]
  }

  return(list(finemeshed.array=smr_0,
              finemeshed.values1=valeurs_var1_smr,
              finemeshed.values2=valeurs_var2_smr,
              finemeshed.values3=valeurs_var3_smr))
}

