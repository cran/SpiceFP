
#' hist_3d
#'
#' This function can be used in order to construct a 3D histogram based on 3
#' variables and relative breaks directly provided as inputs.
#'
#' @param x either a numerical vector to be partitioned or a matrix with 3 numerical 
#' columns to be partitioned.
#' @param y a numerical vector to be partitioned. Not required if x is a matrix.
#' @param z a numerical vector to be partitioned. Not required if x is a matrix
#' @param breaks_x a numerical vector. Contains the breaks related to x for the
#' histogram
#' @param breaks_y a numerical vector. Contains the breaks related to y for the
#' histogram
#' @param breaks_z a numerical vector. Contains the breaks related to z for the
#' histogram
#' @param same.scale logical. Default to FALSE. If TRUE, breaks_x will be used
#' for x, y and z
#' @param na.rm logical. Default to TRUE. Indicates whether missing values
#' should be removed
#' @param FUN function used to summarize bin contents.
#'
#'
#' @details  The default function used for the argument FUN is the function
#' length. When another function is used, it is applied on x or on the first column
#' of x if this is a three-column matrix. The lower limit of each class interval is
#' included in the class and the upper limit is not.
#'
#' @return Using a given set of breaks per each variable, the function returns :
#' \describe{
#' \item{Hist.Values}{a 3 dimensional array. The 1st (respectively 2nd, 3rd) dimension  is 
#' related to the class intervals of x (resp. y, z). Contingency table is returned if
#'  FUN=length}
#' \item{breaks_x, breaks_y, breaks_z }{ same as the inputs of the function}
#' \item{Midpoints.x, Midpoints.y, Midpoints.z }{the midpoints for each bin
#' per variable}
#' \item{nobs.x , nobs.y, nobs.z }{number of observations of x, y and z}
#' \item{n.bins}{ vector of 3 elements containing the number of bins for x, y
#' and z }
#' }
#'
#' @examples
#' set.seed(4)
#' hist_3d(x = rnorm(1000),
#'         y = rnorm( 1000,5,0.1),
#'         z = rnorm( 1000,2,1),
#'         breaks_x = seq(-4, 4, by =1),
#'         breaks_y = seq(2, 8, by =1),
#'         breaks_z = seq(-2, 6, by =1))
#' @export
hist_3d <- function (x,
                     y,
                     z,
                     breaks_x,
                     breaks_y,
                     breaks_z,
                     same.scale = FALSE,
                     na.rm = TRUE,
                     FUN = length){

  # Define x, y and z
  if (is.null(y) + is.null(z) >=1) {
    if (ncol(x) != 3) {
      stop("If y and z are ommitted, x must be a 3 column matrix")
    } else {
      z <- x[, 3]
      y <- x[, 2]
      x <- x[, 1]
    }
  }

  # Management of NAs
  nas <- is.na(x) | is.na(y) | is.na(z)
  if (na.rm) {
    x <- x[!nas]
    y <- y[!nas]
    z <- z[!nas]
  }
  else stop("missing values not permitted if na.rm=FALSE")

  # Use of same scale
  if(same.scale){
    x.cuts <- breaks_x
    y.cuts <- breaks_x
    z.cuts <- breaks_x
  } else {
    x.cuts <- breaks_x
    y.cuts <- breaks_y
    z.cuts <- breaks_z
  }

  index.x <- cut(x, x.cuts, include.lowest = TRUE,right =FALSE,
                 labels = paste0("[",x.cuts[1:length(x.cuts)-1],
                                 ",",x.cuts[2:length(x.cuts)],"["))
  index.y <- cut(y, y.cuts, include.lowest = TRUE,right =FALSE,
                 labels = paste0("[",y.cuts[1:length(y.cuts)-1],
                                 ",",y.cuts[2:length(y.cuts)],"["))
  index.z <- cut(z, z.cuts, include.lowest = TRUE,right =FALSE,
                 labels = paste0("[",z.cuts[1:length(z.cuts)-1],
                                 ",",z.cuts[2:length(z.cuts)],"["))

  # Apply FUN to vector x via tapply. Aside from lenght, note that the other
  # functions are applied on x
  m <- tapply(x, list(index.x, index.y, index.z), FUN)

  # If length is used, replace NA with 0
  if (identical(FUN, base::length))  m[is.na(m)] <- 0

  # xlab, ylab,zlab
  # if (missing(xlab)) xlab <- deparse(substitute(xlab))
  # if (missing(ylab)) ylab <- deparse(substitute(ylab))
  # if (missing(zlab)) zlab <- deparse(substitute(zlab))


  # Point des résultats dans une liste
  midpoints <- function(x) (x[-1] + x[-length(x)])/2
  retval <- list()
  retval$Hist.Values <- m
  retval$breaks_x = x.cuts
  retval$breaks_y = y.cuts
  retval$breaks_z = z.cuts
  retval$Midpoints.x = midpoints(x.cuts)
  retval$Midpoints.y = midpoints(y.cuts)
  retval$Midpoints.z = midpoints(z.cuts)
  retval$nobs.x = length(x)
  retval$nobs.y = length(y)
  retval$nobs.z = length(z)
  retval$n.bins = c(length(x.cuts),length(y.cuts),length(z.cuts))
  retval$call <- match.call()
  retval
}

