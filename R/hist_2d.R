
#' hist_2d
#'
#' This function results from a modification of the
#' \code{\link[gplots:hist2d]{hist2d}} function of the gplots package in order
#' to build the 2D histogram with breaks directly provided as inputs of the
#' new function. 
#'
#' @param x either a numerical vector to be partitioned or a matrix of 2 numerical columns to
#' be partitioned.
#' @param y a numerical vector to be partitioned. Not required if x is a matrix.
#' @param breaks_x a numerical vector. Contains the breaks related to x for the
#' histogram
#' @param breaks_y a numerical vector. Contains the breaks related to y for the
#' histogram
#' @param same.scale logical. Default to FALSE. If TRUE, breaks_x will be used
#' for x and y
#' @param na.rm logical. Default to TRUE. Indicates whether missing values
#' should be removed
#' @param FUN function used to summarize bin contents.
#'
#'
#' @details  The default function used for the argument FUN is the function
#' length. When another function is used, it is applied on x, or on the first column of x
#' if this is a two-column matrix. The lower limit of each class interval is
#' included in the class and the upper limit is not.
#'
#'
#' @return Using a given set of breaks per each variable, the function returns :
#' \describe{
#' \item{Hist.Values}{a matrix with in rows class intervals of x and in columns
#' class intervals of y. Contingency table is returned if FUN=length}
#' \item{breaks_x, breaks_y }{ same as the inputs of the function }
#' \item{Midpoints.x, Midpoints.y }{the midpoints for each bin per variable }
#' \item{nobs.x , nobs.y }{number of observations of x and y  }
#' \item{n.bins }{ vector of 2 elements containing the number of bins for
#' x and y }
#' }
#'
#' @examples
#' set.seed(45)
#' hist_2d(x = rnorm(1000),
#'         y = rnorm( 1000,5,0.1),
#'         breaks_x = seq(-4, 4, by =1),
#'         breaks_y = seq(2, 8, by =1))
#' @export


hist_2d <- function (x,
                     y,
                     breaks_x,
                     breaks_y,
                     same.scale = FALSE,
                     na.rm = TRUE,
                     FUN = base::length)  {
  # Define x and y
  if (is.null(y)) {
    if (ncol(x) != 2)
      stop("If y is ommitted, x must be a 2 column matrix")
    y <- x[, 2]
    x <- x[, 1]
  }

  # Management of NAs
  nas <- is.na(x) | is.na(y)
  if (na.rm) {
    x <- x[!nas]
    y <- y[!nas]
  }
  else stop("missing values not permitted if na.rm=FALSE")

  # Use of same scale
  if(same.scale){
    x.cuts = breaks_x;
    y.cuts = breaks_x;
  }else{
    x.cuts <- breaks_x
    y.cuts <- breaks_y
  }

  index.x <- cut(x, x.cuts, include.lowest = TRUE,right =FALSE,
                 labels = paste0("[",x.cuts[1:length(x.cuts)-1],
                                 ",",x.cuts[2:length(x.cuts)],"["))
  index.y <- cut(y, y.cuts, include.lowest = TRUE,right =FALSE,
                 labels = paste0("[",y.cuts[1:length(y.cuts)-1],
                                 ",",y.cuts[2:length(y.cuts)],"["))
  m <- tapply(x, list(index.x, index.y), FUN)

  # If length is used, replace NA with 0
  if (identical(FUN, base::length))  m[is.na(m)] <- 0
  # xlab, ylab
  # if (missing(xlab)) xlab <- deparse(substitute(xlab))
  # if (missing(ylab)) ylab <- deparse(substitute(ylab))

  # Point des résultats dans une liste
  midpoints <- function(x) (x[-1] + x[-length(x)])/2
  retval = list()
  retval$Hist.Values = m
  retval$breaks_x = x.cuts
  retval$breaks_y = y.cuts
  retval$Midpoints.x = midpoints(x.cuts)
  retval$Midpoints.y = midpoints(y.cuts)
  retval$nobs.x = length(x)
  retval$nobs.y = length(y)
  retval$n.bins = c(length(x.cuts),length(y.cuts))
  retval$call = match.call()
  retval
}

