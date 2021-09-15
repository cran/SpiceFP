

#' getD3dSparse
#'
#' getD3dSparse is a function that helps to construct generalized lasso
#' penalty matrix D when using the \code{\link[genlasso:fusedlasso]{fusedlasso}}
#' function over a 3 dimensional grid
#'
#' @param dim1 positive integer. Based on a 3 dimensional grid, dim1 represents
#'  the number of units represented on the first dimension
#' @param dim2 positive integer which represents the number of units represented
#'  on the second dimension
#' @param dim3 positive integer which represents the number of units represented
#'  on the third dimension
#'
#' @details  The function returns a sparse penalty matrix providing information
#' on the connections between the variables during the implementation of a
#' generalizad fused lasso.
#'
#'
#' @return a matrix with dim1 x dim2 x dim3 columns. Each row represents an
#' edge (a link between 2 variables) and is constructed with the couple (-1, 1),
#'  relative to these 2 variables and 0 for all others. In the context of a
#'  generalized fused lasso, this matrix penalizes only the differences in
#'  coefficients (fusion). To obtain parsimony in addition to the fusion, a
#'  diagonal matrix with the same number of columns must be bound to the
#'  penalty matrix constructed by getD3dSparse. This matrix will contain
#'  diagonally the ratio: parsimony penalty parameter on fusion penalty
#'  parameter. When using \code{\link[genlasso:fusedlasso]{fusedlasso}}
#'  function, this operation is performed when you provide the argument gamma.
#'
#' @importFrom Matrix bandSparse rowSums
#' @importFrom genlasso getGraph
#'
#' @examples
#' library(genlasso)
#' library(Matrix)
#' D<-getD3dSparse(2,3,2)
#' plot(getGraph(D))
#'
#' @export

getD3dSparse <- function (dim1, dim2, dim3) {
  D1 = bandSparse(dim1 * dim2 * dim3, m = dim1 * dim2 * dim3, k = c(0, 1),
                  diagonals = list(rep(-1, dim1 * dim2 * dim3),
                                   rep(1, dim1 * dim2 * dim3 - 1)))
  D2 = bandSparse(dim1 * dim2 * dim3 , m = dim1 * dim2 * dim3,
                  k = c(0, dim1),
                  diagonals = list(rep(-1, dim1 * dim2 * dim3),
                                   rep(1,  dim1 * dim2 * dim3 - 1)))
  D3 = bandSparse(dim1 * dim2 * dim3 , m = dim1 * dim2 * dim3,
                  k = c(0, dim1 * dim2),
                  diagonals = list(rep(-1, dim1 * dim2 * dim3),
                                   rep(1,  dim1 * dim2 * dim3 - 1)))
  D <- rbind(D1, D2, D3)
  D1 <- D[-c(which(rowSums(D) != 0)),]
  return(D1)
}



