
#' logbreaks
#'
#' A function that allows to obtain histogram class limits following a
#' logarithmic scale. It also has a parameter that allows to set the scale at
#' your convenience.
#'
#' @param x either a numeric vector to be partitioned or a numeric vector
#' containing the minimum and maximum of the vector to be partitioned.
#' @param parlist a list of 2 elements. The first one is alpha, a numeric and
#' positive value. It is a parameter affecting the number of breaks closed to
#' the minimum. The second one is J. It is a nonnegative and nonzero integer and
#' represent the selected number of classes.
#' @param round_breaks a nonnegative integer. Equal to 0 by default, it is the
#' number of decimal values of the breaks.
#' @param plot_breaks logical. FALSE by default. If TRUE, the breaks are
#' plotted.
#' @param effect.threshold.begin NA by default. Numeric value between the
#' minimum and maximum of x. If it isn't NA, the first class is created with
#' xmin and effect.threshold.begin.
#' @param effect.threshold.end NA by default. Numeric value between the minimum
#' and maximum of x. If it isn't NA, the last class is created with xmax and
#' effect.threshold.end.
#'
#' @details  The breaks are obtained as follows: 
#'
#' \deqn{L(w) = \min(x) + \frac{e^{\alpha \frac{w-1}{J}} - 1 }{e^{\alpha}-1} 
#' (\max(x) -\min(x)), \  w= 1, \ldots, J+1.}
#' 
#'
#' @return The return is a numeric vector of length J+1 with the breaks
#' obtained following a log scale.
#'
#' @importFrom graphics plot
#'
#' @examples
#' logbreaks(c(10,1000), parlist=list(0.2,5))
#' logbreaks(c(10,1000), parlist=list(0.2,5),plot_breaks=TRUE)
#' @export


logbreaks = function(x,
                     parlist=list(alpha,J),
                     round_breaks = 0,
                     plot_breaks = FALSE,
                     effect.threshold.begin=NA,
                     effect.threshold.end=NA){

  try(if(length(x)<2 & mode(x) != "numeric")
    stop("x must be a numeric vector having at least two different elements"))

  xmin= min(x)
  xmax= max(x)
  alpha= parlist[[1]]
  J= parlist[[2]]

  # Set J according to the effect.threshold
  if(is.na(effect.threshold.begin)+is.na(effect.threshold.end)==2){J=J}
  else if(is.na(effect.threshold.begin)+is.na(effect.threshold.end)==1){J=J-1}
  else if(is.na(effect.threshold.begin)+is.na(effect.threshold.end)==0){J=J-2}

  # Index of breaks
  j=seq(1,J+1,length.out =J+1)

  # Fix effect.threshold if NA
  if(is.na(effect.threshold.begin)){effect.threshold.begin=xmin}
  if(is.na(effect.threshold.end)){effect.threshold.end=xmax}

  # Compute the breaks
  x_j= effect.threshold.begin +
       ((1/alpha)*((1+
       alpha* (effect.threshold.end - effect.threshold.begin))^( (j-1)/J)-1))
  x_j=unique(round(c(xmin,x_j,xmax),digits = round_breaks))

  # Plot
  if (plot_breaks == TRUE) {
    plot(x_j)
  }

  return(x_j)
}
