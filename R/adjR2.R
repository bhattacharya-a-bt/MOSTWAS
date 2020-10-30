#' Compute McNemar's adjusted R2
#'
#' The function computes the adjusted R2 between two numeric vectors
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param method character, correlation method
#'
#' @return adjusted R2 of x and y
#'
#'
#' @export
adjR2 <- function(x,
                  y,
                  method = 'pearson'){

  n = length(x)
  if (var(x) == 0 | var(y) == 0){
    return(0)
  }
  r2 = cor(x,y,method=method)^2
  ar2 = 1 - (1 - r2)*((n-1)/(n-2))
  return(max(0,ar2))

}
