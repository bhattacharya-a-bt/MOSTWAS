#' Simulate trait from latent genetic values
#'
#' The function takes in a vector of latent genetic
#' values and simulates a trait with controlled heritable variance
#' explained
#'
#' @param g vector, latent genetic values
#' @param h2g numeric, heritability of trait
#'
#' @return vector of trait
#'
#' @export
simTrait <- function(g,h2g){

  n = length(g)

  if (h2g > 0){

    s2g = var(g)
    s2e = s2g * ((1/h2g)-1)
    e = rnorm(n,mean = 0, sd = sqrt(s2e))
    y = g + e

  }  else {

    e = rnorm(n)
    y = e

  }

  y = as.numeric(scale(y))
  return(y)

}
