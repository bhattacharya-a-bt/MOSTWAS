#' Parallel implementation of replicate()
#'
#' The function parallelizes replicate(), as in dbframe-R-library
#'
#' @param n number of replications
#' @param expr expression
#' @param simplify character, simplify to array
#'
#' @return array of replications
#'
#'
#' @export
RepParallel <- function(n, expr, simplify = "array") {
  answer <-
    mclapply(integer(n), eval.parent(substitute(function(...) expr)),...)
  if (!identical(simplify, FALSE) && length(answer))
    return(simplify2array(answer, higher = (simplify == "array")))
  else return(answer)
}
