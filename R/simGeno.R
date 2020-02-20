#' Simulate genotypes from reference data
#'
#' The function takes in a reference genotype data
#' and generates simulated genotypes from the LD structure.
#' Modified from Nicholas Mancuso's twas_sim.
#'
#' @param geno data.frame, reference panel SNPs
#' @param n integer, samples size
#'
#' @return matrix of simulated genotypes
#'
#' @export
simGeno <- function(geno,n){

  ## estimate LD
  G = as.matrix(geno[,-1])
  require(pbapply)
  mafs = pbapply(G,1,function(x) mean(x)/2)
  G = pbapply(G,1,scale)
  #LD = (t(G) %*% G) / nrow(G) + diag(ncol(G)) * 0.1
  LD = (t(G) %*% G) / nrow(G) + diag(ncol(G))
  L = chol(LD)
  Z = t(L %*% t(matrix(rnorm(ncol(G) * n),ncol = ncol(G))))
  Z = as.matrix(scale(Z))

  return(Z)

}
