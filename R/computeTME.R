#' Compute total absolute mediation effect
#'
#' The function takes in a SNP-mediator set-gene triplet
#' and computes the total absolute mediation effect
#'
#' @param snp vector, SNP of interest
#' @param expression vector, gene expression of interest
#' @param mediators data frame, mediators of interest
#' @param covs data frame, covariates
#'
#' @return estimate of TAME
#'
#' @export
computeTME <- function(snp,
                       expression,
                       mediators,
                       covs,
                       indices){

  numMed = ncol(mediators)
  snp = c(snp)
  snp = snp[indices]

  snp = t(limma::removeBatchEffect(t(snp),covariates = covs))
  expression = t(limma::removeBatchEffect(t(snp),covariates = covs))
  mediators = t(limma::removeBatchEffect(t(mediators),
                                         covariates = covs))

  beta_M = c(solve(t(mediators) %*% mediators) %*% t(mediators) %*% expression)
  alpha_X = c(solve(t(snp) %*% snp) %*% t(snp) %*% mediators)

  return(as.numeric(alpha_X %*% beta_M))
}
