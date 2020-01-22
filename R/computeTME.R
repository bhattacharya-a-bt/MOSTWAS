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
                        covs){

  numMed = ncol(mediators)
  snp = c(snp)

  snp = snp - c(covs %*% solve(t(covs) %*% covs) %*% t(covs) %*% snp)
  expression = expression - c(covs %*% solve(t(covs) %*% covs) %*% t(covs) %*% expression)
  mediators = t(limma::removeBatchEffect(t(mediators),
                                         covariates = covs))

  beta_M = c(solve(t(mediators) %*% mediators) %*% t(mediators) %*% expression)
  alpha_X = c(solve(t(snp) %*% snp) %*% t(snp) %*% mediators)

  TME = c(alpha_X %*% beta_M)
  return(TME)
}
