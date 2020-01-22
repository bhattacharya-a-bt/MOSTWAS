  #' Compute total absolute mediation effect
#'
#' The function takes in a SNP-mediator set-gene triplet
#' and computes the total absolute mediation effect
#'
#' @param snp vector, SNP of interest
#' @param expression vector, gene expression of interest
#' @param mediators data frame, mediators of interest
#' @param covs data frame, covariates
#' @param numMed integer, number of mediators
#' @param permute logical, permute the SNP genotypes
#'
#' @return estimate of TAME
#'
#' @export
computeTAME <- function(snp,
                        expression,
                        mediators,
                        covs,
                        permute = F){

  numMed = ncol(mediators)
  snp = c(snp)

  if (permute){ snp = sample(snp,replace = T) }

  snp = snp - c(covs %*% solve(t(covs) %*% covs) %*% t(covs) %*% expression)
  expression = expression - c(covs %*% solve(t(covs) %*% covs) %*% t(covs) %*% expression)
  for (i in 1:numMed){
    m = c(mediators[,i])
    mediators[,i] = m - c(covs %*% solve(t(covs) %*% covs) %*% t(covs) %*% m)
  }

  d = as.data.frame(cbind(expression,
                          snp,
                          mediators))
  colnames(d) = c('Gene',
                  'SNP',
                  paste0('Med',1:numMed))

  reg.tot = lm(Gene ~ .,d)
  beta_M = c(coefficients(reg.tot)[-c(1:2)])
  alpha_X = c(solve(t(snp) %*% snp) %*% t(snp) %*% mediators)

  TME = as.numeric(alpha_X %*% beta_M)
  return(TME)


}
