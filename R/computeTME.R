#' Compute total mediation effect
#'
#' The function takes in a SNP-mediator set-gene triplet
#' and computes the total mediation effect with boots structure
#'
#' @param snp vector, SNP of interest
#' @param expression vector, gene expression of interest
#' @param mediators data frame, mediators of interest
#' @param covs data frame, covariates
#' @param indices blank, index for boot package
#'
#' @return estimate of TME
#'
#' @export
computeTME <- function(snp,
                       expression,
                       mediators,
                       covs,
                       indices){
  snp = c(snp)
  snp = snp[indices]
  direct = lm(expression ~ snp + mediators + covs)
  indirect = lm(mediators ~ snp + covs)

  TME =  coef(direct)[3:7] %*% coef(indirect)[2,]

  return(as.numeric(TME))
}
