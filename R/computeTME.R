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

  if (ncol(mediators) == 0){return(0)}

  direct = lm(expression ~ snp + mediators + covs)
  indirect = lm(mediators ~ snp + covs)
  b = coef(direct)[grepl('mediators',names(coef(direct)))]
  if (ncol(mediators) > 1){
    a = coef(indirect)[2,]
  } else {
    a = coef(indirect)[2]
  }

  TME = a %*% b

  return(as.numeric(TME))
}
