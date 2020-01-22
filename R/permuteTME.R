#' Perform permutation test for the total absolute mediation effect
#'
#' The function takes in a SNP-mediator set-gene triplet
#' and computes the total absolute mediation effect and permutation
#' test P-value
#'
#' @param snp vector, SNP of interest
#' @param expression vector, gene expression of interest
#' @param mediators data frame, mediators of interest
#' @param covs data frame, covariates
#' @param nperms integer, number of permutations for the null distribution
#'
#' @return estimate of TAME and the permutation P-value
#'
#' @export
permuteTME = function(snp,
                       expression,
                       mediators,
                       covs,
                       nperms = 1000,
                       cores = 1){

  test.stat = computeTME(snp = snp,
                            expression = expression,
                            mediators = mediators,
                            covs = covs)

  boots = replicate(nperms,
                    sample(snp,replace=T))
  null.dist = apply(boots,
                    MARGIN = 2,
                    FUN = computeTME,
                    expression = expression,
                    mediators = mediators,
                    covs = covs)
  p = mean(abs(test.stat) <= abs(null.dist))

  return(list(test.stat = test.stat,
              p.value = p))

}
