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
#' @param numMed integer, number of mediators
#' @param numCov integer, number of covariates
#' @param nperms integer, number of permutations for the null distribution
#'
#' @return estimate of TAME and the permutation P-value
#'
#' @export
permuteTAME = function(snp,
                       expression,
                       mediators,
                       covs,
                       nperms = 1000,
                       cores = 1){

  test.stat = computeTAME(snp = snp,
                            expression = expression,
                            mediators = mediators,
                            covs = covs,
                            permute = F)

  null.dist = RepParallel(nperms,
                          computeTAME(snp = snp,
                                      expression = expression,
                                      mediators = mediators,
                                      covs = covs,
                                      permute = T),
                          mc.cores = cores)

  p = mean(abs(test.stat) <= abs(null.dist))

  return(list(test.stat = test.stat,
              p.value = p))

}
