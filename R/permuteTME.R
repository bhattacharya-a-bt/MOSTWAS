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
                      parallel = 'no',
                      nc){

  a = boot::boot(data = snp,
                 statistic = computeTME,
                 R = nperms,
                 sim = 'permutation',
                 expression = expression,
                 mediators = mediators,
                 covs = covs,
                 parallel = parallel,
                 ncpus = nc)
  p = (sum(abs(a$t0) >= abs(a$t))+1)/(nperms + 1)
  return(list(test.stat = a$t0,
              p.value = p))

}
