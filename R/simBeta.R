#' Simulate eQTL effect sizes
#'
#' The function takes in a proportion
#' of causal eQTLs and total heritability of eQTLs
#' and generates n_snps eQTL effect sizes
#'
#' @param p.causal numeric, proportion of causal eQTLs
#' @param eqtl_h2 numeric, heritability of expression from eQTLs
#' @param n_snps integer, number of eQTLs
#'
#' @return vector of effect sizes
#'
#' @export
simBeta <- function(p.causal, eqtl_h2, n_snps){

  # number of QTLs
  n_qtls = max(1,floor(p.causal * n_snps))

  # select which SNPs are causal
  c_qtls = sample(1:n_snps,n_qtls)
  b_qtls = rep(0,n_snps)

  # sample effects from normal prior
  b_qtls[c_qtls] = rnorm(n_qtls,
                         mean = 0,
                         sd = sqrt(eqtl_h2/n_qtls))

  return(b_qtls)

}
