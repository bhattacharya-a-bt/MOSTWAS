#' Generate summary statistics for genotype-trait association
#'
#' The function takes in simulated genotypes and trait
#' and outputs summary statistics for associations.
#'
#' @param Z matrix, simulated genotypes
#' @param pheno vector, trait
#'
#' @return data frame of effect sizes, standard errors and P-values of associations
#'
#' @export
regress <- function(Z,pheno){

  extractLM = function(x,pheno){

    a = summary(lm(pheno~x))
    return(coef(a)[2,])

  }

  require(pbapply)
  gwas = t(as.data.frame(apply(Z,2,extractLM,pheno = pheno)))
  colnames(gwas) = c('Beta','SE','t','P')
  gwas = as.data.frame(gwas)
  return(gwas[,-3])

}
