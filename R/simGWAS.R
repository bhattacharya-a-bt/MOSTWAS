#' Simulate eQTL reference panel
#'
#' The function creates a eQTL reference panel for MOSTWAS
#' methods with genotypes, expression, mediator intensities, and QTLs.
#'
#' @param geno.loc data.frame, reference genotypes for local SNPs
#' @param geno.dis data.frame, reference genotypes for distal SNPs
#' @param ngwas integer, sample size
#' @param b_qtls.loc vector, local eQTLs effect sizes
#' @param b_qtls.dis vector, distal eQTLs effect sizes
#' @param var_explained numeric, total variance explained in phenotype
#'
#' @return list of GWAS summary statistics, total eQTL effect size,
#' matrix of simulated genotypes, vector of trait
#'
#' @export
simGWAS <- function(geno.loc,
                    geno.dis,
                    ngwas,
                    b_qtls.loc,
                    b_qtls.dis,
                    var_explained){

  Z_gwas.cis = simGeno(geno.cis, ngwas)
  colnames(Z_gwas.cis) = paste0('Cis',1:ncol(Z_gwas.cis))
  Z_gwas.trans = simGeno(geno.trans, ngwas)
  colnames(Z_gwas.trans) = paste0('Trans',1:ncol(Z_gwas.trans))

  #var_explained only reflects that due to genetics
  gwas_expr = cbind(Z_gwas.cis,Z_gwas.trans) %*% c(b_qtls.cis,b_qtls.trans)

  if (var_explained > 0){
    alpha = rnorm(1)
  }  else {alpha = 0}

  y = simTrait(gwas_expr * alpha, var_explained)
  Z_gwas = cbind(Z_gwas.cis,Z_gwas.trans)
  gwas = regress(Z_gwas,y)
  return(list(gwas = gwas, alpha = alpha, Z_gwas = Z_gwas,y=y))
}
