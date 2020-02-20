#' Simulate eQTL reference panel
#'
#' The function creates a eQTL reference panel for MOSTWAS
#' methods with genotypes, expression, mediator intensities, and QTLs.
#'
#' @param geno.loc data.frame, reference genotypes for local SNPs
#' @param geno.dis data.frame, reference genotypes for distal SNPs
#' @param nqtl integer, sample size
#' @param b_qtls.loc vector, local eQTLs effect sizes
#' @param b_qtls.dis vector, distal eQTLs effect sizes
#' @param eqtl_h2.loc numeric, local heritability of expression
#' @param eqtl_h2.dis numeric, distal heritability of expression
#' @param numMed integer, number of mediators
#'
#' @return list of matrices of local and distal simulated genotypes, expression,
#' mediator intensities, vector of final distal eQTL effect sizes
#'
#' @export
simEQTL = function(geno.loc,
                    geno.dis,
                    nqtl,
                    b_qtls.loc,
                    b_qtls.dis,
                    eqtl_h2.loc,
                    eqtl_h2.dis,
                    numMed){

  ## loc-HERTIABLE PORTION
  Z_qtl.loc = simGeno(geno.loc, nqtl)
  n.loc = nrow(Z_qtl.loc)
  p.loc = ncol(Z_qtl.loc)
  # GRM AND LD
  A.loc = (Z_qtl.loc %*% t(Z_qtl.loc)) / p.loc
  LD_qtl.loc = (t(Z_qtl.loc) %*% Z_qtl.loc) / n.loc

  # sim gene expression
  gexpr.loc = simTrait(Z_qtl.loc %*% b_qtls.loc, eqtl_h2.loc)

  mediator.exp = matrix(nrow = length(gexpr.loc),
                        ncol = numMed)
  Z_qtl.dis = simGeno(geno.dis, nqtl)
  n.dis = nrow(Z_qtl.dis)
  p.dis = ncol(Z_qtl.dis)
  # GRM AND LD
  A.dis = (Z_qtl.dis %*% t(Z_qtl.dis)) / p.dis
  LD_qtl.dis = (t(Z_qtl.dis) %*% Z_qtl.dis) / n.dis
  # sim gene expression
  for (i in 1:numMed){
    mediator.exp[,i] = simTrait(Z_qtl.dis %*% b_qtls.dis[,i],
                                 eqtl_h2.dis)
  }
  w_med = rnorm(numMed,mean = 0,sd = sqrt(eqtl_h2.dis/numMed))
  gexpr.dis = c(mediator.exp %*% w_med) +
    rnorm(nrow(mediator.exp),
          mean = 0,
          sd = 1 - eqtl_h2.dis)
  fin.qtls.dis = b_qtls.dis %*% w_med
  gexpr = gexpr.loc + gexpr.dis

  return(list(Z.loc = Z_qtl.loc,
              Z.dis = Z_qtl.dis,
              exp = gexpr,
              med = mediator.exp,
              fin.qtls.dis = fin.qtls.dis))


}
