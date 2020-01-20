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
#' @param numCov integer, number of covariates
#' @param permute logical, permute the SNP genotypes
#'
#' @return estimate of TAME
#'
#' @export
computeTAME <- function(snp,
                        expression,
                        mediators,
                        covs,
                        numMed,
                        numCov,
                        permute = F){

  snp = as.vector(snp)
  d = as.data.frame(cbind(as.vector(expression),
                          as.vector(snp),
                          mediators,
                          covs))
  colnames(d) = c('Gene',
                  'SNP',
                  paste0('Med',1:numMed),
                  paste0('Cov',1:numCov))

  if (permute){
    d$SNP = sample(d$SNP)
  }

  reg.tot = lm(Gene ~ .,d)

  beta_M = as.vector(reg.tot$coefficients[grepl('Med',
                                                names(reg.tot$coefficients))])
  alpha_X = vector(length = numMed,
                   mode = 'numeric')

  for (j in 1:numMed){

    med.d = d[,c(paste0('Med',j),'SNP',paste0('Cov',1:numCov))]
    colnames(med.d)[1] = 'Mediator'
    med.reg = lm(Mediator~.,med.d)
    alpha_X[j] = med.reg$coefficients[2]

  }

  TAME = as.numeric(abs(alpha_X) %*% abs(beta_M))
  return(TAME)

}
