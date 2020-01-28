#' Compute total absolute mediation effect
#'
#' The function takes in a SNP-mediator set-gene triplet
#' and computes the total absolute mediation effect
#'
#' @param snp vector, SNP of interest
#' @param expression vector, gene expression of interest
#' @param mediators data frame, mediators of interest
#' @param covs data frame, covariates
#'
#' @return estimate of TAME
#'
#' @export
computeTME <- function(snp,
                       expression,
                       mediators,
                       covs,
                       indices){

  numMed = ncol(mediators)
  numCovs = ncol(covs)
  snp = c(snp)
  snp = snp[indices]
  tot = as.data.frame(cbind(expression,snp,mediators,covs))
  colnames(tot) = c('GE',
                    'SNP',
                    paste0('Med',1:numMed),
                    paste0('Cov',1:numCovs))

  tot.reg = lm(GE ~ ., data = tot)
  beta_M = as.numeric((coef(tot.reg)[paste0('Med',1:numMed)]))
  alpha_X = sapply(1:numMed,function(i) {
    cur = tot[,c(paste0('Med',i),'SNP',paste0('Cov',1:numCovs))]
    colnames(cur)[1] = 'Med'
    return(as.numeric((coef(lm(Med~.,data = cur))['SNP'])))
  }
  )

  return(as.numeric(alpha_X %*% beta_M))
}
