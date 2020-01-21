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
                        permute = F,
                        mc.p = F,
                        mc.rep = 20000){

  snp = c(snp)

  if (permute){ snp = sample(snp,replace = T) }

  expression = expression - c(covs %*% solve(t(covs) %*% covs) %*% t(covs) %*% expression)
  expression = expression - c(snp %*% solve(t(snp) %*% snp) %*% t(snp) %*% expression)
  for (i in 1:numMed){

    m = c(mediators[,i])
    mediators[,i] = c(covs %*% solve(t(covs) %*% covs) %*% t(covs) %*% m)

  }

  d = as.data.frame(cbind(expression,
                          snp,
                          mediators))
  colnames(d) = c('Gene',
                  'SNP',
                  paste0('Med',1:numMed))

  reg.tot = lm(Gene ~ .,d)
  beta_M = c(coefficients(reg.tot)[-c(1:2)])
  alpha_X = c(solve(t(snp) %*% snp) %*% t(snp) %*% mediators)

  TME = as.numeric(alpha_X %*% beta_M)
  mcmc.p = NA

  if (mc.p){


    cov_beta = vcov(reg.tot)[paste0('Med',1:numMed),paste0('Med',1:numMed)]/nrow(mediators)
    cov_alpha = (t(mediators) %*%
                   (diag(rep(1,nrow(mediators))) - snp %*%
                      solve(t(snp) %*% snp) %*% t(snp)) %*% mediators)/(nrow(mediators)-1)
    conf = 95
    pest = c(alpha_X,beta_M)
    acov = rbind(cbind(cov_alpha,matrix(rep(0,nrow(cov_alpha) * ncol(cov_beta)),
                                  nrow = nrow(cov_alpha))),
                 cbind(matrix(rep(0,nrow(cov_beta) * ncol(cov_beta)),
                              nrow = nrow(cov_beta)),cov_beta))
    mcmc <- MASS::mvrnorm(mc.rep,pest,acov,empirical=FALSE)
    ab <- vector(mode='numeric',length = mc.rep)
    for (i in 1:numMed){
      ab = ab + mcmc[,i]*mcmc[,i+numMed]
    }
    mcmc.p = mean(abs(TME) > abs(ab))

  }

  ifelse(mc.p,
         return(list(TME = TME, MCMC.P = mcmc.p)),
         return(list(TME = TME)))


}
