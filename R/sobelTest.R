#' Compute total mediation effect
#'
#' The function takes in a SNP-mediator set-gene triplet
#' and computes the total mediation effect and tests
#' the total mediator effect via a Sobel test (1982)
#'
#' @param snp vector, SNP of interest
#' @param expression vector, gene expression of interest
#' @param mediators data frame, mediators of interest
#' @param covs data frame, covariates
#'
#' @return estimate of TME and Sobel test P-value
#'
#' @export
sobelTest <- function(snp,
                       expression,
                       mediators,
                       covs){


  direct = lm(expression ~ snp + mediators + covs)
  indirect = lm(mediators ~ snp + covs)
  b = coef(direct)[3:7]
  a = coef(indirect)[2,]

  TME = a %*% b
  direct.cov = vcov(direct)[3:7,3:7]
  indirect.cov = vcov(indirect)[grepl('snp',rownames(vcov(indirect))),
                                grepl('snp',colnames(vcov(indirect)))]
  zeros = matrix(rep(0,nrow(direct.cov)*ncol(indirect.cov)),
                 nrow = nrow(direct.cov))
  V = rbind(cbind(direct.cov,zeros),
            cbind(t(zeros),indirect.cov))
  cov = sqrt(as.numeric(t(c(b,a)) %*% V %*% c(b,a)))

  return(list(test.stat = TME,
              p.value = pnorm(TME/cov)))

}
