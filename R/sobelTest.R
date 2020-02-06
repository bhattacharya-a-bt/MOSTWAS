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

  if (ncol(mediators) == 0){
    return(list(test.stat = 0,
                p.value = 1))
  }

  direct = lm(expression ~ snp + mediators + covs)
  indirect = lm(mediators ~ snp + covs)
  b = coef(direct)[grepl('mediators',names(coef(direct)))]
  if (ncol(mediators) > 1){
    a = coef(indirect)[2,]
  } else {
    a = coef(indirect)[2]
  }

  TME = a %*% b
  direct.cov = vcov(direct)[grepl('mediators',rownames(vcov(direct))),
                            grepl('mediators',colnames(vcov(direct)))]
  indirect.cov = vcov(indirect)[grepl('snp',rownames(vcov(indirect))),
                                grepl('snp',colnames(vcov(indirect)))]
  if (ncol(mediators) > 1){
    zeros = matrix(rep(0,nrow(direct.cov)*ncol(indirect.cov)),
                 nrow = nrow(direct.cov))
    V = rbind(cbind(direct.cov,zeros),
              cbind(t(zeros),indirect.cov))
  } else {
    V = rbind(c(direct.cov,0),
              c(0,indirect.cov))
  }
  cov = sqrt(as.numeric(t(c(b,a)) %*% V %*% c(b,a)))

  return(list(test.stat = TME,
              p.value = pnorm(abs(TME/cov))))

}
