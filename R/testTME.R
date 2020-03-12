#' Run TME tests through a data.frame of mediators
#'
#' The function runs Sobel or permutation testing
#'
#' @param i integer, row of data.frame
#' @param mediator data.frame, mediator intensities
#' @param qtMed data.frame, mediator QTL results
#' @param transSNPs vector, character vector of distal-eSNPs
#' @param sobel logical, T/F to use Sobel test
#' @param nperms integer, number of permutations
#' @param cores integer, number of parallel ocres
#' @param covariates data.frame, covariates
#' @param parallel character, boots input for parallelization
#'
#' @return list with TME and P-value
#'
#' @export
testTME <- function(i,
                    mediator,
                    qtMed,
                    transSNPs,
                    pheno,
                    sobel = F,
                    nperms = 1000,
                    cores,
                    covariates,
                    parallel = 'no'){

  thisMed = subset(mediator,
                   Mediator %in% qtMed$gene[qtMed$SNP == transSNPs$SNP[i]])

  if (nrow(thisMed) == 0){
    return(list(TME = 0,
                TME.P = 1))
  }

  if (nrow(thisMed) > 0){
    if (sobel == T){
        test = sobelTest(snp = as.numeric(as.vector(transSNPs[i,-1])),
                         expression = pheno,
                         mediators = t(as.matrix(thisMed[,-1])),
                         covs = t(as.matrix(covariates[,-1])))
    }
    if (sobel == F){
      test = permuteTME(snp = as.numeric(as.vector(transSNPs[i,-1])),
                        expression = pheno,
                        mediators = t(as.matrix(thisMed[,-1])),
                        covs = t(as.matrix(covariates[,-1])),
                        nperms = nperms,
                        parallel = parallel,
                        nc = cores)
    }
    return(list(TME = test$test.stat,
               TME.P = test$p.value))
  }
}
