testTME <- function(i,
                    mediator,
                    qtMed,
                    transSNPs,
                    pheno,
                    sobel = F,
                    nperms = 1000,
                    cores,
                    covariates,
                    parallel = 'no',
                    thisMed){

  thisMed = subset(mediator,Mediator %in% qtMed$gene[qtMed$SNP == transSNPs$SNP[i]])
  if (nrow(thisMed) == 0){
    TME = 0
    TME.P = 1
    } else {
      if (sobel){
        test = sobelTest(snp = as.numeric(as.vector(transSNPs[i,-1])),
                         expression = pheno,
                         mediators = t(as.matrix(thisMed[,-1])),
                         covs = t(as.matrix(covariates[,-1])))
        } else {
          test = permuteTME(snp = as.numeric(as.vector(transSNPs[i,-1])),
                            expression = pheno,
                            mediators = t(as.matrix(thisMed[,-1])),
                            covs = t(as.matrix(covariates[,-1])),
                            nperms = nperms,
                            parallel = parallel,
                            nc = cores)
        }

      TME = test$test.stat
      TME.P = test$p.value
      }
  return(list(TME = TME,
              TME.P = TME.P))
}
