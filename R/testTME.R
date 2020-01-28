testTME <- function(i,
                    mediator,
                    qtMed,
                    transSNPs,
                    pheno,
                    nperms = 1000,
                    cores,
                    covariates,
                    thisMed){

  thisMed = subset(mediator,Mediator %in% qtMed$gene[qtMed$SNP == transSNPs$SNP[i]])
  if (nrow(thisMed) == 0){
    TME = 0
    TME.P = 1
  }
  if (nrow(thisMed) > 0){
    test = permuteTME(snp = as.numeric(as.vector(transSNPs[i,-1])),
                      expression = pheno,
                      mediators = t(as.matrix(thisMed[,-1])),
                      covs = t(as.matrix(covariates[,-1])),
                      nperms = 1000,
                      nc = 1)
    TME = test$test.stat
    TME.P = test$p.value
  }
  return(list(TME = TME,
              TME.P = TME.P))
}
