trainSNPPheno <- function(pheno,
                          snpCur,
                          snpList,
                          thisSnp,
                          fileName,
                          prune = T,
                          windowSize = 50,
                          numSNPShift = 5,
                          ldThresh = .5){
  if (prune){
    pruneObj = LDprune(W = t(snpCur),
                       snpList = snpList,
                       snpLocs = thisSNP,
                       fileName = fileName,
                       windowSize = windowSize,
                       numSNPShift = numSNPShift,
                       ldThresh = ldThresh)
    snpCur = t(pruneObj$W)
    snpList = pruneObj$snpList
    thisSNP = pruneObj$onlyThese
    rm(pruneObj)}

  enet = glmnet::cv.glmnet(y = pheno,
                               x = t(snpCur),
                               nfolds = 5)
  blup = rrBLUP::mixed.solve(y = pheno,
                                 Z = t(snpCur),
                                 K = cov(t(snpCur)))
  return(list(enet = tot.enet,
              blup = tot.blup))
}
