trainSNPPheno <- function(pheno,
                          snpCur,
                          snpList,
                          thisSNP,
                          fileName,
                          prune = T,
                          windowSize = 50,
                          numSNPShift = 5,
                          ldThresh = .5,
                          verbose = T){
  if (prune){
    pruneObj = LDprune(W = t(snpCur),
                       snpList = snpList,
                       snpLocs = thisSNP,
                       fileName = fileName,
                       windowSize = windowSize,
                       numSNPShift = numSNPShift,
                       ldThresh = ldThresh,
                       verbose = verbose)
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
  return(list(enet = enet,
              blup = blup))
}