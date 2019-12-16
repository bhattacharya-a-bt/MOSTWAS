trainExpression <- function(geneInt,
                            snps,
                            snpLocs,
                            mediator,
                            medLocs,
                            covariates,
                            qtlFileNames,
                            numMed = 5,
                            seed,
                            k,
                            fileName,
                            cisDist = 1e6,
                            parallel = T,
                            prune = T,
                            windowSize = 50,
                            numSNPShift = 5,
                            ldThresh = .5,
                            cores = 5){

  set.seed(seed)
  parts.train = caret::createFolds(colnames(mediator)[-1],k=k,returnTrain = T)
  set.seed(seed)
  parts.test = caret::createFolds(colnames(mediator)[-1],k=k,returnTrain = F)

  medList = gatherMediators(geneInt,qtlFull,numMed)
  medTrainList = lapply(medList,
                        trainMediator,
                        mediator = mediator,
                        medLocs = medLocs,
                        snps = snps,
                        snpLocs = snpLocs,
                        covariates = covariates,
                        seed = seed,
                        k = k,
                        fileName = fileName,
                        cisDist = 1e6,
                        parallel = T,
                        prune = T,
                        windowSize = 50,
                        numSNPShift = 5,
                        ldThresh = .5,
                        cores = 5)
  names(medTrainList) = medList

  medTrainList = medTrainList[as.numeric(which(sapply(medTrainList,length) == 3))]
  fixedEffects = as.data.frame(matrix(ncol = length(medTrainList),
                                     nrow = ncol(mediator)-1))
  colnames(fixedEffects) = names(medTrainList)
  for (i in 1:ncol(fixedEffects)){
    fixedEffects[,i] = medTrainList[[i]][2]
  }

  fixedEffects$pheno = pheno

  ctrl = caret::trainControl(method = 'cv', number = 10)
  lmCVFit <- caret::train(pheno~.,
                          data = fixedEffects,
                          method = 'lm',
                          trControl = ctrl,
                          metrix = 'Rsquared')



}
