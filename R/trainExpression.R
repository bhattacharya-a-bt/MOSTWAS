trainExpression <- function(geneInt,
                            snps,
                            snpLocs,
                            mediator,
                            medLocs,
                            covariates,
                            qtlFull,
                            numMed = 5,
                            seed,
                            k,
                            fileName,
                            cisDist = 1e6,
                            parallel = T,
                            prune = T,
                            windowSize = 50,
                            numSNPShift = 5,
                            ldThresh = .5){

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
                        cisDist = cisDist,
                        parallel = parallel,
                        prune = parallel,
                        windowSize = windowSize,
                        numSNPShift = numSNPShift,
                        ldThresh = ldThresh,
                        cores = cores)
  names(medTrainList) = medList
  foreach::registerDoSEQ()

  if (length(medTrainList) >= 0){
    medTrainList = medTrainList[as.numeric(which(sapply(medTrainList,
                                                      function(x) x[3]) >= .01))]
    if (length(medTrainList) >= 0){
      fixedEffects = as.data.frame(matrix(ncol = length(medTrainList),
                                          nrow = ncol(mediator)-1))
      colnames(fixedEffects) = names(medTrainList)
      for (i in 1:ncol(fixedEffects)){
        fixedEffects[,i] = medTrainList[[i]][2]
        }
      pheno = as.numeric(mediator[mediator$Mediator == geneInt,-1])
      fixedEffects$pheno = pheno
      ctrl = caret::trainControl(method = 'cv',
                                 number = k,
                                 savePredictions = 'final')
      lmCVFit <- caret::train(pheno~.,
                              data = fixedEffects,
                              method = 'lm',
                              trControl = ctrl,
                              metric = 'Rsquared')

      fe.R2 = adjR2(lmCVFit$pred$pred,lmCVFit$pred$obs)
      pheno = as.numeric(resid(lmCVFit))
    }
  }



  trans.mod.df = as.data.frame(abind::abind(lapply(1:length(medTrainList),
                                     amplifyTrans,
                                     medTrainList = medTrainList,
                                     lmCaretObj = lmCVFit),
                              along = 1))
  trans.mod.df$Effect = as.numeric(as.character(trans.mod.df$Effect))
  trans.mod.df = subset(trans.mod.df,SNP != 'Intercept')
  rownames(trans.mod.df) = NULL

  cisGenoMod = trainMediator(medInt = geneInt,
                             pheno = pheno,
                             mediator = mediator,
                             medLocs = medLocs,
                             snps = snps,
                             snpLocs = snpLocs,
                             covariates = covariates,
                             seed = seed,
                             k = k,
                             fileName = geneInt,
                             cisDist = 1e6,
                             parallel = T,
                             prune = T,
                             windowSize = 50,
                             numSNPShift = 5,
                             ldThresh = .5,
                             cores = 5)

  cisGenoMod$Model$Mediator = 'Cis'
  cisGenoMod$Model = rbind(cisGenoMod$Model,trans.mod.df)
  cisGenoMod$CVR2 = cisGenoMod$CVR2 + fe.R2

  cleanup = list.files('temp/')
  file.remove(paste0('temp/',cleanup))

  return(cisGenoMod)


}

