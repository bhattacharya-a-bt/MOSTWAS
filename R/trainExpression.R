#' Train and predict gene's predictive model with mediators
#'
#' The function trains a predictive model of a given gene using top mediators
#' as fixed effects and assesses in-sample performance with cross-validation.
#'
#' @param geneInt character, identifier for gene of interest
#' @param snps data frame, SNP dosages
#' @param snpLocs data frame, MatrixEQTL locations for SNPs
#' @param mediator data frame, mediator intensities
#' @param medLocs data frame, MatrixEQTL locations for mediators
#' @param covariates data frame, covariates
#' @param qtlFull data frame, all QTLs (cis and trans) between mediators and genes
#' @param integer numMed, number of top mediators to include
#' @param seed integer, random seed for splitting
#' @param k integer, number of training-test splits
#' @param fileName character, throw away name for PLINK files
#' @param parallel logical, TRUE/FALSE to run glmnet in parallel
#' @param prune logical, TRUE/FALSE to LD prune the genotypes
#' @param windowSize integer, window size for PLINK pruning
#' @param numSNPShift integer, shifting window for PLINK pruning
#' @param ldThresh numeric, LD threshold for PLINK pruning
#' @param cores integer, number of parallel cores
#' @param outputAll logical, include mediator information
#'
#' @return final model for gene along with CV R2 and predicted values
#'
#' @importFrom caret createFolds
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom glmnet cv.glmnet
#' @importFrom rrBLUP mixed.solve
#' @importForm parallel mclapply
#' @importFrom abind abind
#'
#' @export
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
                            cisDist = 5e5,
                            parallel = T,
                            prune = T,
                            windowSize = 50,
                            numSNPShift = 5,
                            ldThresh = .5,
                            cores,
                            outputAll = F){

  set.seed(seed)
  medList = gatherMediators(geneInt,qtlFull,numMed)
  if (parallel) {
  medTrainList = parallel::mclapply(medList,
           trainMediator,
           mediator = mediator,
           medLocs = medLocs,
           snps = snps,
           snpLocs = snpLocs,
           covariates = covariates,
           seed = seed,
           k = k,
           cisDist = cisDist,
           prune = prune,
           windowSize = windowSize,
           numSNPShift = numSNPShift,
           ldThresh = ldThresh,
           mc.cores = cores)
  }
  if (!parallel){
    medTrainList = lapply(medList,
                          trainMediator,
                          mediator = mediator,
                          medLocs = medLocs,
                          snps = snps,
                          snpLocs = snpLocs,
                          covariates = covariates,
                          seed = seed,
                          k = k,
                          cisDist = cisDist,
                          prune = prune,
                          windowSize = windowSize,
                          numSNPShift = numSNPShift,
                          ldThresh = ldThresh)
  }
  names(medTrainList) = medList

  pheno = as.numeric(mediator[mediator$Mediator == geneInt,-1])
  fe.R2 = 0
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
      fixedEffects$pheno = pheno
      ctrl = caret::trainControl(method = 'cv',
                                 number = k,
                                 savePredictions = 'final')
      lmCVFit <- caret::train(pheno~.,
                              data = fixedEffects,
                              method = 'lm',
                              trControl = ctrl,
                              metric = 'Rsquared')

      fe.R2 = fe.R2 + adjR2(lmCVFit$pred$pred,lmCVFit$pred$obs)
      pheno = as.numeric(resid(lmCVFit))

      trans.mod.df = as.data.frame(abind::abind(lapply(1:length(medTrainList),
                                                       amplifyTrans,
                                                       medTrainList = medTrainList,
                                                       lmCaretObj = lmCVFit),
                                                along = 1))
      trans.mod.df$Effect = as.numeric(as.character(trans.mod.df$Effect))
      trans.mod.df = subset(trans.mod.df,SNP != 'Intercept')
      rownames(trans.mod.df) = NULL
    }
  }




  cisGenoMod = trainMediator(medInt = geneInt,
                             pheno = pheno,
                             mediator = mediator,
                             medLocs = medLocs,
                             snps = snps,
                             snpLocs = snpLocs,
                             covariates = covariates,
                             seed = seed,
                             k = k,
                             cisDist = cisDist,
                             prune = prune,
                             windowSize = windowSize,
                             numSNPShift = numSNPShift,
                             ldThresh = ldThresh)

  cisGenoMod$Model$Mediator = 'Cis'
  if (exists('trans.mod.df')){
    cisGenoMod$Model = rbind(cisGenoMod$Model,trans.mod.df)
    }
  cisGenoMod$CVR2 = cisGenoMod$CVR2 + fe.R2
  cisGenoMod$CVR2.cis = cisGenoMod$CVR2 - fe.R2

  if (dir.exists('temp')){
    fff = list.files('temp/')
    cleanup = c(medList,geneInt)
    file.remove(paste0('temp/',fff[grepl(paste(cleanup,collapse = '|'),fff)]))
    }

  if (!outputAll){
  return(cisGenoMod)}

  if (outputAll){
    cisGenoMod$medList = medList
    cisGenoMod$medTrainList = medTrainList
    return(cisGenoMod)
  }


}

