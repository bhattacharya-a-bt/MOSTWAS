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
#' @param h2Pcutoff numeric, P-value cutoff for heritability
#' @param numMed integer, number of top mediators to include
#' @param seed integer, random seed for splitting
#' @param k integer, number of training-test splits
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
                            dimNumeric,
                            qtlFull,
                            h2Pcutoff = 0.1,
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
                            verbose = T,
                            LDMS = F,
                            modelDir,
                            ldScrRegion = 200,
                            snpAnnot = NULL){

  set.seed(seed)

  print('GATHERING MEDIATORS')
  pheno = as.numeric(mediator[mediator$Mediator == geneInt,-1])
  medList = gatherMediators(geneInt,qtlFull,numMed)

  print('ESTIMATING HERITABILITY')
  herit = estimateHeritability(biomInt = geneInt,
                               snps = snps,
                               pheno = pheno,
                               snpLocs = snpLocs,
                               mediator = mediator,
                               medLocs = medLocs,
                               covariates = covariates,
                               fileName = geneInt,
                               dimNumeric = dimNumeric,
                               numMed = numMed,
                               cisDist = 5e5,
                               needMed = T,
                               medList = medList,
                               verbose = verbose,
                               windowSize = windowSize,
                               numSNPShift = numSNPShift,
                               ldThresh = ldThresh,
                               ldScrRegion = ldScrRegion,
                               LDMS = LDMS,
                               supplySNP = F,
                               prune = prune,
                               snpAnnot = snpAnnot)


  if (herit$P > h2Pcutoff | herit$h2 <= 0) {
    return(paste(geneInt,
                 'is not germline heritable at P <',
                 h2Pcutoff))
    }

    print('FITTING CIS-GENOTYPES')
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
                               ldThresh = ldThresh,
                               snpAnnot = snpAnnot)
    pheno = pheno - cisGenoMod$Predicted


    print('TRAINING MEDIATORS')

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


  print('FITTING MEDIATORS')
  fe.R2 = 0
  if (length(medTrainList) > 0){
    medTrainList = medTrainList[as.numeric(which(sapply(medTrainList,
                                                      function(x) x[3]) >= .01))]
    if (length(medTrainList) > 0){
      fixedEffects = as.data.frame(matrix(ncol = length(medTrainList),
                                          nrow = ncol(mediator)-1))
      colnames(fixedEffects) = names(medTrainList)
      for (i in 1:ncol(fixedEffects)){
        fixedEffects[,i] = medTrainList[[i]][2]
      }
      fixedEffects = as.data.frame(apply(fixedEffects,2,scale))
      fixedEffects$pheno = pheno
      lmCVFit <- lm(pheno~.,
                    data = fixedEffects)

      fe.R2 = fe.R2 + adjR2(as.numeric(predict(lmCVFit)),pheno)

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




  cisGenoMod$Model$Mediator = 'Cis'
  if (exists('trans.mod.df')){
    cisGenoMod$Model = rbind(cisGenoMod$Model,trans.mod.df)
    }
  cisGenoMod$CVR2 = cisGenoMod$CVR2 + fe.R2
  cisGenoMod$CVR2.cis = cisGenoMod$CVR2 - fe.R2
  cisGenoMod$h2 = herit$h2
  cisGenoMod$h2.P = herit$P

  if (dir.exists('temp')){
    fff = list.files('temp/')
    cleanup = c(medList,geneInt)
    file.remove(paste0('temp/',fff[grepl(paste(cleanup,collapse = '|'),fff)]))
  }


  Model = cisGenoMod$Model
  R2 = cisGenoMod$CVR2
  Predicted = cisGenoMod$Predicted
  Mediators = cisGenoMod$medlist
  CisR2 = cisGenoMod$CVR2.cis
  h2 = abs(herit$h2)
  h2.Pvalue = herit$P
  save(Model,R2,Predicted,Mediators,CisR2,h2,h2.Pvalue,
         file = paste0(modelDir,geneInt,'.wgt.med.RData'))
  }

