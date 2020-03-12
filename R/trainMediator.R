#' Train and predict mediator predictive model
#'
#' The function trains a predictive model of a given mediator and then
#' predicts the genetically regulated intensity of it in the training set via
#' cross-validation.
#'
#' @param medInt character, identifier for mediator of interest
#' @param mediator data frame, mediator intensities
#' @param medLocs data frame, MatrixEQTL locations for mediators
#' @param snps data frame, SNP dosages
#' @param snpLocs data frame, MatrixEQTL locations for SNPs
#' @param covariates data frame, covariates
#' @param seed integer, random seed for splitting
#' @param k integer, number of training-test splits
#' @param parallel logical, TRUE/FALSE to run glmnet in parallel
#' @param prune logical, TRUE/FALSE to LD prune the genotypes
#' @param windowSize integer, window size for PLINK pruning
#' @param numSNPShift integer, shifting window for PLINK pruning
#' @param ldThresh numeric, LD threshold for PLINK pruning
#' @param cores integer, number of parallel cores
#'
#' @return final model for mediator along with CV R2 and predicted values
#'
#' @importFrom caret createFolds
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom glmnet cv.glmnet
#' @importFrom rrBLUP mixed.solve
#' @importFrom parallel mclapply
#' @importFrom abind abind
#'
#' @export
trainMediator <- function(medInt,
                          pheno = NULL,
                          mediator,
                          medLocs,
                          snps,
                          snpLocs,
                          covariates,
                          seed,
                          k,
                          cisDist = 5e5,
                          prune = T,
                          windowSize = 50,
                          numSNPShift = 5,
                          ldThresh = .5,
                          snpAnnot = NULL){

  fileName = medInt
  colnames(mediator)[1] = 'Mediator'
  if (is.null(pheno)){ pheno = as.numeric(mediator[mediator$Mediator == medInt,-1]) }

  if (!is.null(covariates)){
  res = as.data.frame(cbind(pheno,t(covariates[,-1])))
  pheno = as.numeric(resid(lm(pheno~.,data=res)))}

  cisGeno = getCisGenotypes(biomInt = medInt,
                            locs = medLocs,
                            snps = snps,
                            snpLocs = snpLocs,
                            cisDist = cisDist)
  snpCur = cisGeno$snpCur
  snpList = cisGeno$snpList
  thisSNP = cisGeno$thisSNP

  if (prune){
    if (ncol(snpCur) == length(pheno)){
      snpCur = t(snpCur)
    }
    pruneObj = LDprune(W = t(snpCur),
                  snpList = snpList,
                  snpLocs = snpLocs,
                  fileName = fileName,
                  windowSize = windowSize,
                  numSNPShift = numSNPShift,
                  ldThresh = ldThresh,
                  snpAnnot = snpAnnot)
  snpCur = t(pruneObj$W)
  if (length(pheno) == ncol(pruneObj$W)){
    snpCur = t(snpCur)
  }
  snpList = pruneObj$snpList
  thisSNP = pruneObj$onlyThese
  rm(pruneObj)}

  data = as.data.frame(cbind(pheno,t(snpCur)))
  data = data[,!duplicated(colnames(data))]

  set.seed(seed)
  train = caret::createFolds(y = pheno,
                             k=k,
                             returnTrain = T)
  set.seed(seed)
  test = caret::createFolds(y = pheno,
                            k = k,
                            returnTrain = F)

  pred.blup = pred.enet = vector(mode = 'numeric',length = length(pheno))
  for (i in 1:k){

    blup = rrBLUP::mixed.solve(y = pheno[train[[i]]],
                               Z = t(snpCur[,train[[i]]]))
    pred.blup[test[[i]]] = as.numeric(t(snpCur[,test[[i]]]) %*% blup$u)
    enet = glmnet::cv.glmnet(y = pheno[train[[i]]],
                             x = t(snpCur[,train[[i]]]),
                             nfolds = 5)
    pred.enet[test[[i]]] = as.numeric(predict(enet,
                                   newx = t(snpCur[,test[[i]]]),
                                   s = 'lambda.min'))

  }

  model = ifelse(adjR2(pheno,pred.blup) >= adjR2(pheno,pred.enet),
                 'LMM','Elastic net')

  fin.model.enet = glmnet::cv.glmnet(y = pheno,
                                     x = t(snpCur),
                                     nfolds = 5)
  mod.df.enet = data.frame(SNP = c(thisSNP$snpid),
                      Chromosome = c(thisSNP$chr),
                      Position = c(thisSNP$pos),
                      Effect = as.numeric(coef(fin.model.enet,s = 'lambda.min'))[-1])
  mod.df.enet = subset(mod.df.enet,Effect != 0)

  fin.model.blup = rrBLUP::mixed.solve(y = pheno,
                                       Z = t(snpCur))
  mod.df.blup = data.frame(SNP = c(thisSNP$snpid),
                           Chromosome = c(thisSNP$chr),
                           Position = c(thisSNP$pos),
                           Effect = as.numeric(fin.model.blup$u))

  if (model == 'Elastic net' & nrow(mod.df.enet) != 0){

    return(list(Model = mod.df.enet,
                Predicted = pred.enet,
                CVR2 = adjR2(pheno,pred.enet)))
  }

  if (model == 'LMM' | nrow(mod.df.enet) == 0){

    return(list(Model = mod.df.blup,
                Predicted = pred.blup,
                CVR2 = adjR2(pheno,pred.blup)))
  }
  }
