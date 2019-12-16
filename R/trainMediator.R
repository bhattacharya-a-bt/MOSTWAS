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
#' @param fileName character, throw away name for PLINK files
#' @param parallel logical, TRUE/FALSE to run glmnet in parallel
#' @param prune logical, TRUE/FALSE to LD prune the genotypes
#' @param windowSize integer, window size for PLINK pruning
#' @param numSNPShift integer, shifting window for PLINK pruning
#' @param ldThresh numeric, LD threshold for PLINK pruning
#'
#' @return final model for mediator along with CV R2 and predicted values
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
trainMediator <- function(medInt,
                          pheno = NULL,
                          mediator,
                          medLocs,
                          snps,
                          snpLocs,
                          covariates,
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

  colnames(mediator)[1] = 'Mediator'
  if (is.null(pheno)){ pheno = as.numeric(mediator[mediator$Mediator == medInt,-1]) }

  res = as.data.frame(cbind(pheno,t(covariates[,-1])))
  pheno = as.numeric(resid(lm(pheno~.,data=res)))

  cisGeno = getCisGenotypes(biomInt = medInt,
                            locs = medLocs,
                            snps = snps,
                            snpLocs = snpLocs,
                            cisDist = 1e6)
  snpCur = cisGeno$snpCur
  snpList = cisGeno$snpList
  thisSNP = cisGeno$thisSNP

  if (prune){
    pruneObj = LDprune(W = t(snpCur),
                  snpList = snpList,
                  snpLocs = snpLocs,
                  fileName = fileName,
                  windowSize = windowSize,
                  numSNPShift = numSNPShift,
                  ldThresh = ldThresh)
  snpCur = t(pruneObj$W)
  snpList = pruneObj$snpList
  thisSNP = pruneObj$onlyThese
  rm(pruneObj)}

  data = as.data.frame(cbind(pheno,t(snpCur)))

  set.seed(seed)
  control = caret::trainControl(method = "cv",
                                number = 5,
                                savePredictions = 'final')
  suppressWarnings({
    cl <- parallel::makePSOCKcluster(cores)
    doParallel::registerDoParallel(cl)
    model.enet = caret::train(pheno~.,
                       data = data,
                       method = 'glmnet',
                       trControl=control,
                       tuneLength = 5,
                       metric = 'Rsquared')
    ParallelLogger::stopCluster(cl)
    })
  best.model = model.enet$finalModel
  best.lambda = model.enet$results$lambda[which.max(model.enet$results$Rsquared)]
  best.alpha = model.enet$results$alpha[which.max(model.enet$results$Rsquared)]
  pred = model.enet$pred
  pred = pred[order(pred$rowIndex),]
  r2 = max(model.enet$results$Rsquared)

  mod.df = data.frame(SNP = c(thisSNP$snpid),
                      Chromosome = c(thisSNP$chr),
                      Position = c(thisSNP$pos),
                      Effect = as.numeric(coef(best.model,s = best.lambda))[-1])
  mod.df = subset(mod.df,Effect!=0)
  return(list(Model = mod.df,
              Predicted = pred$pred,
              CVR2 = r2))

  }
