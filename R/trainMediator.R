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

  if (parallel){
    doParallel::registerDoParallel(5)
  }

  set.seed(seed)
  parts.train = caret::createFolds(colnames(mediator)[-1],k=k,returnTrain = T)
  set.seed(seed)
  parts.test = caret::createFolds(colnames(mediator)[-1],k=k,returnTrain = F)

  colnames(mediator)[1] = 'Mediator'
  pheno = as.numeric(mediator[mediator$Mediator == medInt,-1])

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

  if (prune == T){
    prune = LDprune(W = t(snpCur),
                  snpList = snpList,
                  snpLocs = snpLocs,
                  fileName = fileName,
                  windowSize = windowSize,
                  numSNPShift = numSNPShift,
                  ldThresh = ldThresh)
  snpCur = t(prune$W)
  snpList = prune$snpList
  thisSNP = prune$onlyThese
  rm(prune)}

  if (parallel){
    out = abind::abind(parallel::mclapply(1:k,
                                          CVFit,
                                          parts.train=parts.train,
                                          parts.test=parts.test,
                                          pheno=pheno,
                                          snpCur=snpCur,
                                          parallel = F,
                                          mc.cores = cores),
                       along = 1)
  }
  if (!parallel){
    out = abind::abind(lapply(1:k,
                              CVFit,
                              parts.train=parts.train,
                              parts.test=parts.test,
                              pheno=pheno,
                              snpCur=snpCur,
                              parallel = F),
                       along = 1)
  }

  out = as.data.frame(out)
  colnames(out) = c('id','enet','lasso','lmm')
  out = out[order(out$id),]

  r2.lmm = adjR2(pheno,out$lmm)
  r2.enet = adjR2(pheno,out$enet)
  r2.lasso = adjR2(pheno,out$lasso)

  if (max(r2.lmm,r2.enet,r2.lasso) <= 0.01){return('Mediator is not well-predicted')}

  if (r2.lmm >= max(r2.lasso,r2.enet)){

    tot.mod = rrBLUP::mixed.solve(y = pheno,
                                  Z = t(snpCur))
    mod.df = data.frame(SNP = c('Intercept',thisSNP$snpid),
                        Chromosome = c('-',thisSNP$chr),
                        Position = c('-',thisSNP$pos),
                        Effect = c(tot.mod$beta,tot.mod$u))
    return(list(Model = mod.df,
                Predicted = out$lmm,
                CVR2 = r2.lmm))

  }

  if (r2.enet >= max(r2.lasso,r2.lmm)){

    tot.mod = glmnet::cv.glmnet(y = pheno,
                                x = t(snpCur),
                                nfolds = 5,
                                alpha = 0.5,
                                parallel = T,
                                intercept = T)
    mod.df = data.frame(SNP = c('Intercept',thisSNP$snpid),
                        Chromosome = c('-',thisSNP$chr),
                        Position = c('-',thisSNP$pos),
                        Effect = as.numeric(coef(tot.mod,s='lambda.min')))
    return(list(Model = subset(mod.df,Effect != 0),
                Predicted = out$enet,
                CVR2 = r2.enet))

  }

  if (r2.lasso >= max(r2.enet,r2.lmm)){

    tot.mod = glmnet::cv.glmnet(y = pheno,
                                x = t(snpCur),
                                nfolds = 5,
                                alpha = 1,
                                parallel = T,
                                intercept = T)
    mod.df = data.frame(SNP = c('Intercept',thisSNP$snpid),
                        Chromosome = c('-',thisSNP$chr),
                        Position = c('-',thisSNP$pos),
                        Effect = as.numeric(coef(tot.mod,s='lambda.min')))
    return(list(Model = subset(mod.df,Effect != 0),
                Predicted = out$lasso,
                CVR2 = r2.lasso))

  }

}
