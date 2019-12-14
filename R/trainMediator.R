#' Train and predict mediator predictive model
#' 
#' The function trains a predictive model of a given mediator and then
#' predicts the genetical regulated intensity of it in the training set via
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
                          parallel = T,
                          prune = T,
                          windowSize = 50,
                          numSNPShift = 5,
                          ldThresh = .5){
  
  set.seed(seed)
  parts.train = caret::createFolds(colnames(mediator)[-1],k=k,returnTrain = T)
  set.seed(seed)
  parts.test = caret::createFolds(colnames(mediator)[-1],k=k,returnTrain = F)
  
  pheno = as.numeric(mediator[grepl(medInt,mediator$Mediator),-1])
  
  res = as.data.frame(cbind(pheno,t(covariates[,-1])))
  pheno = as.numeric(resid(lm(pheno~.,data=res)))
  
  medLocs = subset(medLocs, grepl(medInt,snpid))
  
  thisSNP = subset(snpLocs, chr == paste0('chr',medLocs$chr[1]) & pos < medLocs$pos+1e6 & pos > medLocs$pos - 1e6)
  snpCur = subset(snps,SNP %in% thisSNP$snpid)
  snpList = snpCur$SNP
  snpCur = as.matrix(snpCur[,-1])
  rownames(snpCur) = snpList
  
  if (parallel){
    doParallel::registerDoParallel(5)
  }
  
  pheno.pred.enet = pheno.pred.lasso = pheno.pred.lmm = vector(length = length(pheno))
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
  
  for (i in 1:k){
    
    training = parts.train[[i]]
    test = parts.test[[i]]
    y.train = pheno[training]
    
    W.train = snpCur[,training]
    W.test = snpCur[,test]
    
    mod.enet = glmnet::cv.glmnet(y = y.train,
                                 x = t(W.train),
                                 nfolds = 5,
                                 alpha = 0.5,
                                 parallel = T)
    pheno.pred.enet[test] = predict(mod.enet,
                                    newx = t(W.test),
                                    s = 'lambda.min')
    
    mod.lasso = glmnet::cv.glmnet(y = y.train,
                                  x = t(W.train),
                                  nfold = 5,
                                  alpha = 1,
                                  parallel = T)
    pheno.pred.lasso[test] = predict(mod.lasso,
                                     newx = t(W.test),
                                     s = 'lambda.min')
    
    mod.lmm = rrBLUP::mixed.solve(y = y.train,
                                  Z = t(W.train))
    pheno.pred.lmm[test] = rep(mod.lmm$beta,ncol(W.test)) + 
      as.numeric(mod.lmm$u %*% W.test)
    
  }
  r2.lmm = adjR2(pheno,pheno.pred.lmm)
  r2.enet = adjR2(pheno,pheno.pred.enet)
  r2.lasso = adjR2(pheno,pheno.pred.lasso)
  
  if (max(r2.lmm,r2.enet,r2.lasso) <= 0.01){return('Mediator is not well-predicted')}
  
  if (r2.lmm >= max(r2.lasso,r2.enet)){
    
    tot.mod = rrBLUP::mixed.solve(y = pheno,
                                  Z = t(snpCur))
    mod.df = data.frame(SNP = c('Intercept',thisSNP$snpid),
                        Chromosome = c('-',thisSNP$chr),
                        Position = c('-',thisSNP$pos),
                        Effect = c(tot.mod$beta,tot.mod$u))
    return(list(Model = mod.df,
                Predicted = pheno.pred.lmm,
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
                Predicted = pheno.pred.enet,
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
                Predicted = pheno.pred.lasso,
                CVR2 = r2.lasso))
    
  }
  
}