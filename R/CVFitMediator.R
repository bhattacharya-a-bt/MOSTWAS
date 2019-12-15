#' Fit a single split in cross-validation for a given biomarker
#'
#' The function trains a predictive model on the given training split and then
#' predicts the genetically regulated intensity of it in the test split.
#'
#' @param foldNum integer, which fold of the k splits
#' @param parts.train list, list of training indices
#' @param parts.test list, list of test indices
#' @param pheno vector, vector of standardized intensities of biomarker
#' @param snpCur matrix, matrix of cis-genotypes to biomarker
#' @param parallel logical, run in parallel?
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom rrBLUP mixed.solve
#'
#' @export

CVFitMediator <- function(foldNum,
                  parts.train,
                  parts.test,
                  pheno,
                  snpCur,
                  parallel = T,
                  fixedEffects = NULL){

  training = parts.train[[foldNum]]
  test = parts.test[[foldNum]]
  y.train = pheno[training]

  W.train = snpCur[,training]
  W.test = snpCur[,test]

  mod.enet = glmnet::cv.glmnet(y = y.train,
                               x = t(W.train),
                               nfolds = 5,
                               alpha = 0.5,
                               parallel = parallel)

  mod.lasso = glmnet::cv.glmnet(y = y.train,
                                x = t(W.train),
                                nfold = 5,
                                alpha = 1,
                                parallel = parallel)

  mod.lmm = rrBLUP::mixed.solve(y = y.train,
                                Z = t(W.train))
  df = data.frame(id = test,
                  enet = predict(mod.enet,
                                 newx = t(W.test),
                                 s = 'lambda.min'),
                  lasso = predict(mod.lasso,
                                  newx = t(W.test),
                                  s = 'lambda.min'),
                  lmm = rep(mod.lmm$beta,ncol(W.test)) +
                    as.numeric(mod.lmm$u %*% W.test))

  return(df)

}
