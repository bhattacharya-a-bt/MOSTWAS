#' Train a local-only model for a phenotype
#'
#' The function runs FUSION-type local-only modeling for a phenotype
#'
#' @param pheno vector, phenotypes
#' @param snpCur matrix, SNP dosages
#' @param snpList vector, vector of SNP identifiers
#' @param thisSNP data.frame, SNP locations
#' @param fileName character, prefix for LDprune
#' @param prune logical, T/F to LD prune
#' @param windowSize numeric, PLINK LD prune window size
#' @param numSNPShift integer, PLINK LD shifting interval
#' @param ldThresh numeric, LD threshold
#' @param verbose logical, T/F for verbosity
#' @param snpAnnot data.frame, SNP annotations
#'
#' @return list with TME and P-value
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom rrBLUP mixed.solve
#'
#' @export
trainSNPPheno <- function(pheno,
                          snpCur,
                          snpList,
                          thisSNP,
                          fileName,
                          prune = T,
                          windowSize = 50,
                          numSNPShift = 5,
                          ldThresh = .5,
                          verbose = T,
                          snpAnnot = NULL){
  if (prune){
    if (length(snpList) == nrow(snpCur)){
      snpCur = t(snpCur)
    }
    pruneObj = LDprune(W = t(snpCur),
                       snpList = snpList,
                       snpLocs = thisSNP,
                       fileName = fileName,
                       windowSize = windowSize,
                       numSNPShift = numSNPShift,
                       ldThresh = ldThresh,
                       verbose = verbose,
                       snpAnnot = snpAnnot)
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
