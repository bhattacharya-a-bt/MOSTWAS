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
trainDeP <- function(geneInt,
                     snps,
                     snpLocs,
                     mediator,
                     medLocs,
                     covariates,
                     qtlTra,
                     qtMed,
                     h2Pcutoff,
                     seed,
                     k,
                     parallel,
                     prune,
                     windowSize = 50,
                     numSNPShift = 5,
                     ldThresh = .5,
                     cores,
                     outputAll){

  set.seed(seed)
  pheno = as.numeric(mediator[mediator$Mediator == geneInt,-1])

  print('GATHERING GENOTYPES')
  cisGeno = getCisGenotypes(biomInt = geneInt,
                            locs = medLocs,
                            snps = snps,
                            snpLocs = snpLocs,
                            cisDist = cisDist)
  tra.eSNP = qtlTra$SNP[qtlTra$gene == geneInt]
  snpsThis = subset(snps,
                    SNP %in% c(cisGeno$snpList,tra.eSNP))

  print('ESTIMATING HERITABILITY')
  herit = estimateHeritability(biomInt = geneInt,
                               snps = snpsThis,
                               pheno = pheno,
                               snpLocs = snpLocs,
                               mediator = mediator,
                               medLocs = medLocs,
                               covariates = covariates,
                               fileName = geneInt,
                               dimNumeric = dimNumeric,
                               numMed = 5,
                               cisDist = 5e5,
                               needMed = F,
                               medList = medList,
                               verbose = verbose,
                               windowSize = windowSize,
                               numSNPShift = numSNPShift,
                               ldThresh = ldThresh,
                               ldScrRegion = ldScrRegion,
                               LDMS = F,
                               supplySNP = T)

  if (herit$P < h2Pcutoff) {return(paste0(geneInt,' is not heritable at P < ',h2Pcutoff))}

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
                             ldThresh = ldThresh)

  cisR2 = cisGenoMod$CVR2
  pheno = pheno - cisGenoMod$Predicted

  print('WEIGHTS FOR DISTAL-GENOTYPES')
  qtMed = subset(qtMed,SNP %in% tra.eSNP)
  if (nrow(qtMed) != 0){
    transSNPs = subset(snps,SNP %in% qtMed$SNP)
    TME = TME.P = vector(length = nrow(transSNPs),
                         mode = 'numeric')
    for (i in 1:length(TME)){

      thisMed = subset(mediator,Mediator %in% qtMed$gene[qtMed$SNP == transSNPs$S[i]])
      if (nrow(thisMed) == 0){
        TME[i] = 0
        TME.P[i] = 1
      }
      if (nrow(thisMed) > 0){
        test = permuteTME(snp = as.numeric(as.vector(transSNPs[i,-1])),
                          expression = pheno,
                          mediators = t(as.matrix(thisMed[,-1])),
                          covs = t(as.matrix(covariates[,-1])),
                          nperms = 1000,
                          nc = 1)
        TME[i] = test$test.stat
        TME.P[i] = test$p.value
        }

    }
  }

  pen = IHW::ihw(TME.P ~ rowMeans(transSNPs[,-1])/2,alpha = .1)
  set.seed(seed)
  train = caret::createFolds(y = pheno,
                             k=k,
                             returnTrain = T)
  set.seed(seed)
  test = caret::createFolds(y = pheno,
                            k = k,
                            returnTrain = F)

  pred.wenet = pred.enet = vector(mode = 'numeric',length = length(pheno))
  transSNPMat = as.matrix(transSNPs[,-1])
  rownames(transSNPMat) = transSNPs$SNP
  wenet = glmnet::cv.glmnet(y = pheno[train[[i]]],
                            x = t(transSNPMat[,train[[i]]]),
                            penalty.factor = TME,
                            nfolds = 10,
                            type.measure = 'deviance')
  transLocs = subset(snpLocs,snpid %in% transSNPs$SNP)
  transLocs = transLocs[match(rownames(transSNPMat),transLocs$snpid),]
  transMod = data.frame(SNP = transLocs$snpid,
                        Chromosome = transLocs$chr,
                        Position = transLocs$pos,
                        Effect = as.vector(coef(wenet,s = 'lambda.min'))[-1])

}
