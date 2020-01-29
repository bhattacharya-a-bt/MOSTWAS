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
                     cisDist = 5e5,
                     qtlTra,
                     qtMed,
                     h2Pcutoff,
                     dimNumeric,
                     verbose,
                     seed,
                     nperms = 1000,
                     k,
                     parallel,
                     parType = 'no',
                     prune,
                     windowSize = 50,
                     numSNPShift = 5,
                     ldThresh = .5,
                     cores,
                     qtlTra_parts,
                     qtMed_parts,
                     modelDir){

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

  if (herit$P > h2Pcutoff) {return(paste0(geneInt,' is not heritable at P < ',h2Pcutoff))}

  print('MEDIATION ANALYSIS ON DISTAL SNPS')
  qtMed = subset(qtMed,SNP %in% tra.eSNP)
  if (nrow(qtMed) != 0){
    transSNPs = subset(snps,SNP %in% qtMed$SNP)
    if (parallel){
      medTest = parallel::mclapply(1:nrow(transSNPs),
                                   testTME,
                                   mediator = mediator,
                                   qtMed = qtMed,
                                   nperms = nperms,
                                   transSNPs = transSNPs,
                                   pheno = pheno,
                                   covariates = covariates,
                                   parallel = parType,
                                   cores = cores,
                                   mc.cores = ceiling(cores/2))
    }
    if (!parallel){
      medTest = pbapply::pblapply(1:nrow(transSNPs),
                                  testTME,
                                  mediator = mediator,
                                  nperms = nperms,
                                  qtMed = qtMed,
                                  transSNPs = transSNPs,
                                  pheno = pheno,
                                  covariates = covariates,
                                  cores = cores)
    }
    TME = sapply(medTest,function(x) x[[1]])
    TME.P = sapply(medTest,function(x) x[[2]])
    p_weights = IHW::ihw(TME.P ~ rowMeans(transSNPs[,-1])/2,
                         alpha = .05)
    include.trans = transSNPs$SNP[IHW::adj_pvalues(p_weights) < 0.05]
    snpCur = subset(snps,
                      SNP %in% c(cisGeno$snpList,include.trans))
    snpList = snpCur$SNP
    thisSNP = subset(snpLocs,
                     snpid %in% snpList)
    snpCur = as.matrix(snpCur[,-1])
    rownames(snpCur) = snpList
  }
    print('FITTING FULL MODEL')
    tot.mods = trainSNPPheno(pheno,
                             snpCur,
                             snpList,
                             thisSNP,
                             fileName = geneInt,
                             prune = prune)
    rm(snpCur)
    totSNP = subset(thisSNP,
                    snpid %in% rownames(coef(tot.mods$enet,s='lambda.min'))[-1])
    totSNP = totSNP[match(rownames(coef(tot.mods$enet,s='lambda.min'))[-1],
                          totSNP$snpid),]

    print('CROSS-VALIDATION ON DISTAL SNPS')
    set.seed(seed)
    train = caret::createFolds(1:length(pheno),
                               k = 3,
                               returnTrain = T)
    pred.enet = pred.blup = vector(mode = 'numeric',
                                   length = length(pheno))

    for (i in 1:k){

      qtlTraP = data.table::fread(qtlTra_parts[i])
      qtMedP = data.table::fread(qtMed_parts[i])
      tra.eSNP = qtlTra$SNP[qtlTraP$gene == geneInt]
      rm(qtlTraP)
      qtMedP = subset(qtMedP,SNP %in% tra.eSNP & FDR < 0.05)

      if (nrow(qtMedP) != 0){
        transSNPs = subset(snps,SNP %in% qtMedP$SNP)
        if (parallel){
          medTest = parallel::mclapply(1:nrow(transSNPs),
                                       testTME,
                                       mediator = mediator,
                                       qtMed = qtMedP,
                                       nperms = nperms,
                                       transSNPs = transSNPs,
                                       covariates = covariates,
                                       pheno = pheno,
                                       cores = cores,
                                       mc.cores = ceiling(cores/2))
        }
        if (!parallel){
          medTest = pbapply::pblapply(1:nrow(transSNPs),
                                      testTME,
                                      mediator = mediator,
                                      nperms = nperms,
                                      qtMed = qtMedP,
                                      pheno = pheno,
                                      transSNPs = transSNPs,
                                      covariates = covariates,
                                      cores = cores)
        }
        TME = sapply(medTest,function(x) x[[1]])
        TME.P = sapply(medTest,function(x) x[[2]])

        p_weights = IHW::ihw(TME.P ~ rowMeans(transSNPs[,-1])/2,
                             alpha = .10)
        include.trans = transSNPs$SNP[IHW::adj_pvalues(p_weights) < 0.10]
        snpCur = subset(snps,
                        SNP %in% c(cisGeno$snpList,include.trans))
        snpList = snpCur$SNP
        thisSNP = subset(snpLocs,
                         snpid %in% snpList)
        snpCur = as.matrix(snpCur[,-1])

        if (prune){
          pruneObj = LDprune(W = t(snpCur),
                             snpList = snpList,
                             snpLocs = thisSNP,
                             fileName = geneInt,
                             windowSize = windowSize,
                             numSNPShift = numSNPShift,
                             ldThresh = ldThresh)
          snpCur = t(pruneObj$W)
          snpList = pruneObj$snpList
          thisSNP = pruneObj$onlyThese
          rm(pruneObj)}

        thisMod = trainSNPPheno(pheno = pheno[train[[i]]],
                                snpCur[,train[[i]]],
                                snpList,
                                thisSnp = thisSNP,
                                fileName = geneInt,
                                prune = F)

        pred.enet[-train[[i]]] = as.numeric(predict(thisMod$enet,
                                                    newx = t(snpCur[,-train[[i]]]),
                                                    s = 'lambda.min'))
        pred.blup[-train[[i]]] = as.numeric(t(snpCur[,-train[[i]]]) %*% thisMod$blup$u)


      }
    }

    if (adjR2(pheno,pred.blup) < adjR2(pheno,pred.enet) | mean(coef(tot.mods$enet,s='lambda.min')[-1,] == 0) != 1){
      Model = data.frame(SNP = c(totSNP$snpid),
                         Chromosome = c(totSNP$chr),
                         Position = c(totSNP$pos),
                         Effect = as.numeric(coef(tot.mods$enet,
                                                  s='lambda.min'))[-1])
      Predicted = pred.enet

      } else {
        Model = data.frame(SNP = c(totSNP$snpid),
                           Chromosome = c(totSNP$chr),
                           Position = c(totSNP$pos),
                           Effect = as.numeric(tot.mods$blup$u))
        Predicted = pred.blup
        }


    Model = subset(Model,Effect != 0)
    R2 = max(adjR2(pheno,pred.blup),adjR2(pheno,pred.enet))
    h2 = herit$h2
    h2.Pvalue = herit$P
    save(Model,R2,Predicted,h2,h2.Pvalue,
       file = paste0(modelDir,geneInt,'.wgt.med.RData'))
}
