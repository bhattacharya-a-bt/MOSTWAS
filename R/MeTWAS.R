#' Train and predict gene's predictive model with mediators using MeTWAS
#'
#' The function trains a predictive model of a given gene using top mediators
#' as fixed effects and assesses in-sample performance with cross-validation.
#'
#' @param geneInt character, identifier for gene of interest
#' @param snpObj binsnp object, SNP dosages
#' @param mediator data frame, mediator intensities
#' @param medLocs data frame, MatrixEQTL locations for mediators
#' @param covariates data frame, covariates
#' @param dimNumeric numeric, number of numeric covariates
#' @param qtlFull data frame, all QTLs (cis and trans) between mediators and genes
#' @param h2Pcutoff numeric, P-value cutoff for heritability
#' @param numMed integer, number of top mediators to include
#' @param seed integer, random seed for splitting
#' @param k integer, number of training-test splits
#' @param parallel logical, TRUE/FALSE to run glmnet in parallel
#' @param prune logical, TRUE/FALSE to LD prune the genotypes
#' @param ldThresh numeric, LD threshold for PLINK pruning
#' @param cores integer, number of parallel cores
#' @param verbose logical, output everything
#' @param R2Cutoff numeric, cutoff for model R2
#' @param modelDir character, directory for saving models
#' @param tempFolder character, directory of saving snp backing files
#'
#' @return final model for gene along with CV R2 and predicted values
#'
#' @importFrom caret createFolds
#' @importFrom bigsnpr snp_attach
#' @importFrom bigsnpr snp_clumping
#' @importFrom bigsnpr snp_writeBed
#' @importFrom parallel mclapply
#' @importFrom abind abind
#'
#' @export
MeTWAS <- function(geneInt,
                   snpObj,
                   mediator,
                   medLocs,
                   covariates,
                   dimNumeric,
                   qtlFull,
                   h2Pcutoff = .1,
                   numMed = 5,
                   seed = 1218,
                   k = 5,
                   cisDist = 1e6,
                   parallel = T,
                   prune = F,
                   ldThresh = .5,
                   cores,
                   verbose = T,
                   R2Cutoff = 0.01,
                   modelDir,
                   tempFolder,
                   gctaFolder = NULL){

  set.seed(seed)
  if (!dir.exists(modelDir)){
    dir.create(modelDir)
  }

  if (is.null(gctaFolder)){
    gctaFolder = ''
  }

  print('GATHERING MEDIATORS')
  pheno = as.numeric(mediator[mediator$Mediator == geneInt,-1])
  pheno = (pheno - mean(pheno))/sd(pheno)
  medList = gatherMediators(geneInt,qtlFull,numMed)

  print('ESTIMATING HERITABILITY')
  lll = c(geneInt,medList)
  lll = lll[!is.na(lll)]
  w = c()
  for (i in 1:length(lll)){

    ml = subset(medLocs,geneid == lll[i])
    w = c(w,which(snpObj$map$chromosome == ml$chr[1] &
                    snpObj$map$physical.pos < ml$right[1] + cisDist &
                    snpObj$map$physical.pos > ml$left[1] - cisDist))

  }
  if (length(w) == 0){
    return('SNPs not found')
  }
  midSNPfile = subset(snpObj,ind.col = w)
  midSNP = bigsnpr::snp_attach(midSNPfile)
  midSNP$fam$affection = pheno

  if (prune){
    print('RUNNING LD CLUMPING')
    keep = bigsnpr::snp_clumping(
      midSNP$genotypes,
      infos.chr = midSNP$map$chromosome,
      ind.row = rows_along(midSNP$genotypes),
      S = NULL,
      thr.r2 = ldThresh,
      size = 100/ldThresh,
      exclude = NULL,
      ncores = ifelse(parallel,cores,1)
    )
    midSNP = bigsnpr::snp_attach(subset(midSNP,ind.col=keep))

  }

  tmpBed = paste0(tempFolder,"MeTWAS_",geneInt,".bed")
  bigsnpr::snp_writeBed(midSNP,tmpBed)

  system(paste(paste0(gctaFolder,'gcta64'),
               '--bfile',strsplit(tmpBed,'.bed')[[1]][1],
               '--autosome --make-grm',
               '--out',strsplit(tmpBed,'.bed')[[1]][1]),
         intern = !verbose)

  phenFile = paste0(strsplit(tmpBed,'.bed')[[1]][1],'.phen')
  covarFile = paste0(strsplit(tmpBed,'.bed')[[1]][1],'.qcovar')
  write.table(midSNP$fam[,c(1,2,6)],phenFile,row.names = F,
              col.names = F, quote = F)
  write.table(cbind(midSNP$fam[,1:2],
                    t(covariates[1:dimNumeric,-1])),
              covarFile,row.names = F,col.names = F,quote=F)
  system(paste(paste0(gctaFolder,'gcta64'),
               '--reml --reml-no-constrain --reml-bendV',
               '--grm',strsplit(tmpBed,'.bed')[[1]][1],
               '--pheno',phenFile,
               '--qcovar',covarFile,
               '--out',paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi')),
         intern = !verbose)

  if (file.exists(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'))){
    a = data.table::fread(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'),fill=T)
    herit = list(h2 = a$Variance[a$Source == 'V(G)/Vp'],
                 P = a$Variance[a$Source == 'Pval'])
  } else {
    system(paste(paste0(gctaFolder,'gcta64'),
                 '--reml',
                 '--grm',strsplit(tmpBed,'.bed')[[1]][1],
                 '--pheno',phenFile,
                 '--qcovar',covarFile,
                 '--out',paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi')),
           intern = !verbose)
    if (file.exists(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'))){
      a = data.table::fread(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'),fill=T)
      herit = list(h2 = a$Variance[a$Source == 'V(G)/Vp'],
                   P = a$Variance[a$Source == 'Pval'])
    } else {
    herit = list(h2 = 0,P=1)
  }}

  if (herit$P > h2Pcutoff) {
    return(paste(geneInt,
                 'is not germline heritable at P <',
                 h2Pcutoff))
  }


  print('TRAINING MEDIATORS')

  if (parallel) {
    medTrainList = parallel::mclapply(medList,
                                      trainLocalModel,
                                      midSNP = midSNP,
                                      mediator = mediator,
                                      medLocs = medLocs,
                                      covariates = covariates,
                                      cisDist = cisDist,
                                      seed = seed,
                                      nfolds = k,
                                      mc.cores = cores)
  }
  if (!parallel){
    medTrainList = lapply(medList,
                          trainLocalModel,
                          midSNP = midSNP,
                          mediator = mediator,
                          medLocs = medLocs,
                          covariates = covariates,
                          cisDist = cisDist,
                          seed = seed,
                          nfolds = k)
  }
  names(medTrainList) = medList
  medTraitList = medTrainList[sapply(medTrainList,length) > 0]

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
      pheno = pheno - as.numeric(predict(lmCVFit))
    } else {
      pheno = pheno
      fixedEffects = NULL}
  }


  print('FITTING LOCAL-GENOTYPES')

  cisGenoMod = trainLocalModel(phenoInt = geneInt,
                               midSNP = midSNP,
                               mediator = mediator,
                               medLocs = medLocs,
                               covariates = covariates,
                               seed = seed,
                               nfolds = k,
                               cisDist = cisDist)
  if (is.null(cisGenoMod$Model)){
    cisGenoMod$Model = as.data.frame(matrix(nrow = 1,ncol = 4))
    colnames(cisGenoMod$Model) = c('SNP',
                                   'Chromosome',
                                   'Position',
                                   'Effect')
    cisGenoMod$Model[1,] = 0
  }


  cisGenoMod$Model$Mediator = 'Cis'
  if (exists('trans.mod.df')){
    cisGenoMod$Model = rbind(cisGenoMod$Model,trans.mod.df)
  }
  cisGenoMod$Model = subset(cisGenoMod$Model,Effect!=0)
  cisGenoMod$CVR2 = cisGenoMod$CVR2 + fe.R2
  cisGenoMod$CVR2.cis = cisGenoMod$CVR2 - fe.R2
  cisGenoMod$h2 = herit$h2
  cisGenoMod$h2.P = herit$P



  Model = cisGenoMod$Model
  R2 = cisGenoMod$CVR2
  Predicted = cisGenoMod$Predicted
  Mediators = cisGenoMod$medlist
  CisR2 = cisGenoMod$CVR2.cis
  h2 = abs(herit$h2)
  h2.Pvalue = herit$P
  if (R2 < R2Cutoff){ return('CV R2 < 0.01.') }
  ## REMOVE THE NEXT LINE
  CorMat = cbind(Predicted,fixedEffects)
  save(Model,R2,Predicted,Mediators,CisR2,h2,h2.Pvalue,CorMat,
       file = paste0(modelDir,geneInt,'.wgt.med.RData'))
  rm(midSNP)
  file.remove(midSNPfile)

}
