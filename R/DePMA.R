#' Train and predict gene's predictive model with mediators using DePMA
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
#' @param sobel logical, Sobel asymptotic test T/F
#' @param nperms numeric, number of permutations
#' @param parType character, parallelization type for boots
#' @param qtlTra_parts character vector, files for SNP to gene distal-eqtls
#' @param qtMed_parts character vector, files for SNP to mediators local-eqtls
#'
#' @return final model for gene along with CV R2 and predicted values
#'
#' @importFrom caret createFolds
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom parallel mclapply
#'
#' @export
DePMA <- function(geneInt,
                  snpObj,
                  mediator,
                  medLocs,
                  covariates,
                  cisDist = 1e6,
                  qtlTra,
                  qtMed,
                  h2Pcutoff,
                  dimNumeric,
                  verbose,
                  seed,
                  sobel = F,
                  nperms = 1000,
                  k,
                  parallel,
                  parType = 'no',
                  prune,
                  ldThresh = .5,
                  cores,
                  qtlTra_parts,
                  qtMed_parts,
                  modelDir,
                  tempFolder,
                  R2Cutoff){


  if (!dir.exists(modelDir)){
    dir.create(modelDir)
  }

  set.seed(seed)
  colnames(mediator)[1] = 'Mediator'
  pheno = as.numeric(mediator[mediator$Mediator == geneInt,-1])
  if (var(pheno) == 0){
    return('Gene expression has variance 0')
  }
  pheno = (pheno - mean(pheno))/sd(pheno)

  print('GATHERING GENOTYPES')
  ml = subset(medLocs, geneid == geneInt)
  w = which(snpObj$map$chromosome == ml$chr[1] &
                  snpObj$map$physical.pos < ml$right[1] + cisDist &
                  snpObj$map$physical.pos > ml$left[1] - cisDist)
  tra.eSNP = qtlTra$SNP[qtlTra$gene == geneInt]
  tra.w = which(snpObj$map$marker.ID %in% tra.eSNP)
  tot.w = c(w,tra.w)
  if (length(tot.w) == 0) {return('SNPs not found')}
  midSNPfile = subset(snpObj,ind.col = tot.w)
  midSNP = bigsnpr::snp_attach(midSNPfile)

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

  print('ESTIMATING HERITABILITY')
  midSNP$fam$affection = pheno
  tmpBed = tempfile(fileext = ".bed")
  bigsnpr::snp_writeBed(midSNP,tmpBed)

  system(paste('gcta64',
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
  system(paste('gcta64',
               '--reml',
               '--grm',strsplit(tmpBed,'.bed')[[1]][1],
               '--pheno',phenFile,
               '--qcovar',covarFile,
               '--reml-no-constrain --reml-bendV',
               '--out',paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi')),
         intern = !verbose)

  if (file.exists(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'))){
  a = fread(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'),fill=T)
  herit = list(h2 = a$Variance[4],
               P = a$Variance[9])} else {
                 herit = list(h2 = 0,P = 1)
               }

  if (herit$P > h2Pcutoff) {
    return(paste0(geneInt,' is not heritable at P < ',h2Pcutoff))
  }


  print('MEDIATION ANALYSIS ON DISTAL SNPS')
  qtMed = subset(qtMed,SNP %in% tra.eSNP)
  if (nrow(qtMed) != 0){
    snps = as.data.frame(cbind(midSNP$map$marker.ID,
                 t(midSNP$genotypes[])))
    colnames(snps) = c('SNP',midSNP$fam$family.ID)
    transSNPs = subset(snps,SNP %in% qtMed$SNP)
    if (nrow(transSNPs) != 0){
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
                                     mc.cores = ceiling(cores/2),
                                     sobel = sobel)
      }
      if (!parallel){
        medTest = lapply(1:nrow(transSNPs),
                         testTME,
                         mediator = mediator,
                         nperms = nperms,
                         qtMed = qtMed,
                         transSNPs = transSNPs,
                         pheno = pheno,
                         covariates = covariates,
                         cores = cores,
                         sobel = sobel)
      }
      TME = sapply(medTest,function(x) x[[1]])
      TME.P = sapply(medTest,function(x) x[[2]])
      p_weights = p.adjust(TME.P,'BH')
      include.trans = transSNPs$SNP[p_weights < 0.10]
    } else {include.trans = NULL}
  } else {include.trans = NULL}
  w = c(which(midSNP$map$chromosome == ml$chr[1] &
              midSNP$map$physical.pos < ml$right[1] + cisDist &
              midSNP$map$physical.pos > ml$left[1] - cisDist),
        which(midSNP$map$marker.ID %in% include.trans))
  midSNPfile = subset(snpObj,ind.col = w)
  midSNP = bigsnpr::snp_attach(midSNPfile)

  print('FITTING FULL MODEL')
  geno.var = which(apply(as.matrix(midSNP$genotypes[]),2,var) != 0)
  if (length(geno.var) == 0){return('Final model has no SNPs')}
  midSNP = bigsnpr::snp_attach(subset(midSNP,ind.col = geno.var))
  fin.model.enet = bigstatsr::big_spLinReg(midSNP$genotypes,
                                           pheno,
                                           alphas = .5,
                                           warn = FALSE)
  mod.df.enet = data.frame(SNP =midSNP$map$marker.ID[attr(fin.model.enet,'ind.col')],
                           Chromosome = midSNP$map$chromosome[attr(fin.model.enet,'ind.col')],
                           Position = midSNP$map$physical.pos[attr(fin.model.enet,'ind.col')],
                           Effect = summary(fin.model.enet)$beta)
  colnames(mod.df.enet) = c('SNP','Chromosome','Position','Effect')
  mod.df.enet = subset(mod.df.enet,Effect != 0)


  fin.model.blup = rrBLUP::mixed.solve(y = pheno,
                                       Z = midSNP$genotypes[])
  mod.df.blup = data.frame(SNP = midSNP$map$marker.ID,
                           Chromosome = midSNP$map$chromosome,
                           Position = midSNP$map$physical.pos,
                           Effect = as.numeric(fin.model.blup$u))


  print('CROSS-VALIDATION ON DISTAL SNPS')
  set.seed(seed)
  train = caret::createFolds(1:length(pheno),
                             k = length(qtlTra_parts),
                             returnTrain = T)
  set.seed(seed)
  test = caret::createFolds(1:length(pheno),
                             k = length(qtlTra_parts),
                             returnTrain = F)
  pred.enet = pred.blup = vector(mode = 'numeric',
                                 length = length(pheno))
  k = length(qtlTra_parts)
  for (i in 1:k){
    if (var(pheno[train[[i]]]) == 0){
      pred.enet[test[[i]]] = pred.blup[test[[i]]] = 0
    } else {
    qtlTraP = data.table::fread(qtlTra_parts[i])
    qtMedP = data.table::fread(qtMed_parts[i])
    tra.eSNP = qtlTra$SNP[qtlTraP$gene == geneInt]
    rm(qtlTraP)
    qtMedP = subset(qtMedP,SNP %in% tra.eSNP & FDR < 0.05)

    if (nrow(qtMedP) != 0){
      snps = as.data.frame(cbind(midSNP$map$marker.ID,
                                 t(midSNP$genotypes[])))
      colnames(snps) = c('SNP',midSNP$fam$family.ID)
      transSNPs = subset(snps,SNP %in% qtMed$SNP)
      if (nrow(transSNPs) != 0){
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
                                       sobel = sobel,
                                       mc.cores = ceiling(cores/2))
        }
        if (!parallel){
          medTest = lapply(1:nrow(transSNPs),
                           testTME,
                           mediator = mediator,
                           nperms = nperms,
                           qtMed = qtMedP,
                           pheno = pheno,
                           transSNPs = transSNPs,
                           covariates = covariates,
                           cores = cores,
                           sobel = sobel)
        }
        TME = sapply(medTest,function(x) x[[1]])
        TME.P = sapply(medTest,function(x) x[[2]])
        p_weights = p.adjust(TME.P,'BH')
        include.trans = transSNPs$SNP[p_weights < 0.10]
      } else {include.trans = NULL} }
    w = which(snpObj$map$chromosome == ml$chr[1] &
                snpObj$map$physical.pos < ml$right[1] + cisDist &
                snpObj$map$physical.pos > ml$left[1] - cisDist)
    current.w = c(w,which(snpObj$map$marker.ID %in% include.trans))
    if (length(current.w) == 0){
      pred.enet[test[[i]]] = pred.blup[test[[i]]] = 0
    } else {
    snpCur = bigsnpr::snp_attach(subset(snpObj,ind.col = current.w))

    if (prune){

      keep = bigsnpr::snp_clumping(
        snpCur$genotypes,
        infos.chr = snpCur$map$chromosome,
        ind.row = rows_along(snpCur$genotypes),
        S = NULL,
        thr.r2 = ldThresh,
        size = 100/ldThresh,
        exclude = NULL,
        ncores = ifelse(parallel,cores,1)
      )
      snpCur = bigsnpr::snp_attach(subset(snpCur,ind.col=keep))

    }

    geno.var = which(apply(as.matrix(midSNP$genotypes[]),2,var) != 0)
    if (length(geno.var) == 0){
      pred.enet[test[[i]]] = pred.blup[test[[i]]] = 0
    } else {
    mod.enet <- bigstatsr::big_spLinReg(snpCur$genotypes,
                                        pheno[train[[i]]],
                                        ind.train = train[[i]],
                                        alphas = .5,K = 5, warn = FALSE)
    pred.enet[test[[i]]] = predict(mod.enet,
                                   snpCur$genotypes,
                                   ind.row = test[[i]])

    mod.blup = rrBLUP::mixed.solve(y = pheno[train[[i]]],
                                   Z = snpCur$genotypes[train[[i]],])
    pred.blup[test[[i]]] = as.numeric(snpCur$genotypes[test[[i]],] %*%
                                        mod.blup$u)
    }}}
  }

  r2.blup = adjR2(pheno,pred.blup)
  r2.enet = adjR2(pheno,pred.enet)
  if (max(c(r2.blup,r2.enet)) < R2Cutoff){
    return(paste0('CV R2 < ',R2Cutoff))
  }

  if (r2.blup <= r2.enet & nrow(mod.df.enet) != 0){
    Model = mod.df.enet
    R2 = r2.enet
    h2 = herit$h2
    h2.Pvalue = herit$P
    Predicted = pred.enet
    save(Model,R2,Predicted,h2,h2.Pvalue,
         file = paste0(modelDir,geneInt,'.wgt.med.RData'))

  } else {
    Model = mod.df.blup
    R2 = r2.blup
    h2 = herit$h2
    h2.Pvalue = herit$P
    Predicted = pred.blup
    save(Model,R2,Predicted,h2,h2.Pvalue,
         file = paste0(modelDir,geneInt,'.wgt.med.RData'))
  }


  if (dir.exists('temp')){
    fff = list.files('temp/')
    cleanup = geneInt
    file.remove(paste0('temp/',fff[grepl(geneInt,fff)]))
  }
  rm(midSNP)
  file.remove(midSNPfile)

}

