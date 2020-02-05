#' Estimate germline heritability of a biomarker
#'
#' The function invokves GCTA in the system and estimates the germline
#' heritability of a biomarker given cis-regions around the biomarker
#' and any associated mediators of interest.
#'
#' @param biomInt character, identifier for biomarker of interest
#' @param snps data frame, SNP dosages
#' @param snpLocs data frame, MatrixEQTL locations for SNPs
#' @param mediator data frame, mediator intensities
#' @param medLocs data frame, MatrixEQTL locations for mediators
#' @param covariates data frame, covariates
#' @param numMed integer, number of top mediators to include
#' @param fileName character, throw away name for PLINK files
#' @param dimNumeric integer, number of continuous covariates
#' @param numMed integer, maximum number of mediators to find
#' @param cisDist numeric, base pair window of cis definition
#' @param needMed logical, T/F if the biomarker of interest has associated mediators
#' @param medList character vector,  vector of mediator names
#' @param verbose logical, T/F for detailed output
#' @param windowSize integer, window size for PLINK pruning
#' @param numSNPShift integer, shifting window from PLINK pruning
#' @param ldThresh numeric, maximum threshold for linkage disequilibrium
#' @param ldScrRegion numeric, number of SNPs per bin in GCTA-LDMS
#' @param LDMS logical, T/F to run GCTA-LDMS
#' @param supplySNP logical, T/F to run GCTA on the full provided snps
#'
#'
#' @return heritability estimate and associated LRT P-value
#'
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom abind abind
#'
#' @export
estimateHeritability <- function(biomInt,
                                   snps,
                                   pheno,
                                   snpLocs,
                                   mediator,
                                   medLocs,
                                   covariates,
                                   fileName,
                                   dimNumeric,
                                   numMed = 5,
                                   cisDist = 5e5,
                                   needMed = T,
                                   medList,
                                   verbose = F,
                                   windowSize,
                                   numSNPShift,
                                   ldThresh,
                                   ldScrRegion = 200,
                                   LDMS = F,
                                   supplySNP = F){

    if (!supplySNP){
    if (needMed) {
      cisGeno = lapply(c(biomInt,medList),
                       getCisGenotypes,
                       locs = medLocs,
                       snps = snps,
                       snpLocs = snpLocs,
                       cisDist = cisDist)
      W = abind::abind(lapply(cisGeno,function(x) x[[1]]),
                       along=1)
      snpList = unlist(lapply(cisGeno, function(x) x[[2]]))
      thisSNP = as.data.frame(abind::abind(lapply(cisGeno,function(x) x[[3]]),
                                           along=1))
    }

    if (!needMed) {
      cisGeno = getCisGenotypes(biomInt,
                               medLocs,
                               snps,
                               snpLocs,
                               cisDist)
      W = cisGeno$snpCur
      snpList = cisGeno$snpList
      thisSNP = cisGeno$thisSNP
    }}

    if (supplySNP){
      snps = snps[!duplicated(snps$SNP),]
      snpList = snps$SNP
      W = subset(snps,SNP %in% snpList)[,-1]
      thisSNP = subset(snpLocs,snpid %in% snps$SNP)
    }

    fileName = paste0('h2_',biomInt)
    tempDir = paste0(biomInt,'_h2_temp/')
    if (!dir.exists(tempDir)){ dir.create(tempDir) }
    genfile = paste0(tempDir,fileName,'.gen')
    samplefile = paste0(tempDir,fileName,'.sample')
    bedfile = paste0(tempDir,fileName)
    outFile = paste0(tempDir,fileName)
    phenFile = paste0(tempDir,fileName,'.phen')
    covarFile = paste0(tempDir,fileName,'.qcovar')

    df = as.data.frame(matrix(nrow = 1,ncol =2))
    ids = colnames(W)
    geno = as.data.frame(matrix(ncol=ncol(W)+4,nrow = nrow(W)))
    colnames(geno) <- c('SNP','Pos','A1','A2',ids)
    geno[,5:ncol(geno)] = W
    W = W[match(snpList,rownames(W)),]
    geno$SNP = snpList
    onlyThese <- thisSNP[thisSNP$snpid %in% geno$SNP,]
    geno <- geno[geno$SNP %in% snpLocs$snpid,]
    onlyThese <- onlyThese[match(onlyThese$snpid,geno$SNP),]
    chr <- onlyThese$chr
    geno$Pos <- onlyThese$pos
    geno$A1 <- unlist(lapply(strsplit(geno$SNP,':'),function(x) as.character(x[3])))
    geno$A2 <- unlist(lapply(strsplit(geno$SNP,':'),function(x) as.character(x[4])))
    chr_dosage <- cbind(chr,geno)
    rm(geno)

    new.levels <- c('1 0 0','0 1 0','0 0 1')
    matrix.alleles <- as.matrix(chr_dosage[,6:ncol(chr_dosage)] + 1)
    impute2.format <- matrix(new.levels[matrix.alleles],ncol=ncol(matrix.alleles))
    gen <- cbind(chr_dosage[,1:5],impute2.format)
    gen[is.na(gen)] <- '<NA>'

    require(data.table)
    data.table::fwrite(gen,genfile,row.names=FALSE,
                       col.names = FALSE, quote = FALSE, sep='\t')

    sample <- as.data.frame(matrix(nrow=(ncol(chr_dosage)-4),ncol=5))
    colnames(sample) <- c('ID_1','ID_2','missing','gender','pheno')
    sample$ID_1[2:nrow(sample)] <- paste('1A',2:nrow(sample),sep='')
    sample$ID_2[2:nrow(sample)] <- colnames(chr_dosage)[6:ncol(chr_dosage)]
    sample$missing[2:nrow(sample)] <- 0
    sample$gender[2:nrow(sample)] <- 2
    sample$pheno[2:nrow(sample)] <- pheno
    sample[1,] <- c(0,0,0,'D','P')
    write.table(sample,samplefile,row.names=FALSE,
                col.names = TRUE, quote = FALSE)
    write.table(sample[-1,c(1,2,5)],phenFile,row.names = F,
                col.names = F, quote = F)
    write.table(cbind(sample[-1,1:2],t(covariates[1:dimNumeric,-1])),
                covarFile,row.names = F,col.names = F,quote=F)

    if (!verbose){
      sink('/dev/null')
    }

  system(paste('plink',
               '--gen',genfile,
               '--sample',samplefile,
               '--make-bed',
               '--allow-no-sex',
               '--out',bedfile),
         intern = !verbose)

  a = data.table::fread(paste0(bedfile,'.fam'))
  a$V5 = 2
  data.table::fwrite(a,paste0(bedfile,'.fam'),col.names=F,row.names=F,quote=F,sep='\t')

  system(paste('plink','--bfile',bedfile,
               '--indep-pairwise',windowSize,numSNPShift,ldThresh,
               '--out',bedfile),
         intern = !verbose)

  system(paste('plink','--bfile',bedfile,
               '--extract',paste0(bedfile,'.prune.in'),
               '--make-bed',
               '--out',paste0(bedfile,'_prune')),
         intern = !verbose)

  bedfile = paste0(bedfile,'_prune')

  if (LDMS){

  system(paste('gcta64',
               '--bfile',bedfile,
               '--ld-score-region',ldScrRegion,
               '--out',bedfile),
         intern = !verbose)


  lds_seg = read.table(paste0(bedfile,".score.ld"),
                              header=T,
                       colClasses=c("character",rep("numeric",8)))
  quartiles=summary(lds_seg$ldscore_SNP)

  lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
  lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
  lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
  lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

  lb1_snp = lds_seg$SNP[lb1]
  lb2_snp = lds_seg$SNP[lb2]
  lb3_snp = lds_seg$SNP[lb3]
  lb4_snp = lds_seg$SNP[lb4]

  write.table(lb1_snp,
              paste0(bedfile,"snp_group1.txt"),
              row.names=F,
              quote=F,
              col.names=F)
  write.table(lb2_snp,
              paste0(bedfile,"snp_group2.txt"),
              row.names=F,
              quote=F,
              col.names=F)
  write.table(lb3_snp,
              paste0(bedfile,"snp_group3.txt"),
              row.names=F,
              quote=F,
              col.names=F)
  write.table(lb4_snp,
              paste0(bedfile,"snp_group4.txt"),
              row.names=F,
              quote=F,
              col.names=F)

  system(paste('gcta64',
               '--bfile',bedfile,
               '--extract',paste0(bedfile,"snp_group1.txt"),
               '--make-grm',
               '--out',paste0(bedfile,'_group1')),
         intern = !verbose)
  system(paste('gcta64',
               '--bfile',bedfile,
               '--extract',paste0(bedfile,"snp_group2.txt"),
               '--make-grm',
               '--out',paste0(bedfile,'_group2')),
         intern = !verbose)
  system(paste('gcta64',
               '--bfile',bedfile,
               '--extract',paste0(bedfile,"snp_group3.txt"),
               '--make-grm',
               '--out',paste0(bedfile,'_group3')),
         intern = !verbose)
  system(paste('gcta64',
               '--bfile',bedfile,
               '--extract',paste0(bedfile,"snp_group4.txt"),
               '--make-grm',
               '--out',paste0(bedfile,'_group4')),
         inter = !verbose)

  write.table(paste0(bedfile,'_group',1:4),
              paste0(bedfile,'_multiGRM.txt'),
              col.names = F,
              row.names = F,
              quote=F)

  system(paste('gcta64',
               '--reml',
               '--mgrm',paste0(bedfile,'_multiGRM.txt'),
               '--pheno',phenFile,
               '--qcovar',covarFile,
               '--reml-no-constrain',
               '--out',paste0(bedfile,'_multi')),
         intern = !verbose)

  if (file.exists(paste0(bedfile,'_multi.hsq'))) {
    hsq = data.table::fread(paste0(bedfile,'_multi.hsq'),
                            fill = T)
    h2 = hsq$Variance[10]
    P = hsq$Variance[17]
  }

  if (!file.exists(paste0(bedfile,'_multi.hsq'))){
    h2 = 0
    P = 1
  }

  system(paste('rm -r',
                tempDir))
  }


  if (!LDMS) {


    system(paste('gcta64',
                 '--bfile',bedfile,
                 '--autosome --make-grm',
                 '--out',bedfile),
           intern = !verbose)


    system(paste('gcta64',
                 '--reml',
                 '--grm',bedfile,
                 '--pheno',phenFile,
                 '--qcovar',covarFile,
                 '--reml-no-constrain',
                 '--out',paste0(bedfile,'_multi')),
           intern = !verbose)

    if (file.exists(paste0(bedfile,'_multi.hsq'))){
    hsq = data.table::fread(paste0(bedfile,'_multi.hsq'),
                            fill = T)
    h2 = hsq$Variance[4]
    P = hsq$Variance[9]
    }

    if (!file.exists(paste0(bedfile,'_multi.hsq'))){
      h2 = 0
      P = 1
    }


    system(paste('rm -r',
                 tempDir))

  }


    if (!verbose){
      sink()
    }

  return(list(h2 = h2, P = P))


}
