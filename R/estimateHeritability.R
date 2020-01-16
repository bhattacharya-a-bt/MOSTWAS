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
#' @param cisDist numeric, base pair window of cis definition
#' @param needMed logical, T/F if the biomarker of interest has associated mediators
#' @param medList character vector,  vector of mediator names
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
                                 geneInt,
                                 snps,
                                 snpLocs,
                                 mediator,
                                 medLocs,
                                 covariates,
                                 fileName,
                                 dimNumeric,
                                 numMed = 5,
                                 cisDist = 5e5,
                                 needMed = F,
                                 medList){

  cisGeno = as.data.frame(ifelse(needMed,
                  abind::abind(lapply(c(biomInt,medList),
                                      getCisGenotypes
                                      locs = medLocs,
                                      snps = snps,
                                      snpLocs = snpLocs,
                                      cisDist = cisDist),
                               along = 1),
                  getCisGenotype(biomInt,
                                 medLocs,
                                 snps,
                                 snpLocs,
                                 cisDist)))

  W = t(cisGeno$snpCur)
  snpList = cisGeno$snpList

  if (!dir.exists('temp/')){ dir.create('temp/') }
  genfile = paste0('temp/',fileName,'.gen')
  samplefile = paste0('temp/',fileName,'.sample')
  bedfile = paste0('temp/',fileName)
  outFile = paste0('temp/',fileName)
  phenFile = paste0('temp/',fileName,'.phen')
  covarFile = paste0('temp/',fileName,'.qcovar')

  df = as.data.frame(matrix(nrow = 1,ncol =2))
  ids = row.names(W)
  geno = as.data.frame(matrix(ncol=nrow(W)+4,nrow = ncol(W)))
  colnames(geno) <- c('SNP','Pos','A1','A2',ids)
  geno[,5:ncol(geno)] = t(W)
  geno$SNP = snpList
  onlyThese <- snpLocs[snpLocs$snpid %in% geno$SNP,]
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
  sample$pheno[2:nrow(sample)] <- rnorm(nrow(sample)-1)
  sample[1,] <- c(0,0,0,'D','P')
  write.table(sample,samplefile,row.names=FALSE,
              col.names = TRUE, quote = FALSE)
  write.table(sample[,c(1,2,6)],phenFile,row.names = F,
              col.names = F, quote = F)
  write.table(cbind(sample[,1:2],t(covariates[1:dimNumeric,-1])),
              covarFile,row.names = F,col.names = F,quote=F)

  system(paste('plink','--gen',genfile,'--sample',samplefile,'--make-bed','--out',bedfile))

  a = data.table::fread(paste0(bedfile,'.fam'))
  a$V5 = 2
  data.table::fwrite(a,paste0(bedfile,'.fam'),col.names=F,row.names=F,quote=F,sep='\t')

  system(paste('gcta64',
               '--bfile',bedfile,
               '--autosome',
               '--make-grm',
               '--out',bedfile))
  system(paste('gcta64',
               '--reml',
               '--grm',bedfile,
               '--pheno',phenFile,
               '--reml-pred-rand',
               '-qcovar',covFile,
               '--out',bedfile))

  hsq = data.table::fread(paste0(bedfile,'.hsq'),
                          fill = T)

  h2 = hsq$Variance[4]
  P = hsq$Variance[9]

  return(list(h2 = h2, P = P))


}
