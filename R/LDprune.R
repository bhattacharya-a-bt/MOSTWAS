#' LD prune a genotype matrix
#'
#' The function prunes SNPs that are in LD using PLINK
#'
#' @param W matrix, dosage of SNPs
#' @param snpList character vector, vector of SNPs in W
#' @param snpLocs data frame, SNP locations in MatrixEQTL form
#' @param fileName character, name for PLINK intermediate files
#' @param windowSize integer, window size for PLINK pruning
#' @param numSNPShift integer, shifting window for PLINK pruning
#' @param ldThresh numeric, LD threshold for PLINK pruning
#'
#' @return list of pruned genotype matrix, vector of SNP names, and location data frame
#'
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#'
#' @export
LDprune <- function(W,
                    snpList,
                    snpLocs,
                    fileName,
                    windowSize,
                    numSNPShift,
                    ldThresh,
                    verbose = F){

  if (!dir.exists('temp/')){ dir.create('temp/') }
  genfile = paste0('temp/',fileName,'.gen')
  samplefile = paste0('temp/',fileName,'.sample')
  bedfile = paste0('temp/',fileName)
  outFile = paste0('temp/',fileName)

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

  if (!verbose) {sink('/dev/null')}
  system(paste('plink',
               '--gen',genfile,
               '--sample',samplefile,
               '--make-bed --allow-no-sex',
               '--out',bedfile),
         intern = !verbose)

  a = data.table::fread(paste0(bedfile,'.fam'))
  a$V5 = 2
  data.table::fwrite(a,paste0(bedfile,'.fam'),
                     col.names=F,
                     row.names=F,
                     quote=F,
                     sep='\t')

  system(paste('plink','--bfile',bedfile,
               '--indep-pairwise',windowSize,numSNPShift,ldThresh,
               '--out',outFile),
         intern = !verbose)

  if (!verbose) {sink()}

  file.remove(genfile,samplefile)

  s <- as.character(data.table::fread(paste0(outFile,
                                             '.prune.in'),
                                      header=F)$V1)

  W <- W[,which(snpList %in% s)]
  q <- snpList[snpList %in% s]
  snpList <- q[match(s,q)]
  onlyThese <- subset(onlyThese,snpid %in% snpList)
  snpList = snpList[snpList %in% onlyThese$snpid]

  return(list(W = W,
              snpList = snpList,
              onlyThese = onlyThese))


}
