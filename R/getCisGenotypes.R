#' Gather the cis-genotypes for a given biomarker
#'
#' The function takes in a biomarker of interest (mediator or gene)
#' and gather the cis-genotypes in a user-defined window around
#' the biomarker.
#'
#' @param biomInt character, identifier for biomarker of interest
#' @param locs data frame, location file for the biomarkers
#' @param snps data frame, SNP dosages
#' @param snpLocs data frame, MatrixEQTL locations for SNPs
#' @param cisDist numeric, window for cis distance
#'
#' @return list with matrix of cis-genotypes, vector of SNP names, and data frame of SNP locations
#'
#' @export
getCisGenotypes <- function(biomInt,
                            locs,
                            snps,
                            snpLocs,
                            cisDist = 5e5){

  colnames(locs)[1] = 'snpid'
  colnames(snpLocs)[1] = 'snpid'

  locs = subset(locs, grepl(biomInt,snpid))
  thisSNP = subset(snpLocs, chr == locs$chr[1] &
                     pos < locs$pos[1]+cisDist & pos > locs$pos[1]-cisDist)
  thisSNP = thisSNP[!duplicated(thisSNP$snpid),]
  snpCur = subset(snps,SNP %in% thisSNP$snpid)
  snpCur = snpCur[!duplicated(snpCur$SNP),]
  snpList = snpCur$SNP
  snpCur = as.matrix(snpCur[,-1])
  rownames(snpCur) = snpList
  return(list(snpCur = snpCur,
              snpList = snpList,
              thisSNP = thisSNP))

}
