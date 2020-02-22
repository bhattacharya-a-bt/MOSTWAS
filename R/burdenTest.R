#' Compute weighted burden test
#'
#' The function takes in a gene expression model in MOSTWAS form
#' and GWAS summary statistics and
#' carries out the weighted burden Z-test for a trait
#'
#' @param wgt character, name of expression model .RData file
#' @param snps data.frame, reference panel SNPs
#' @param sumStats data frame, GWAS summary statistics
#' @param snpAnnot data.frame, SNP REF/ALT annotations
#' @param onlyCis logical, T/F to consider only cis-component
#' @param beta character, colnames in sumStats that keeps the effect sizes
#' @param se character, colnames in sumStats that keeps the standard errors
#' @param chr character, colnames in sumStats that keeps the chromosome
#' @param pos character, colnames in sumStats that keeps the position
#' @param ref character, colnames in sumStats that keeps the reference allele
#' @param R2cutoff numeric, predictive R2 cutoff
#'
#' @return estimate of TME
#'
#' @export
burdenTest <- function(wgt,
                       snps,
                       sumStats,
                       snpAnnot = NULL,
                       onlyCis = F,
                       beta,
                       se,
                       chr,
                       pos,
                       ref,
                       R2cutoff,
                       alpha,
                       nperms = 1e3){


  bigZ = abs(qnorm(alpha))
  load(wgt)

  pieces = strsplit(wgt,'/')
  fn = pieces[[1]][length(pieces[[1]])]
  geneInt = strsplit(fn,'[.]')[[1]][1]

  if (R2 <= R2cutoff){
    return(paste0(geneInt,
                  ' is not predicted at R2 > ',
                  R2cutoff))
  }

  if (!'GenPos' %in% colnames(sumStats)){
    sumStats$GenPos = paste(sumStats$Chromosome,sumStats$Position,sep = ':')
  }


  if (onlyCis){

    Model = subset(Model,Mediator == 'Cis')
    if (CisR2 <= R2cutoff){
      return(paste0(geneInt,
                    ' is not locally predicted at R2 > ',
                    R2cutoff))
    }

  }

  require(dplyr)
  if ('Mediator' %in% colnames(Model)){
    Model =
      as.data.frame(Model %>%
                      dplyr::group_by(SNP,Chromosome,Position) %>%
                      dplyr::summarize(sum(Effect)))
    colnames(Model) = c('SNP','Chromosome','Position','Effect')
  } else {
  Model =
    as.data.frame(Model %>%
                    dplyr::group_by(SNP,Chromosome,Position) %>%
                    dplyr::summarize(sum(Effect)))
  colnames(Model) = c('SNP','Chromosome','Position','Effect')}
  Model$GenPos = paste(Model$Chromosome,Model$Position,sep = ':')

  sumS = subset(sumStats,GenPos %in% Model$GenPos)
  colnames(sumS)[which(colnames(sumS) == beta)] = 'Beta'
  colnames(sumS)[which(colnames(sumS) == se)] = 'SE'
  colnames(sumS)[which(colnames(sumS) == chr)] = 'Chromosome'
  colnames(sumS)[which(colnames(sumS) == pos)] = 'Position'
  colnames(sumS)[which(colnames(sumS) == ref)] = 'REF'
  Model = subset(Model,GenPos %in% sumS$GenPos)

  sumS = sumS[match(Model$GenPos,sumS$GenPos),]

  if (nrow(sumS) == 0){
    return('SNPs not found.')
  }

  if (is.null(snpAnnot)){

    snpAnnot = as.data.frame(matrix(nrow = nrow(Model),
                                    ncol = 3))
    colnames(snpAnnot) = c('SNP','REF','ALT')
    snpAnnot$SNP = Model$SNP
    snpAnnot$REF = sumS$REF
    snpAnnot$ALT = sapply(strsplit(as.character(snpAnnot$SNP),':'),
                          function(x) x[4])
    }

  annot = subset(snpAnnot, SNP %in% Model$SNP)
  Model = subset(Model, SNP %in% annot$SNP)
  sumS = subset(sumS, GenPos %in% Model$GenPos)

  annot = merge(annot,
                Model[,c('SNP','GenPos')],
                by = 'SNP')

  flip = function(i,
                  snpAnnot,
                  sumS){

    search = sumS[i,]
    ref = annot[annot$GenPos == search$GenPos,]
    if (ref$REF[1] != search$REF[1]){
      return(-1*search$Beta[1])
    } else {return(search$Beta)}

  }

  sumS$Flip = sapply(1:nrow(sumS),
                     flip,
                     snpAnnot = annot,
                     sumS = sumS)

  Z = as.numeric(sumS$Flip)/as.numeric(sumS$SE)
  snpCur = subset(snps, SNP %in% Model$SNP)
  twasZ = as.numeric(Model$Effect %*% Z)
  genos = as.matrix(snpCur[,-1])
  LD = genos %*% t(genos) / (ncol(genos)-1)
  twasr2pred = as.numeric(Model$Effect %*% LD %*% Model$Effect)

  if (twasr2pred <= 0){
    return(paste0(geneInt,' has zero predictive accuracy. Try a different model. There
    are plenty these days.'))
  } else {
    twas = as.numeric(twasZ/sqrt(twasr2pred))
    permute.p = 1
    if (abs(twas) >= bigZ){
      perms = replicate(nperms,sample(Model$Effect,
                                      replace = F))
      nullZ = t(perms) %*% Z
      permute.p = mean(as.numeric(abs(twas)) <= as.vector(nullZ))
      }
    return(list(Gene = geneInt,
                Z = twas,
                P = 2*pnorm(-abs(twas)),
                permute.P = permute.p))
  }


}
