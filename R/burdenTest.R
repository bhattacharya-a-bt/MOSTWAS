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
                       R2cutoff){


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

  }

  Model = with(Model,aggregate(Effect,
                               list(SNP,
                                    Chromosome,
                                    Position),
                               sum))
  colnames(Model) = c('SNP','Chromosome','Position','Effect')
  Model$GenPos = paste(Model$Chromosome,Model$Position,sep = ':')

  sumS = subset(sumStats,GenPos %in% Model$GenPos)
  colnames(sumS)[which(colnames(sumS) == beta)] = 'Beta'
  colnames(sumS)[which(colnames(sumS) == se)] = 'SE'
  colnames(sumS)[which(colnames(sumS) == chr)] = 'Chromosome'
  colnames(sumS)[which(colnames(sumS) == pos)] = 'Position'
  colnames(sumS)[which(colnames(sumS) == ref)] = 'REF'
  colnames(sumS)[which(colnames(sumS) == alt)] = 'ALT'
  Model = subset(Model,GenPos %in% sumS$GenPos)

  sumS = sumS[match(Model$GenPos,sumS$GenPos),]

  if (is.null(snpAnnot)){

    snpAnnot = as.data.frame(matrix(nrow = nrow(Model),
                                    ncol = 3))
    colnames(snpAnnot) = c('SNP','REF','ALT')
    snpAnnot$SNP = Model$SNP
    snpAnnot$REF = sapply(strsplit(snpAnnot$SNP,':'),
                          function(x) x[3])
    snpAnnot$ALT = sapply(strsplit(snpAnnot$SNP,':'),
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

  Z = sumS$Flip/sumS$SE
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
    return(list(Gene = geneInt,
                Z = twas,
                P = 2*pnorm(-abs(twas))))
  }


}
