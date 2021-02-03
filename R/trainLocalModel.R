#' Train a local-only model for a phenotype
#'
#' The function runs FUSION-type local-only modeling for a phenotype
#'
#' @param phenoInt character, phenotype name
#' @param midSNP snpObject, snp Object in the bigsnpr format
#' @param mediator data frame, mediator intensities
#' @param medLocs data frame, mediator locations
#' @param covariates data frame, covariates in MatrixEQTL format
#' @param cisDist numeric, local window, default .5 Mb
#' @param nfolds integer, number of folds
#' @param seed numeric, random seed
#' @param verbose logical, print gene name
#'
#' @return list with Model, predicted values, and CV-R2
#'
#' @importFrom bigsnpr snp_attach
#' @importFrom glmnet cv.glmnet
#'
#' @export
trainLocalModel <- function(phenoInt,
                           midSNP,
                           mediator,
                           medLocs,
                           covariates,
                           cisDist = 5e5,
                           nfolds = 5,
                           seed = 1218,
                           verbose = T){

  if (verbose){
    print(phenoInt)
  }
  colnames(mediator)[1] = 'Mediator'
  pheno = as.numeric(mediator[mediator$Mediator == phenoInt,-1])
  ml = subset(medLocs,geneid == phenoInt)
  w = which(midSNP$map$chromosome == ml$chr[1] &
              midSNP$map$physical.pos < ml$right[1] + cisDist &
              midSNP$map$physical.pos > ml$left[1] - cisDist)
  if (length(w) == 0){
    return(list(Model = data.frame(SNP = NA,
                                   Chromosome = NA,
                                   Position = NA,
                                   Effect = NA),
                Predicted = rep(NA,length = length(pheno)),
                CVR2 = 0))
  }
  midSNP.cur = bigsnpr::snp_attach(subset(midSNP,ind.col = w))

  set.seed(seed)
  train = caret::createFolds(y = pheno,
                             k=nfolds,
                             returnTrain = T)
  set.seed(seed)
  test = caret::createFolds(y = pheno,
                            k = nfolds,
                            returnTrain = F)

  pred.blup = pred.enet = vector(mode = 'numeric',
                                 length = length(pheno))
  df = cbind(pheno,t(as.matrix(covariates[,-1])))
  colnames(df)[1] = 'pheno'
  pheno = as.numeric(resid(lm(pheno ~ .,data = as.data.frame(df))))

  for (i in 1:5){
    mod.enet <- bigstatsr::big_spLinReg(midSNP.cur$genotypes,
                         pheno[train[[i]]],
                         ind.train = train[[i]],
                         alphas = .5,K = 5, warn = FALSE)
    pred.enet[test[[i]]] = predict(mod.enet,
                                   midSNP.cur$genotypes,
                        ind.row = test[[i]])

    mod.blup = rrBLUP::mixed.solve(y = pheno[train[[i]]],
                               Z = midSNP.cur$genotypes[train[[i]],])
    pred.blup[test[[i]]] = as.numeric(midSNP.cur$genotypes[test[[i]],] %*%
                                        mod.blup$u)
  }
  enet.R2 = adjR2(pheno,pred.enet)
  blup.R2 = adjR2(pheno,pred.blup)
  model = ifelse(blup.R2 >= enet.R2,
                 'LMM','Elastic net')

  if (model == 'Elastic net'){

    fin.model.enet = bigstatsr::big_spLinReg(midSNP.cur$genotypes,
                                             pheno,
                                             alphas = .5, K = 5, warn = FALSE)
    mod.df.enet = data.frame(SNP = midSNP.cur$map$marker.ID[attr(fin.model.enet,'ind.col')],
                             Chromosome = midSNP.cur$map$chromosome[attr(fin.model.enet,'ind.col')],
                             Position = midSNP.cur$map$physical.pos[attr(fin.model.enet,'ind.col')],
                             Effect = summary(fin.model.enet)$beta)
    colnames(mod.df.enet) = c('SNP','Chromosome','Position','Effect')
    mod.df.enet = subset(mod.df.enet,Effect != 0)
    if (nrow(mod.df.enet) == 0){
      model = 'LMM'
    } else {
      return(list(Model = mod.df.enet,
                  Predicted = pred.enet,
                  CVR2 = adjR2(pheno,pred.enet)))
    }
  }
}
