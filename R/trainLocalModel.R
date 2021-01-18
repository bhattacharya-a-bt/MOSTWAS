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
    
        df = cbind(pheno,t(as.matrix(covariates[,-1])))
    colnames(df)[1] = 'pheno'
    pheno = as.numeric(resid(lm(pheno ~ .,data = as.data.frame(df))))
    
    set.seed(seed)
    fin.model.enet = glmnet::cv.glmnet(midSNP.cur$genotypes[],
                                       pheno,
                                       alpha = .5,
                                       K = 5,
                                       keep = T)
    fin.model.blup = glmnet::cv.glmnet(midSNP.cur$genotypes[],
                                       pheno,
                                       alpha = 0,
                                       K = 5,
                                       keep = T)
    pred.enet = fin.model.enet$fit.preval[,which.min(fin.model.enet$cvm)]
    pred.blup = fin.model.blup$fit.preval[,which.min(fin.model.blup$cvm)]
    enet.R2 = adjR2(pheno,pred.enet)
    blup.R2 = adjR2(pheno,pred.blup)
    model = ifelse(blup.R2 >= enet.R2,
                   'LMM','Elastic net')
    
    if (model == 'Elastic net'){
        
       mod.df.enet = data.frame(SNP = midSNP.cur$map$marker.ID,
                                 Chromosome = midSNP.cur$map$chromosome,
                                 Position = midSNP.cur$map$physical.pos,
                                 Effect = as.numeric(coef(fin.model.enet,s='lambda.min'))[-1])
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
    
    if (model == 'LMM') {
        mod.df.blup = data.frame(SNP = midSNP.cur$map$marker.ID,
                                 Chromosome = midSNP.cur$map$chromosome,
                                 Position = midSNP.cur$map$physical.pos,
                                 Effect = as.numeric(coef(fin.model.blup,s='lambda.min'))[-1])
        return(list(Model = mod.df.blup,
                    Predicted = pred.blup,
                    CVR2 = adjR2(pheno,pred.blup)))
    }
}
