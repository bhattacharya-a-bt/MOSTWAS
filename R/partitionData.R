#' Create training-test splits on training set for MOSTWAS
#' 
#' The function splits SNP dosage, expression, mediator, and 
#' covariate datasets into a user-defined number of
#' training-test splits.
#' 
#' @param snpFile character, filename for SNP dosage file
#' @param expFile character, filename for expression file
#' @param mediatorfile character, filename for mediator file
#' @param covsFile character, filename for covsFile
#' @param k numeric, number of training-test splits
#' @param seed numeric, seed for splitting
#' @param modifyFile character, modification on filensames for each split
#' 
#' @return write out k files for each input files for each split
#' 
#' @example
#'  
#' partitionData('snps.txt','exp.txt','med.txt','covs.txt',k = 5,seed=1218,modifyFile = 'part')
#'
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom caret createFolds
#'
#' @export
partitionData <- function(snpFile,
                          expFile,
                          mediatorFile,
                          covsFile = NULL,
                          k,
                          seed = 1218,
                          modifyFile = 'partition'){
  
  interCols = colnames(data.table::fread(expFile,nrow=1))
  
  set.seed(seed)
  training_split = caret::createFolds(y = interCols[-1],
                                      k=k,
                                      returnTrain=T)
  
  extractPartition <- function(part,train,df,filename,modify){
    
    cur = df[,c(1,train[[part]]+1),with=F]
    data.table::fwrite(cur,paste(modify,k,filename,sep='_'),
           quote = F, col.names = T, 
           row.names = F, sep='\t',verbose = F)
    
  }
  
  snp = data.table::fread(snpFile)
  lapply(1:k,extractPartition,train = training_split,
         filename = snpFile, modify = modifyFile,df = snp)
  exp = data.table::fread(expFile)
  lapply(1:k,extractPartition,train = training_split,
         filename = expFile, modify = modifyFile,df = exp)
  mediator = data.table::fread(mediatorFile)
  lapply(1:k,extractPartition,train = training_split,
         filename = mediatorFile, modify = modifyFile, df = mediator)
  if (!is.null(covsFile)){
    cvrt = data.table::fread(covsFile)
  lapply(1:k,extractPartition,train = training_split,
         filename = covsFile, modify = modifyFile,df = cvrt)}
  
}