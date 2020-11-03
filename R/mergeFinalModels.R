#' Output final MOSTWAS models
#'
#' The function computes the adjusted R2 between two numeric vectors
#'
#' @param MeTWASFolder character, folder name for MeTWAS models
#' @param DePMAFolder character, folder name for DePMA models
#' @param outFolder character, folder name for final models folder
#'
#' @return no return object, creates a folder with final MOSTWAS models
#'
#'
#' @export
mergeFinalModels <- function(MeTWASFolder,
                             DePMAFolder,
                             outFolder,
                             sumFile,
                             R2Cutoff = .01,
                             h2PCutoff = 0.05){

  if (!dir.exists(outFolder)){
    dir.create(outFolder,recursive = T)
  }

  metwas.mods = list.files(MeTWASFolder)
  depma.mods = list.files(DePMAFolder)

  df = as.data.frame(matrix(nrow=0,ncol = 3))
  colnames(df) = c('Gene','R2','Method')

  for (f in metwas.mods){

    load(file.path(MeTWASFolder,f))
    df.cur = data.frame(Gene = strsplit(f,'.wgt')[[1]][1],
                        R2 = R2,
                        h2 = h2,
                        h2.PValue = h2.Pvalue,
                        Method = 'MeTWAS')
    df = rbind(df,
               df.cur)

  }


  for (f in depma.mods){

    load(file.path(DePMAFolder,f))
    df.cur = data.frame(Gene = strsplit(f,'.wgt')[[1]][1],
                        R2 = R2,
                        h2 = h2,
                        h2.PValue = h2.Pvalue,
                        Method = 'DePMA')
    df = rbind(df,df.cur)

  }


  df = subset(df,R2 > R2Cutoff & h2.PValue < h2PCutoff)

  df = df[order(df$Gene,df$R2,decreasing = T),]
  df = df[!duplicated(df$Gene),]

  finalMods = c(file.path(MeTWASFolder,
                          paste0(df$Gene[df$Method == 'MeTWAS'],
                                 '.wgt.med.RData')),
                file.path(DePMAFolder,
                          paste0(df$Gene[df$Method == 'DePMA'],
                                 '.wgt.med.RData')))

  file.copy(finalMods,outFolder)
  write.table(df,sumFile,col.names = T,row.names = F, quote = F)



}
