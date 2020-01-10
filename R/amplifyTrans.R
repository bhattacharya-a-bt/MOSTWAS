#' Write out mediator-based trans-eSNP effects
#'
#' The function takes in the SNP-weights for mediators and
#' writes out the amplified SNP-weights based on
#'
#' @param geneInt character, identifier for gene of interest
#' @param qtlFull data frame, mediator-gene qtl results
#' @param numMed numeric, user-defined maximum number of correlated mediators
#'
#' @return a vector of correlated biomarkers to be treated as mediators
#'
#'
#' @export
amplifyTrans <- function(whichMed,
                         medTrainList,
                         lmCaretObj){
  cur.mod = medTrainList[[whichMed]]$Model
  cur.mod$Effect[-1] = cur.mod$Effect[-1] * as.numeric(coef(summary(lmCaretObj))[whichMed+1,1])
  cur.mod$Mediator = names(medTrainList)[whichMed]
  cur.mod$SNP = as.character(cur.mod$SNP)
  return(cur.mod)
}
