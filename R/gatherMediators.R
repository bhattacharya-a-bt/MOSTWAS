#' Extract mediating biomarkers for a given gene
#' 
#' The function runs aggregates top correlated biomarkers
#' to a given gene of interest and exports the list
#' 
#' @param geneInt character, identifier for gene of interest
#' @param qtlFull data frame, mediator-gene qtl results
#' @param numMed numeric, user-defined maximum number of correlated mediators
#' 
#' @return a vector of correlated biomarkers to be treated as mediators
#'
#' 
#' @export
gatherMediators <- function(geneInt,
                            qtlFull,
                            numMed){
  
  qtl = subset(qtlFull,gene == geneInt)
  qtl = qtl[order(qtl$FDR),]
  qtl = qtl[(1:(min(numMed,nrow(qtl)))),]
  return(qtl$SNP)
    
}