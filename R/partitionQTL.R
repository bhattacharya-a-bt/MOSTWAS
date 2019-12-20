#' Run QTL analysis of mediators v. mRNA expression
#'
#' The function runs QTL analysis between mediators and
#' gene expression for a single training-test split
#' using MatrixEQTL (Shabalin 2012)
#'
#' @param mediator_file_name character, filename for mediators file
#' @param mediator_location_file_name character, filename for genomic positions of mediators
#' @param expression_file_name character, filename for gene expression file
#' @param gene_location_file_name character, filename for genomic positions of genes
#' @param covariates_file_name character, filename for covariates file
#' @param output_file_name_cis character, output file for cis-QTLs
#' @param output_file_name_tra character, output file for trans-QTLs
#' @param useModel constant, model for QTL analysis for Matrix_eQTL_engine
#' @param cisP numeric, P-value cutoff for cis-QTLs
#' @param traP numeric, P-value cutoff for tra-QTLs
#' @param errorCovariance matrix, error covariance matrix, defaults to identity
#' @param cisDist numeric, cis-QTL window
#'
#' @return writes out cis- and trans-QTLs between mediators and gene expression
#'
#' @import MatrixEQTL
#'
#' @export
partitionQTL <- function(mediator_file_name,
                         mediator_location_file_name,
                         expression_file_name,
                         gene_location_file_name,
                         covariates_file_name,
                         output_file_name_cis,
                         output_file_name_tra,
                         useModel = modelLINEAR,
                         cisP = 1e-6,
                         traP = 1e-6,
                         errorCovariance = numeric(),
                         cisDist = 1e6){

  if (!require(MatrixEQTL)) install.packages('MatrixEQTL')
  require(MatrixEQTL)

  ## Load genotype data

  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(mediator_file_name);

  ## Load gene expression data

  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);

  ## Load covariates

  if (!is.null(covariates_file_name)){
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }}

  ## Run the analysis
  snpspos = read.table(mediator_location_file_name,
                       header = TRUE,
                       stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name,
                       header = TRUE,
                       stringsAsFactors = FALSE);

  if (!is.null(covariates_file_name)){
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name      = output_file_name_tra,
    pvOutputThreshold     = traP,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis  = output_file_name_cis,
    pvOutputThreshold.cis = cisP,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);}

  if (is.null(covariates_file_name)){
    me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      output_file_name      = output_file_name_tra,
      pvOutputThreshold     = traP,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      output_file_name.cis  = output_file_name_cis,
      pvOutputThreshold.cis = cisP,
      snpspos = snpspos,
      genepos = genepos,
      cisDist = cisDist,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);}



  a = fread(output_file_name_cis)
  a = subset(a, SNP != gene)
  fwrite(a,output_file_name_cis,quote=F,col.names = T,row.names = F,sep='\t')

  a = fread(output_file_name_tra)
  a = subset(a, SNP != gene)
  fwrite(a,output_file_name_tra,quote=F,col.names = T,row.names = F,sep='\t')

}
