<<<<<<< HEAD
#' LD prune a genotype matrix
#'
#' The function prunes SNPs that are in LD using PLINK
#'
#' @param SNP_file_name character, file name for SNP dosages in MatrixEQTL format
#' @param snps_location_file_name character, file name for SNP locations in MatrixEQTL format
#' @param SNP_file_name_dis character, file name for SNP distal dosages
#' @param expression_file_name character, file name for gene expressions
#' @param gene_location_file_name character, file name for gene locations
#' @param mediators_file_name character, file name for mediator intensities
#' @param meds_location_file_name character, file name for mediator locations as genes
#' @param covariates_files_name character, file name for covariaties
#' @param output_file_name_loc_qtl character, file name for first local QTLs
#' @param output_file_name_dis_qtl character, file name for first distal QTLs
#' @param output_file_name_loc_med character, file name for second local QTLs
#' @param output_file_name_dis_med character, file name for second distal QTLs
#' @param p_loc_qtl numeric, P-value threshold for first QTLs for DePMA
#' @param p_dis_qtl numeric, P-value threshold for first QTLs for DePMA
#' @param FDRcut numeric, FDR-adjusted P-value threshold for second QTLs for DePMA
#' @param useModel MatrixEQTL model type
#' @param DePMA logical, T/F if this analysis is for DePMA
#'
#' @return list of pruned genotype matrix, vector of SNP names, and location data frame
#'
#' @import MatrixEQTL
#'
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#'
#' @export
eQTL_MOSTWAS <- function(SNP_file_name,
                         snps_location_file_name,
                         SNP_file_name_dis = NULL,
                         expression_file_name,
                         gene_location_file_name,
                         mediators_file_name,
                         meds_location_file_name,
                         covariates_file_name,
                         output_file_name_loc_qtl,
                         output_file_name_dis_qtl,
                         output_file_name_loc_med = NULL,
                         output_file_name_dis_med = NULL,
                         p_loc_qtl = 1e-6,
                         p_dis_qtl = 1e-6,
                         FDRcut = 0.05,
                         useModel = modelLINEAR,
                         DePMA = F
                       ){

  require(MatrixEQTL)

  errorCovariance = numeric()
  locDist = 3e6;

  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);

  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);

  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }

  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name      = output_file_name_dis_qtl,
    pvOutputThreshold     = p_dis_qtl,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis  = output_file_name_loc_qtl,
    pvOutputThreshold.cis = p_loc_qtl,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = locDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

  if (DePMA){

  a = data.table::fread(output_file_name_dis_qtl)
  snp_dis = a$SNP
  snp = data.table::fread(SNP_file_name)
  snp = subset(snp,SNP %in% snp_dis)
  data.table::fwrite(snp,
                     SNP_file_name_dis,
                     quote = F,
                     col.names = T,
                     row.names = F,
                     sep = '\t')


  pvOutputThreshold_loc = 1e-3;
  pvOutputThreshold_dis = 1e-30;
  locDist = 3e6;

  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name_dis);


  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(mediators_file_name);


  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }

  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(meds_location_file_name, header = TRUE, stringsAsFactors = FALSE);

  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name      = output_file_name_dis_med,
    pvOutputThreshold     = pvOutputThreshold_dis,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis  = output_file_name_loc_med,
    pvOutputThreshold.cis = pvOutputThreshold_loc,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = locDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

  a = fread(output_file_name_loc_med)
  a = subset(a,FDR < FDRcut)
  fwrite(a,
         output_file_name_loc_med,
         quote=F,
         col.names=T,
         row.names=F,
         sep='\t')
  }
}
=======
#' LD prune a genotype matrix
#'
#' The function prunes SNPs that are in LD using PLINK
#'
#' @param SNP_file_name character, file name for SNP dosages in MatrixEQTL format
#' @param snps_location_file_name character, file name for SNP locations in MatrixEQTL format
#' @param SNP_file_name_dis character, file name for SNP distal dosages
#' @param expression_file_name character, file name for gene expressions
#' @param gene_location_file_name character, file name for gene locations
#' @param mediators_file_name character, file name for mediator intensities
#' @param meds_location_file_name character, file name for mediator locations as genes
#' @param covariates_files_name character, file name for covariaties
#' @param output_file_name_loc_qtl character, file name for first local QTLs
#' @param output_file_name_dis_qtl character, file name for first distal QTLs
#' @param output_file_name_loc_med character, file name for second local QTLs
#' @param output_file_name_dis_med character, file name for second distal QTLs
#' @param p_loc_qtl numeric, P-value threshold for first QTLs for DePMA
#' @param p_dis_qtl numeric, P-value threshold for first QTLs for DePMA
#' @param FDRcut numeric, FDR-adjusted P-value threshold for second QTLs for DePMA
#' @param useModel MatrixEQTL model type
#' @param DePMA logical, T/F if this analysis is for DePMA
#'
#' @return list of pruned genotype matrix, vector of SNP names, and location data frame
#'
#' @import MatrixEQTL
#'
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#'
#' @export
eQTL_MOSTWAS <- function(SNP_file_name,
                         snps_location_file_name,
                         SNP_file_name_dis = NULL,
                         expression_file_name,
                         gene_location_file_name,
                         mediators_file_name,
                         meds_location_file_name,
                         covariates_file_name,
                         output_file_name_loc_qtl,
                         output_file_name_dis_qtl,
                         output_file_name_loc_med = NULL,
                         output_file_name_dis_med = NULL,
                         p_loc_qtl = 1e-6,
                         p_dis_qtl = 1e-6,
                         FDRcut = 0.05,
                         useModel = modelLINEAR,
                         DePMA = F
                       ){

  require(MatrixEQTL)

  errorCovariance = numeric()
  locDist = 3e6;

  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);

  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);

  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }

  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name      = output_file_name_dis_qtl,
    pvOutputThreshold     = p_dis_qtl,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis  = output_file_name_loc_qtl,
    pvOutputThreshold.cis = p_loc_qtl,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = locDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

  if (DePMA){

  a = data.table::fread(output_file_name_dis_qtl)
  snp_dis = a$SNP
  snp = data.table::fread(SNP_file_name)
  snp = subset(snp,SNP %in% snp_dis)
  data.table::fwrite(snp,
                     SNP_file_name_dis,
                     quote = F,
                     col.names = T,
                     row.names = F,
                     sep = '\t')


  pvOutputThreshold_loc = 1e-3;
  pvOutputThreshold_dis = 1e-30;
  locDist = 3e6;

  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name_dis);


  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(mediators_file_name);


  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }

  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(meds_location_file_name, header = TRUE, stringsAsFactors = FALSE);

  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name      = output_file_name_dis_med,
    pvOutputThreshold     = pvOutputThreshold_dis,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis  = output_file_name_loc_med,
    pvOutputThreshold.cis = pvOutputThreshold_loc,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = locDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

  a = fread(output_file_name_loc_med)
  a = subset(a,FDR < FDRcut)
  fwrite(a,
         output_file_name_loc_med,
         quote=F,
         col.names=T,
         row.names=F,
         sep='\t')
  }
}
