#' Train and predict gene's predictive model with mediators
#'
#' The function trains a predictive model of a given gene using top mediators
#' as fixed effects and assesses in-sample performance with cross-validation.
#'
#' @param format character, 'plink' or 'vcf'
#' @param fileName character, prefix only for plink or vcf files
#' @param doseFile character, output file name for dosage file
#' @param locsFile character, output file name for SNP location file
#' @param annotFile character, output file name for SNP annotation file
#' @param append logical, T/F if we are appending to existing files for large plink/vcf files
#'
#' @return final model for gene along with CV R2 and predicted values
#'
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#'
#' @export
convertToDosage <- function(format = 'plink',
                            fileName,
                            doseFile,
                            locsFile,
                            annotFile,
                            append = F){

  if (format == 'vcf'){

    system(paste('plink --vcf',
                 paste0(fileName,'.vcf'),
                 ' --make-bed --out',
                 fileName))

  }

  system(paste('plink --bfile',
               fileName,
               '--recode oxford --out',
               fileName))

  ### MAKE LOCATION FILES
  locGen = data.table::fread(paste0(fileName,'.gen'),
                             select = 1:3)
  locs = data.frame(snpid = locGen$V2,
                    chr = locGen$V1,
                    pos = locGen$V3)
  data.table::fwrite(locs,
                     locsFile,
                     append = append,
                     row.names = F,
                     quote = F,
                     sep = '\t')


  ### MAKE ANNOTATION FILES
  annotGen = data.table::fread(paste0(fileName,'.gen'),
                             select = c(2,4:5))
  colnames(annotGen) = c('SNP','REF','ALT')
  data.table::fwrite(annotGen,
                     annotFile,
                     append = append,
                     row.names = F,
                     quote = F,
                     sep = '\t')


  ### MAKE DOSAGE FILES
  tot = data.table::fread(paste0(fileName,'.gen'))
  tot = tot[,c(2,6:ncol(tot)),with=F]
  new.tot = data.frame(matrix(nrow = nrow(tot),
                              ncol = (ncol(tot)-1)/3 + 1))
  colnames(new.tot) = c('SNP',
                        data.table::fread(paste0(fileName,
                                                 '.fam'))$V2)

  genToDose = function(v){
    v = as.numeric(v)
    d = vector(mode = 'numeric',
               length = length(v)/3)
    for (i in 1:(length(v)/3)){
      this = v[(3*i-2):(3*i)]
      d[i] = which(this == 1)-1
    }
    return(d)
  }
  new.tot$SNP = tot$V2
  new.tot[,2:ncol(new.tot)] = apply(tot[,-1],1,genToDose)
  data.table::fwrite(new.tot,
                     doseFile,
                     append = append,
                     row.names = F,
                     quote = F,
                     sep = '\t')

}
