% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simGWAS.R
\name{simGWAS}
\alias{simGWAS}
\title{Simulate eQTL reference panel}
\usage{
simGWAS(geno.loc, geno.dis, ngwas, b_qtls.loc, b_qtls.dis, var_explained)
}
\arguments{
\item{geno.loc}{data.frame, reference genotypes for local SNPs}

\item{geno.dis}{data.frame, reference genotypes for distal SNPs}

\item{ngwas}{integer, sample size}

\item{b_qtls.loc}{vector, local eQTLs effect sizes}

\item{b_qtls.dis}{vector, distal eQTLs effect sizes}

\item{var_explained}{numeric, total variance explained in phenotype}
}
\value{
list of GWAS summary statistics, total eQTL effect size,
matrix of simulated genotypes, vector of trait
}
\description{
The function creates a eQTL reference panel for MOSTWAS
methods with genotypes, expression, mediator intensities, and QTLs.
}
