% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permuteTME.R
\name{permuteTME}
\alias{permuteTME}
\title{Perform permutation test for the total mediation effect}
\usage{
permuteTME(
  snp,
  expression,
  mediators,
  covs,
  nperms = 1000,
  parallel = "no",
  nc
)
}
\arguments{
\item{snp}{vector, SNP of interest}

\item{expression}{vector, gene expression of interest}

\item{mediators}{data frame, mediators of interest}

\item{covs}{data frame, covariates}

\item{nperms}{integer, number of permutations for the null distribution}

\item{parallel}{character, boot parallel input}

\item{nc}{integer, number of cores}
}
\value{
estimate of TME and the permutation P-value
}
\description{
The function takes in a SNP-mediator set-gene triplet
and computes the total absolute mediation effect and permutation
test P-value
}
