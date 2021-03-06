% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputePIContinuousTraits.R
\name{imputePIContinousTraits}
\alias{imputePIContinousTraits}
\title{Imputation of missing data for continuous traits, mandatory phylogenetic information}
\usage{
imputeContinous(missingData, tree)
}
\arguments{
\item{missingData}{data.frame of 1 or more numeric columns containing NAs}

\item{tree}{phylo object}
}
\value{
a data.frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the
imputation
}
\description{
This function imputes missing data for continuous traits applying the phylopars approach.
}
