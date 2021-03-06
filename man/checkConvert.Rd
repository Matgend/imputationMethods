% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilsPack.R
\name{checkConvert}
\alias{checkConvert}
\title{Check + conversion columns class.}
\usage{
checkConvert(missingData)
}
\arguments{
\item{missingData}{array or dataframe of data with missing values}
}
\value{
return the missing dataset with the columns containing integers or characters as factors and columns
containing float as numeric.
}
\description{
This function checks if character and integer columns are factors and columns composed of float are
numeric. If not the case the columns are converted.
}
