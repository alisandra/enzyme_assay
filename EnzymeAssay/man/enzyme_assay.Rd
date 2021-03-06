% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/enzyme_assayClass.R
\name{enzyme_assay}
\alias{enzyme_assay}
\title{Getting Maximal Slopes From Enzyme Assay}
\usage{
enzyme_assay(raw_data, x, seconds = FALSE, decreasing = FALSE,
  windowsize = 8, method = c("individual", "together"))
}
\arguments{
\item{raw_data}{A matrix of measurements (rows) and wells (columns).}

\item{x}{The time measurements corresponding to each row, either in seconds, or formatted "H:MM:SS"}

\item{seconds}{Whether or not x is given in seconds. Defaults to FALSE}

\item{decreasing}{Whether or not assay is decreasing. Defaults to FALSE}

\item{windowsize}{Length of ranges to consider.}

\item{method}{Method for calculating range of maximum slope, "individual" for well by well, "together" for assay mean}
}
\description{
This functions generates an object of Enzyme_Assay class from raw data
}
\examples{
enzyme_assay(enzyme_assay.data, enzyme_assay.time)
}

