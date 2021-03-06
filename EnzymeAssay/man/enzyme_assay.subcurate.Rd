% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/enzyme_assayClass.R
\name{enzyme_assay.subcurate}
\alias{enzyme_assay.subcurate}
\title{Curate One Column}
\usage{
enzyme_assay.subcurate(object, ...)
}
\arguments{
\item{object}{An object of the Enzyme_Assay class}

\item{column}{The curent column}

\item{status}{'curated' or 'auto', whether the initial values are taken from slot 'auto_range' or 'curated_range'}

\item{at_once}{The number of columns/wells to plot at once}
}
\description{
This internal function handles curation of one well (column) of an enzyme assay
}
\examples{
ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
eaCurated <- enzyme_assay.curate(ea)
}

