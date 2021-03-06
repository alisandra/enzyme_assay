% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/enzyme_assayClass.R
\name{enzyme_assay.getinput}
\alias{enzyme_assay.getinput}
\title{Handling User Input}
\usage{
enzyme_assay.getinput(object, ...)
}
\arguments{
\item{object}{An object of the Enzyme_Assay class}

\item{oldvalue}{Previous start or stop value}
}
\description{
This internal function asks the user to curate start or stop measurement for determining slope in an enzyme assay
}
\examples{
ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
eaCurated <- enzyme_assay.curate(ea)
}

