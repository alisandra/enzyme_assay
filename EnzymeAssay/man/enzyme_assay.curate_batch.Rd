% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/enzyme_assayClass.R
\name{enzyme_assay.curate_batch}
\alias{enzyme_assay.curate_batch}
\title{Batch Curation}
\usage{
enzyme_assay.curate_batch(object, new_range)
}
\arguments{
\item{object}{An object of the Enzyme_Assay class}

\item{new_Range}{Matrix with start values for each column in slot 'data' in row one and end values in row two}
}
\description{
This function replaces slot 'curated_range' with user supplied ranges, and recalculates 'curated_slope'
}
\examples{
ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
new_range <- ea@auto_range
new_range[1,] <- 1
eaCurated <- enzyme_assay.curate_batch(ea, new_range)
}

