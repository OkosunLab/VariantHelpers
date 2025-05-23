% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadingVCF.R
\name{split_vep}
\alias{split_vep}
\title{split_vep}
\usage{
split_vep(VCF, header = FALSE, pick = TRUE, verbose = FALSE, ...)
}
\arguments{
\item{VCF}{the VCF file as a dataframe}

\item{header}{Vector of the names for the new annotations columns. Must be the same length as the number of annotation columns.}
}
\value{
A dataframe of the calls with the VEP annotations separated into columns. Multiple transcripts will be split into different columns.
}
\description{
A function to split VEP annotations into columns
}
\examples{

split_vep(VCF, header = FALSE)
}
\keyword{VCF}
