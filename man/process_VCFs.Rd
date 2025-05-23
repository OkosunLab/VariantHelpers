% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadingVCF.R
\name{process_VCFs}
\alias{process_VCFs}
\title{process_VCFs}
\usage{
process_VCFs(VCFs, JoinAll = TRUE, ...)
}
\arguments{
\item{VCFs}{a vector of paths to VCF files}

\item{JoinAll}{bool of whether to join all the files together (default TRUE)}
}
\value{
A dataframe of the VCF skipping the commented lines
}
\description{
A function to process a vector of VCF files
}
\examples{

process_VCFs(file)
}
\keyword{VCF}
