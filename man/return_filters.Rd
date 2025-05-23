% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Filtering.R
\name{return_filters}
\alias{return_filters}
\title{return_filters}
\usage{
return_filters(
  biotype = c("protein_coding"),
  impact = c("HIGH", "MODERATE"),
  existing = TRUE,
  population = 0.01,
  vaf = 0,
  depth = 1,
  alt.depth = 1,
  returnDF = FALSE,
  ...
)
}
\arguments{
\item{biotype}{A vector of biotypes to select (default: c("protein_coding"))}

\item{impact}{A vector of biotypes to select (default: c("HIGH", "MODERATE"))}

\item{existing}{bool of wether to filter by dbSNP/COSMIC (default: TRUE)}

\item{population}{numeric value to filter 1000g/gnomAD (default: 0.01)}

\item{vaf}{numeric value to filter vaf (default: 0)}

\item{depth}{numeric value to filter depth (default: 1)}

\item{alt.depth}{numeric value to filter alt allele depth (default: 1)}

\item{returnDF}{bool return dataframe rather than vector (default: FALSE)}
}
\value{
A list of the selected filters
}
\description{
A function to return the chosen filters
}
\examples{

return_filters()
return_filters(vaf = 0.10, depth = 50, alt.depth = 5)
return_filters(impact = c("LOW", "MODIFIER", "HIGH", "MODERATE"))

}
\keyword{VCF}
