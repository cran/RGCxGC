% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_loadings.R
\name{make_loadings}
\alias{make_loadings}
\title{Import foreign model loadings}
\usage{
make_loadings(floadings, time, mod_time, acq_rate)
}
\arguments{
\item{floadings}{a numeric matrix with foreign loadings. With variables in
columns and eigenvalues as rows.}

\item{time}{a vector of length two with the time range of the
chromatographic run}

\item{mod_time}{the modulation time of the second dimension.}

\item{acq_rate}{the acquisition rate of the mass analyzer. 
The acquisition rate  is printed when read_chrom function is performed}
}
\description{
`make_loading` method takes the loading matrix obtained
 by a mixOmixs package and fold them into two-dimensional matrix
}
\details{
We strongly recommend to use the plsda function in the mixOmics
 package to perform partial least squares-discriminant analysis. The result
 of this model is a list containing a loading matrix.
 The method retrieves a matrix A with \emph{m} and \emph{n} dimensions.
 Where \emph{m} is the eigenvalues and \emph{n}
 is the number of loadings which the model returns.
}
