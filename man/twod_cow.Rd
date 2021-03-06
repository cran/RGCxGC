% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twod_cow.R
\name{twod_cow}
\alias{twod_cow}
\title{Two-dimensional correlation optimized warping alignment}
\usage{
twod_cow(sample_chrom, ref_chrom, segments, max_warp)
}
\arguments{
\item{sample_chrom}{A GCxGC class chromatogram imported by read_chrom 
function or a preprocessed chromatogram.}

\item{ref_chrom}{A representative GCxGC chromatogram chosen to be the 
template which sample_chrom will be aligned.}

\item{segments}{A two integer vector with number of segments
which the first and second dimension will be divided, respectively.}

\item{max_warp}{A two integer vector with the maximum warping parameter.}
}
\description{
This is an adaptation of two-dimesional COW alignment, first implemented 
in MATLAB \insertCite{Tomasi2004}{RGCxGC}. 
This functions takes a sample chromatogram to be aligned 
against a reference. The argument [segment] will be used to split the whole
chromatogram in \emph{n} and \emph{m} parts the first and the second
dimension, respectively. The [max_warp] argument provides de maximum
tolerance of the signal transformation for the first and the second dimension
\insertCite{DabaoZhang2008}{RGCxGC}.
}
\examples{

GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
GB08 <- read_chrom(GB08_fl, 5L)
GB09 <- read_chrom(GB09_fl, 5L)
\donttest{
GB09_al <- twod_cow(sample_chrom = GB09, ref_chrom = GB08,
                    segments = c(20, 40), max_warp = c(2, 8))
}

}
\references{
\insertAllCited{}
}
