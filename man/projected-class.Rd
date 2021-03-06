% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/A_classes.R
\docType{class}
\name{projected-class}
\alias{projected-class}
\title{Class projected}
\description{
The \emph{projected} class defines the superclass for projection methods,
specially for multiway principal component analysis and discriminant
analysis based on partial least squares. The class represents the 
convergence of in-package results (m_prcomp) and the foreing 
building model (PLS-DA) procedure.
}
\section{Slots}{

\describe{
\item{\code{loadings}}{The eigenvectors of each principal component.}

\item{\code{time}}{The time range of chromatographic run}

\item{\code{mod_time}}{modulation time of the second dimension}
}}

