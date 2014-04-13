\name{as.matrix.plsres}
\alias{as.matrix.plsres}
\title{as.matrix method for PLS results}
\usage{
\method{as.matrix}{plsres}(x, ncomp = NULL, ny = 1, ...)
}
\arguments{
  \item{x}{PLS results (object of class \code{plsres})}

  \item{ncomp}{number of components to calculate the
  statistics for}

  \item{ny}{for which response variable calculate the
  statistics for}

  \item{...}{other arguments}
}
\description{
Returns a matrix with model performance statistics for PLS
results
}

