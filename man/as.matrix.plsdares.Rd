\name{as.matrix.plsdares}
\alias{as.matrix.plsdares}
\title{as.matrix method for PLS-DA results}
\usage{
\method{as.matrix}{plsdares}(x, ncomp = NULL, nc = NULL, ...)
}
\arguments{
  \item{x}{PLS-DA results (object of class
  \code{plsdares})}

  \item{ncomp}{number of components to calculate the
  statistics for}

  \item{nc}{for which class to calculate the statistics
  for}

  \item{...}{other arguments}
}
\description{
Returns a matrix with model performance statistics for
PLS-DA results
}

