\name{as.matrix.regcoeffs}
\alias{as.matrix.regcoeffs}
\title{as.matrix method for regression coefficients class}
\usage{
\method{as.matrix}{regcoeffs}(x, ncomp = 1, ny = 1, ...)
}
\arguments{
  \item{x}{regression coefficients object (class
  \code{regcoeffs})}

  \item{ncomp}{number of components to return the
  coefficients for}

  \item{ny}{number of response variable to return the
  coefficients for}

  \item{...}{other arguments}
}
\description{
returns matrix with regression coeffocoents for given
response number and amount of components
}

