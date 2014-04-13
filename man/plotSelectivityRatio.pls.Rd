\name{plotSelectivityRatio.pls}
\alias{plotSelectivityRatio.pls}
\title{Selectivity ratio plot for PLS model}
\usage{
\method{plotSelectivityRatio}{pls}(obj, ncomp = NULL, ny = 1, type = "l",
  main = NULL, xlab = "Variables", ylab = "", ...)
}
\arguments{
  \item{obj}{a PLS model (object of class \code{pls})}

  \item{ncomp}{number of components to get the values for
  (if NULL user selected as optimal will be used)}

  \item{ny}{which response to get the values for (if y is
  multivariate)}

  \item{type}{type of the plot}

  \item{main}{main title for the plot}

  \item{xlab}{label for x axies}

  \item{ylab}{label for y axis}

  \item{...}{other plot parameters (see \code{mdaplot} for
  details)}
}
\description{
Shows a plot with selectivity ratio values for given number
of components and response variable
}
\references{
[1] Tarja Rajalahti et al. Chemometrics and Laboratory
Systems, 95 (2009), pp. 35-48.
}

