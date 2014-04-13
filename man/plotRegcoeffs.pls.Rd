\name{plotRegcoeffs.pls}
\alias{plotRegcoeffs.pls}
\title{Regression coefficient plot for PLS}
\usage{
\method{plotRegcoeffs}{pls}(obj, ncomp = NULL, ny = 1, main = NULL,
  ylab = NULL, ...)
}
\arguments{
  \item{obj}{a PLS model (object of class \code{pls})}

  \item{ncomp}{how many components to use (if NULL - user
  selected optimal value will be used)}

  \item{ny}{number of response variable to make the plot
  for (if y is multivariate)}

  \item{main}{main plot title}

  \item{ylab}{label for y axis}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows plot with regression coefficient values for selected
components.
}
\details{
See examples in help for \code{\link{pls}} function.
}

