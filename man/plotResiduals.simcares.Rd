\name{plotResiduals.simcares}
\alias{plotResiduals.simcares}
\title{Residuals plot for SIMCA results}
\usage{
\method{plotResiduals}{simcares}(obj, ncomp = NULL, show.limits = T,
  type = "p", main = NULL, xlab = "T2", ylab = "Q2", legend = NULL,
  ...)
}
\arguments{
  \item{obj}{SIMCA results (object of class
  \code{simcares})}

  \item{ncomp}{which principal components to show the plot
  for}

  \item{show.limits}{logical, show or not lines with
  statistical limits for the residuals}

  \item{type}{type of the plot}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{legend}{vector with legend items}

  \item{...}{other plot parameters (see \code{mdaplot} for
  details)}
}
\description{
Shows a plot with Q2 vs. T2 residuals for SIMCA results
}
\details{
See examples in help for \code{\link{simcares}} function.
}

