\name{plotYCumVariance.pls}
\alias{plotYCumVariance.pls}
\title{Cumulative explained Y variance plot for PLS}
\usage{
\method{plotYCumVariance}{pls}(obj, type = "b",
  main = "Y cumulative variance", xlab = "Components",
  ylab = "Explained variance, \%", show.legend = T, ...)
}
\arguments{
  \item{obj}{a PLS model (object of class \code{pls})}

  \item{type}{type of the plot('b', 'l' or 'h')}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{show.legend}{logical, show or not legend for the
  plot}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows plot with cumulative explained Y variance vs. number
of components.
}
\details{
See examples in help for \code{\link{pls}} function.
}

