\name{plotVariance.pls}
\alias{plotVariance.pls}
\title{Variance plot for PLS}
\usage{
\method{plotVariance}{pls}(obj, decomp = "xdecomp", variance = "expvar",
  type = "b", main = "X variance", xlab = "Components",
  ylab = "Explained variance, \%", show.labels = F, show.legend = T, ...)
}
\arguments{
  \item{obj}{a PLS model (object of class \code{pls})}

  \item{decomp}{which decomposition to use ('xdecomp' for x
  or 'ydecomp' for y)}

  \item{variance}{which variance to use ('expvar',
  'cumexpvar)}

  \item{type}{type of the plot('b', 'l' or 'h')}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{show.labels}{logical, show or not labels for the
  plot elements}

  \item{show.legend}{logical, show or not legend for the
  plot}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows plot with variance values vs. number of components.
}
\details{
See examples in help for \code{\link{pls}} function.
}

