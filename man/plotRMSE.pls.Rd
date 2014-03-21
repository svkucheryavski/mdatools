\name{plotRMSE.pls}
\alias{plotRMSE.pls}
\title{RMSE plot for PLS}
\usage{
\method{plotRMSE}{pls}(obj, ny = 1, type = "b", main = "RMSE",
  xlab = "Components", ylab = NULL, show.legend = T, show.labels = F,
  ...)
}
\arguments{
  \item{obj}{a PLS model (object of class \code{pls})}

  \item{ny}{number of response variable to make the plot
  for (if y is multivariate)}

  \item{type}{type of the plot('b', 'l' or 'h')}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{show.legend}{logical, show or not legend for the
  plot}

  \item{show.labels}{logical, show or not labels for the
  plot elements}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows plot with root mean squared error values vs. number
of components for PLS model.
}
\details{
See examples in help for \code{\link{pls}} function.
}

