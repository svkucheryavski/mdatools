\name{plotXLoadings.pls}
\alias{plotXLoadings.pls}
\title{X loadings plot for PLS}
\usage{
\method{plotXLoadings}{pls}(obj, comp = c(1, 2), type = "p",
  main = "X loadings", ylab = NULL, xlab = NULL, show.axes = T, ...)
}
\arguments{
  \item{obj}{a PLS model (object of class \code{pls})}

  \item{comp}{which components to show the plot for (one or
  vector with several values)}

  \item{type}{type of the plot}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{show.axes}{logical, show or not a axes lines
  crossing origin (0,0)}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows plot with X loading values for selected components.
}
\details{
See examples in help for \code{\link{pls}} function.
}

