\name{plotXYLoadings.pls}
\alias{plotXYLoadings.pls}
\title{XY loadings plot for PLS}
\usage{
\method{plotXYLoadings}{pls}(obj, comp = c(1, 2), main = "XY loadings",
  ylab = NULL, xlab = NULL, show.axes = F, ...)
}
\arguments{
  \item{obj}{a PLS model (object of class \code{pls})}

  \item{comp}{which components to show the plot for (one or
  vector with several values)}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{show.axes}{logical, show or not a axes lines
  crossing origin (0,0)}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows plot with X and Y loading values for selected
components.
}
\details{
See examples in help for \code{\link{pls}} function.
}

