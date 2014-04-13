\name{plotCooman.simcamres}
\alias{plotCooman.simcamres}
\title{Cooman's plot for SIMCAM results}
\usage{
\method{plotCooman}{simcamres}(obj, nc = c(1, 2), type = "p",
  main = "Cooman's plot", xlab = NULL, ylab = NULL, show.limits = T,
  legend = NULL, ...)
}
\arguments{
  \item{obj}{SIMCAM results (object of class
  \code{simcamres})}

  \item{nc}{vector with two values - classes (SIMCA models)
  to show the plot for}

  \item{type}{type of the plot}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{show.limits}{logical, show or not lines with
  statistical limits for the residuals}

  \item{legend}{vector with legend items}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows a Cooman's plot for a pair of SIMCA models
}
\details{
See examples in help for \code{\link{simcamres}} function.
}

