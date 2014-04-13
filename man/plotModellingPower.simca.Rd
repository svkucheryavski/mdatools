\name{plotModellingPower.simca}
\alias{plotModellingPower.simca}
\title{Modelling power plot for SIMCA model}
\usage{
\method{plotModellingPower}{simca}(obj, ncomp = NULL, type = "h",
  main = NULL, xlab = "Variables", ylab = "", ...)
}
\arguments{
  \item{obj}{a SIMCA model (object of class \code{simca})}

  \item{ncomp}{number of components to show the values for}

  \item{type}{type of the plot}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows a plot with modelling power values for each predictor
}

