\name{plotModelDistance.simcam}
\alias{plotModelDistance.simcam}
\title{Modelling distance plot for SIMCAM model}
\usage{
\method{plotModelDistance}{simcam}(obj, nc = 1, type = "h", main = NULL,
  xlab = "Models", ylab = "", xticklabels = NULL, ...)
}
\arguments{
  \item{obj}{a SIMCAM model (object of class
  \code{simcam})}

  \item{nc}{for which class (SIMCA model) to show the plot
  for}

  \item{type}{type of the plot}

  \item{main}{main plot title}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{xticklabels}{labels for x axis ticks}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows a plot with distance from data objects to a SIMCA
model
}
\details{
See examples in help for \code{\link{simcam}} function.
}

