\name{plotModellingPower.simcam}
\alias{plotModellingPower.simcam}

\title{
Modelling power plot for SIMCA multiclass model
}

\description{
\code{plotModellingPower.simcam} is used to show plot with modelling power of variables for individual SIMCA model. 
}
\usage{
plotModellingPower(model, nc = 1, ...)
}

\arguments{
  \item{model}{an object of class \code{simcam} (a SIMCA multiclass model).}
  \item{nc}{number of individual SIMCA model (class) to make the plot for.}
  \item{\dots}{most of \code{\link{mdaplot}} parameters can be used.}
}

\details{
See examples for \code{\link{simcam}} class.
}

\seealso{
\code{\link{simcam}} to make SIMCAM model, \code{\link{simcamres}} to store and explore SIMCA multiclass classification results.
}

