\name{plotDiscriminationPower.simcam}
\alias{plotDiscriminationPower.simcam}

\title{
Discrimination power plot for SIMCA multiclass model   
}

\description{
\code{plotDiscriminationPower.simcam} is used to show plot with discrimination power of variables for two SIMCA models. 
}
\usage{
plotDiscriminationPower(model, nc = c(1, 2), ...)
}

\arguments{
  \item{model}{an object of class \code{simcam} (a SIMCA multiclass model).}
  \item{nc}{vector with two values - number of models to make the plot for.}  
  \item{\dots}{most of \code{\link{mdaplot}} parameters can be used.}
}

\details{
See examples for \code{\link{simcam}} class.
}

\seealso{
\code{\link{simcam}} to make SIMCAM model, \code{\link{simcamres}} to store and explore SIMCA multiclass classification results.
}

