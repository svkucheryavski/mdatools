\name{plotResiduals.simcamres}
\alias{plotResiduals.simcamres}

\title{
Residuals plot for SIMCA multiclass classification results
}

\description{
\code{plotResiduals.simcamres} is used to show plot with T2 vs Q2 residuals for SIMCA multiclass classification results. 
}

\usage{
plotResiduals(obj, nc = 1, ...)
}

\arguments{
  \item{obj}{an object of class \code{simcamres} (SIMCA multiclass classification results).}
  \item{nc}{number of model (class) to make the plot for.}
  \item{\dots}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
See examples for \code{\link{simcamres}} class.
}

\seealso{
\code{\link{simcam}} to make SIMCAM model, \code{\link{simcamres}} to store and explore SIMCA multiclass classification results.
}

