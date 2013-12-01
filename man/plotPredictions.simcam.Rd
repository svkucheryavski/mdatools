\name{plotPredictions.simcam}
\alias{plotPredictions.simcam}

\title{
Predictions plot for SIMCA multiclass model
}

\description{
\code{plotPredictions.simcam} is used to show plot with predicted values for SIMCA multiclass model.
}
\usage{
plotPredictions(model, ...)
}

\arguments{
  \item{model}{an object of class \code{simcam} (a SIMCA multiclass model).}
  \item{\dots}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
The plot is made for calibration data used to build individual SIMCA models.

See examples for \code{\link{simcam}} class.
}

\seealso{
\code{\link{simcam}} to make SIMCAM model, \code{\link{simcamres}} to store and explore SIMCA multiclass classification results.
}

