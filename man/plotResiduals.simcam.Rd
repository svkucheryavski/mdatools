\name{plotResiduals.simcam}
\alias{plotResiduals.simcam}

\title{
Residuals plot for SIMCA multiclass model
}

\description{
\code{plotResiduals.simcam} is used to show plot with T2 vs Q2 residuals for SIMCA multiclass model. 
}

\usage{
plotResiduals(model, nc = 1, ...)
}

\arguments{
  \item{model}{an object of class \code{simcam} (a SIMCA multiclass model).}
  \item{nc}{number of model (class) to make the plot for.}
  \item{\dots}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
The plot is made for calibration data used to build individual SIMCA models.

See examples for \code{\link{simcam}} class.
}

\seealso{
\code{\link{simcam}} to make SIMCAM model, \code{\link{simcamres}} to store and explore SIMCA multiclass classification results.
}

