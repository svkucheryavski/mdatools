\name{plotCooman.simcam}
\alias{plotCooman.simcam}

\title{
Cooman's plot for SIMCA multiclass models
}

\description{
\code{plotCooman.simcam} is used to show Cooman's plot for two SIMCA models. 
}

\usage{
plotCooman(model, nc = c(1, 2), ...)
}

\arguments{
  \item{model}{an object of class \code{simcam} (a SIMCA mutliclass model).}
  \item{nc}{number of individual models (classes) to make the plot for.}
  \item{\dots}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
The plot is made for calibration data used to build individual SIMCA models.

See examples for \code{\link{simcam}} class.
}

\seealso{
\code{\link{simcam}} to make SIMCAM model, \code{\link{simcamres}} to store and explore SIMCA multiclass classification results.
}

\examples{
}
