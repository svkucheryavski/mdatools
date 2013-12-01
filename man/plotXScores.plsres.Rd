\name{plotXScores.plsres}
\alias{plotXScores.plsres}

\title{
X scores plot for PLS results
}

\description{
\code{plotXScores.plsres} is used to make a plot with scores for predictor variables decomposition. 
}

\usage{
plotXScores(obj, comp = c(1, 2), show.axes = F, ...)
}

\arguments{
  \item{obj }{PLS results (object of class \code{plsres}.)}
  \item{comp}{either one value or vector with two values - number of components to make the plot for. }
  \item{show.axes}{logical, show or not coordinate axes, that cross the origin (0, 0).}
  \item{...}{most of \code{\link{mdaplot}} parameters can be used.}
}

\details{
The plot shows predictor scores (T) for objects from calibration set and test set, if the last was provided. 
If one value is provided for \code{comp} parameters, the plot will be made for scores vs. objects number.
}

\examples{
## Examples of using X scores plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 8, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 2))
plotXScores(model$calres, show.labels = T)
plotXScores(model$calres, c(1, 3), show.axes = T, cgroup = conc.c[, 1])
plotXScores(model$testres)
plotXScores(model$testres, 2, type = 'h')
par(mfrow = c(1, 1))

}
