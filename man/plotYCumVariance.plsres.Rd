\name{plotYCumVariance.plsres}
\alias{plotYCumVariance.plsres}

\title{
Cumulative Y variance plot for PLS results
}

\description{
\code{plotYCumVariance.plsres} is used to show cumulative variance of response variables, explained by each component in PLS model.
}

\usage{
plotYCumVariance(obj, type = 'b', ...)
}

\arguments{
  \item{obj}{PLS results (object of class \code{plsres}).}
  \item{type}{type of the plot.}
  \item{...}{most of the \code{\link{mdaplot}} parameters can be used.}
}

\seealso{
\code{\link{plotYVariance.plsres}} - makes variance plot for responses.
}

\examples{
## Examples of using cumulative Y variance plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 1))
plotYCumVariance(model$calres, main = 'Y cumvariance, calibration', show.labels = T)
plotYCumVariance(model$testres, main = 'Y cumvariance, test', show.labels = T)
par(mfrow = c(1, 1))

}
