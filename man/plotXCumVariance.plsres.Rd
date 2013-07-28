\name{plotXCumVariance.plsres}
\alias{plotXCumVariance.plsres}

\title{
Cumulative X variance plot for PLS results
}

\description{
\code{plotXCumVariance.plsres} is used to show cumulative variance of predictor variables, explained by each component in PLS model.
}

\usage{
plotXCumVariance(obj, type = 'b', ...)
}

\arguments{
  \item{obj}{PLS results (object of class \code{plsres}).}
  \item{type}{type of the plot.}
  \item{...}{most of the \code{\link{mdaplot}} parameters can be used.}
}

\seealso{
\code{\link{plotXVariance.plsres}} - makes variance plot for predictors.
}

\examples{
## Examples of using cumulative X variance plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 1))
plotXCumVariance(model$calres, main = 'X cumvariance, calibration', show.labels = T)
plotXCumVariance(model$testres, main = 'X cumvariance, test', show.labels = T)
par(mfrow = c(1, 1))

}
