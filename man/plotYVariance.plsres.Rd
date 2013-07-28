\name{plotYVariance.plsres}
\alias{plotYVariance.plsres}

\title{
Y variance plot for PLS results
}

\description{
\code{plotYVariance.plsres} is used to show variance of response variables, explained by each component in PLS model.
}

\usage{
plotYVariance(obj, type = 'b', ...)
}

\arguments{
  \item{obj}{PLS results (object of class \code{plsres}).}
  \item{type}{type of the plot.}
  \item{...}{most of the \code{\link{mdaplot}} parameters can be used.}
}

\seealso{
\code{\link{plotYCumVariance.plsres}} - makes cumulative variance plot for predictors.
}

\examples{
## Examples of using Y variance plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 1))
plotYVariance(model$calres, main = 'Y variance, calibration', show.labels = T)
plotYVariance(model$testres, main = 'Y variance, test', show.labels = T)
par(mfrow = c(1, 1))

}
