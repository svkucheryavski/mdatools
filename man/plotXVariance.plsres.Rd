\name{plotXVariance.plsres}
\alias{plotXVariance.plsres}

\title{
X variance plot for PLS results
}

\description{
\code{plotXVariance.plsres} is used to show variance of predictor variables, explained by each component in PLS model.
}

\usage{
plotXVariance(obj, type = 'b', ...)
}

\arguments{
  \item{obj}{PLS results (object of class \code{plsres}).}
  \item{type}{type of the plot.}
  \item{...}{most of the \code{\link{mdaplot}} parameters can be used.}
}

\seealso{
\code{\link{plotXCumVariance.plsres}} - makes cumulative variance plot for predictors.
}

\examples{
## Examples of using X variance plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 1))
plotXVariance(model$calres, main = 'X variance, calibration', show.labels = T)
plotXVariance(model$testres, main = 'X variance, test', show.labels = T)
par(mfrow = c(1, 1))

}
