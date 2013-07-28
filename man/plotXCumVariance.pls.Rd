\name{plotXCumVariance.pls}
\alias{plotXCumVariance.pls}

\title{
Cumulative X variance plot for PLS model
}

\description{
\code{plotXCumVariance.pls} is used to show cumulative variance of predictor variables, explained by components in PLS model.
}

\usage{
plotXCumVariance(model, type = 'b', ...)
}

\arguments{
  \item{model}{a PLS model (object of class \code{pls}).}
  \item{type}{type of the plot.}
  \item{...}{most of the \code{\link{mdaplotg}} parameters can be used.}
}

\details{
Shows variance both for calibration and validation (cross-validation, test) results, if validation was used.
}

\seealso{
\code{\link{plotXVariance.pls}} - makes variance plot for X.
}

\examples{
## Examples of using cumulative X variance plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 1))
plotXCumVariance(model, legend.position = 'topleft')
plotXCumVariance(model, type = 'h', show.labels = T)
par(mfrow = c(1, 1))

}
