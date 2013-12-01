\name{plotYCumVariance.pls}
\alias{plotYCumVariance.pls}

\title{
Cumulative Y variance plot for PLS model
}

\description{
\code{plotYCumVariance.pls} is used to show cumulative variance of response variables, explained by each component in PLS model.
}

\usage{
plotYCumVariance(model, type = 'b', ...)
}

\arguments{
  \item{model}{a PLS model (object of class \code{pls}).}
  \item{type}{type of the plot.}
  \item{...}{most of the \code{\link{mdaplotg}} parameters can be used.}
}

\details{
Shows cumulative variance both for calibration and validation (cross-validation, test) results, if validation was used.
}

\seealso{
\code{\link{plotYVariance.pls}} - makes variance plot for Y.
}

\examples{
## Examples of using cumulative Y variance plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 1))
plotYCumVariance(model, legend.position = 'bottomright')
plotYCumVariance(model, type = 'h', show.labels = T, ylim = c(0, 120))
par(mfrow = c(1, 1))

}
