\name{plotYVariance.pls}
\alias{plotYVariance.pls}

\title{
Y variance plot for PLS model
}

\description{
\code{plotYVariance.pls} is used to show variance of response variables, explained by each component in PLS model.
}

\usage{
plotYVariance(model, type = 'b', ...)
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
\code{\link{plotYCumVariance.pls}} - makes cumulative variance plot for Y.
}

\examples{
## Examples of using Y variance plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 1))
plotYVariance(model)
plotYVariance(model, type = 'h', show.labels = T)
par(mfrow = c(1, 1))

}
