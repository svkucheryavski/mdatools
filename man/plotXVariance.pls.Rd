\name{plotXVariance.pls}
\alias{plotXVariance.pls}

\title{
X variance plot for PLS model
}

\description{
\code{plotXVariance.pls} is used to show variance of predictor variables, explained by each component in PLS model.
}

\usage{
plotXVariance(model, type = 'b', ...)
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
\code{\link{plotXCumVariance.pls}} - makes cumulative variance plot for predictors.
}

\examples{
## Examples of using X variance plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 1))
plotXVariance(model)
plotXVariance(model, type = 'h', show.labels = T)
par(mfrow = c(1, 1))

}
