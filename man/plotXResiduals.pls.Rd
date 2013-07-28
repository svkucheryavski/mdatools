\name{plotXResiduals.pls}
\alias{plotXResiduals.pls}

\title{
X residuals plot for PLS model
}

\description{
\code{plotXResiduals.pls} is used to show plot for Hotelling T2 values vs. Q2 residuals for X decomposition in PLS model.
}
\usage{
plotXResiduals.pls(model, ncomp = NULL, ...)
}

\arguments{
  \item{model}{a PLS model (object of class \code{pls})}
  \item{ncomp}{a number of components to show the residuals for (default: \code{ncomp.selected}).} 
  \item{...}{most of the \code{\link{mdaplotg}} parameters can be used.}
}

\examples{
## Examples of using X residuals plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)
model = selectCompNum(model, 3)

par(mfrow = c(2, 1))
plotXResiduals(model)
plotXResiduals(model, ncomp = 2, show.labels = T)
par(mfrow = c(1, 1))

}
