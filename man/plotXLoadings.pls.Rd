\name{plotXLoadings.pls}
\alias{plotXLoadings.pls}

\title{
X loadings plot for PLS model
}

\description{
\code{plotXLoadings.pls} is used to make a plot with loadings for predictor values decomposition.
}
\usage{
plotXLoadings(model, comp = c(1, 2), type = 'p', ...)
}

\arguments{
  \item{model }{a PLS model (object of \code{pls} class.)}
  \item{comp }{either one value or vector with two values - number of components to make the plot for. }  
  \item{show.axes }{logical, show or not coordinate axes, that cross the origin (0, 0).}
  \item{type }{type of the plot.}
  \item{...}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
}

\examples{
## Examples of using X loadings plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c[, 2], ncomp = 8)

par(mfrow = c(2, 2))
plotXLoadings(model)
plotXLoadings(model, c(1, 3), show.axes = T, show.labels = T)
plotXLoadings(model, c(1, 3), type = 'l')
plotXLoadings(model, 1:4, type = 'l')
par(mfrow = c(1, 1))
}
