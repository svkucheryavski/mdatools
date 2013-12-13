\name{plotXYScores.pls}
\alias{plotXYScores.pls}

\title{
XY scores plot for PLS model
}

\description{
\code{plotXYScores.pls} is used to make a plot with scores values for PLS decomposition of predictors and responses.
}

\usage{
plotXYScores(model, ncomp = 1, show.axes = T, ...)
}

\arguments{
  \item{model }{a PLS model (object of class \code{pls}.)}
  \item{comp}{for which component to make the plot for. }
  \item{show.axes}{logical, show or not coordinate axes, that cross the origin (0, 0).}
  \item{...}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
Also know as T-U plot. Useful to see importance of components and to detect outliers.
}

\examples{
## Examples of using XY scores plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 8, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 2))
plotXYScores(model)
plotXYScores(model, 2, legend.position = 'topleft')
plotXYScores(model, 3)
plotXYScores(model, colmap = c('orange', 'green'), pch = c(18, 17), main = 'T-U plot')
par(mfrow = c(1, 1))
}

