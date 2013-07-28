\name{plotXYScores.plsres}
\alias{plotXYScores.plsres}

\title{
XY scores plot for PLS results
}

\description{
\code{plotXYScores.plsres} is used to make a plot with scores values for PLS decomposition of predictors and responses.
}

\usage{
plotXYScores(obj, ncomp = 1, show.axes = T, ...)
}

\arguments{
  \item{obj }{PLS results (object of class \code{plsres}.)}
  \item{comp}{for which component to make the plot for. }
  \item{show.axes}{logical, show or not coordinate axes, that cross the origin (0, 0).}
  \item{...}{most of \code{\link{mdaplot}} parameters can be used.}
}

\details{
Also know as T-U plot. Useful to see importance of components and to detect outliers. The plot is not available for cross-validation results, since to scores are calculated.
}

\examples{
## Examples of using XY scores plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 8, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 2))
plotXYScores(model$calres)
plotXYScores(model$calres, 2, cgroup = conc.c[, 2])
plotXYScores(model$testres, show.axes = T)
plotXYScores(model$testres, 2, cgroup = conc.t[, 2])
par(mfrow = c(1, 1))
}
