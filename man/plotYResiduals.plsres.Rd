\name{plotYResiduals.plsres}
\alias{plotYResiduals.plsres}

\title{
Y residuals plot for PLS results
}

\description{
\code{plotYResiduals.plsres} is used to show plot for Y residuals in PLS results.
}

\usage{
plotYResiduals.pls(obj, ncomp = NULL, ny = 1, ...)
}

\details{
Residuals are shown versus object number. Use color grouping (\code{cgroup}) to show relationship between residuals and actual y values.
}

\arguments{
  \item{obj}{PLS results (object of class \code{plsres})}
  \item{ncomp}{a number of components to show the residuals for.} 
  \item{ny }{if response is multivariate, which variable to make the plot for.}
  \item{...}{most of the \code{\link{mdaplot}} parameters can be used.}
}

\examples{
## Examples of using Y residuals plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)
model = selectCompNum(model, 3)

par(mfrow = c(2, 2))
plotYResiduals(model$calres)
plotYResiduals(model$calres, cgroup = conc.c[, 2], ncomp = 1, ny = 2)
plotYResiduals(model$testres)
plotYResiduals(model$testres, cgroup = conc.c[, 2], ncomp = 1, ny = 2)
par(mfrow = c(1, 1))

}
