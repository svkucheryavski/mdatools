\name{plotRegcoeffs.pls}
\alias{plotRegcoeffs.pls}

\title{
Regression coefficients plot for PLS model
}

\description{
\code{plotRegcoeffs.pls} is used to make a plot with regression coefficients for each predictor variable.
}

\usage{
plotRegcoeffs(model, ncomp = NULL, ny = 1, ...)
}

\arguments{
  \item{model }{a PLS model (object of \code{pls} class.)}
  \item{ncomp }{number of components to make the plot for (default: \code{ncomp.selected}).}
  \item{ny }{if response is multivariate, which variable to make the plot for.}
  \item{\dots}{most of the \code{\link{mdaplot}} parameters can be used.}
}

\examples{
## Examples of using regressions plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 8, cv = 1)
model = selectCompNum(model, 3)

par(mfrow = c(2, 2))
plotRegcoeffs(model)
plotRegcoeffs(model, ncomp = 4)
plotRegcoeffs(model, ny = 2)
plotRegcoeffs(model, ny = 2, ncomp = 4)
par(mfrow = c(1, 1))

}
