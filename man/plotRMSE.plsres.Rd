\name{plotRMSE.plsres}
\alias{plotRMSE.plsres}

\title{
RMSE plot for PLS results
}

\description{
\code{plotRMSE.plsres} is used to make a plot with root mean squared error values for calibration and validation vs. number of PLS components.
}

\usage{
plotRMSE(obj, ny = 1, ...)
}

\arguments{
  \item{obj }{PLS results (object of class \code{plsres}.)}
  \item{ny }{if response is multivariate, which variable to make the plot for.}
  \item{\dots}{most of the \code{\link{mdaplot}} parameters can be used.}
}

\examples{
## Examples of using RMSE plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 2))
plotRMSE(model$calres)
plotRMSE(model$calres, ny = 2, type = 'h', col = 'green')
plotRMSE(model$testres, main = 'RMSE, test set')
plotRMSE(model$testres, ny = 2, type = 'h', col = 'red')
par(mfrow = c(1, 1))

}
