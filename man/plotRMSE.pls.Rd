\name{plotRMSE.pls}
\alias{plotRMSE.pls}

\title{
RMSE plot for PLS model
}

\description{
\code{plotRMSE.pls} is used to make a plot with root mean squared error values for calibration and validation vs. number of PLS components.
}

\usage{
plotRMSE(model, ny = 1, ...)
}

\arguments{
  \item{model }{a PLS model (object of \code{pls} class.)}
  \item{ny }{if response is multivariate, which variable to make the plot for.}
  \item{\dots}{most of the \code{\link{mdaplotg}} parameters can be used.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\examples{
## Examples of using RMSE plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 2))
plotRMSE(model)
plotRMSE(model, ny = 2, colmap = c('black', 'orange'))
plotRMSE(model, ny = 3, type = 'h', main = 'RMSE for C3', ylab = 'RMSE')
plotRMSE(model, ny = 1, type = 'h', show.labels = T)
par(mfrow = c(1, 1))

}
