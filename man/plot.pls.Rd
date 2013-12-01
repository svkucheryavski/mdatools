\name{plot.pls}
\alias{plot.pls}

\title{
Plot overview of PLS model
}

\description{
\code{plot.pls} is used to show a set of plots (X residuals, regression coefficients, RMSE and predictions) for inspecting of a PLS model.
}

\usage{
plot(model, ncomp = NULL, ny = 1, show.labels = F)
}

\arguments{
  \item{model }{a PLS model (object of class \code{pls}).}
  \item{ncomp }{number of components to show the plot for (default: \code{ncomp.selected}).}
  \item{ny }{if response is multivariate, which variable to show the plot for.}
  \item{show.labels}{logical, show or not object labels for the data points.}
}

\examples{
## Examples of using plot function for PLS model

library(mdatools)
data(Simdata)
model = pls(spectra.c, conc.c, ncomp = 5, cv = 4, x.test = spectra.t, y.test = conc.t)
model = selectCompNum(model, 3)

# default
plot(model)

# for 2 components
plot(model, ncomp = 2, show.labels = T)

# for 5 components and second response variable
plot(model, ny = 2, ncomp = 5)

}
