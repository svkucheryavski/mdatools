\name{summary.pls}
\alias{summary.pls}

\title{
Summary for PLS model
}

\description{
\code{summary.pls} is used to show performance statistics for calibration and validation of a PLS model.
}

\usage{
summary(model, ncomp = NULL, ny = NULL)
}

\arguments{
  \item{model }{a PLS model (object of class \code{pls}).}
  \item{ncomp }{number of components to show the summary for (default: \code{ncomp.selected}).}
  \item{ny }{if response is multivariate, which variable to show the summary for.}
}

\details{
If \code{ny} is not specified, summary is shown for all response variables.
}

\examples{
## Examples of using summary method for PLS model

library(mdatools)
data(Simdata)
model = pls(spectra.c, conc.c, ncomp = 5, cv = 4, x.test = spectra.t, y.test = conc.t)
model = selectCompNum(model, 3)

summary(model)
summary(model, ncomp = 2)
summary(model, ny = 2, ncomp = 4)
}
