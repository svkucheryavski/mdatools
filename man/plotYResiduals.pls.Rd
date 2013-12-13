\name{plotYResiduals.pls}
\alias{plotYResiduals.pls}

\title{
Y residuals plot for PLS model
}

\description{
\code{plotYResiduals.pls} is used to show plot for Y residuals in PLS model.
}
\usage{
plotYResiduals.pls(model, ncomp = NULL, ny = 1, ...)
}

\details{
By default scatter plot is used and residuals are shown versus reference response values. If line or bar plots are used, the residuals shown versus object number.
}

\arguments{
  \item{model}{a PLS model (object of class \code{pls})}
  \item{ncomp}{a number of components to show the residuals for.} 
  \item{ny }{if response is multivariate, which variable to make the plot for.}
  \item{...}{most of the \code{\link{mdaplotg}} parameters can be used.}
}

\examples{
## Examples of using Y residuals plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, x.test = spectra.t, y.test = conc.t)
model = selectCompNum(model, 3)

par(mfrow = c(2, 2))
plotYResiduals(model)
plotYResiduals(model, ny = 2, ncomp = 2, show.labels = T)
plotYResiduals(model, ny = 2, ncomp = 2, type = 'l')
plotYResiduals(model, type = 'h', show.labels = T)
par(mfrow = c(1, 1))

}
