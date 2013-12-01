\name{selectCompNum.pls}
\alias{selectCompNum.pls}

\title{
Set number of optimal components for PLS model
}

\description{
\code{selectCompNum.pls} is used to select number of optimal components in a PLS model.
}

\usage{
selectCompNum(model, ncomp)
}

\arguments{
  \item{model}{a PLS model (object of class \code{pls}).}
  \item{ncomp}{number of components to select.}
}

\details{
}

\value{
Returns the same model as provided as an argument, but where the parameter \code{ncomp.selected} is set to \code{ncomp} and all distances and limits are recalculated.
}

\examples{

library(mdatools)
data(Simdata)
model = pls(spectra.c, conc.c[, 1], ncomp = 8, cv = 1)
model1 = selectCompNum(model, 3)
model2 = selectCompNum(model, 5)

par(mfrow = c(2, 2))
plotPredictions(model1)
plotRegcoeffs(model1)
plotPredictions(model2)
plotRegcoeffs(model2)
par(mfrow = c(1, 1))

}

