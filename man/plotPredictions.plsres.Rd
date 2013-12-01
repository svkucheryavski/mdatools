\name{plotPredictions.plsres}
\alias{plotPredictions.plsres}

\title{
Predictions plot for PLS results
}

\description{
\code{plotPredictions.plsres} is used to make a plot with predicted vs. reference (measured) response values.
}

\usage{
plotPredictions(obj, ncomp = NULL, ny = 1, show.line = T, ...)
}

\arguments{
  \item{obj }{PLS results (object of \code{plsres} class.)}
  \item{ncomp }{number of components to make the plot for (default: \code{ncomp.selected}).}
  \item{ny }{if response is multivariate, which variable to make the plot for.}
  \item{show.line }{logical, show or not a target line for predicted vs. measured points}
  \item{\dots}{most of the \code{\link{mdaplot}} parameters can be used.}
}


\examples{
## Examples of using predictions plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 8, cv = 4)
model = selectCompNum(model, 3)

par(mfrow = c(1, 2))
plotPredictions(model$cvres, ny = 2)
plotPredictions(model$cvres, ny = 2, ncomp = 5, col = 'red', pch = 18)
par(mfrow = c(1, 1))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
