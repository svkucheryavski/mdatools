\name{plotPredictions.pls}
\alias{plotPredictions.pls}

\title{
Predictions plot for PLS model
}

\description{
\code{plotPredictions.pls} is used to make a plot with predicted vs. reference response values.
}

\usage{
plotPredictions(model, ncomp = NULL, ny = 1, show.lines = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model }{a PLS model (object of \code{pls} class.)}
  \item{ncomp }{number of components to make the plot for (default: \code{ncomp.selected}).}
  \item{ny }{if response is multivariate, which variable to make the plot for.}
  \item{show.lines }{logical, show or not the target lines for predicted vs. measured points}
  \item{\dots}{most of the \code{\link{mdaplotg}} parameters can be used.}
}

\details{
In fact only \code{model} parameter is needed to make the plot, other parameters are optional. 
}


\examples{
## Examples of using predictions plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 8, cv = 1)
model = selectCompNum(model, 3)

par(mfrow = c(1, 2))
plotPredictions(model)
plotPredictions(model, ny = 2, ncomp = 5, colmap = c('orange', 'green'))
par(mfrow = c(1, 1))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
