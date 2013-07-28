\name{plotXScores.pls}
\alias{plotXScores.pls}

\title{
X scores plot for PLS model
}

\description{
\code{plotXScores.pls} is used to make a plot with scores for predictor variables decomposition. 
}

\usage{
plotXScores(model, comp = c(1, 2), show.axes = F, ...)
}

\arguments{
  \item{model }{a PLS model (object of \code{pls} class.)}
  \item{comp}{either one value or vector with two values - number of components to make the plot for. }
  \item{show.axes}{logical, show or not coordinate axes, that cross the origin (0, 0).}
  \item{...}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
The plot shows x scores (T) for objects from calibration set and test set, if the last was provided. 
If one value is provided for \code{comp} parameters, the plot will be made for scores vs. objects number.
}

\examples{
## Examples of using X scores plot

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 8, x.test = spectra.t, y.test = conc.t)

par(mfrow = c(2, 2))
plotXScores(model)
plotXScores(model, c(1, 3), show.axes = T)
plotXScores(model, 2, colmap = c('orange', 'green'))
plotXScores(model, 1, type = 'l', ylab = 'LV1', xlab = 'Samples')
par(mfrow = c(1, 1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
