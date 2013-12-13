\name{predict.pls}
\alias{predict.pls}

\title{
Predict method for PLS model
}

\description{
predict.pls is used to apply a PLS model to a new data set.
}

\usage{
predict(model, x, y.ref = NULL)
}

\arguments{
  \item{model}{a PLS model (object of class \code{pls})}
  \item{x}{a matrix with predictor values.}
  \item{y.ref}{a matrix with reference values for response (optional).}
}

\details{
The method applies a PLS model to a new data. If response values are not provided, no prediction statistics will be available.
}

\value{Returns an object of \code{\link{plsres}} class.}


\seealso{
}

\examples{
## make a PLS model a apply it to a new data

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c[, 1], ncomp = 8)
model = selectCompNum(model, 3)

res = predict(model, spectra.t, conc.t[, 1])
summary(res)
plot(res)

}

