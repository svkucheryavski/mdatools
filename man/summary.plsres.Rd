\name{summary.plsres}
\alias{summary.plsres}

\title{
Summary for PLS results
}

\description{
\code{summary.plsres} is used to show summary statistics for PLS results.
}

\usage{
summary(obj, ncomp = NULL, ny = NULL)
}

\arguments{
  \item{obj }{PLS results (object of class \code{plsres}).}
  \item{ncomp }{which components to show the summary for (default: \code{ncomp.selected}), can be a vector.}
  \item{ny }{if response is multivariate, which variable to show the summary for.}
}

\details{
If \code{ny} is not specified, summary is shown for all response variables.
}

\examples{
## Examples of using summary method for PLS results

library(mdatools)
data(Simdata)

model = pls(spectra.c, conc.c, ncomp = 5, cv = 4)
model = selectCompNum(model, 3)

# show summary for cross-validated results
summary(model$cvres)
summary(model$cvres, ncomp = 1:4)
summary(model$cvres, ny = 2, ncomp = 3)
}
