\name{plotPredictions.classres}
\alias{plotPredictions.classres}

\title{
Plot with predicted class values
}

\description{
\code{plotPredictions.classres} is used to show plot with predicted class values for classification results
}

\usage{
plotPredictions(obj, nc = 1, ncomp = NULL, ...)
}
\arguments{
  \item{obj}{classification results (object of \code{classress} class or similar).}
  \item{nc}{if results include several classes, for which class to show the plot for.}
  \item{ncomp}{number of components to show the plot for.}
  \item{\dots}{most of \code{\link{mdaplot}} parameters can be used.}
}

\details{
If a classification method returns continuous values used to make decision about class membership (e.g. probability or PLS predictions) the values are shown instead of +1 or -1 for classification.

See examples for \code{\link{simcares}}, \code{\link{simcamres}} and \code{\link{plsdares}}.
}

\seealso{
\code{\link{showPredictions.classres}} to show predicted class values as table.
}
