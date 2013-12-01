\name{plotPredictions.classres}
\alias{plotPredictions.classres}

\title{
Plot with predicted class values
}

\description{
\code{plotPredictions.classres} is used to show plot with predicted class values for classification results
}

\usage{
plotPredictions(obj, nc = NULL, ncomp = NULL, ...)
}
\arguments{
  \item{obj}{classification results (object of \code{classress} class or similar).}
  \item{nc}{if results include several classes, for which class to show the plot for.}
  \item{ncomp}{number of components to show the plot for.}
  \item{\dots}{most of \code{\link{mdaplot}} parameters can be used.}
}

\details{
If class number is not provided, the results will be shown for all available classes. If number of
components is not provided, the results will be shown for selected (optimal) number.

See examples for \code{\link{simcares}}, \code{\link{simcamres}} and \code{\link{plsdares}}.
}

\seealso{
\code{\link{showPredictions.classres}} to show predicted class values as table.
}
