\name{plotSensitivity.classres}
\alias{plotSensitivity.classres}

\title{
Sensitivity plot
}

\description{
\code{plotSensitivity.classres} is used to show a plot with sensitivity values vs. model complexity for classification results.
}

\usage{
plotSensitivity(obj, nc = 1, ...)
}

\arguments{
  \item{obj}{classification results (object of \code{classress} class or similar).}
  \item{nc}{if results include several classes, which class to show the plot for.}
  \item{\dots}{most of \code{\link{mdaplot}} parameters can be used.}
}

\details{
See examples for \code{\link{simcares}}, \code{\link{simcamres}} and \code{\link{plsdares}}.
}

\seealso{
\code{\link{plotSpecificity.classres}} to show specificity plot, \code{\link{plotPerformance.classres}} to show both specificity and sensitivity vs. complexity of classification model on the same plot.
}
