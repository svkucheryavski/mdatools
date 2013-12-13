\name{plotMisclassified.classres}
\alias{plotMisclassified.classres}

\title{
Misclassified ratio plot
}

\description{
\code{plotMisclassified.classres} is used to show a plot with misclassified ratio values for different number of components.
}

\usage{
plotMisclassified(obj, nc = 1, ...)
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
\code{\link{plotSensitivity.classres}} to show sensitivity plot, \code{\link{plotSpecificity.classres}} to show specificity plot, 
\code{\link{plotPerformance.classres}} to show all performance statistics vs. complexity of
classification model on the same plot.
}
