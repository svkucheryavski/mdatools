\name{plotPerformance.classres}
\alias{plotPerformance.classres}

\title{
Classification performance plot
}

\description{
\code{plotPerformance.classres} is used to show a plot with specificity and sensitivity values vs. model complexity for classification results.
}

\usage{
plotPerformance(obj, nc = 1, ...)
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
\code{\link{plotSensitivity.classres}} to show sensitivity plot, \code{\link{plotSpecificity.classres}} to show specificity plot.
}
