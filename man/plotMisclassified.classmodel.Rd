\name{plotMisclassified.classmodel}
\alias{plotMisclassified.classmodel}
\title{Misclassified ratio plot for classification model}
\usage{
\method{plotMisclassified}{classmodel}(obj, nc = NULL, ...)
}
\arguments{
  \item{obj}{classification model (object of class
  \code{plsda}, \code{simca}, etc.).}

  \item{nc}{if there are several classes, which class to
  make the plot for (NULL - summary for all classes).}

  \item{...}{most of the graphical parameters from
  \code{\link{mdaplotg}} function can be used.}
}
\description{
Makes a plot with misclassified ratio values vs. model
complexity (e.g. number of components)
}
\details{
See examples in description of \code{\link{plsda}},
\code{\link{simca}} or \code{\link{simcam}}.
}

