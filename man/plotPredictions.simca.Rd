\name{plotPredictions.simca}
\alias{plotPredictions.simca}

\title{
Plot with predicted class values for SIMCA model
}

\description{
\code{plotPredictions.simca} is used to show plot with class values for calibration and cross-validation predictions.
}

\usage{
plotPredictions(model, ncomp = NULL, ...)
}

\arguments{
  \item{model}{an object of class \code{simca} (a SIMCA model).}
  \item{ncomp}{number of components to show the plot for.}
  \item{\dots}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
See examples for \code{\link{simca}} class.
}

\seealso{
\code{\link{simca}} to make SIMCA model, \code{\link{simcares}} to store and explore SIMCA one-class classification results.
}
