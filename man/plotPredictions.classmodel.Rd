\name{plotPredictions.classmodel}
\alias{plotPredictions.classmodel}

\title{
Plot with predicted class values for a classification model
}

\description{
\code{plotPredictions.classmodel} is used to show plot with predicted class values for calibration validation predictions.
}

\usage{
plotPredictions(model, result = 'calres', ncomp = NULL, nc = NULL, ...)
}

\arguments{
  \item{model}{a classification model (object of class \code{simca}, \code{simcam}, \code{plsda}, etc).}
  \item{result}{which result to show predictions for (possible values: \code{"calres"}, \code{"cvres"}, \code{"testres"}.}
  \item{ncomp}{number of components to show the predictions for (if NULL, selected (optimal) number will be used).}
  \item{nc}{which class to make the plot for (if NULL predictions for all classes in the model will be shown).}
  \item{\dots}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
See examples for \code{\link{simca}}, \code{\link{simcam}} or \code{\link{plsda}} class.
}

