\name{plotMisclassified.classmodel}
\alias{plotMisclassified.classmodel}

\title{
Misclassified ratio plot for classification model
}

\description{
\code{plotMisclassified.classmodel} is used to show plot with misclassified ratio values for calibration and validation predictions for any classification model.
}

\usage{
plotMisclassified(model, nc = NULL, ...)
}

\arguments{
  \item{model}{a classification model (object of class \code{simca}, \code{simcam}, \code{plsda}, etc).}
  \item{nc}{which class to make the plot for (if NULL summary values for all classes from the model will be used).} 
  \item{\dots}{most of \code{\link{mdaplotg}} parameters can be used.}
}

\details{
See examples for \code{\link{simca}}, \code{\link{simcam}} or \code{\link{plsda}} class.
}

\seealso{
\code{\link{plotSensitivity.classmodel}},  \code{\link{plotSpecificity.classmodel}}.
}
