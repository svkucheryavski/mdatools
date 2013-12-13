\name{plotSensitivity.classmodel}
\alias{plotSensitivity.classmodel}

\title{
Sensitivity plot for classification model
}

\description{
\code{plotSensitivity.classmodel} is used to show plot with sensitivity values for calibration and validation predictions for any classification model.
}

\usage{
plotSensitivity(model, nc = NULL, ...)
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
\code{\link{plotSpecificity.classmodel}},  \code{\link{plotMisclassified.classmodel}}.
}
