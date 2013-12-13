\name{plotSpecificity.classmodel}
\alias{plotSpecificity.classmodel}

\title{
Specificity plot for classification model
}

\description{
\code{plotSpecificity.classmodel} is used to show plot with specificity values for calibration and validation predictions for any classification model.
}

\usage{
plotSpecificity(model, nc = NULL, ...)
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
\code{\link{plotSensitivity.classmodel}},  \code{\link{plotMisclassified.classmodel}}.
}
