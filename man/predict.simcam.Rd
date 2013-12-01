\name{predict.simcam}
\alias{predict.simcam}

\title{
Predict method for SIMCA multiclass model
}

\description{
\code{predict.simcam} is used to classify a new data with SICMA multiclass model.
}

\usage{
predict(model, data, c.ref = NULL)
}

\arguments{
  \item{model}{an object of class \code{simcam} (a SIMCA multiclass model).}
  \item{data}{a numerical matrix with data.}
  \item{c.ref}{a vector with class membership for the data objects.}
}

\details{
The \code{c.ref} argument must be either a vector with class names or class numbers for each data object. In case if
names are used, the names should correspond to the class names of individual models. If numbers are used, it will 
be considered that 1 is for all objects from first individual model of \code{simcam} object, 2 - objects from the
second model and so on.

For details see examples for \code{\link{simcam}} object.
}

\value{
Returns object of class \code{simcamres} with results of SIMCA multiclass classification 
of the data. See \code{\link{simcamres}} examples for details.
}

\seealso{
\code{\link{simcam}} to make SIMCA multiclass model, \code{\link{simcamres}} to store and explore SIMCA multiclass classification results.
}

