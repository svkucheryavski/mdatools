\name{showPredictions.classres}
\alias{showPredictions.classres}

\title{
Predicted class values
}

\description{
\code{showPredictions.classres} is used to show predicted class values for classification results.
}

\usage{
showPredictions(obj, ncomp = NULL)
}

\arguments{
  \item{obj}{classification results (object of class \code{classres} or similar).}
  \item{ncomp}{number of components to show the predictions for.}
}

\details{
If \code{ncomp} is not provided, the results are shown for selected (optimal) number of components. For classification methods, which do not use latent variables, \code{ncomp} is equal to one or corresponds to another model complexity parameter if any.  

See examples for \code{\link{simcares}}, \code{\link{simcamres}} and \code{\link{plsdares}}.
}

