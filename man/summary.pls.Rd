\name{summary.pls}
\alias{summary.pls}
\title{Summary method for PLS model object}
\usage{
\method{summary}{pls}(object, ncomp = NULL, ny = NULL, ...)
}
\arguments{
  \item{object}{a PLS model (object of class \code{pls})}

  \item{ncomp}{how many components to use (if NULL - user
  selected optimal value will be used)}

  \item{ny}{which y variable to show the summary for (if
  NULL, will be shown for all)}

  \item{...}{other arguments}
}
\description{
Shows performance statistics for the model.
}

