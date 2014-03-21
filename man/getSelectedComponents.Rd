\name{getSelectedComponents}
\alias{getSelectedComponents}
\title{Get selected components}
\usage{
getSelectedComponents(obj, ncomp = NULL)
}
\arguments{
  \item{obj}{a model object (e.g. PCA or PLS model)}

  \item{ncomp}{user selection of components}
}
\description{
Service function for user selection of components
}
\details{
the function is used with different plots, when user
selection is expected. It returns optimal number if user
did not specify the \code{ncomp} value and checks if the
specified value is correct.
}

