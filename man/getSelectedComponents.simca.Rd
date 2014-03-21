\name{getSelectedComponents.simca}
\alias{getSelectedComponents.simca}
\title{Get selected components}
\usage{
\method{getSelectedComponents}{simca}(obj, ncomp = NULL)
}
\arguments{
  \item{obj}{a SIMCA model (object of class \code{simca})}

  \item{ncomp}{number of components to select, provided by
  user}
}
\description{
returns number of components depending on a user choice
}
\details{
Depedning on a user choice it returns optimal number of
component for the model (if use did not provide any value)
or check the user choice for correctness and returns it
back
}

