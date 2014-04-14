\name{predict.simca}
\alias{predict.simca}
\title{SIMCA predictions}
\usage{
\method{predict}{simca}(object, x, c.ref = NULL, cv = F, ...)
}
\arguments{
  \item{object}{a SIMCA model (object of class
  \code{simca})}

  \item{x}{a matrix with x values (predictors)}

  \item{c.ref}{a vector with reference class names (same as
  class names for models)}

  \item{cv}{logical, are predictions for cross-validation
  or not}

  \item{...}{other arguments}
}
\value{
SIMCA results (an object of class \code{simcares})
}
\description{
Applies SIMCA model to a new data set
}
\details{
See examples in help for \code{\link{simca}} function.
}

