\name{predict.plsda}
\alias{predict.plsda}
\title{PLS-DA predictions}
\usage{
\method{predict}{plsda}(object, x, c = NULL, cv = F, ...)
}
\arguments{
  \item{object}{a PLS-DA model (object of class
  \code{plsda})}

  \item{x}{a matrix with x values (predictors)}

  \item{c}{a vector with reference class values}

  \item{cv}{logical, are predictions for cross-validation
  or not}

  \item{...}{other arguments}
}
\value{
PLS-DA results (an object of class \code{plsdares})
}
\description{
Applies PLS-DA model to a new data set
}
\details{
See examples in help for \code{\link{plsda}} function.
}

