\name{predict.pls}
\alias{predict.pls}
\title{PLS predictions}
\usage{
\method{predict}{pls}(object, x, y.ref = NULL, cv = F, ...)
}
\arguments{
  \item{object}{a PLS model (object of class \code{pls})}

  \item{x}{a matrix with x values (predictors)}

  \item{y.ref}{a matrix with reference y values
  (responses)}

  \item{cv}{logical, are predictions for cross-validation
  or not}

  \item{...}{other arguments}
}
\value{
PLS results (an object of class \code{plsres})
}
\description{
Applies PLS model to a new data set
}
\details{
See examples in help for \code{\link{pls}} function.
}

