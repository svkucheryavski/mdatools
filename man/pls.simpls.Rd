\name{pls.simpls}
\alias{pls.simpls}
\title{SIMPLS algorithm}
\usage{
pls.simpls(x, y, ncomp)
}
\arguments{
  \item{x}{a matrix with x values (predictors)}

  \item{y}{a matrix with y values (responses)}

  \item{ncomp}{number of components to calculate}
}
\value{
a list with computed regression coefficients, loadings and
scores for x and y matrices, and weights.
}
\description{
SIMPLS algorithm for calibration of PLS model
}
\references{
[1]. S. de Jong. SIMPLS: An Alternative approach to partial
least squares regression. Chemometrics and Intelligent
Laboratory Systems, 18, 1993 (251-263).
}

