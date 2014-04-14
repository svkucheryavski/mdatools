\name{regcoeffs}
\alias{regcoeffs}
\title{Regression coefficients}
\usage{
regcoeffs(coeffs, ci.coeffs = NULL, ci.alpha = 0.1)
}
\arguments{
  \item{coeffs}{vector or matrix with regression
  coefficients}

  \item{ci.coeffs}{array (nobj x ncomp x ny x cv) with
  regression coefficients for computing confidence
  intervals (e.g. from jack-knifing)}

  \item{ci.alpha}{significance level for computing of the
  confidence intervals}
}
\value{
a list (object of \code{regcoeffs} class) with fields,
including: \tabular{ll}{ \code{values} \tab an array (nvar
x ncomp x ny) with regression coefficients \cr \code{ci}
\tab an array (nvar x ncomp x ny) with confidence intervals
for coefficients\cr \code{p.values} \tab an array (nvar x
ncomp x ny) with p-values for coefficients \cr } last two
fields are available if proper values for calculation of
the statistics were provided.
}
\description{
class for storing and visualisation of regression
coefficients for regression models
}
