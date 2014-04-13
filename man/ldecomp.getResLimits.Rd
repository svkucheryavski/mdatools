\name{ldecomp.getResLimits}
\alias{ldecomp.getResLimits}
\title{Statistical limits for Q2 and T2 residuals}
\usage{
ldecomp.getResLimits(eigenvals, nobj, ncomp, alpha = 0.05)
}
\arguments{
  \item{eigenvals}{vector with eigenvalues}

  \item{nobj}{number of objects in data}

  \item{ncomp}{number of calculated components}

  \item{alpha}{significance level}
}
\value{
Returns a list with two vectors: \code{T2lim} and
\code{Q2lim}.
}
\description{
Computes statisticsl limits for Q2 and T2 residuals
}
\details{
T2 limits are calculated using Hotelling statistics.
}

