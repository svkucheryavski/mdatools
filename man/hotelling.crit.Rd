% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ldecomp.R
\name{hotelling.crit}
\alias{hotelling.crit}
\title{Calculate critical limits for distance values using Hotelling T2 distribution}
\usage{
hotelling.crit(nobj, ncomp, alpha = 0.05, gamma = 0.01)
}
\arguments{
\item{nobj}{number of objects in calibration set}

\item{ncomp}{number of components}

\item{alpha}{significance level for extreme objects}

\item{gamma}{significance level for outliers}
}
\value{
vector with four values: critical limits for given alpha and gamma, mean distance and DoF.
}
\description{
Calculate critical limits for distance values using Hotelling T2 distribution
}
