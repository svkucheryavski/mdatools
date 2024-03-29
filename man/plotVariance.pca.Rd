% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.R
\name{plotVariance.pca}
\alias{plotVariance.pca}
\title{Explained variance plot for PCA model}
\usage{
\method{plotVariance}{pca}(
  obj,
  type = "b",
  labels = "values",
  variance = "expvar",
  xticks = seq_len(obj$ncomp),
  res = obj$res,
  ylab = "Explained variance, \%",
  ...
)
}
\arguments{
\item{obj}{a PCA model (object of class \code{pca})}

\item{type}{type of the plot ("b", "l", "h")}

\item{labels}{what to use as labels (if \code{show.labels = TRUE})}

\item{variance}{which variance to show}

\item{xticks}{vector with ticks for x-axis}

\item{res}{list with result objects to show the variance for}

\item{ylab}{label for y-axis}

\item{...}{other plot parameters (see \code{mdaplotg} for details)}
}
\description{
Shows a plot with explained variance or cumulative explained variance for components.
}
\details{
See examples in help for \code{\link{pca}} function.
}
