\name{plotCumVariance.pca}
\alias{plotCumVariance.pca}
\title{Cumulative explained variance plot for PCA}
\usage{
\method{plotCumVariance}{pca}(obj, xlab = "Components",
  ylab = "Explained variance, \%", main = "Cumulative variance", ...)
}
\arguments{
  \item{obj}{a PCA model (object of class \code{pca})}

  \item{main}{main title for the plot}

  \item{xlab}{label for x axis}

  \item{ylab}{label for y axis}

  \item{...}{other plot parameters (see \code{mdaplotg} for
  details)}
}
\description{
Shows a plot with cumulative explained variance for
components.
}
\details{
See examples in help for \code{\link{pca}} function.
}

