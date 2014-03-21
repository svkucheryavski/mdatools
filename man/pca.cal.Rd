\name{pca.cal}
\alias{pca.cal}
\title{PCA model calibration}
\usage{
pca.cal(x, ncomp, center = T, scale = F, method = "svd")
}
\arguments{
  \item{x}{matrix with data values}

  \item{ncomp}{number of principal components to calculate}

  \item{center}{logical, do mean centering or not}

  \item{scale}{logical, do standardization or not}

  \item{method}{algorithm for compiting PC space (only
  'svd' is supported so far)}
}
\value{
an object with calibrated PCA model
}
\description{
Calibrates (builds) a PCA model for given data and
parameters
}

