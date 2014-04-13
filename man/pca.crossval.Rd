\name{pca.crossval}
\alias{pca.crossval}
\title{Cross-validation of a PCA model}
\usage{
pca.crossval(model, x, cv, center = T, scale = F)
}
\arguments{
  \item{model}{a PCA model (object of class \code{pca})}

  \item{x}{a matrix with data values (calibration set)}

  \item{cv}{number of segments (if cv = 1, full
  cross-validation will be used)}

  \item{center}{logical, do mean centering or not}

  \item{scale}{logical, do standardization or not}
}
\value{
object of class \code{pcares} with results of
cross-validation
}
\description{
Does the cross-validation of a PCA model
}

