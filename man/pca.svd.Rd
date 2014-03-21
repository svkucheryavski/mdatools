\name{pca.svd}
\alias{pca.svd}
\title{Singular Values Decomposition based PCA algorithm}
\usage{
pca.svd(x, ncomp = NULL)
}
\arguments{
  \item{x}{a matrix with data values (preprocessed)}

  \item{ncomp}{number of components to calculate}
}
\value{
a list with scores, loadings and eigencalues for the
components
}
\description{
Computes principal component space using Singular Values
Decomposition
}

