\name{pca.nipals}
\alias{pca.nipals}
\title{NIPALS based PCA algorithm}
\usage{
pca.nipals(x, ncomp)
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
Calculates principal component space using non-linear
iterative partial least squares algorithm (NIPALS)
}
\references{
[1] Geladi, Paul; Kowalski, Bruce (1986), "Partial Least
Squares Regression:A Tutorial", Analytica Chimica Acta 185:
1-17
}

