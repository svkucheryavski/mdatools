\name{plotCumVariance.pca}
\alias{plotCumVariance.pca}

\title{
Cumulative explained variance plot for PCA model
}

\description{
\code{plotCumVariance.pca} is used to make a plot with cumulative variance, explained by each component in PCA model.
}

\usage{
plotCumVariance(obj, show.labels = F, show.legend = T)
}
\arguments{
  \item{obj}{an object of \code{pca} class (PCA model).}
  \item{show.labels}{logical, show or not labels for the plot points.}
  \item{show.legend}{logical, show or not a legend for the plot.}
}
\details{
If test set or cross-validation was used, the plot shows cumulative explained variance for the validation results as well. 
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
}

\keyword{ ~pca }
