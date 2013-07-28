\name{plotCumVariance.pca}
\alias{plotCumVariance.pca}

\title{
Cumulative explained variance plot for PCA model
}

\description{
\code{plotCumVariance.pca} is used to make a plot with cumulative variance, explained by each component in PCA model.
}

\usage{
plotCumVariance(model, type = 'b')
}
\arguments{
  \item{obj}{a PCA model (object of class \code{pca}).}
  \item{type}{type of the plot.}
  \item{...}{most of the \code{\link{mdaplotg}} parameters can be used.}
}

\details{
If test set or cross-validation was used, the plot shows cumulative explained variance for the validation results as well. 
}

\examples{
library(mdatools)
data(People)

model = pca(people, scale = T, cv = 1)

par(mfrow = c(2, 1))
plotCumVariance(model, legend.position = 'topleft')
plotCumVariance(model, show.labels = T, type = 'h')
par(mfrow = c(1, 1))
}

\keyword{ ~pca }
