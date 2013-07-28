\name{plotVariance.pca}
\alias{plotVariance.pca}

\title{
Explained variance plot for PCA model
}

\description{
\code{plotVariance.pca} is used to make a plot with variance, explained by each component in PCA model.
}

\usage{
plotVariance(obj, type = 'b', ...)
}

\arguments{
  \item{obj}{an object of \code{pca} class (PCA model).}
  \item{type}{type of the plot.}
  \item{...}{most of the \code{\link{mdaplotg}} parameters can be used.}
}

\details{
If test set or cross-validation was used, the plot shows explained variance for the validation results as well. 
}

\examples{
## make a PCA model for People data with full cross-validation and show variance plots

library(mdatools)
data(People)

model = pca(people, scale = T, cv = 1)

par(mfrow = c(2, 1))
plotVariance(model)
plotVariance(model, show.labels = T, type = 'h')
par(mfrow = c(1, 1))
}

\keyword{ ~pca }
