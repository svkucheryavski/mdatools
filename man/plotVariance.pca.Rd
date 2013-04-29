\name{plotVariance.pca}
\alias{plotVariance.pca}

\title{
Explained variance plot for PCA model
}

\description{
\code{plotVariance.pca} is used to make a plot with variance, explained by each component in PCA model.
}

\usage{
plotVariance(obj, show.labels = F, show.legend = T)
}

\arguments{
  \item{obj}{an object of \code{pca} class (PCA model).}
  \item{show.labels}{logical, show or not labels for the plot points.}
  \item{show.legend}{logical, show or not a legend for the plot.}
}

\details{
If test set or cross-validation was used, the plot shows explained variance for the validation results as well. 
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
## make a PCA model for People data with full cross-validation and show variance plots

library(mdatools)
data(People)

model = pca(people, scale = T, cv = 1)

par(mfrow = c(1, 2))
plotVariance(model)
plotVariance(model, show.labels = T, show.legend = F)
par(mfrow = c(1, 1))
}

\keyword{ ~pca }
