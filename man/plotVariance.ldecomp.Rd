\name{plotVariance.ldecomp}
\alias{plotVariance.ldecomp}

\title{
Explained variance plot for data decomposition
}

\description{
\code{plotVariance} is used to make a plot with data variance, explained by each component.
}
\usage{
plotVariance(obj, type = 'b', ...)
}
\arguments{
  \item{obj}{an object of \code{ldecomp} class (e.g. PCA results).}
  \item{type}{type of the plot.}  
  \item{...}{most of the \code{\link{mdaplot}} parameters can be used.}
}

\details{
The variance values are in percent of total variance.
}

\seealso{
\code{\link{plotCumVariance.ldecomp}} makes cumulative explained variance plot.
}

\examples{
## make a PCA model and show explained variance plot for calibration set

library(mdatools)
data(People)

model = pca(people, scale = T)

par(mfrow = c(2, 1))
plotVariance(model$calres)
plotVariance(model$calres, show.labels = T)
par(mfrow = c(1, 1))

}

\keyword{ ~kwd1 }
