\name{plotCumVariance.ldecomp}
\alias{plotCumVariance.ldecomp}

\title{
Cumulative explained variance plot for data decomposition
}
\description{
\code{plotVariance.ldecomp} is used to make a plot with cumulative data variance, explained by components.
}
\usage{
plotCumVariance(obj, type = 'b', ...)
}

\arguments{
  \item{obj }{an object of \code{ldecomp} class (e.g. PCA results).}
  \item{type}{type of the plot.}  
  \item{...}{most of \code{\link{mdaplot}} parameters can be used.}
}

\details{
Cumulative explained variance for component N is a sum of explained variance values for components from 1 to N.
}

\seealso{
\code{\link{plotVariance.ldecomp}} makes explained variance plot.
}

\examples{
## make a PCA model and show cumulative explained variance plot for calibration set

library(mdatools)
data(People)

model = pca(people, scale = T)

par(mfrow = c(2, 1))
plotCumVariance(model$calres)
plotCumVariance(model$calres, show.labels = T)
par(mfrow = c(1, 1))

}

\keyword{ ~kwd1 }
