\name{plotVariance.ldecomp}
\alias{plotVariance.ldecomp}

\title{
Explained variance plot
}

\description{
\code{plotVariance.ldecomp} is used to make a plot with data variance, explained by each component.
}
\usage{
plotVariance(obj, show.labels = F)
}
\arguments{
  \item{obj}{an object of \code{ldecomp} class (e.g. PCA results).}
  \item{show.labels}{logical, show or not labels with variance values for plot points.}
}

\details{
The variance values are in percent of total variance.
}

\seealso{
\code{\link{plotCumVariance.ldecomp}} makes cumulative explained variance plot.
}
\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
## make a PCA model and show explained variance plot for calibration set

library(mdatools)
data(People)

model = pca(people, scale = T)

par(mfrow = c(1, 2))
plotVariance(model$calres)
plotVariance(model$calres, show.labels = T)
par(mfrow = c(1, 1))

}

\keyword{ ~kwd1 }
