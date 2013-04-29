\name{plotCumVariance.ldecomp}
\alias{plotCumVariance.ldecomp}

\title{
Cumulative explained variance plot
}
\description{
\code{plotVariance.ldecomp} is used to make a plot with cumulative data variance, explained by components.
}
\usage{
plotCumVariance(obj, show.labels = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj }{an object of \code{ldecomp} class (e.g. PCA results).}
  \item{show.labels }{logical, show or not labels with variance values for plot points.}
}

\details{
Cumulative explained variance for component N is a sum of explained variance values for components from 1 to N.
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\seealso{
\code{\link{plotVariance.ldecomp}} makes explained variance plot.
}


\examples{
## make a PCA model and show cumulative explained variance plot for calibration set

library(mdatools)
data(People)

model = pca(people, scale = T)

par(mfrow = c(1, 2))
plotCumVariance(model$calres)
plotCumVariance(model$calres, show.labels = T)
par(mfrow = c(1, 1))

}

\keyword{ ~kwd1 }
