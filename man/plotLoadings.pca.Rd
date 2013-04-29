\name{plotLoadings.pca}
\alias{plotLoadings.pca}

\title{
Loadings plot for PCA model
}
\description{
\code{plotLoadings.pca} is used to make loadings plot (scatter or line) for a PCA model.
}

\usage{
plotLoadings(obj, comp = c(1, 2), show.labels = T, show.legend = T, type = "p")
}
\arguments{
  \item{obj}{an object of \code{pca} class (PCA model).}
  \item{comp}{number of components to make the plot for. If two values is used a scatter plot will be shown (otherwise line plot). }
  \item{show.labels}{logical, show or not labels for the plot points.}
  \item{show.legend}{logical, show or not a legend if a line plot for several components is shown.}
  \item{type}{if \code{comp} has two values \code{type} allows to force line plot ('l') instead of scatter 'p'. Also possible to use line-scatter plot (type = 'b').}
}
\details{
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
## make a PCA model for People data and show loading plots

library(mdatools)
data(People)

model = pca(people, scale = T)

par(mfrow = c(2, 2))
plotLoadings(model)
plotLoadings(model, c(1, 3))
plotLoadings(model, c(1, 3), type = 'l')
plotLoadings(model, 1:3, type = 'b')
par(mfrow = c(1, 1))
}

\keyword{ ~pca }

