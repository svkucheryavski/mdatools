\name{plotLoadings.pca}
\alias{plotLoadings.pca}

\title{
Loadings plot for PCA model
}

\description{
\code{plotLoadings.pca} is used to make loadings plot (scatter or line) for a PCA model.
}

\usage{
plotLoadings(model, comp = c(1, 2), show.axes = T, ...)
}
\arguments{
  \item{model}{a PCA model (object of class \code{pca}).}
  \item{comp}{number of components to make the plot for (vector with 1 or more values).} 
  \item{show.axes}{logical, show or not coordinate axes, that cross the origin (0, 0).}
  \item{...}{most of the \code{\link{mdaplotg}} parameters can be used.}
}

\details{
If \code{comp} is a vector with two values and plot type is scatter, the loadings plot is made for one component vs. the second one. Otherwise the plot shows loading values vs. variables for each component.
}

\examples{
## make a PCA model for People data and show loading plots

library(mdatools)
data(People)

model = pca(people, scale = T)

par(mfrow = c(2, 2))
plotLoadings(model)
plotLoadings(model, c(1, 3, 4))
plotLoadings(model, c(1, 3), type = 'h')
plotLoadings(model, 1:3, type = 'l')
par(mfrow = c(1, 1))
}

\keyword{ ~pca }

