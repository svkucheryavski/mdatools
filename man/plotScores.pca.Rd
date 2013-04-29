\name{plotScores.pca}
\alias{plotScores.pca}

\title{
Scores plot for PCA model
}

\description{
\code{scoresPlot.pca} is used to make scores plot for a PCA model.
}

\usage{
plotScores(obj, comp = c(1, 2), show.labels = F, show.legend = T, show.axes = T)
}

\arguments{
  \item{obj}{an object of \code{pca} class (PCA model).}
  \item{comp}{either one value or vector with two values - number of components 
  to make the plot for. If one value is provided, the plot will be made for scores vs. objects.}
  \item{show.labels}{logical, show or not labels for the plot points.}
  \item{show.axes}{logical, show or not coordinate axes, that cross the origin (0, 0)ю}
  \item{show.legend}{logical, show or not a legend if test results are also presented in the plot.}
}

\details{
The plot shows scores for objects from calibration set and test set, if the last was provided. Color grouping (as in scores plot for PCA results) is not available here, because color is used to separate calibration and test objects. 
}


\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
## make a PCA model with test set validation and show score plots

library(mdatools)
data(People)

tdata = people[seq(1, nrow(people), 4), ]
cdata = people[-seq(1, nrow(people), 4), ]

model = pca(cdata, scale = T, test.data = tdata)

par(mfrow = c(2, 2))
plotScores(model)
plotScores(model, c(1, 3), show.labels = T)
plotScores(model, c(1, 3), show.labels = T, show.legend = F)
plotScores(model, 2, show.labels = T)
par(mfrow = c(1, 1))
}

\keyword{ ~pca }
