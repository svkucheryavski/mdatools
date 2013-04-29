\name{plotScores.ldecomp}
\alias{plotScores.ldecomp}

\title{
Scores plot
}
\description{
\code{plotScores.pca} is used to make a scores plot for \code{ldecomp} objects (PCA or PLS results, etc) 
}

\usage{
plotScores(obj, comp = c(1, 2), cgroup = NULL, show.labels = F, 
show.colorbar = T, show.axes = T)
}
\arguments{
  \item{obj}{an object of \code{ldecomp} class (e.g. PCA results).}
  \item{comp}{either one value or vector with two values - number of components 
  to make the plot for. If one value is provided, the plot will be made for scores vs. objects.}
  \item{cgroup}{vector with values used for color grouping of plot points (number of values should 
  be the same as number of objects).}
  \item{show.labels}{logical, show or not labels for plot points.}
  \item{show.colorbar}{logical, show or not a colorbar legend if \code{cgroup} is provided.}
  \item{show.axes}{logical, show or not coordinate axes, that cross the origin (0, 0)ю}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
## make a PCA model and show various score plots for calibration set

library(mdatools)
data(People)

model = pca(people, scale = T)

par(mfrow = c(2, 2))
plotScores(model$calres, show.labels = T)
plotScores(model$calres, c(1, 3), show.labels = T, cgroup = people[, 'Income'])
plotScores(model$calres, 2, cgroup = people[, 'Beer'])
plotScores(model$calres, 2, cgroup = people[, 'Beer'], show.colorbar = F)
par(mfrow = c(1, 1))

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
