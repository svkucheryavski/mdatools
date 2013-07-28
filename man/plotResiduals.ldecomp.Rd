\name{plotResiduals.ldecomp}
\alias{plotResiduals.ldecomp}

\title{
Residuals plot for data decomposition
}

\description{
\code{ploResiduals} is used to make a plot with Hotelling T2 values vs. Q2 residuals. 
}

\usage{
plotResiduals(obj, ncomp = NULL, ...)
}

\arguments{
  \item{obj}{an object of \code{ldecomp} class (e.g. PCA results).}
  \item{ncomp}{number of "selected" components, the distances are calculated for.}
  \item{...}{most of the \code{\link{mdaplot}} parameters can be used.}
}

\details{
By default number of components, selected in a particular model, is used. But since the T2 and Q2 values in \code{ldecomp} are calculated for all components, it is possible to look at the plot for any number. 
}


\examples{
## make a PCA model and show residual plots for calibration set

library(mdatools)
data(People)

model = pca(people, scale = T)
model = selectCompNum(model, 5)

par(mfrow = c(2, 2))
plotResiduals(model$calres)
plotResiduals(model$calres, 2, show.labels = T)
plotResiduals(model$calres, cgroup = people[, 'Height'])
plotResiduals(model$calres, 2, cgroup = people[, 'Height'], show.colorbar = F)
par(mfrow = c(1, 1))
}

\keyword{ ~kwd1 }
