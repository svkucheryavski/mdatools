\name{pcacvres}
\alias{pcacvres}
\title{
Cross-validation results for PCA model
}

\description{
\code{pcacvres} is used to store results for cross-validation of PCA model.
}
\usage{
pcacvres(T2, Q2, fullvar, ncom.selected = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T2 }{matrix with T2 distances calculated using cross-validation.}
  \item{Q2 }{matrix with Q2 residuals calculated using cross-validation.}
  \item{fullvar }{full variance of the data matrix.}
  \item{ncomp.selected }{selected number of components.}
}

\details{
There is no need to create a pcacvres object manually, it is created automatically when build a pca model (see \code{\link{pca}}). The object can be used to show summary and plots for the results.
}

\value{
Returns an object (list) of \code{pcacvres} and \code{ldecomp} classes with following fields:
\item{T2 }{matrix with T2 distances (nobj x ncomp).}
\item{Q2 }{matrix with Q2 distances (nobj x ncomp).}
\item{ncomp.selected }{selected number of components.}
\item{expvar }{explained variance for each component.}
\item{cumexpvar }{cumulative explained variance.}
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\note{
}

\seealso{
Methods for \code{pcacvres} objects:
  \tabular{ll}{
      \code{print.pcacvres} \tab shows information about the object.\cr
      \code{summary.pcacvres} \tab shows some statistics for the object.\cr
   }

Methods, inherited from \code{\link{ldecomp}} class:   
  \tabular{ll}{
      \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
      \code{\link{plotCumVariance.ldecomp}} \tab makes commulative explained variance plot.\cr
      \code{\link{plotResiduals.ldecomp}} \tab makes Q2 vs. T2 residuals plot.\cr
   }
   Check also \code{\link{pca}} and \code{\link{pcares}}.    
}
\examples{
## make PCA model with full cross-validation and show cv results
library(mdatools)
data(People)

model = pca(people, scale = T, cv = 1)
model = pca.selectCompNum(model, 5)

print(model$cvres)
summary(model$cvres)

par(mfrow = c(1, 2))
plotResiduals(model$cvres, show.labels = T)
plotVariance(model$cvres)
par(mfrow = c(1, 1))
}

\keyword{ ~pca }
