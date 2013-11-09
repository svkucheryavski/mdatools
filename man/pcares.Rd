\name{pcares}
\alias{pcares}

\title{
Results of PCA decomposition
}

\description{
\code{pcares} is used to store results for PCA decomposition of data.
}

\usage{
pcares(scores, loadings, residuals, totvar, tnorm = NULL, ncomp.selected = NULL, 
   T2 = NULL, Q2 = NULL)
}
\arguments{
  \item{scores}{matrix with score values (nobj x ncomp).}
  \item{loadings}{matrix with loading values (nvar x ncomp).}
  \item{residuals }{matrix with data residuals.}
  \item{totvar }{total variance for data matrix.}
  \item{tnorm }{vector with singular values for scores normalization.}
  \item{ncomp.selected }{selected number of components for decomposition.}
  \item{T2 }{matrix with Hotelling T2 values for the data objects.}
  \item{Q2 }{matrix with Q2 residuals for the data objects.}
}

\details{
In fact \code{pcares} is a wrapper for \code{\link{ldecomp}} - general class for storing
results for linear decomposition X = TP' + E. So, most of the methods, arguments and 
returned values are inherited from \code{ldecomp}.

There is no need to create a \code{pcares} object manually, it is created automatically when build a PCA model (see \code{\link{pca}}) or apply the model to a new data (see \code{\link{predict.pca}}). The object can be used to show summary and plots for the results.
}
\value{
Returns an object (list) of class \code{pcares} and \code{ldecomp} with following fields:
\item{scores }{matrix with score values (nobj x ncomp).}
\item{residuals }{matrix with data residuals (nobj x nvar).}
\item{T2 }{matrix with T2 distances (nobj x ncomp).}
\item{Q2 }{matrix with Q2 distances (nobj x ncomp).}
\item{tnorm }{vector with singular values used for scores normalization.}
\item{ncomp.selected }{selected number of components.}
\item{expvar }{explained variance for each component.}
\item{cumexpvar }{cumulative explained variance.}
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}


\seealso{
Methods for \code{pcares} objects:
  \tabular{ll}{
      \code{print.pcares} \tab shows information about the object.\cr
      \code{summary.pcares} \tab shows statistics for the PCA results.\cr
   }

Methods, inherited from \code{\link{ldecomp}} class:   
  \tabular{ll}{
      \code{\link{plotScores.ldecomp}} \tab makes scores plot.\cr
      \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
      \code{\link{plotCumVariance.ldecomp}} \tab makes cumulative explained variance plot.\cr
      \code{\link{plotResiduals.ldecomp}} \tab makes Q2 vs. T2 residuals plot.\cr
   }
   Check also \code{\link{pca}} and \code{\link{ldecomp}}.    
}

\examples{
## make PCA model and show results for calibration set
library(mdatools)
data(People)

model = pca(people, scale = T)
model = selectCompNum(model, 5)

print(model$calres)
summary(model$calres)

par(mfrow = c(1, 2))
plotScores(model$calres, show.labels = T)
plotResiduals(model$calres, show.labels = T)
par(mfrow = c(1, 1))
}

\keyword{ ~pca }
