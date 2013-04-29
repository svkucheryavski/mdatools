\name{ldecomp}
\alias{ldecomp}

\title{
Linear decomposition results.
}

\description{
\code{ldecomp} is used to store and visualize results for linear decomposition of data matrix X in form X = TP' + E (e.g. for PCA, MCR, PLS, etc). 
}

\usage{
ldecomp(scores, loadings, residuals, fullvar, tnorm = NULL, ncomp.selected = NULL)
}

\arguments{
  \item{scores }{a matrix with score values (nobj x ncomp).}
  \item{loadings }{a matrix with loading values (nvar x ncomp).}
  \item{residuals }{matrix with data residuals.}
  \item{fullvar }{full variance for data matrix.}
  \item{tnorm }{vector with singular values for scores normalization.}
  \item{ncomp.selected }{selected number of components for decomposition.}
}

\details{
\code{ldecomp} is a general class for decomposition X = TP' + E. It is used, for example, for PCA results (\code{\link{pcares}}), in PLS and other methods. The class also includes methods for calculation of residuals, variances, and so on.

There is no need to use the \code{ldecomp} manually. For example, when build PCA model with \code{\link{pca}}, the results will automatically inherit all methods of \code{ldecomp}.
}

\value{
Returns an object (list) of \code{ldecomp} class with following fields:
\item{scores }{matrix with score values (nobj x ncomp).}
\item{residuals }{matrix with data residuals (nobj x nvar).}
\item{T2 }{matrix with T2 distances (nobj x ncomp).}
\item{Q2 }{matrix with Q2 distances (nobj x ncomp).}
\item{ncomp.selected }{selected number of components.}
\item{expvar }{explained variance for each component.}
\item{cumexpvar }{cumulative explained variance.}
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\seealso{
Methods for \code{ldecomp} objects:
  \tabular{ll}{
      \code{print.ldecomp} \tab shows information about the object.\cr
      \code{summary.ldecomp} \tab shows statistics for the decomposition.\cr
      \code{\link{plotScores.ldecomp}} \tab makes scores plot.\cr
      \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
      \code{\link{plotCumVariance.ldecomp}} \tab makes commulative explained variance plot.\cr
      \code{\link{plotResiduals.ldecomp}} \tab makes Q2 vs. T2 residuals plot.\cr
   }

Other methods, implemented in the class:   
  \tabular{ll}{
      \code{\link{ldecomp.getDistances}} \tab calculates T2 and Q2 distances.\cr
      \code{\link{ldecomp.getVariances}} \tab calculates explained variance.\cr
      \code{\link{ldecomp.getResLimits}} \tab computes statistical limits for T2 and Q2 distances.\cr
   }
}

\keyword{ ~kwd1 }
