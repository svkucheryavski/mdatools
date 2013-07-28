\name{ldecomp}
\alias{ldecomp}

\title{
Linear decomposition of data
}

\description{
\code{ldecomp} is used to store and visualize results for linear decomposition of data matrix X in form X = TP' + E (e.g. for PCA, MCR, PLS, etc). 
}

\usage{
ldecomp(scores, loadings, residuals, totvar, tnorm = NULL, ncomp.selected = NULL, 
   T2 = NULL, Q2 = NULL)
}

\arguments{
  \item{scores }{a matrix with score values (nobj x ncomp).}
  \item{loadings }{a matrix with loading values (nvar x ncomp).}
  \item{residuals }{matrix with data residuals.}
  \item{totvar }{total variance for data matrix.}
  \item{tnorm }{vector with singular values for scores normalization.}
  \item{ncomp.selected }{selected number of components for decomposition.}
  \item{T2 }{matrix with Hotelling T2 values for the data objects.}
  \item{Q2 }{matrix with Q2 residuals for the data objects.}
}

\details{
\code{ldecomp} is a general class for decomposition X = TP' + E. Here, X is a data matrix, T - matrix with scores, P - matrix with loadings and E - matrix with residuals. It is used, for example, for PCA results (\code{\link{pcares}}), in PLS and other methods. The class also includes methods for calculation and plotting residuals, variances, and so on.

There is no need to use the \code{ldecomp} manually. For example, when build PCA model with \code{\link{pca}} or apply it to a new data, the results will automatically inherit all methods of \code{ldecomp}.
}

\value{
Returns an object (list) of \code{ldecomp} class with following fields:
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
There is no need to use \code{ldecomp} methods manually. However, some of them might be useful:
  \tabular{ll}{
      \code{\link{ldecomp.getDistances}} \tab calculates T2 and Q2 distances.\cr
      \code{\link{ldecomp.getVariances}} \tab calculates explained variance.\cr
   }
}

\keyword{ ~kwd1 }
