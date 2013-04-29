\name{ldecomp.getDistances}
\alias{ldecomp.getDistances}

\title{
Calculate T2 and Q2 distances
}

\description{
\code{ldecomp.getDistances} is used to calculate T2 and Q2 distances for every object in the
\code{ldecomp} decomposition.
}

\usage{
ldecomp.getDistances(scores, loadings, residuals, tnorm = NULL)
}

\arguments{
  \item{scores}{a matrix with score values (nobj x ncomp).}
  \item{loadings}{a matrix with loading values (nvar x ncomp).}
  \item{residuals}{a matrix with data residuals.}
  \item{tnorm}{a vector with singular values for scores, used for normalization.}
}

\details{
The T2 values (Hotelling T2 values) are calculated for normalized scores, so in fact it is a squared Mahalanobis distance between the projection of an object to latent variable space and origin of the space. 

The Q2 values are squared distances from an original object to its projection. They show individual object variance, which is not explained by latent variable space.

The distances are calculated automatically when a \code{ldecomp} based model (e.g. PCA) is calibrated and/or applied.
}

\value{
Returns a list with following fields:
\item{Q2 }{a matrix with calculated Q2 values (nobj x ncomp)}
\item{T2 }{a matrix with calculated T2 values (nobj x ncomp)}
\item{tnorm}{a vector with singular values for scores normalization}
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\keyword{ ~kwd1 }
