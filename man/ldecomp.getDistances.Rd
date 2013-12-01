\name{ldecomp.getDistances}
\alias{ldecomp.getDistances}

\title{
Calculate T2 and Q2 distances for data decomposition
}

\description{
\code{ldecomp.getDistances} is used to calculate Hotelling T2 values and Q2 residuals for every object in the \code{ldecomp} decomposition.
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

If \code{tnorm} values are not provided they calculated from \code{scores}. 
}

\value{
Returns a list with following fields:
\item{Q2 }{a matrix with calculated Q2 values (nobj x ncomp)}
\item{T2 }{a matrix with calculated T2 values (nobj x ncomp)}
\item{tnorm}{a vector with singular values for scores normalization}
}

\examples{
## In the example we use an R function princomp to get scores and loadings and 
## then calculate Q2 and T2 distances using ldecomp.getDistances() function

library(mdatools)
data(People)
people = prep.autoscale(people, center = T, scale = T)
pcares = princomp(people)

residuals = people - pcares$scores \%*\% t(pcares$loadings)
dist = ldecomp.getDistances(pcares$scores, pcares$loadings, residuals)

# show residuals plot for 5 components
nlv = 5
plot(dist$T2[, nlv], dist$Q2[, nlv], type = "p", col = 'blue', 
   pch = 16, xlab = 'T2', ylab = 'Q2')
text(dist$T2[, nlv], dist$Q2[, nlv], rownames(people), pos = 3)
}

\keyword{ ~kwd1 }
