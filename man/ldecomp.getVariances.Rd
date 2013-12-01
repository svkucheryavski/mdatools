\name{ldecomp.getVariances}
\alias{ldecomp.getVariances}

\title{
Explained and cumulative explained variances for data decomposition
}

\description{
\code{ldecomp.getVariances} is used to calculate explained variance for every component in the \code{ldecomp} decomposition as well as cumulative explained variance.
}

\usage{
ldecomp.getVariances(Q2, totvar)
}

\arguments{
   \item{Q2 }{a matrix Q2 residuals (nobj x ncomp).}
   \item{totvar}{a total variance of the decomposed data.}
}

\details{
The explained variance for component N shows how much of data variance (in percent) is captured by this component. The cumulative explained variance shows how much of data variance is captured by components from 1 to N together. 
}

\value{
Returns a list with following fields:
\item{expvar }{a vector with explained variance values}
\item{cumexpvar }{a vector with cumulative explained variance values}
}

\examples{
## In the example we use an R function princomp to get scores and loadings and 
## then calculate Q2 and T2 distances and explained variance using ldecomp functions

library(mdatools)
data(People)
people = prep.autoscale(people, center = T, scale = T)
pcares = princomp(people)

residuals = people - pcares$scores \%*\% t(pcares$loadings)
dist = ldecomp.getDistances(pcares$scores, pcares$loadings, residuals)
var = ldecomp.getVariances(dist$Q2, sum(people^2))

show(var$cumexpvar)
}

\keyword{ ~kwd1 }
