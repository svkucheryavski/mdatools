\name{randtest}
\alias{randtest}
\title{Randomization test for PLS regression}

\usage{
randtest(x, y, ncomp = 15, center = T, scale = F, nperm = 1000, 
   sig.level = 0.05, silent = TRUE)
}

\description{
\code{randtest} is used to carry out randomization/permutation test for a PLS regression model 
}

\arguments{
  \item{x }{matrix with predictors.}
  \item{y }{vector or one-column matrix with response.}
  \item{ncomp }{maximum number of components to test.}
  \item{center }{logical, center or not predictors and response values.}
  \item{scale }{logical, scale (standardize) or not predictors and response values.}
  \item{nperm }{number of permutations.}
  \item{sig.level }{significance level.}
  \item{silent}{logical, show or not test progress.}
}

\value{
Returns an object of \code{randtest} class with following fields:
   \item{nperm }{number of permutations used for the test.} 
   \item{stat }{statistic values calculated for each component.} 
   \item{alpha }{alpha values calculated for each component.} 
   \item{statperm }{matrix with statistic values for each permutation.} 
   \item{corrperm }{matrix with correlation between predicted and reference y-vales for each permutation.} 
   \item{ncomp.selected }{suggested number of components.} 
}

\details{
The class implements a method for selection of optimal number of components in PLS1 regression based on the randomization test [1]. The basic idea is that for each component from 1 to \code{ncomp} a statistic T, which is a covariance between t-score (X score, derived from a PLS model) and the reference Y values, is calculated. By repeating this for randomly permuted Y-values a distribution of the statistic is obtained. A parameter \code{alpha} is computed to show how often the statistic T, calculated for permuted Y-values, is the same or higher than the same statistic, calculated for original data without permutations.

If a component is important, then the covariance for unpermuted data should be larger than the covariance for permuted data and therefore the value for \code{alpha} will be quie small (there is still a small chance to get similar covariance). This makes \code{alpha} very similar to p-value in a statistical test.

The \code{randtest} procedure calculates alpha for each component, the values can be observed using \code{summary} or \code{plot} functions. There are also several function, allowing e.g. to show distribution of statistics and the critical value for each component.

}

\references{
1. S. Wiklund et al. Journal of Chemometrics 21 (2007) 427-439.
}

\seealso{
Methods for \code{randtest} objects:
  \tabular{ll}{
      \code{print.randtest} \tab prints information about a \code{randtest} object.\cr
      \code{\link{summary.randtest}} \tab shows summary statistics for the test.\cr
      \code{\link{plot.randtest}} \tab shows bar plot for alpha values.\cr
      \code{\link{plotHist.randtest}} \tab shows distribution of statistic plot.\cr
      \code{\link{plotCorr.randtest}} \tab shows determination coefficient plot.\cr
   }
   
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
### Examples of using the test

## Get the spectral data from Simdata set and apply SNV transformation

data(simdata)

y = simdata$conc.c[, 3]
x = simdata$spectra.c
x = prep.snv(x)

## Run the test and show summary 
## (normally use higher nperm values > 1000)
r = randtest(x, y, ncomp = 4, nperm = 200, silent = FALSE)
summary(r)

## Show plots

par( mfrow = c(3, 2))
plot(r)
plotHist(r, comp = 3)
plotHist(r, comp = 4)
plotCorr(r, 3)
plotCorr(r, 4)
par( mfrow = c(1, 1))

}

