\name{pls}
\alias{pls}
\title{Partial Least Squares regression}

\usage{
  pls(x, y, ncomp = 20, center = T, scale = F, cv = NULL,
    x.test = NULL, y.test = NULL, method = 'simpls', info = '')
}

\description{
\code{pls} is used to calibrate, validate and use of partial least squares (PLS) regression model.
}

\arguments{
  \item{x}{matrix with predictors.}
  \item{y}{matrix with responses.}
  \item{ncomp }{maximum number of components to calculate.}
  \item{center }{logical, center or not predictors and response values.}
  \item{scale }{logical, scale (standardize) or not predictors and response values.}
  \item{cv}{number of segments for cross-validation (if cv = 1, full cross-validation will be used).}
  \item{x.test }{matrix with predictors for test set.}
  \item{y.test}{matrix with responses for test set.}
  \item{method }{method for calculating PLS model.}
  \item{info }{short text with information about the model.}
}

\value{
Returns an object of \code{pls} class with following fields:
   \item{ncomp }{number of components included to the model.} 
   \item{ncomp.selected }{selected (optimal) number of components.} 
   \item{xloadings }{matrix with loading values for x decomposition.} 
   \item{yloadings }{matrix with loading values for y decomposition.} 
   \item{weights }{matrix with PLS weights.} 
   \item{coeffs }{matrix with regression coefficients calculated for each component.}   
   \item{info }{information about the model, provided by user when build the model.} 
   \item{calres }{an object of class \code{\link{plsres}} with PLS results for a calibration data.} 
   \item{testres }{an object of class \code{\link{plsres}} with PLS results for a test data, if it was provided.} 
   \item{cvres }{an object of class \code{\link{plsres}} with PLS results for cross-validation, if this option was chosen.} 
}

\details{
So far only SIMPLS method [1] is available, more coming soon. Implementation works both with one and multiple response variables.

Like in \code{\link{pca}}, \code{pls} uses number of components (\code{ncomp}) as a minimum of number of objects - 1, number of x variables and the default or provided value. Regression coefficients, predictions and other results are calculated for each set of components from 1 to \code{ncomp}: 1, 1:2, 1:3, etc. Besides that, there is also a  function (\code{\link{selectCompNum.pls}}) for selecting an optimal number of components in a model (\code{ncomp.selected}). The selected optimal number of components is used for all default operations — predictions, plots, etc. 

}

\references{
1. S. de Jong, SIMPLS: an alternative approach to partial least squares regression, Chemometrics and Intelligent Laboratory Systems 18 (1993) 251–263.
}

\seealso{
Methods for \code{pls} objects:
  \tabular{ll}{
      \code{print} \tab prints information about a \code{pls} object.\cr
      \code{\link{summary.pls}} \tab shows performance statistics for the model.\cr
      \code{\link{plot.pls}} \tab shows plot overview of the model.\cr
      \code{\link{pls.simpls}} \tab implementation of SIMPLS algorithm.\cr
      \code{\link{predict.pls}} \tab applies PLS model to a new data.\cr
      \code{\link{selectCompNum.pls}} \tab set number of optimal components in the model.\cr
      \code{\link{plotPredictions.pls}} \tab shows predicted vs. measured plot.\cr
      \code{\link{plotRegcoeffs.pls}} \tab shows regression coefficients plot.\cr      
      \code{\link{plotXScores.pls}} \tab shows scores plot for x decomposition.\cr
      \code{\link{plotXYScores.pls}} \tab shows scores plot for x and y decomposition.\cr
      \code{\link{plotXLoadings.pls}} \tab shows loadings plot for x decomposition.\cr
      \code{\link{plotXYLoadings.pls}} \tab shows loadings plot for x and y decomposition.\cr
      \code{\link{plotRMSE.pls}} \tab shows RMSE plot.\cr
      \code{\link{plotXVariance.pls}} \tab shows explained variance plot for x decomposition.\cr
      \code{\link{plotYVariance.pls}} \tab shows explained variance plot for y decomposition.\cr
      \code{\link{plotXCumVariance.pls}} \tab shows cumulative explained variance plot for y decomposition.\cr
      \code{\link{plotYCumVariance.pls}} \tab shows cumulative explained variance plot for y decomposition.\cr
      \code{\link{plotXResiduals.pls}} \tab shows T2 vs. Q2 plot for x decomposition.\cr
      \code{\link{plotYResiduals.pls}} \tab shows residuals plot for y values.\cr
   }
   
   Most of the methods for plotting data (except loadings and regression coefficients) are also available for PLS results 
   (\code{\link{plsres}}) objects.    
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\examples{
## Make a PLS model with full cross-validation and show model overview

data(Simdata)
model = pls(spectra.c, conc.c[, 1], ncomp = 8, cv = 1)
model = selectCompNum(model, 3)
summary(model)
plot(model)

}

