\name{plsres}
\alias{plsres}

\title{
PLS results
}

\description{
\code{plsres} is used to store and visualize results of applying a PLS model to a new data.
}

\usage{
plsres(y.pred, y.ref = NULL, ncomp.selected = NULL, xdecomp = NULL, ydecomp = NULL,
   info = '')
}

\arguments{
  \item{y.pred}{predicted y values. }
  \item{y.ref}{reference (measured) y values. }
  \item{ncomp.selected}{selected (optimal) number of components. }
  \item{xdecomp}{PLS decomposition of X data (object of class \code{ldecomp}).}
  \item{ydecomp}{PLS decomposition of Y data (object of class \code{ldecomp}).}
  \item{info}{information about the object. }
}

\details{
Do not use \code{plsres} manually, the object is created automatically when one applies a PLS model to a new data set, e.g. when calibrate and validate a PLS model (all calibration and validation results in PLS model are stored as objects of \code{plsres} class) or use function \code{\link{predict.pls}}.

The object gives access to all PLS results as well as to the plotting methods for visualisation of the results. The \code{plsres} class also inherits all properties and methods of \code{regres} - general class for regression results.  

If no reference values provided, regression statistics will not be calculated and most of the plots not available. The class is also used for cross-validation results, in this case some of the values and methods are not available (e.g. scores and scores plot, etc.).

All plots are based on \code{\link{mdaplot}} function, so most of its options can be used (e.g. color grouping, etc.).

RPD is ratio of standard deviation of response values to standard error of prediction (SDy/SEP).
}

\value{
Returns an object of \code{plsres} class with following fields:
   \item{ncomp }{number of components included to the model.} 
   \item{ncomp.selected }{selected (optimal) number of components.} 
   \item{y.ref }{a matrix with reference values for responses.} 
   \item{y.pred }{a matrix with predicted values for responses.} 
   \item{rmse }{a matrix with root mean squared error values for each response and component.} 
   \item{slope }{a matrix with slope values for each response and component.} 
   \item{r2 }{a matrix with determination coefficients for each response and component.} 
   \item{bias }{a matrix with bias values for each response and component.} 
   \item{sep }{a matrix with standard error values for each response and component.} 
   \item{rpd }{a matrix with RPD values for each response and component.} 
   \item{xdecomp }{decomposition of predictors (object of class \code{ldecomp}).} 
   \item{ydecomp }{decomposition of responses (object of class \code{ldecomp}).} 
   \item{info }{information about the object.} 
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\seealso{
Methods for \code{plsres} objects:
  \tabular{ll}{
      \code{print} \tab prints information about a \code{plsres} object.\cr
      \code{\link{summary.plsres}} \tab shows performance statistics for the results.\cr
      \code{\link{plot.plsres}} \tab shows plot overview of the results.\cr
      \code{\link{plotPredictions.plsres}} \tab shows predicted vs. measured plot.\cr
      \code{\link{plotXScores.plsres}} \tab shows scores plot for x decomposition.\cr
      \code{\link{plotXYScores.plsres}} \tab shows scores plot for x and y decomposition.\cr
      \code{\link{plotRMSE.plsres}} \tab shows RMSE plot.\cr
      \code{\link{plotXVariance.plsres}} \tab shows explained variance plot for x decomposition.\cr
      \code{\link{plotYVariance.plsres}} \tab shows explained variance plot for y decomposition.\cr
      \code{\link{plotXCumVariance.plsres}} \tab shows cumulative explained variance plot for y decomposition.\cr
      \code{\link{plotYCumVariance.plsres}} \tab shows cumulative explained variance plot for y decomposition.\cr
      \code{\link{plotXResiduals.plsres}} \tab shows T2 vs. Q2 plot for x decomposition.\cr
      \code{\link{plotYResiduals.plsres}} \tab shows residuals plot for y values.\cr
   }
See also \code{\link{pls}} - a class for PLS models.
}

\examples{
## Examples of using PLS result methods

library(mdatools)
data(Simdata)
model = pls(spectra.c, conc.c, ncomp = 5, cv = 4)
model = selectCompNum(model, 3)

# get and explore cross-validated results from the model
calres = model$calres
plot(calres)
summary(calres)
summary(calres, ncomp = 1:5, ny = 2)

# apply model to a new data and show the results
newres = predict(model, spectra.t, conc.t)
plot(newres)
summary(newres)
summary(newres, ncomp = 1:5, ny = 2)


}
