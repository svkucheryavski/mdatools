\name{plsdares}
\alias{plsdares}

\title{
PLS-DA results
}

\description{
\code{plsdares} is used to store and visualize results of applying a PLS-DA model to a new data.
}

\usage{
plsdares(plsres, cres)
}

\arguments{
  \item{plsres}{PLS results for the data. }
  \item{cres}{Classification results for the data. }
}

\details{
Do not use \code{plsdares} manually, the object is created automatically when one applies a PLS-DA model to a new data set, 
e.g. when calibrate and validate a PLS-DA model (all calibration and validation results in PLS-DA model are stored as objects 
of \code{plsdares} class) or use function \code{\link{predict.plsda}}.

The object gives access to all PLS-DA results as well as to the plotting methods for visualisation of the results. 
The \code{plsidares} class also inherits all properties and methods of \code{classres} and \code{plsres} classes.

If no reference values provided, classification statistics will not be calculated and performance plots will not be available. 
}

\value{
Returns an object of \code{plsdares} class with fields, inherited from \code{\link{classres}} and \code{\link{plsres}}.
}

\author{
Sergey Kucheryavskiy (svkucheryavski@gmail.com)
}

\seealso{
Methods for \code{plsda} objects:
  \tabular{ll}{
      \code{print.plsda} \tab shows information about the object.\cr
      \code{summary.plsda} \tab shows statistics for results of classification.\cr
      \code{plot.plsda} \tab shows plots for overview of the results.\cr
   }

Methods, inherited from \code{\link{classres}} class:
  \tabular{ll}{
      \code{\link{showPredictions.classres}} \tab show table with predicted values.\cr
      \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
      \code{\link{plotSensitivity.classres}} \tab makes plot with sensitivity vs. components values.\cr
      \code{\link{plotSpecificity.classres}} \tab makes plot with specificity vs. components values.\cr
      \code{\link{plotPerformance.classres}} \tab makes plot with both specificity and sensitivity values.\cr
   }

Methods, inherited from \code{plsres} class:
  \tabular{ll}{
      \code{\link{plotPredictions.plsres}} \tab shows predicted vs. measured plot.\cr
      \code{\link{plotXScores.plsres}} \tab shows scores plot for x decomposition.\cr
      \code{\link{plotXYScores.plsres}} \tab shows scores plot for x and y decomposition.\cr
      \code{\link{plotRMSE.plsres}} \tab shows RMSE plot.\cr
      \code{\link{plotXVariance.plsres}} \tab shows explained variance plot for x decomposition.\cr
      \code{\link{plotYVariance.plsres}} \tab shows explained variance plot for y decomposition.\cr
      \code{\link{plotXCumVariance.plsres}} \tab shows cumulative explained variance plot for y decomposition.\cr
      \code{\link{plotYCumVariance.plsres}} \tab shows cumulative explained variance plot for y decomposition.\cr
      \code{\link{plotXResiduals.plsres}} \tab shows T2 vs. Q2 plot for x decomposition.\cr
      \code{\link{plotYResiduals.regres}} \tab shows residuals plot for y values.\cr
   }

See also \code{\link{plsda}} - a class for PLS-DA models, \code{\link{predict.plsda}} applying
PLS-DA model for a new dataset.
}

\examples{
### Examples for PLS-DA results class

library(mdatools)

## 1. Make a PLS-DA model with full cross-validation, get
## calibration results and show overview

# make a calibration set from iris data (3 classes)
# use names of classes as class vector
x.cal = iris[seq(1, nrow(iris), 2), 1:4] 
c.cal = iris[seq(1, nrow(iris), 2), 5]

model = plsda(x.cal, c.cal, ncomp = 3, cv = 1, info = 'IRIS data example')
model = selectCompNum(model, 1)

res = model$calres

# show summary and basic plots for calibration results
summary(res)
plot(res)

## 2. Apply the calibrated PLS-DA model to a new dataset

# make a new data
x.new = iris[seq(2, nrow(iris), 2), 1:4] 
c.new = iris[seq(2, nrow(iris), 2), 5]

res = predict(model, x.new, c.new)
summary(res)
plot(res)

## 3. Show performance plots for the results
par(mfrow = c(2, 2))
plotSpecificity(res)
plotSensitivity(res)
plotMisclassified(res)
plotMisclassified(res, nc = 2)
par(mfrow = c(1, 1))

## 3. Show both class and y values predictions
par(mfrow = c(2, 2))
plotPredictions(res)
plotPredictions(res, ncomp = 2, nc = 2)
plotPredictions.plsres(res)
plotPredictions.plsres(res, ncomp = 2, ny = 2)
par(mfrow = c(1, 1))

## 4. All plots from ordinary PLS results can be used, e.g.:
par(mfrow = c(2, 2))
plotXYScores(res)
plotYVariance(res, type = 'h')
plotXVariance(res, type = 'h')
plotXResiduals(res)
par(mfrow = c(1, 1))

}


