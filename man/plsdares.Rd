\name{plsdares}
\alias{plsdares}

\title{
PLS-DA results
}

\description{
\code{plsdares} is used to store and visualize results of applying a PLS-DA model to a new data.
}

\usage{
plsdares()
}

\arguments{
  \item{pres}{PLS results for the data. }
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
      \code{\link{plotYResiduals.plsres}} \tab shows residuals plot for y values.\cr
   }

See also \code{\link{plsda}} - a class for PLS-DA models. Examples can be found in help for \code{\link{predict.plsda}} function.
}

