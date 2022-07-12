#' PLS-DA results
#' @description
#' \code{plsdares} is used to store and visualize results of applying a PLS-DA model to a new data.
#'
#' @param plsres
#' PLS results for the data.
#' @param cres
#' Classification results for the data.
#'
#' @details
#' Do not use \code{plsdares} manually, the object is created automatically when one applies a
#' PLS-DA model to a new data set, e.g. when calibrate and validate a PLS-DA model (all calibration
#' and validation results in PLS-DA model are stored as objects of \code{plsdares} class) or use
#' function \code{\link{predict.plsda}}.
#'
#' The object gives access to all PLS-DA results as well as to the plotting methods for
#' visualisation of the results. The \code{plsidares} class also inherits all properties and methods
#' of \code{classres} and \code{plsres} classes.
#'
#' If no reference values provided, classification statistics will not be calculated and
#' performance plots will not be available.
#'
#' @return
#' Returns an object of \code{plsdares} class with fields, inherited from \code{\link{classres}}
#' and \code{\link{plsres}}.
#'
#' @seealso
#' Methods for \code{plsda} objects:
#' \tabular{ll}{
#'  \code{print.plsda} \tab shows information about the object.\cr
#'  \code{summary.plsda} \tab shows statistics for results of classification.\cr
#'  \code{plot.plsda} \tab shows plots for overview of the results.\cr
#' }
#'
#' Methods, inherited from \code{\link{classres}} class:
#' \tabular{ll}{
#'  \code{\link{showPredictions.classres}} \tab show table with predicted values.\cr
#'  \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
#'  \code{\link{plotSensitivity.classres}} \tab makes plot with sensitivity vs. components
#' values.\cr
#'  \code{\link{plotSpecificity.classres}} \tab makes plot with specificity vs. components
#' values.\cr
#'  \code{\link{plotPerformance.classres}} \tab makes plot with both specificity and sensitivity
#'   values.\cr
#' }
#'
#' Methods for \code{plsres} objects:
#' \tabular{ll}{
#'    \code{print} \tab prints information about a \code{plsres} object.\cr
#'    \code{\link{summary.plsres}} \tab shows performance statistics for the results.\cr
#'    \code{\link{plot.plsres}} \tab shows plot overview of the results.\cr
#'    \code{\link{plotXScores.plsres}} \tab shows scores plot for x decomposition.\cr
#'    \code{\link{plotXYScores.plsres}} \tab shows scores plot for x and y decomposition.\cr
#'    \code{\link{plotXVariance.plsres}} \tab shows explained variance plot for x decomposition.\cr
#'    \code{\link{plotYVariance.plsres}} \tab shows explained variance plot for y decomposition.\cr
#'    \code{\link{plotXCumVariance.plsres}} \tab shows cumulative explained variance plot for y
#'    decomposition.\cr
#'    \code{\link{plotYCumVariance.plsres}} \tab shows cumulative explained variance plot for y
#'    decomposition.\cr
#'    \code{\link{plotXResiduals.plsres}} \tab shows T2 vs. Q plot for x decomposition.\cr
#'    \code{\link{plotYResiduals.plsres}} \tab shows residuals plot for y values.\cr
#' }
#'
#' Methods inherited from \code{regres} class (parent class for \code{plsres}):
#' \tabular{ll}{
#'    \code{\link{plotPredictions.regres}} \tab shows predicted vs. measured plot.\cr
#'    \code{\link{plotRMSE.regres}} \tab shows RMSE plot.\cr
#' }
#'
#' See also \code{\link{plsda}} - a class for PLS-DA models, \code{\link{predict.plsda}} applying
#' PLS-DA model for a new dataset.
#'
#' @examples
#' ### Examples for PLS-DA results class
#'
#' library(mdatools)
#'
#' ## 1. Make a PLS-DA model with full cross-validation, get
#' ## calibration results and show overview
#'
#' # make a calibration set from iris data (3 classes)
#' # use names of classes as class vector
#' x.cal = iris[seq(1, nrow(iris), 2), 1:4]
#' c.cal = iris[seq(1, nrow(iris), 2), 5]
#'
#' model = plsda(x.cal, c.cal, ncomp = 3, cv = 1, info = 'IRIS data example')
#' model = selectCompNum(model, 1)
#'
#' res = model$calres
#'
#' # show summary and basic plots for calibration results
#' summary(res)
#' plot(res)
#'
#' ## 2. Apply the calibrated PLS-DA model to a new dataset
#'
#' # make a new data
#' x.new = iris[seq(2, nrow(iris), 2), 1:4]
#' c.new = iris[seq(2, nrow(iris), 2), 5]
#'
#' res = predict(model, x.new, c.new)
#' summary(res)
#' plot(res)
#'
#' ## 3. Show performance plots for the results
#' par(mfrow = c(2, 2))
#' plotSpecificity(res)
#' plotSensitivity(res)
#' plotMisclassified(res)
#' plotMisclassified(res, nc = 2)
#' par(mfrow = c(1, 1))
#'
#' ## 3. Show both class and y values predictions
#' par(mfrow = c(2, 2))
#' plotPredictions(res)
#' plotPredictions(res, ncomp = 2, nc = 2)
#' plotPredictions(structure(res, class = "regres"))
#' plotPredictions(structure(res, class = "regres"), ncomp = 2, ny = 2)
#' par(mfrow = c(1, 1))
#'
#' ## 4. All plots from ordinary PLS results can be used, e.g.:
#' par(mfrow = c(2, 2))
#' plotXYScores(res)
#' plotYVariance(res, type = 'h')
#' plotXVariance(res, type = 'h')
#' plotXResiduals(res)
#' par(mfrow = c(1, 1))
#'
#' @export
plsdares <- function(plsres, cres) {
   obj <- c(plsres, cres)
   class(obj) <- c("plsdares", "classres", "plsres", "regres")
   obj$call <- match.call()

   return(obj)
}

#' Overview plot for PLS-DA results
#'
#' @description
#' Shows a set of plots (x residuals, y variance, classification performance and predictions)
#' for PLS-DA results.
#'
#' @param x
#' PLS-DA results (object of class \code{plsdares})
#' @param ncomp
#' how many components to use
#' @param nc
#' which class to show the plot for
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plot.plsdares <- function(x, nc = 1, ncomp = x$ncomp.selected, show.labels = FALSE, ...) {

   if (is.null(x$c.ref)) {
      par(mfrow = c(2, 2))
      plotXResiduals.plsres(x, ncomp = ncomp, show.labels = show.labels)
      plotYVariance.plsres(x, ncomp = ncomp, show.labels = show.labels)
      plotPredictions.classres(x, ncomp = ncomp, show.labels = show.labels)
      par(mfrow = c(1, 1))
      return()
   }

   par(mfrow = c(2, 2))
   plotXResiduals.plsres(x, ncomp = ncomp, show.labels = show.labels)
   plotYVariance.plsres(x, ncomp = ncomp, show.labels = show.labels)
   plotPerformance.classres(x, nc = nc, show.labels = show.labels)
   plotPredictions.classres(x, ncomp = ncomp, show.labels = show.labels)
   par(mfrow = c(1, 1))
}

#' as.matrix method for PLS-DA results
#'
#' @description
#' Returns a matrix with model performance statistics for PLS-DA results
#'
#' @param x
#' PLS-DA results (object of class \code{plsdares})
#' @param ncomp
#' number of components to calculate the statistics for (if NULL gets for all components)
#' @param nc
#' for which class to calculate the statistics for
#' @param ...
#' other arguments
#'
#' @export
as.matrix.plsdares <- function(x, ncomp = NULL, nc = 1, ...) {
   return(
      cbind(
         as.matrix.plsres(x, ncomp = ncomp, ny = nc)[, 1:4, drop = FALSE],
         as.matrix.classres(x, ncomp = ncomp, nc = nc)
      )
   )
}

#' Summary method for PLS-DA results object
#'
#' @description
#' Shows performance statistics for the results.
#'
#' @param object
#' PLS-DA results (object of class \code{plsdares})
#' @param nc
#' which class to show the summary for (if NULL, will be shown for all)
#' @param ...
#' other arguments
#'
#' @export
summary.plsdares <- function(object, nc = seq_len(object$nclasses), ...) {
   cat("\nPLS-DA results (class plsdares) summary:\n")
   fprintf("Number of selected components: %.0f\n", object$ncomp.selected)

   if (is.null(object$c.ref)) {
      cat("No reference data available.\n")
      return()
   }

   for (n in nc) {
      fprintf("\nClass #%.0f (%s):\n", n, object$classnames[n])
      out <- as.matrix(object, nc = n)
      if (!any(is.na(out[, 1:4]))) out[, 1:4] <- round(out[, 1:4], 3)
      print(out)
      cat("\n")
   }
}

#' Print method for PLS-DA results object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' PLS-DA results (object of class \code{plsdares})
#' @param ...
#' other arguments
#'
#' @export
print.plsdares <- function(x, ...) {

   cat("\nPLS-DA results (class plsdares)\n")
   cat("\nCall:\n")
   print(x$call)

   cat("\nMajor fields:\n")
   cat("$ncomp.selected - number of selected components\n")
   cat("$c.pred - array with predicted class values\n")

   if (!is.null(x$c.ref)) {
      cat("$c.ref - vector with reference class values\n")
      cat("$tp - number of true positives\n")
      cat("$fp - number of false positives\n")
      cat("$fn - number of false negatives\n")
      cat("$specificity - specificity of predictions\n")
      cat("$sensitivity - sensitivity of predictions\n")
      cat("$misclassified - misclassification ratio for predictions\n")
      cat("$ydecomp - decomposition of y values (ldecomp object)\n")
   }

   cat("$xdecomp - decomposition of x values (ldecomp object)\n")
}
