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
#' If no reference values provided, classification statistics will not be calculated and performance 
#' plots will not be available. 
#'
#' @return 
#' Returns an object of \code{plsdares} class with fields, inherited from \code{\link{classres}} and
#' \code{\link{plsres}}.
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
#'  \code{\link{plotSensitivity.classres}} \tab makes plot with sensitivity vs. components values.\cr
#'  \code{\link{plotSpecificity.classres}} \tab makes plot with specificity vs. components values.\cr
#'  \code{\link{plotPerformance.classres}} \tab makes plot with both specificity and sensitivity
#'   values.\cr
#' }
#' 
#' Methods, inherited from \code{plsres} class:
#' \tabular{ll}{
#'  \code{\link{plotPredictions.plsres}} \tab shows predicted vs. measured plot.\cr
#'  \code{\link{plotXScores.plsres}} \tab shows scores plot for x decomposition.\cr
#'  \code{\link{plotXYScores.plsres}} \tab shows scores plot for x and y decomposition.\cr
#'  \code{\link{plotRMSE.plsres}} \tab shows RMSE plot.\cr
#'  \code{\link{plotXVariance.plsres}} \tab shows explained variance plot for x decomposition.\cr
#'  \code{\link{plotYVariance.plsres}} \tab shows explained variance plot for y decomposition.\cr
#'  \code{\link{plotXCumVariance.plsres}} \tab shows cumulative explained variance plot for y 
#'  decomposition.\cr
#'  \code{\link{plotYCumVariance.plsres}} \tab shows cumulative explained variance plot for y 
#'  decomposition.\cr
#'  \code{\link{plotXResiduals.plsres}} \tab shows T2 vs. Q plot for x decomposition.\cr
#'  \code{\link{plotYResiduals.regres}} \tab shows residuals plot for y values.\cr
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
#' plotPredictions(structure(res, class = "plsres"))
#' plotPredictions(structure(res, class = "plsres"), ncomp = 2, ny = 2)
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
plsdares = function(plsres, cres)
{
   obj = c(plsres, cres)
   class(obj) = c('plsdares', 'classres', 'plsres')   
   
   obj$call = match.call()   
   
   obj
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
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param nc
#' which class to show the summary for (if NULL, will be shown for all)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.line
#' logical, show or not target line on predictions plot
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export 
plot.plsdares = function(x, nc = NULL, ncomp = NULL, show.labels = F, show.line = T, ...)
{
   obj = x
   
   par(mfrow = c(2, 2))
   plotXResiduals.plsres(obj, ncomp = ncomp, show.labels = show.labels)
   plotYVariance.plsres(obj, ncomp = ncomp, show.labels = show.labels)
   plotPerformance(obj, nc = nc, ncomp = ncomp, show.labels = show.labels)
   plotPredictions(obj, show.line = show.line, nc = nc, ncomp = ncomp, show.labels = show.labels)
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
#' number of components to calculate the statistics for
#' @param nc
#' for which class to calculate the statistics for
#' @param ...
#' other arguments
#' 
#' @export
as.matrix.plsdares = function(x, ncomp = NULL, nc = NULL, ...)
{
   obj = x
   
   plsmat = as.matrix.plsres(obj, ncomp = ncomp, ny = nc)
   classmat = as.matrix.classres(obj, ncomp = ncomp, nc = nc)
   mat = cbind(plsmat[, 1:4, drop = F], classmat)
   
   mat
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
summary.plsdares = function(object, nc = NULL, ...)
{
   obj = object
   
   if (is.null(nc))
      nc = 1:obj$nclasses
   cat('\nPLS-DA results (class plsdares) summary:\n');
   cat(sprintf('Number of selected components: %d\n', obj$ncomp.selected))
   for (n in nc)
   {
      cat(sprintf('\nClass #%d (%s)\n', n, obj$classnames[n]))
      
      mat = as.matrix(obj, nc = n)
      mat[, 1:4] = round(mat[, 1:4], 2)
      print(mat)
      cat('\n')
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
print.plsdares = function(x, ...)
{
   obj = x
   
   cat('\nPLS-DA results (class plsdares)\n')
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$ncomp.selected - number of selected components\n')
   cat('$c.pred - array with predicted class values\n')
   if (!is.null(obj$c.ref))
   {   
      cat('$c.ref - vector with reference class values\n')
      cat('$tp - number of true positives\n')
      cat('$fp - number of false positives\n')
      cat('$fn - number of false negatives\n')
      cat('$specificity - specificity of predictions\n')
      cat('$sensitivity - sensitivity of predictions\n')
      cat('$misclassified - misclassification ratio for predictions\n')
      cat('$ydecomp - decomposition of y values (ldecomp object)\n')
   }
   cat('$xdecomp - decomposition of x values (ldecomp object)\n')
}
