#' PLS results
#'
#' @description 
#' \code{plsres} is used to store and visualize results of applying a PLS model to a new data.
#'
#' @param y.pred
#' predicted y values.
#' @param y.ref
#' reference (measured) y values.
#' @param ncomp.selected
#' selected (optimal) number of components.
#' @param xdecomp
#' PLS decomposition of X data (object of class \code{ldecomp}).
#' @param ydecomp
#' PLS decomposition of Y data (object of class \code{ldecomp}).
#' @param info
#' information about the object.
#'
#' @details 
#' Do not use \code{plsres} manually, the object is created automatically when one applies a PLS 
#' model to a new data set, e.g. when calibrate and validate a PLS model (all calibration and 
#' validation results in PLS model are stored as objects of \code{plsres} class) or use function 
#' \code{\link{predict.pls}}.
#' 
#' The object gives access to all PLS results as well as to the plotting methods for visualisation 
#' of the results. The \code{plsres} class also inherits all properties and methods of \code{regres}
#'  - general class for regression results.  
#' 
#' If no reference values provided, regression statistics will not be calculated and most of the 
#' plots not available. The class is also used for cross-validation results, in this case some of 
#' the values and methods are not available (e.g. scores and scores plot, etc.).
#' 
#' All plots are based on \code{\link{mdaplot}} function, so most of its options can be used (e.g. 
#' color grouping, etc.).
#' 
#' RPD is ratio of standard deviation of response values to standard error of prediction (SDy/SEP).
#'
#' @return 
#' Returns an object of \code{plsres} class with following fields:
#' \item{ncomp }{number of components included to the model.} 
#' \item{ncomp.selected }{selected (optimal) number of components.} 
#' \item{y.ref }{a matrix with reference values for responses.} 
#' \item{y.pred }{a matrix with predicted values for responses.} 
#' \item{rmse }{a matrix with root mean squared error values for each response and component.} 
#' \item{slope }{a matrix with slope values for each response and component.} 
#' \item{r2 }{a matrix with determination coefficients for each response and component.} 
#' \item{bias }{a matrix with bias values for each response and component.} 
#' \item{sep }{a matrix with standard error values for each response and component.} 
#' \item{rpd }{a matrix with RPD values for each response and component.} 
#' \item{xdecomp }{decomposition of predictors (object of class \code{ldecomp}).} 
#' \item{ydecomp }{decomposition of responses (object of class \code{ldecomp}).} 
#' \item{info }{information about the object.} 
#'
#' @seealso      
#' Methods for \code{plsres} objects:
#' \tabular{ll}{
#'    \code{print} \tab prints information about a \code{plsres} object.\cr
#'    \code{\link{summary.plsres}} \tab shows performance statistics for the results.\cr
#'    \code{\link{plot.plsres}} \tab shows plot overview of the results.\cr
#'    \code{\link{plotPredictions.plsres}} \tab shows predicted vs. measured plot.\cr
#'    \code{\link{plotXScores.plsres}} \tab shows scores plot for x decomposition.\cr
#'    \code{\link{plotXYScores.plsres}} \tab shows scores plot for x and y decomposition.\cr
#'    \code{\link{plotRMSE.regres}} \tab shows RMSE plot.\cr
#'    \code{\link{plotXVariance.plsres}} \tab shows explained variance plot for x decomposition.\cr
#'    \code{\link{plotYVariance.plsres}} \tab shows explained variance plot for y decomposition.\cr
#'    \code{\link{plotXCumVariance.plsres}} \tab shows cumulative explained variance plot for y 
#'    decomposition.\cr
#'    \code{\link{plotYCumVariance.plsres}} \tab shows cumulative explained variance plot for y 
#'    decomposition.\cr
#'    \code{\link{plotXResiduals.plsres}} \tab shows T2 vs. Q plot for x decomposition.\cr
#'    \code{\link{plotYResiduals.regres}} \tab shows residuals plot for y values.\cr
#' }
#' See also \code{\link{pls}} - a class for PLS models.
#' 
#' @examples 
#' ### Examples of using PLS result class
#' library(mdatools)

#' ## 1. Make a PLS model for concentration of first component 
#' ## using full-cross validation and get calibration results
#' 
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#' 
#' model = pls(x, y, ncomp = 8, cv = 1)
#' model = selectCompNum(model, 2)
#' res = model$calres
#' 
#' summary(res)
#' plot(res)
#' 
#' ## 2. Make a PLS model for concentration of first component 
#' ## and apply model to a new dataset
#' 
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#' 
#' model = pls(x, y, ncomp = 6, cv = 1)
#' model = selectCompNum(model, 2)
#' 
#' x.new = simdata$spectra.t
#' y.new = simdata$conc.t[, 1]
#' res = predict(model, x.new, y.new)
#' 
#' summary(res)
#' plot(res)
#' 
#' ## 3. Show variance and error plots for PLS results
#' par(mfrow = c(2, 2))
#' plotXCumVariance(res, type = 'h')
#' plotYCumVariance(res, type = 'b', show.labels = TRUE, legend.position = 'bottomright')
#' plotRMSE(res)
#' plotRMSE(res, type = 'h', show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' ## 4. Show scores plots for PLS results
#' ## (for results plot we can use color grouping)
#' par(mfrow = c(2, 2))
#' plotXScores(res)
#' plotXScores(res, show.labels = TRUE, cgroup = y.new)
#' plotXYScores(res)
#' plotXYScores(res, comp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' ## 5. Show predictions and residuals plots for PLS results
#' par(mfrow = c(2, 2))
#' plotXResiduals(res, show.label = TRUE, cgroup = y.new)
#' plotYResiduals(res, show.label = TRUE)
#' plotPredictions(res)
#' plotPredictions(res, ncomp = 4, xlab = 'C, reference', ylab = 'C, predictions')
#' par(mfrow = c(1, 1))
#'
#' @export
plsres = function(y.pred, y.ref = NULL, ncomp.selected = NULL, xdecomp = NULL, ydecomp = NULL, info = '') {
   obj = regres(y.pred, y.ref = y.ref, ncomp.selected = ncomp.selected)
   obj$ncomp = ncol(y.pred)
   obj$xdecomp = xdecomp
   obj$ydecomp = ydecomp
   obj$info = info
   
   if (is.null(ncomp.selected))
      obj$ncomp.selected = obj$ncomp
   else
      obj$ncomp.selected = ncomp.selected
   
   obj$call = match.call()   
   class(obj) = c("plsres", "regres")
   
   obj
}

#' X scores plot for PLS results
#' 
#' @description
#' Shows plot with scores values for PLS decomposition of x data.
#' 
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
#' @export
plotXScores.plsres = function(obj, comp = c(1, 2), main = 'X scores', ...) {
   if (!is.null(obj$xdecomp$scores))
      plotScores.ldecomp(obj$xdecomp, comp = comp, main = main, ...)
   else
      warning('Scores values are not available.')
}

#' XY scores plot for PLS results
#' 
#' @description
#' Shows plot with X vs. Y scores values for PLS results.
#' 
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' which component to show the plot for
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
#' @export
plotXYScores.plsres = function(obj, ncomp = 1, ...) {
   
   if (is.null(obj$xdecomp$scores) || is.null(obj$ydecomp$scores)) {   
      warning('X or Y scores are not available.')
      return()
   }

   if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp)
      stop('Wrong value for ncomp argument!')
   
      
   data = cbind(obj$xdecomp$scores[, ncomp, drop = F], obj$ydecomp$scores[, ncomp, drop = F])
   data = mda.setattr(data, mda.getattr(obj$xdecomp$scores))
   attr(data, 'name') = sprintf('XY scores')
   rownames(data) = rownames(obj$xdecomp$scores)
   colnames(data) = c(
      sprintf('X scores (Comp %d, %.2f%%)', ncomp, obj$xdecomp$expvar[ncomp]), 
      sprintf('Y scores (Comp %d, %.2f%%)', ncomp, obj$ydecomp$expvar[ncomp])
   )

   mdaplot(data, type = 'p', ...)
}

#' X residuals plot for PLS results
#' 
#' @description
#' Shows a plot with Q residuals vs. Hotelling T2 values for PLS decomposition of x data.
#' 
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param main
#' main title for the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
#' @export
plotXResiduals.plsres = function(obj, ncomp = NULL, main = NULL, ...) {   
   if (is.null(main)) {   
      if (is.null(ncomp))
         main = 'X residuals'
      else   
         main = sprintf('X residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected

   if (length(ncomp) != 1 || ncomp < 0 || ncomp > obj$ncomp)
      stop('Wrong value for ncomp argument!')
   
   plotResiduals.ldecomp(obj$xdecomp, ncomp = ncomp, main = main, ...)
}

#' Explained X variance plot for PLS results
#' 
#' @description
#' Shows plot with explained X variance vs. number of components.
#' 
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
#' @export
plotXVariance.plsres = function(obj, main = 'X variance', ...) {   
   plotVariance.ldecomp(obj$xdecomp, main = main, ...)
}

#' Explained Y variance plot for PLS results
#' 
#' @description
#' Shows plot with explained Y variance vs. number of components.
#' 
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
#' @export
plotYVariance.plsres = function(obj, main = 'Y variance', ...) {   
   plotVariance.ldecomp(obj$ydecomp, main = main, ...)
}

#' Explained cumulative X variance plot for PLS results
#' 
#' @description
#' Shows plot with cumulative explained X variance vs. number of components.
#' 
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
#' @export
plotXCumVariance.plsres = function(obj, main = 'X cumulative variance', ...) {   
   plotCumVariance.ldecomp(obj$xdecomp, main = main, ...)
}

#' Explained cumulative Y variance plot for PLS results
#' 
#' @description
#' Shows plot with cumulative explained Y variance vs. number of components.
#' 
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
#' @export
plotYCumVariance.plsres = function(obj, main = 'Y cumulative variance', ...) {   
   plotCumVariance.ldecomp(obj$ydecomp, main = main, ...)
}

#' Predictions plot for PLS results
#' 
#' @description
#' Shows plot with predicted vs. reference (measured) y values for selected components.
#' 
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' number of response variable to make the plot for (if y is multivariate)
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
#' @seealso
#' \code{\link{plotPredictions.regres}} - prediction plot for regression results.
#' 
#' @export
plotPredictions.plsres = function(obj, ny = 1, ncomp = NULL, main = NULL, ...) {
   if (is.null(main)) {
      if (is.null(ncomp))
         main = 'Predictions'
      else
         main = sprintf('Predictions (ncomp = %d)', ncomp)
   }   
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   else if (ncomp < 0 || ncomp > obj$ncomp)
      stop('Wrong number of components!')
   
   plotPredictions.regres(obj, ny = ny, ncomp = ncomp, main = main, ...)
}

#' Overview plot for PLS results
#' 
#' @description
#' Shows a set of plots for PLS results.
#' 
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' which y variable to show the summary for (if NULL, will be shown for all)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
#' @export
plot.plsres = function(x, ncomp = NULL, ny = 1, show.labels = F, ...) {
   obj = x
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected 
   
   if (is.null(obj$y.ref)) {
      par(mfrow = c(1, 2))
      plotXResiduals(obj, ...)
      plotPredictions.plsres(obj, ncomp = ncomp, ny = ny, ...)
      par(mfrow = c(1, 1))      
   } else {   
      par(mfrow = c(2, 2))
      plotXResiduals(obj, ncomp = ncomp, ...) 
      plotYVariance(obj, ...)   
      plotRMSE(obj, ny = ny, ...)
      plotPredictions.plsres(obj, ncomp = ncomp, ny = ny, ...)
      par(mfrow = c(1, 1))
   }
}   

#' as.matrix method for PLS results
#' 
#' @description
#' Returns a matrix with model performance statistics for PLS results
#' 
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' number of components to calculate the statistics for
#' @param ny
#' for which response variable calculate the statistics for
#' @param ...
#' other arguments
#' 
#' @export
as.matrix.plsres = function(x, ncomp = NULL, ny = 1, ...) {
   obj = x
   if (is.null(ncomp)) {   
      res = cbind(
         obj$xdecomp$expvar,
         obj$xdecomp$cumexpvar,
         obj$ydecomp$expvar,
         obj$ydecomp$cumexpvar,
         as.matrix.regres(obj, ny = ny)
      )
      rownames(res) = colnames(obj$xdecomp$Q)
   } else {
      res = cbind(
         obj$xdecomp$expvar[ncomp],
         obj$xdecomp$cumexpvar[ncomp],
         obj$ydecomp$expvar[ncomp],
         obj$ydecomp$cumexpvar[ncomp],
         as.matrix.regres(obj, ncomp = ncomp, ny = ny)
      )
      rownames(res) = colnames(obj$y.pred)[ncomp]
   }
   colnames(res)[1:4] = c('X expvar', 'X cumexpvar', 'Y expvar', 'Y cumexpvar')
   res = res[, -6, drop = FALSE] 
   res
}


#' summary method for PLS results object
#' 
#' @description
#' Shows performance statistics for the results.
#' 
#' @param object
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' for which response variable show the summary for
#' @param ...
#' other arguments
#' 
#' @export
summary.plsres = function(object, ny = NULL, ncomp = NULL, ...) {
   obj = object
   
   cat('\nPLS regression results (class plsres) summary\n')
   if (!is.null(obj$y.ref)) {         
      #if (is.null(ncomp))
      #   ncomp = obj$ncomp.selected
      
      if (is.null(ny))
         ny = 1:ncol(obj$y.ref)
      
      if (length(ncomp) == 1)
         cat(sprintf('\nNumber of selected components: %d\n', ncomp))
      
      for (i in ny) {   
         cat(sprintf('\nResponse variable %s:\n', colnames(obj$y.ref)[i]))
         res = as.matrix.plsres(obj, ny = i, ncomp = ncomp)
         res[, 1:4] = round(res[, 1:4], 3)      
         res[, 5:6] = round(res[, 5:6], 3)
         res[, 7] = round(res[, 7], 4)      
         res[, 8] = round(res[, 8], 1)      
         print(res)
      }
      
   } else {
      cat('No reference data provided to calculate prediction performance.')
   }   
}   

#' print method for PLS results object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ...
#' other arguments
#' 
#' @export
print.plsres = function(x, ...) {
   obj = x
   
   cat('\nPLS results (class plsres)\n')
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$ncomp.selected - number of selected components\n')
   cat('$yp - array with predicted y values\n')
   if (!is.null(obj$y.ref)) {   
      cat('$y - matrix with reference y values\n')
      cat('$rmse - root mean squared error\n')
      cat('$r2 - coefficient of determination\n')
      cat('$slope - slope for predicted vs. measured values\n')
      cat('$bias - bias for prediction vs. measured values\n')
      cat('$ydecomp - decomposition of y values (ldecomp object)\n')
   }
   cat('$xdecomp - decomposition of x values (ldecomp object)\n')
}   

