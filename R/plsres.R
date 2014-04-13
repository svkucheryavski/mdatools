plsres = function(y.pred, y.ref = NULL, ncomp.selected = NULL, xdecomp = NULL, ydecomp = NULL, info = '')
{
   # Class for storing and viasualising PLS results:
   #
   # Arguments:
   #   y.pred: vector or matrix with predicted values
   #   y.ref: vector with reference (measured) values
   #   ncomp.selected: if yp is a matrix, wich column to use
   #   xdecomp: PLS decomposition of X data (class "ldecomp")
   #   ydecomp: PLS decomposition of Y data (class "ldecomp")
   #   info: information about the object
   
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

#' RMSE plot for PLS results
#' 
#' @description
#' Shows plot with root mean squared error values vs. number of components for PLS results.
#' 
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param xlab
#' label for x axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
plotRMSE.plsres = function(obj, xlab = 'Components', ...)
{
   plotRMSE.regres(obj, xlab = xlab, ...)
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
plotXScores.plsres = function(obj, comp = c(1, 2), main = NULL, ...)
{
   if (is.null(main))
      main = 'X scores'
   
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
#' @param comp
#' which component to show the plot for
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#' 
plotXYScores.plsres = function(obj, comp = 1, type = 'p', main = 'XY scores', 
                               xlab = 'X scores', ylab = 'Y scores', ...)
{
   if (is.null(obj$xdecomp$scores) || is.null(obj$ydecomp$scores))
   {   
      warning('X or Y scores are not available.')
   }
   else
   {   
      if (comp < 0 || comp > obj$ncomp)
         stop('Wrong value for ncomp argument!')
      
      data = cbind(obj$xdecomp$scores[, comp, drop = F],
                   obj$ydecomp$scores[, comp, drop = F])
      colnames(data) = c(xlab, ylab)
      mdaplot(data, type = type, main = sprintf('%s (ncomp = %d)', main, comp), ...)
   }   
}

#' X residuals plot for PLS results
#' 
#' @description
#' Shows a plot with Q2 residuals vs. Hotelling T2 values for PLS decomposition of x data.
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
plotXResiduals.plsres = function(obj, ncomp = NULL, main = NULL, ...)
{   
   if (is.null(main))
   {   
      if (is.null(ncomp))
         main = 'X residuals'
      else   
         main = sprintf('X residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected

   if (ncomp < 0 || ncomp > obj$ncomp)
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
plotXVariance.plsres = function(obj, main = 'X variance', ...)
{   
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
plotYVariance.plsres = function(obj, main = 'Y variance', ...)
{   
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
plotXCumVariance.plsres = function(obj, main = 'X cumulative variance', ...)
{   
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
plotYCumVariance.plsres = function(obj, main = 'Y cumulative variance', ...)
{   
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
plotPredictions.plsres = function(obj, ny = 1, ncomp = NULL, main = NULL, ...)
{
   if (is.null(main))
   {
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
plot.plsres = function(x, ncomp = NULL, ny = 1, show.labels = F, ...)
{
   obj = x
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected 
   
   if (is.null(obj$y.ref))
   {
      par(mfrow = c(1, 2))
      plotXResiduals(obj, ...)
      plotPredictions.plsres(obj, ncomp = ncomp, ny = ny, ...)
      par(mfrow = c(1, 1))      
   }  
   else
   {   
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
#' @method as.matrix plsres
#' @S3method as.matrix plsres
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
as.matrix.plsres = function(x, ncomp = NULL, ny = 1, ...)
{
   obj = x
   
   if (is.null(ncomp))
   {   
      res = cbind(
         obj$xdecomp$expvar,
         obj$xdecomp$cumexpvar,
         obj$ydecomp$expvar,
         obj$ydecomp$cumexpvar,
         as.matrix.regres(obj, ny = ny)
      )
   }
   else
   {
      res = cbind(
         obj$xdecomp$expvar[ncomp],
         obj$xdecomp$cumexpvar[ncomp],
         obj$ydecomp$expvar[ncomp],
         obj$ydecomp$cumexpvar[ncomp],
         as.matrix.regres(obj, ncomp = ncomp, ny = ny)
      )
   }   
   colnames(res)[1:4] = c('X expvar', 'X cumexpvar', 'Y expvar', 'Y cumexpvar')
   
   res
}


#' summary method for PLS results object
#' 
#' @method summary plsres
#' @S3method summary plsres
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
summary.plsres = function(object, ny = NULL, ncomp = NULL, ...)
{
   obj = object
   
   cat('\nPLS regression results (class plsres) summary\n')
   if (!is.null(obj$y.ref))
   {         
      if (is.null(ncomp))
         ncomp = obj$ncomp.selected
      
      if (is.null(ny))
         ny = 1:ncol(obj$y.ref)
      
      if (length(ncomp) == 1)
         cat(sprintf('\nNumber of selected components: %d\n', ncomp))
      
      for (i in ny)
      {   
         cat(sprintf('\nResponse variable %s:\n', colnames(obj$y.ref)[i]))
         res = as.matrix.plsres(obj, ny = i, ncomp = ncomp)
         res[, 1:4] = round(res[, 1:4], 3)      
         res[, 6:7] = round(res[, 6:7], 3)  
         res[, 5] = mdaplot.formatValues(res[, 5], round.only = T)
         res[, 8] = round(res[, 8], 4)      
         res[, 9] = round(res[, 9], 1)      
         rownames(res) = colnames(obj$y.pred)[ncomp]
         print(res)
      }
      
   }      
   else
   {
      cat('No reference data provided to calculate prediction performance.')
   }   
}   

#' print method for PLS results object
#' 
#' @method print plsres
#' @S3method print plsres
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ...
#' other arguments
#' 
print.plsres = function(x, ...)
{
   obj = x
   
   cat('\nPLS results (class plsres)\n')
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$ncomp.selected - number of selected components\n')
   cat('$yp - array with predicted y values\n')
   if (!is.null(obj$y.ref))
   {   
      cat('$y - matrix with reference y values\n')
      cat('$rmse - root mean squared error\n')
      cat('$r2 - coefficient of determination\n')
      cat('$slope - slope for predicted vs. measured values\n')
      cat('$bias - bias for prediction vs. measured values\n')
      cat('$ydecomp - decomposition of y values (ldecomp object)\n')
   }
   cat('$xdecomp - decomposition of x values (ldecomp object)\n')
   
}   

