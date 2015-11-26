#' Regression results
#' 
#' @description
#' Class for storing and visualisation of regression predictions 
#'
#' @param y.pred
#' vector or matrix with y predicted values
#' @param y.ref
#' vector with reference (measured) y values
#' @param ncomp.selected
#' if y.pred calculated for different components, which to use as default
#' 
#' @return
#' a list (object of \code{regres} class) with fields, including:
#' \tabular{ll}{
#'    \code{y.pred} \tab a matrix with predicted values \cr
#'    \code{y.pred} \tab a matrix with predicted values \cr
#'    \code{y.ref} \tab a vector with reference (measured) values \cr
#'    \code{ncomp.selected} \tab selected column/number of components for predictions \cr
#'    \code{rmse} \tab root mean squared error for predicted vs measured values \cr
#'    \code{slope} \tab slope for predicted vs measured values \cr
#'    \code{r2} \tab coefficient of determination for predicted vs measured values \cr
#'    \code{bias} \tab bias for predicted vs measured values \cr
#'    \code{rpd} \tab RPD values \cr
#' }
#' 
#' @export
regres = function(y.pred, y.ref = NULL, ncomp.selected = 1)
{   
   obj = list()
   obj$y.pred = y.pred
   obj$ncomp.selected = ncomp.selected
   
   if (!is.null(y.ref))
   {
      y.ref = as.matrix(y.ref)      
      obj$y.ref = y.ref
      obj$rmse = regres.rmse(y.ref, y.pred)
      obj$slope = regres.slope(y.ref, y.pred)
      obj$r2 = regres.r2(y.ref, y.pred)
      obj$bias = regres.bias(y.ref, y.pred)
      obj$sep = sqrt(obj$rmse^2 - obj$bias^2)
      obj$rpd = apply(y.ref, 2, sd)/obj$sep
   }
         
   obj$call = match.call()   
   class(obj) = "regres"
   
   obj
}

#' Determination coefficient
#' 
#' @description
#' Calculates matrix with coeffient of determination for every response and components 
#'
#' @param y.ref
#' vector with reference values
#' @param y.pred
#' matrix with predicted values
#'
regres.r2 = function(y.ref, y.pred)
{
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   r2 = matrix(0, nrow = nresp, ncol = ncomp)
   
   for (i in 1:nresp)
      r2[i, ] = as.vector(cor(y.ref[, i], y.pred[, , i])^2)   

   rownames(r2) = colnames(y.ref)
   colnames(r2) = dimnames(y.pred)[[2]]
   
   r2
}  

#' Prediction bias 
#' 
#' @description
#' Calculates matrix with bias (average prediction error) for every response and components 
#'
#' @param y.ref
#' vector with reference values
#' @param y.pred
#' matrix with predicted values
#'
regres.bias = function(y.ref, y.pred)
{
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   bias = matrix(0, nrow = nresp, ncol = ncomp)
   
   for (i in 1:nresp)
      bias[i, ] = as.vector(apply(y.ref[, i] - y.pred[, , i, drop = F], 2, mean))

   rownames(bias) = colnames(y.ref)
   colnames(bias) = dimnames(y.pred)[[2]]
   
   bias
}  

#' RMSE
#' 
#' @description
#' Calculates matrix with root mean squared error of prediction for every response and components.
#'
#' @param y.ref
#' vector with reference values
#' @param y.pred
#' matrix with predicted values
#'
regres.rmse = function(y.ref, y.pred)
{
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   rmse = matrix(0, nrow = nresp, ncol = ncomp)
   
   for (i in 1:nresp)
      rmse[i, ] = sqrt(colSums((y.ref[, i] - y.pred[, , i, drop = F])^2)/length(y.ref[, i]))      
   
   rownames(rmse) = colnames(y.ref)
   colnames(rmse) = dimnames(y.pred)[[2]]
   
   rmse
} 

#' Slope 
#' 
#' @description
#' Calculates matrix with slope of predicted and measured values for every response and components.
#'
#' @param y.ref
#' vector with reference values
#' @param y.pred
#' matrix with predicted values
#'
regres.slope = function(y.ref, y.pred)
{
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   slope = matrix(0, nrow = nresp, ncol = ncomp)
   
   for (i in 1:nresp)
   {   
      for (a in 1:ncomp)
      {   
         m = lm(y.pred[, a, i] ~ y.ref[, i])
         slope[i, a] = m$coefficients[[2]]
      }
   }

   rownames(slope) = colnames(y.ref)
   colnames(slope) = dimnames(y.pred)[[2]]
   
   slope
}   

#' RMSE plot for regression results
#' 
#' @description
#' Shows plot with RMSE values vs. model complexity (e.g. number of components).
#'
#' @param obj
#' regression results (object of class \code{regres})
#' @param ny
#' number of predictor to show the plot for (if y is multivariate)
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @export
plotRMSE.regres = function(obj, ny = 1, type = 'b', main = 'RMSE', 
                           xlab = 'Complexity', ylab = NULL, ...)
{
   if (!is.null(obj$rmse))
   {   
      if (is.null(ylab))
         if (ncol(obj$y.ref) > 1 && !is.null(colnames(obj$y.ref)))
            ylab = sprintf('RMSE (%s)', colnames(obj$y.ref)[ny])
      else
         ylab = 'RMSE'
      
      data = cbind(1:ncol(obj$y.pred[, , ny]), obj$rmse[ny, ])
      colnames(data) = c(xlab, ylab)
      rownames(data) = mdaplot.formatValues(obj$rmse[ny, ])
      mdaplot(data, type = type, main = main, ...)
   }
   else
   {
      warning('RMSE values are not available.')
   }   
}

#' Predictions plot for regression results
#' 
#' @description
#' Shows plot with predicted y values.
#'
#' @param obj
#' regression results (object of class \code{regres})
#' @param ny
#' number of predictor to show the plot for (if y is multivariate)
#' @param ncomp
#' complexity of model (e.g. number of components) to show the plot for
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.line
#' logical, show or not line fit for the plot points
#' @param colmap
#' a colormap to use for coloring the plot items
#' @param col
#' a vector with color values for plot items
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' If reference values are available, the function shows a scatter plot with predicted vs. 
#' reference values, otherwise predicted values are shown vs. object numbers.
#' 
#' @export
plotPredictions.regres = function(obj, ny = 1, ncomp = NULL, main = 'Predictions', 
                                  xlab = NULL, ylab = NULL, 
                                  show.line = T, colmap = 'default', col = NULL, ...)
{
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   
   if (is.null(ylab))
   {   
      if (!is.null(dimnames(obj$y.pred)) && !is.null(dimnames(obj$y.pred)[3]))
         ylab = sprintf('%s, predicted', dimnames(obj$y.pred)[[3]][ny])
      else
         ylab = 'y, predicted'
   }
   
   if (is.null(obj$y.ref))
   {   
      if (is.null(xlab))
         xlab = 'Objects'
      
      data = cbind(1:nrow(obj$y.pred[, , ny]), obj$y.pred[, ncomp, ny, drop = F])     
      mdaplot(data, type = 'p', main = main, colmap = colmap, xlab = xlab, ylab = ylab, col = col, ...)
   }
   else
   {      
      if (is.null(xlab))
      {   
         if (ncol(obj$y.ref) > 1)
            xlab = sprintf('%s, measured', colnames(obj$y.ref)[ny])
         else
            xlab = 'y, measured'
      }
      data = cbind(obj$y.ref[, ny], obj$y.pred[, ncomp, ny, drop = F])
      mdaplot(data, type = 'p', main = main, xlab = xlab, ylab = ylab, colmap = colmap, col = col, ...)
      
      if (show.line == T)
         mdaplot.showRegressionLine(data, colmap = colmap, col = col)
   }
}

#' Residuals plot for regression results
#'
#' @description
#' Shows plot with Y residuals (difference between predicted and reference values) for selected 
#' response variable and complexity (number of components). 
#'
#' @param obj
#' regression results (object of class \code{regres})
#' @param ny
#' number of predictor to show the plot for (if y is multivariate)
#' @param ncomp
#' complexity of model (e.g. number of components) to show the plot for
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.line
#' logical, show or zero line on the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @export
plotYResiduals.regres = function(obj, ny = 1, ncomp = NULL, type = 'p', main = NULL,
                                 xlab = 'Objects', ylab = NULL, show.line = T, ...)
{
   if (is.null(obj$y.ref))
   {   
      warning('Y residuals can not be plotted without reference values.')
   }
   else
   {
      if (is.null(ncomp))
         ncomp = obj$ncomp.selected
      else if (ncomp < 1 || ncomp > ncol(obj$y.pred))
         stop('Wrong number of components!')
      
      if (show.line == T)
         show.line = c(NA, 0)
      
      if (is.null(main))
      {   
         if (is.null(ncomp))
            main = 'Y residuals'
         else
            main = sprintf('Y residuals (ncomp = %d)', ncomp)
      }
      
      if (is.null(ylab))
      {   
         if (ncol(obj$y.ref) > 1)
            ylab = sprintf('Residuals (%s)', colnames(obj$y.ref)[ny])
         else
            ylab = 'Residuals'
      }
      data = cbind(1:nrow(obj$y.pred[, , ny]), obj$y.ref[, ny] - obj$y.pred[, ncomp, ny, drop = F])
      colnames(data) = c(xlab, ylab)
      mdaplot(data, type = type, main = main, show.lines = show.line, ...)
   }
}

#' plot method for regression results
#' 
#' @details
#' Shows prediction plot for the results (the same as \code{plotPredictions.regres})
#' 
#' @param x
#' regression results (object of class \code{regres})
#' @param ny
#' which response to show the plot for (if y is multivariate)
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @export
plot.regres = function(x, ny = 1, ...)
{
   obj = x
   
   plotPredictions.regres(obj, ny = ny, ...)
}   

#' as.matrix method for regression results
#' 
#' @description
#' Returns a matrix with model performance statistics for regression results
#' 
#' @param x
#' regression results (object of class \code{regres})
#' @param ncomp
#' model complexity (number of components) to calculate the statistics for
#' @param ny
#' for which response variable calculate the statistics for
#' @param ...
#' other arguments
#' 
#' @export
as.matrix.regres = function(x, ncomp = NULL, ny = 1, ...)
{
   obj = x
   
   if (!is.null(obj$y.ref))
   {  
      if (is.null(ncomp))
         res = cbind(obj$rmse[ny, ], obj$r2[ny, ], obj$slope[ny, ], obj$bias[ny, ], obj$rpd[ny, ])   
      else
         res = cbind(obj$rmse[ny, ncomp], obj$r2[ny, ncomp], obj$slope[ny, ncomp], 
                     obj$bias[ny, ncomp], obj$rpd[ny, ncomp])   
      
      colnames(res) = c('RMSE', 'R^2', 'Slope', 'Bias', 'RPD')
   }
   else
   {
      res = NULL
   }   
   
   res
}

#' summary method for regression results object
#' 
#' @description
#' Shows performance statistics for the regression results.
#' 
#' @param object
#' regression results (object of class \code{regres})
#' @param ncomp
#' model complexity to show the summary for
#' @param ny
#' for which response variable show the summary for
#' @param ...
#' other arguments
#' 
#' @export
summary.regres = function(object, ncomp = NULL, ny = NULL, ...)
{
   obj = object
   
   cat('\nRegression results (class regres) summary\n')
   if (!is.null(obj$y.ref))
   {         
      if (is.null(ncomp))
         ncomp = obj$ncomp.selected
      
      if (is.null(ny))
         ny = 1:ncol(obj$y.ref)
      
      if (!is.null(ncomp))
         cat(sprintf('\nNumber of selected components: %d\n\n', ncomp))
         
      for (i in ny)
      {   
         cat(sprintf('\nResponse variable %s:\n', colnames(obj$y.ref)[i]))
         res = as.matrix.regres(obj, ny = i, ncomp = ncomp)
         rownames(res) = ncomp
         print(res)
      }
      
   }      
   else
   {
      cat('No reference data provided to calculate prediction performance.')
   }   
}   

#' print method for regression results object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' regression results (object of class \code{regres})
#' @param ...
#' other arguments
#' 
#' @export
print.regres = function(x, ...)
{
   obj = x
   
   cat('\nRegression results (class regres)\n')
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$y.pred - matrix or vector with predicted y values\n')
   if (!is.null(obj$y.ref))
   {   
      cat('$y.ref - vector with reference y values\n')
      cat('$rmse - root mean squared error\n')
      cat('$r2 - coefficient of determination\n')
      cat('$slope - slope for predicted vs. measured values\n')
      cat('$bias - bias for prediction vs. measured values\n')
   }
   
   if (ncol(obj$y.pred) > 1)   
      cat('$ncomp.selected - number of selected components for PCR or PLS\n')
}   

