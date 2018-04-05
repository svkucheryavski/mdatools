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
regres = function(y.pred, y.ref = NULL, ncomp.selected = 1) {   
   obj = list()
   obj$y.pred = y.pred
   obj$ncomp.selected = ncomp.selected
   
   if (!is.null(y.ref)) {
      attrs = mda.getattr(y.ref)
      y.ref = as.matrix(y.ref)      
      y.ref = mda.setattr(y.ref, attrs)
      obj$y.ref = y.ref
      
      # remove excluded rows so they are not counted
      # when calculating statistics
      attrs = mda.getattr(y.pred)
      if (length(attrs$exclrows) > 0) {
         y.pred = y.pred[-attrs$exclrows, , , drop = F]
         if (nrow(y.ref) > nrow(y.pred))
            y.ref = y.ref[-attrs$exclrows, , drop = F]
      }
      
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
regres.r2 = function(y.ref, y.pred) {
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   r2 = matrix(0, nrow = nresp, ncol = ncomp)
  
   ytot = colSums((y.ref - mean(y.ref))^2)
   for (i in 1:nresp){
      yp = y.pred[, , i, drop = F]
      dim(yp) = dim(y.pred)[1:2]
      r2[i, ] = (1 - colSums(apply(yp, 2, '-', y.ref[, i])^2)/ytot[i]) * 100   
   }
   
   rownames(r2) = colnames(y.ref)
   colnames(r2) = dimnames(y.pred)[[2]]
   attr(r2, 'name') = 'Coefficient of determination'
   attr(r2, 'xaxis.name') = 'Components'
   attr(r2, 'yaxis.name') = 'Predictors'
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
regres.bias = function(y.ref, y.pred) {
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   bias = matrix(0, nrow = nresp, ncol = ncomp)
   
   for (i in 1:nresp)
      bias[i, ] = as.vector(apply(y.ref[, i] - y.pred[, , i, drop = F], 2, mean))

   rownames(bias) = colnames(y.ref)
   colnames(bias) = dimnames(y.pred)[[2]]
   attr(bias, 'name') = 'Bias'
   attr(bias, 'xaxis.name') = 'Components'
   attr(bias, 'yaxis.name') = 'Predictors'
   
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
regres.rmse = function(y.ref, y.pred) {
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   rmse = matrix(0, nrow = nresp, ncol = ncomp)
   
   for (i in 1:nresp)
      rmse[i, ] = sqrt(colSums((y.ref[, i] - y.pred[, , i, drop = F])^2)/length(y.ref[, i]))      
   
   rownames(rmse) = colnames(y.ref)
   colnames(rmse) = dimnames(y.pred)[[2]]
   attr(rmse, 'name') = 'RMSE'
   attr(rmse, 'xaxis.name') = 'Components'
   attr(rmse, 'yaxis.name') = 'Predictors'
   
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
regres.slope = function(y.ref, y.pred) {
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
   attr(slope, 'name') = 'Slope'
   attr(slope, 'xaxis.name') = 'Components'
   attr(slope, 'yaxis.name') = 'Predictors'
   
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
#' @param labels
#' what to show as labels for plot objects. 
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @export
plotRMSE.regres = function(obj, ny = 1, type = 'b', labels = 'values', ...) {
   
   if (is.null(obj$rmse)) {
      warning('RMSE values are not available.')
      return()
   }
 
   data = mda.subset(obj$rmse, ny)
   attr(data, 'xaxis.name') = 'Components'
   attr(data, 'yaxis.name') = 'RMSE'
   attr(data, 'name') = 'RMSE'
   if (length(ny) == 1)
      mdaplot(data, type = type, labels = labels, ...)
   else
      mdaplotg(data, type = type, labels = labels, ...)
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
#' @param show.line
#' logical, show or not line fit for the plot points
#' @param show.stat 
#' logical, show or not legend with statistics on the plot
#' @param stat.col 
#' color of text in legend with statistics
#' @param stat.cex 
#' size of text in legend with statistics
#' @param axes.equal
#' logical, make limits for x and y axes equal or not
#' @param col
#' color for the plot objects.
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' If reference values are available, the function shows a scatter plot with predicted vs. 
#' reference values, otherwise predicted values are shown vs. object numbers.
#' 
#' @export
plotPredictions.regres = function(obj, ny = 1, ncomp = NULL, show.line = T, 
                                  show.stat = F, stat.col = '#606060', stat.cex = 0.85,
                                  axes.equal = T, col = mdaplot.getColors(1), ...) {
   
   if (length(ny) != 1)
      stop('You can show prediction plot only for one selected response variable!')

   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   else if (length(ncomp) != 1 || ncomp < 1 || ncomp > ncol(obj$y.pred))
      stop('Wrong number of components!')
   
   if (is.null(dimnames(obj$y.pred)) || is.null(dimnames(obj$y.pred)[[3]])) 
      yaxis.name = 'y, predicted'
   else
      yaxis.name = sprintf('%s, predicted', dimnames(obj$y.pred)[[3]][ny])
   
   if (!is.null(dimnames(obj$y.pred)) && !is.null(dimnames(obj$y.pred)[[3]])) 
      xaxis.name = sprintf('%s, reference', dimnames(obj$y.pred)[[3]][ny])
   else
      xaxis.name = 'y, reference'

   attrs = mda.getattr(obj$y.pred)
   if (is.null(obj$y.ref)) {   
      data = matrix(obj$y.pred[, ncomp, ny], ncol = 1)
      xaxis.name = NULL
   } else {      
      data = cbind(obj$y.ref[, ny], obj$y.pred[, ncomp, ny])
   }
   
   lim = c(min(data), max(data))
   data = mda.setattr(data, attrs)
   colnames(data) = c(xaxis.name, yaxis.name)
   rownames(data) = rownames(obj$y.pred)
   attr(data, 'name') = 'Predictions'
   
   if (axes.equal && !is.null(lim))
      p = mdaplot(data, type = 'p', xlim = lim, ylim = lim, ...)
   else
      p = mdaplot(data, type = 'p', ...)
   
   if (show.stat && !is.null(obj$y.ref)) {
      dl = (lim[2] - lim[1])/20
      stat.text = paste(
         'nLV = ', ncomp, '\n',
         'RMSE = ', format(obj$rmse[ncomp], digits = 3), '\n',
         'R2 = ', round(obj$ydecomp$cumexpvar[ncomp]/100, 3),
         sep = ''
      )
      
      text(lim[1] + dl, lim[2], stat.text, adj = c(0, 1), col = stat.col, cex = stat.cex)
   }
   if (show.line == T && ncol(data) == 2)
      mdaplot.showRegressionLine(data, colmap = 'default', col = col)
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
#' @param show.line
#' logical, show or not zero line on the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @export
plotYResiduals.regres = function(obj, ny = 1, ncomp = NULL, show.line = T, ...) {

   if (is.null(obj$y.ref)) {   
      warning('Y residuals can not be plotted without reference values.')
   }

   if (length(ny) != 1)
      stop('You can show prediction plot only for one selected response variable!')
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   else if (length(ncomp) != 1 || ncomp < 1 || ncomp > ncol(obj$y.pred))
      stop('Wrong number of components!')
   
   if (show.line == T)
      show.line = c(NA, 0)
      
   if (is.null(ncomp))
      name = 'Y residuals'
   else
      name = sprintf('Y residuals (ncomp = %d)', ncomp)

   if (!is.null(dimnames(obj$y.pred)) && !is.null(dimnames(obj$y.pred)[[3]])) 
      xaxis.name = sprintf('%s, reference', dimnames(obj$y.pred)[[3]][ny])
   else
      xaxis.name = 'y, reference'
   
   attr = mda.getattr(obj$y.pred)
   data = cbind(obj$y.ref[, ny], obj$y.ref[, ny] - obj$y.pred[, ncomp, ny])
   data = mda.setattr(data, attr)
   colnames(data) = c(xaxis.name, 'Residuals')
   attr(data, 'name') = name

   mdaplot(data, type = 'p', show.lines = show.line, ...)
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
as.matrix.regres = function(x, ncomp = NULL, ny = 1, ...) {
   obj = x
   
   if (!is.null(obj$y.ref)) {  
      if (is.null(ncomp))
         res = cbind(obj$rmse[ny, ], obj$r2[ny, ], obj$slope[ny, ], obj$bias[ny, ], obj$rpd[ny, ])   
      else
         res = cbind(obj$rmse[ny, ncomp], obj$r2[ny, ncomp], obj$slope[ny, ncomp], 
                     obj$bias[ny, ncomp], obj$rpd[ny, ncomp])   
      
      colnames(res) = c('RMSE', 'R^2', 'Slope', 'Bias', 'RPD')
   } else {
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
summary.regres = function(object, ncomp = NULL, ny = NULL, ...) {
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
print.regres = function(x, ...){
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

