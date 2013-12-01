regres = function(y.pred, y.ref = NULL, ncomp.selected = 1)
{
   # Class for storing and visualising of regression predictions 
   #
   # Arguments:
   #   y.pred: vector or matrix with y predicted values
   #   y.ref: vector with reference (measured) y values
   #   ncomp.selected: if y.pred calculated for different components, which to use as default
   #
   # Returns:
   # a list (object of "regres" class) with following fields
   #   y.pred: a matrix with predicted values
   #   y.ref: a vector with reference (measured) values
   #   ncomp.selected: selected column/number of components for predictions
   #   rmse: root mean squared error for predicted vs measured values
   #   slope: slope for predicted vs measured values
   #   r2: coefficient of determination for predicted vs measured values
   #   bias: bias for predicted vs measured values
   #   rpd: RPD values
   
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

regres.r2 = function(y.ref, y.pred)
{
   # Calculates determination coefficients (R2) 
   #
   # Arguments:
   #   y.ref: vector with reference values
   #   y.pred: matrix with predicted values
   #
   # Returns:
   #   r2: matrix with R2 values for each component and response variable
   
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   r2 = matrix(0, nrow = nresp, ncol = ncomp)
   
   for (i in 1:nresp)
      r2[i, ] = as.vector(cor(y.ref[, i], y.pred[, , i])^2)   

   rownames(r2) = colnames(y.ref)
   colnames(r2) = dimnames(y.pred)[[2]]
   
   r2
}  

regres.bias = function(y.ref, y.pred)
{
   # Calculates bias (average prediction error) 
   #
   # Arguments:
   #   y.ref: vector with reference values
   #   y.pred: matrix with predicted values
   #
   # Returns:
   #   bias: matrix with bias values for each component and response variable
   
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   bias = matrix(0, nrow = nresp, ncol = ncomp)
   
   for (i in 1:nresp)
      bias[i, ] = as.vector(apply(y.ref[, i] - y.pred[, , i], 2, mean))

   rownames(bias) = colnames(y.ref)
   colnames(bias) = dimnames(y.pred)[[2]]
   
   bias
}  

regres.rmse = function(y.ref, y.pred)
{
   # Calculates root mean squared error 
   #
   # Arguments:
   #   y.ref: vector with reference values
   #   y.pred: matrix with predicted values
   #
   # Returns:
   #   rmse: matrix with RMSE values for each component and response variable
   
   nresp = ncol(y.ref)
   ncomp = ncol(y.pred)
   rmse = matrix(0, nrow = nresp, ncol = ncomp)
   
   for (i in 1:nresp)
      rmse[i, ] = sqrt(colSums((y.ref[, i] - y.pred[, , i, drop = F])^2)/length(y.ref[, i]))      
   
   rownames(rmse) = colnames(y.ref)
   colnames(rmse) = dimnames(y.pred)[[2]]
   
   rmse
} 

regres.slope = function(y.ref, y.pred)
{
   # Calculates a slope for predicted vs. measured valus 
   #
   # Arguments:
   #   y.ref: vector with reference values
   #   y.pred: matrix with predicted values
   #
   # Returns:
   #   slope: matrix with slope values for each component and response variable
   
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

plotRMSE.regres = function(obj, ny = 1, type = 'b', main = 'RMSE', 
                           xlab = 'Complexity', ylab = NULL, ...)
{
   # Makes plot with RMSE values for each column of y.pred 
   #
   # Arguments:
   #   obj: object of "regres" class
   #   ny: number of response variable y to make the plot for
   #   type: type of the plot 
   #   main: main title for the plot
   #   xlab: text for x axis label
   #   ylab: text for y axis label
   
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

plotPredictions.regres = function(obj, ny = 1, ncomp = NULL, main = 'Predictions', 
                                  xlab = NULL, ylab = NULL, 
                                  show.line = T, colmap = 'default', col = NULL, ...)
{
   # Makes plot with predicted vs. measured values or predicted values vs. object number 
   # for selected column of y.pred 
   #
   # Arguments:
   #   obj: object of "regres" class
   #   ny: number of response variable y to make the plot for
   #   ncomp: which column of y.pred to use
   #   main: main title of the plot
   #   xlab: text for x axis label
   #   ylab: text for y axis label
   #   show.line: logical, show or not target for the points
   #   colmap: colormap for plot and target line
   #   col: color for plot and target line
   
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

plotYResiduals.regres = function(obj, ny = 1, ncomp = NULL, main = NULL, type = 'p',
                                 xlab = 'Objects', ylab = NULL, show.lines = T, ...)
{
   # Makes plot with Y residuals (difference between predicted and reference values) 
   # for selected column of y.pred 
   #
   # Arguments:
   #   obj: object of "regres" class
   #   ny: number of response variable y to make the plot for   
   #   ncomp: which column of y.pred to use
   #   main: main title of the plot
   #   xlab: text for x axis label
   #   ylab: text for y axis label   
   #   show.lines: logical, show or not zero line
   
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
      
      if (show.lines == T)
         show.lines = c(NA, 0)
      
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
      mdaplot(data, type = type, main = main, show.lines = show.lines, ...)
   }
}

plot.regres = function(obj, ny = 1, ...)
{
   # Plot method for "regres" objects
   #
   # Arguments:
   #   obj: object of "regres" class
   #   ny: number of response variable y to make the plot for   
   
   plotPredictions.regres(obj, ny = ny, ...)
}   

as.matrix.regres = function(obj, ncomp = NULL, ny = 1)
{
   # as.matrix method for "regres" objects
   #
   # Arguments:
   #   obj: object of "regres" class
   #   ny: number of response variable y to make the plot for   

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

print.regres = function(obj, ...)
{
   # Print method for "regres" object
   #
   # Arguments:
   #   obj: object of "regres" class
   
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

summary.regres = function(obj, ncomp = NULL, ny = NULL, ...)
{
   # Summary method for "regres" object
   #
   # Arguments:
   #   obj: object of "regres" class
   #   ncomp: which column of y.pred to use
   
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