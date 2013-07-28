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

plotRMSE.plsres = function(obj, xlab = 'Components', ...)
{
   # Makes a plot for RMSE vs number of components
   #
   # Arguments:
   #   obj: object of class "plsres"
   #   xlab: label for X axis
   
   plotRMSE.regres(obj, xlab = xlab, ...)
}

plotXScores.plsres = function(obj, comp = c(1, 2), main = NULL, ...)
{
   # Makes a scores plot for X decomposition
   #
   # Arguments:
   #   obj: object of class "plsres"
   #   ...: possible arguments for plotScores.ldecomp
   
   if (is.null(main))
      main = 'X scores'
   
   if (!is.null(obj$xdecomp$scores))
      plotScores.ldecomp(obj$xdecomp, comp = comp, main = main, ...)
   else
      warning('Scores values are not available.')
}

plotXYScores.plsres = function(obj, ncomp = 1, type = 'p', main = 'XY scores', 
                               xlab = 'X scores', ylab = 'Y scores', ...)
{
   # Makes a scores plot for X and Y decompositions
   #
   # Arguments:
   #   obj: object of class "plsres"
   #   ncomp: for which component makes the plot
   #   type: type of the plot
   #   main: text for main title of the plot
   #   xlab: text for x axis label
   #   ylab: text for y axis label
   
   if (is.null(obj$xdecomp$scores) || is.null(obj$ydecomp$scores))
   {   
      warning('X or Y scores are not available.')
   }
   else
   {   
      if (ncomp < 0 || ncomp > obj$ncomp)
         stop('Wrong value for ncomp argument!')
      
      data = cbind(obj$xdecomp$scores[, ncomp, drop = F],
                   obj$ydecomp$scores[, ncomp, drop = F])
      colnames(data) = c(xlab, ylab)
      mdaplot(data, type = type, main = sprintf('%s (ncomp = %d)', main, ncomp), ...)
   }   
}

plotYResiduals.plsres = function(obj, ncomp = NULL, ny = 1, ...)
{
   plotYResiduals.regres(obj, ncomp = ncomp, ny = ny, ...)   
}  

plotXResiduals.plsres = function(obj, ncomp = NULL, main = NULL, ...)
{   
   # Makes a residuals plot for X decomposition
   #
   # Arguments:
   #   obj: object of class "plsres"
   #   ncomp: number of components to show the plot for
   #   main: main title for the plot
      
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

plotXVariance.plsres = function(obj, main = 'X variance', ...)
{   
   # Makes a plot for X explained variance for each component
   #
   # Arguments:
   #   obj: object of class "plsres"
   #   main: main title for the plot
   
   plotVariance.ldecomp(obj$xdecomp, main = main, ...)
}

plotYVariance.plsres = function(obj, main = 'Y variance', ...)
{   
   # Makes a plot for Y explained variance for each components
   #
   # Arguments:
   #   obj: object of class "plsres"
   #   main: main title for the plot
   
   plotVariance.ldecomp(obj$ydecomp, main = main, ...)
}

plotXCumVariance.plsres = function(obj, main = 'X cumulative variance', ...)
{   
   # Makes a plot for X cumulative explained variance vs number of components
   #
   # Arguments:
   #   obj: object of class "plsres"
   #   main: main title for the plot
   
   plotCumVariance.ldecomp(obj$xdecomp, main = main, ...)
}

plotYCumVariance.plsres = function(obj, main = 'Y cumulative variance', ...)
{   
   # Makes a plot for Y cumulative explained variance vs number of components
   #
   # Arguments:
   #   obj: object of class "plsres"
   #   main: main title for the plot
   
   plotCumVariance.ldecomp(obj$ydecomp, main = main, ...)
}

plotPredictions.plsres = function(obj, ny = 1, ncomp = NULL, main = NULL, ...)
{
   # Makes a plot for predicted y values (vs measured or object number)
   #
   # Arguments:
   #   obj: object of class "plsres"
   #   ny: number of response variable to make the plot for
   #   ncomp: number of components to make the plot for
   #   main: main title for the plot
   
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

plot.plsres = function(obj, ncomp = NULL, ny = 1, show.labels = F, ...)
{
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected 
   
   if (is.null(obj$y.ref))
   {
      par(mfrow = c(1, 2))
      plotXScores(obj, ...)
      plotPredictions(obj, ncomp = ncomp, ny = ny, ...)
      par(mfrow = c(1, 1))      
   }  
   else
   {   
      par(mfrow = c(2, 2))
      if (!is.null(obj$xdecomp$scores))
         plotXScores(obj, ...)      
      else
         plotXVariance(obj, ...)   
      
      plotYVariance(obj, ...)   
      plotRMSE(obj, ny = ny, ...)
      plotPredictions(obj, ncomp = ncomp, ny = ny, ...)
      par(mfrow = c(1, 1))
   }
}   

as.matrix.plsres = function(obj, ncomp = NULL, ny = 1)
{
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

print.plsres = function(obj, ...)
{
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

summary.plsres = function(obj, ny = NULL, ncomp = NULL, ...)
{
   # Summary method for "plsres" object
   #
   # Arguments:
   #   obj: object of "plsres" class
   #   ncomp: which column of yp to use
   
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