plsresult = function(yp, y, ...) UseMethod("plsresult")

plsresult.default = function(yp, y = NULL, X = NULL, nlv = 0, xscores = NULL, yscores = NULL,
                             xloadings = NULL, yloadings = NULL)
{
   yp = as.matrix(yp)
   
   if (!is.null(xloadings))
   {   
      xdist = ldecomp.getDistances(xscores, xloadings, X)
      xvar = ldecomp.getVariances(xdist$T2, sum(X^2))
      Q2 = xdist$Q2
      T2 = xdist$T2
      xexpvar = xvar$expvar
   }
   else
   {
      xexpvar = NULL
      Q2 = NULL
      T2 = NULL
   }   
   
   if (!is.null(yloadings) && !is.null(y))
   {   
      ydist = ldecomp.getDistances(yscores, yloadings, y)
      yvar = ldecomp.getVariances(ydist$T2, sum(y^2))
      yexpvar = yvar$expvar
   }
   else
   {
      yexpvar = NULL
   }   
   
   if (!is.null(y))
   {   
      y = as.numeric(y)
      
      plsresult = list(
         yp = yp,
         y = y,
         xscores = xscores,
         yscores = yscores,
         T2 = T2,
         Q2 = Q2,
         xexpvar = xexpvar,
         yexpvar = yexpvar,
         rmse = plsresult.rmse(y, yp),
         slope = plsresult.slope(y, yp),
         r2 = cor(y, yp)^2,
         bias = apply(y - yp, 2, mean)
      )      
   }
   else
   {
      plsresult = list(
         yp = yp,
         y = y,
         xscores = xscores,
         yscores = yscores,
         T2 = T2,
         Q2 = Q2,
         xexpvar = xexpvar,
         yexpvar = yexpvar,
         rmse = NULL,
         slope = NULL,
         r2 = NULL,
         bias = NULL
      )
   }   
   
   if (nlv == 0) { nlv = ncol(yp) }   
   
   plsresult$nlvselected = nlv
   plsresult$call = match.call()   
   class(plsresult) = "plsresult"
   
   plsresult   
}

# calculates root mean squared error for prediction
plsresult.rmse = function(y, yp)
{
   return (sqrt(colSums((y - yp)^2)/length(y)))
} 

# calculates slope of a line fit for (y, yp) values
plsresult.slope = function(y, yp)
{
   nlv = dim(yp)[2]
   slope = 1:nlv
   for (a in 1:nlv)
   {   
      m = lm(yp[, a] ~ y)
      slope[a] = m$coefficients[2]
   }
   
   return (slope)
}   

plsresult.plot_rmse = function(result, col = 'blue', type = 'b', main = 'RMSE', xlab = 'LVs', 
                               ylab = 'RMSE', pch = 16, xlim = NULL, ylim = NULL)
{
   if (!is.null(result$rmse))
   {
      nlv = length(result$rmse)
      
      if (is.null(xlim) || is.null(ylim))
      {   
         plot(1:nlv, result$rmse, type = type, col = col, main = main, pch = pch,
              xlab = xlab, ylab = ylab)
      } 
      else
      {
         plot(1:nlv, result$rmse, type = type, col = col, main = main, pch = pch,
              xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)         
      }   
      points(result$nlvselected, result$rmse[result$nlvselected], col = 'red')
      grid()   
   }   
}

plsresult.points_rmse = function(result, col = 'blue', type = 'b',  pch = 16)
{
   if (!is.null(result$rmse))
   {
      nlv = length(result$rmse)      
      lines(1:nlv, result$rmse, type = type, col = col, pch = pch)
      points(result$nlvselected, result$rmse[result$nlvselected], col = 'red')
      grid()   
   }   
}

plsresult.plot_predictions = function(result, nlv = 0, col = 'blue', 
                                      labels = F,
                                      main = 'Predictions', 
                                       xlab = NULL, ylab = 'y, predicted', 
                                      pch = 16, ...)
{
   if (nlv == 0) { nlv = result$nlvselected }
   
   if (!is.null(result$y))
   {  
      if (is.null(xlab)) { xlab = 'y, measured' }   
      
      # show predicted versus measured values plot
      plot(result$y, result$yp[, nlv], col = col, pch = pch,
           main = main, xlab = xlab, ylab = ylab, cex.axis = 0.85, ...)
      abline(lm(result$yp[, nlv] ~ result$y), col = col)      
      if (labels == T)
      {   
         text(result$y, result$yp[,nlv], rownames(result$yp), cex = 0.6, pos = 3, col = 'gray')
      }
   }
   else
   {
      if (is.null(xlab)) { xlab = 'Objects' }   
      
      # show predicted y values for every object
      plot(1:length(result$yp[, nlv]), result$yp[, nlv], 
           col = col, pch = pch,
           main = main,
           xlab = xlab,
           ylab = ylab,
           cex.axis = 0.85,           
           ...
      )
      if (labels == T)
      {   
         text(1:length(result$yp[, nlv]), result$yp[,nlv], rownames(result$yp), cex = 0.6, pos = 3, col = 'gray')
      }
   }   
   grid()            
}

plsresult.plot_xscores = function(result, nlv = c(1, 2), col = 'blue', main = 'X scores', 
                                  pch = 16, labels = F, ...)
{
   plot(result$xscores[, nlv[1]], result$xscores[, nlv[2]], col = col, pch = pch,
        xlab = colnames(result$xscores)[nlv[1]],
        ylab = colnames(result$xscores)[nlv[2]],
        main = main, ...)

   if (labels == T)
   {   
      text(result$xscores[, nlv[1]], result$xscores[, nlv[2]], rownames(result$xscores), 
           cex = 0.6, pos = 3, col = 'gray')
   }
}

plsresult.points_xscores = function(result, nlv = c(1, 2), col = 'blue', labels = F, pch = 16, ...)
{
   points(result$xscores[, nlv[1]], result$xscores[, nlv[2]], col = col, pch = pch, ...)
   if (labels == T)
   {   
      text(result$xscores[, nlv[1]], result$xscores[, nlv[2]], rownames(result$xscores), 
           cex = 0.6, pos = 3, col = 'gray')
   }
}

plsresult.plot_xyscores = function(result, nlv = 1)
{
   plot(result$xscores[, nlv], result$yscores[, nlv])
}

plsresult.plot_xresiduals = function()
{
   
}

plsresult.plot_yresiduals = function()
{
   
}

plsresult.points_predictions = function(result, nlv = 0, labels = F, col = 'blue', pch = 16, ...)
{
   if (nlv == 0) { nlv = result$nlvselected }
   
   if (!is.null(result$y))
   {        
      # show predicted versus measured values plot
      points(result$y, result$yp[, nlv], col = col, pch = pch, ... )
      abline(lm(result$yp[, nlv] ~ result$y), col = col)
      if (labels == T)
      {   
         text(result$y, result$yp[,nlv], rownames(result$yp), cex = 0.6, pos = 3, col = 'gray')
      }
   }
   else
   {
      # show predicted y values for every object
      points(1:length(result$yp[, nlv]), result$yp[, nlv], col = col, pch = pch, ... )
      if (labels == T)
      {   
         text(1:length(result$yp[, nlv]), result$yp[,nlv], rownames(result$yp), cex = 0.6, pos = 3, col = 'gray')
      }
   }      
}

plot.plsresult = function(result)
{
   plsresult.plot_predictions(result)
}   
  
as.matrix.plsresult = function(result)
{
   nlv = result$nlvselected
   
   if (!is.null(result$y))
   {   
      res = matrix(c(result$rmse, result$r2, result$slope, result$bias), ncol = 4)   
      colnames(res) = c('RMSE', 'R^2', 'Slope', 'Bias')
   }
   else
   {
      res = NULL
   }   
   return (res)
}

print.plsresult = function(result, ...)
{
   cat('\nPLS results (class plsresult)\n')
   cat('\nPredictions:\n')
}   

summary.plsresult = function(result, ...)
{
   cat('\nPLS results (class plsresult) summary\n')
   if (!is.null(result$y))
   {   
      nlv = result$nlvselected
      res = as.matrix(result)[nlv, , drop = FALSE]
      
      cat('\nPerformance:')
      cat('\n')
      cat(sprintf('\n\t%s: %d', 'Number of LVs', nlv))
      cat(sprintf('\n\t%s: %.4f', colnames(res), res))
   }      
   else
   {
      cat('No reference data provided to calculate prediction performance.')
   }   
}   


