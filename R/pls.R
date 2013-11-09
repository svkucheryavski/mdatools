# class and methods for Partial Least Squares regression #

pls = function(x, y, ncomp = 20, center = T, scale = F, cv = NULL, 
               x.test = NULL, y.test = NULL, method = 'simpls', alpha = 0.05, info = '')
{
   # Calibrate and validate a PLS model.
   #
   # Arguments:
   #   x: a matrix with predictor values    
   #   y: a vector with response values
   #   ncomp: maximum number of components to calculate
   #   center: logical, center or not x and y data
   #   scale: logical, standardize or not x data
   #   cv: number of segments for cross-validation (1 - full CV)
   #   x.test: a matrix with predictor values for test set validation
   #   y.test: a vector with response values for test set validation
   #   method: a method to calculate PLS model
   #   alpha: a sigificance limit for Q2 values
   #   info: a short string with information about the model
   #
   # Returns:
   #   model: a PLS model (object of class pls) 
   
   x = as.matrix(x)
   y = as.matrix(y)
   
   if (is.null(colnames(y)))
      colnames(y) = paste('y', 1:ncol(y), sep = '')
   
   # correct maximum number of components
   ncomp = min(ncol(x), nrow(x) - 1, ncomp)
   
   # build a model and apply to calibration set
   model = pls.cal(x, y, ncomp, center = center, scale = scale, method = method)
   model$ncomp.selected = ncomp
   model$alpha = alpha   
   model$calres = predict.pls(model, x, y)
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = pls.crossval(model, x, y, cv, center = center, scale = scale)    
   
   # do test set validation if provided
   if (!is.null(x.test) && !is.null(y.test))
   {
      x.test = as.matrix(x.test)
      y.test = as.matrix(y.test)
      
      if (is.null(colnames(y.test)))
         colnames(y.test) = paste('y', 1:ncol(y), sep = '')
      
      if (ncol(x.test) != ncol(x) || ncol(y.test) != ncol(y))
         stop('Calibration and data set should have the same number of variables!')
      
      if (nrow(x.test) != nrow(y.test))
         stop('Number of rows in x.test and y.test should be the same!')
      
      model$testres = predict.pls(model, x.test, y.test)
   }
   
   model$call = match.call()
   model$info = info
   
   class(model) = "pls"
   
   model
}

pls.cal = function(x, y, ncomp, center, scale, method = 'simpls', cv = F)
{   
   # Calibrates a PLS model.
   #
   # Arguments:
   #   x: a matrix with predictor values
   #   y: a vector with response values
   #   ncomp: number of components to calculate
   #   center: logical, center or not x and y data
   #   scale: logical, standardize or not y data
   #   method: a method to calculate PLS model
   #   cv: logical, calibrate model for cross-validation or not
   #
   # Returns:
   #   model: calibrated PLS model (list) 
   
   # center and scale data according to arguments
   x = prep.autoscale(as.matrix(x), center = center, scale = scale)
   y = prep.autoscale(as.matrix(y), center = center, scale = scale)
   
   # do PLS
   if (method == 'simpls')
      res = pls.simpls(x, y, ncomp)
   else
      stop('Method with this name is not supported!')

   # return a list with model parameters
   model = list(
      xloadings = res$xloadings,
      yloadings = res$yloadings,
      weights = res$weights,
      coeffs = regcoeffs(res$coeffs),
      method = method,
      xtnorm = sqrt(colSums(res$xscores ^ 2)/(nrow(res$xscores) - 1)),   
      ytnorm = sqrt(colSums(res$yscores ^ 2)/(nrow(res$yscores) - 1)),   
      xcenter = attr(x, 'prep:center'),
      xscale = attr(x, 'prep:scale'),
      ycenter = attr(y, 'prep:center'),
      yscale = attr(y, 'prep:scale'),
      ncomp = ncomp
   )
   
   model
}

pls.simpls = function(x, y, ncomp)
{
   # SIMPLS algorithm for calibrating PLS model.
   #
   # Arguments:
   #   x: a matrix with predictor values
   #   y: a matrix with response values
   #   ncomp: number of components to calculate
   #
   # Returns:
   #   res: a list with PLS calibration results 
   
   x = as.matrix(x)
   y = as.matrix(y)
   
   # get names for objects, variables and components
   objnames = rownames(x);
   prednames = colnames(x);
   respnames = colnames(y);
   compnames = paste('Comp', 1:ncomp)
   
   nobj = nrow(x)
   npred = ncol(x)
   nresp = ncol(y)
   
   # initial estimation
   A = t(x) %*% y
   M = t(x) %*% x
   C = diag(npred)
   
   # prepare space for results
   B = array(0, dim = c(npred, ncomp, nresp))
   W = matrix(0, nrow = npred, ncol = ncomp)
   P = matrix(0, nrow = npred, ncol = ncomp)
   Q = matrix(0, nrow = nresp, ncol = ncomp)
   
   # loop for each components
   for (n in 1:ncomp)
   {
      # get the dominate eigenvector of A'A
      e = eigen(t(A) %*% A)
      q = e$vectors[1:nresp]
      
      # calculate and store weights
      w = A %*% q
      c = t(w) %*% M %*% w
      w = w/sqrt(as.numeric(c))
      W[, n] = w
      
      # calculate and store x loadings
      p = M %*% w
      P[, n] = p
      
      # calculate and store y loadings
      q = t(A) %*% w
      Q[, n] = q
      
      v = C %*% p
      v = v/sqrt(as.numeric(t(v) %*% v))
      
      # calculate and store regression coefficients
      B[, n, ] = W[, 1:n, drop = FALSE] %*% t(Q[, 1:n, drop = FALSE])
      
      # recalculate matrices for the next compnonent
      C = C - v %*% t(v)
      M = M - p %*% t(p)
      A = C %*% A      
   }
   
   # calculate x and y scores
   U = y %*% Q 
   TT = x %*% (W %*% solve(t(P) %*% W))  
   
   # set names for all results
   colnames(Q) = colnames(B) = colnames(P) = colnames(W) = compnames      
   colnames(TT) = colnames(U) = compnames      
   rownames(P) = rownames(B) = rownames(W) = prednames
   rownames(TT) = rownames(U) = objnames
   rownames(Q) = dimnames(B)[[3]] = respnames
   
   res = list(
      coeffs = B,
      weights = W,
      xloadings = P,
      xscores = TT,
      yloadings = Q,
      yscores = U
   )   
}  

pls.crossval = function(model, x, y, cv, center = T, scale = F)
{
   # Cross-validation of a PLS model.
   #
   # Arguments:
   #   model: a PLS model (object of class pls) 
   #   x: a matrix with predictor values
   #   y: a vector with response values
   #   cv: number of segments for cross-validation (1 - full CV)
   #   center: logical, to center x and y data or not
   #   scale: logical, to standardize x data or not
   
   # Returns:
   #   res: a list with cross-validation results 
   
   x = as.matrix(x)
   y = as.matrix(y)
   
   ncomp = model$ncomp
   nobj = nrow(x)
   nvar = ncol(x)
   nresp = ncol(y)
   
   # get matrix with indices for cv segments
   idx = crossval(nobj, cv)
   
   seglen = ncol(idx);
   
   yp = array(0, dim = c(nobj, ncomp, nresp))
   Q2x = matrix(0, ncol = ncomp, nrow = nobj)   
   T2x = matrix(0, ncol = ncomp, nrow = nobj)   
   Q2y = matrix(0, ncol = ncomp, nrow = nobj)   
   T2y = matrix(0, ncol = ncomp, nrow = nobj)   
   
   # loop over segments
   for (i in 1:nrow(idx))
   {
      
      ind = na.exclude(idx[i,])
      if (length(ind) > 0)
      {   
         xc = x[-ind, , drop = F]
         yc = y[-ind, , drop = F]
         xt = x[ind, , drop = F]
         yt = y[ind, , drop = F]
         
         m = pls.cal(xc, yc, ncomp, center = center, scale = scale)               
         res = predict.pls(m, xt, yt, cv = T)
         
         xdist = ldecomp.getDistances(res$xscores, m$xloadings, res$xresiduals, 
                                      model$calres$xdecomp$tnorm)
         ydist = ldecomp.getDistances(res$xscores, m$yloadings, res$yresiduals, 
                                      model$calres$xdecomp$tnorm)
         
         yp[ind, , ] = res$yp
         Q2x[ind, ] = xdist$Q2        
         T2x[ind, ] = xdist$T2        
         Q2y[ind, ] = ydist$Q2        
         T2y[ind, ] = ydist$T2        
      }
      
   }  
   
   dimnames(yp) = list(rownames(x), colnames(model$coeffs$values), colnames(model$calres$y.ref))         
   
   res = plsres(yp, y.ref = y, ncomp.selected = model$ncomp.selected,
                xdecomp = ldecomp(totvar = model$calres$xdecomp$totvar,
                                  tnorm = model$calres$xdecomp$tnorm,
                                  ncomp.selected = model$ncomp.selected,
                                  Q2 = Q2x, T2 = T2x),
                ydecomp = ldecomp(totvar = model$calres$ydecomp$totvar,
                                  tnorm = model$calres$ydecomp$tnorm,
                                  ncomp.selected = model$ncomp.selected,
                                  Q2 = Q2y, T2 = T2y)   
   )
}

selectCompNum.pls = function(model, ncomp)
{
   # Set user defined value for optimal number of components.
   #
   # Arguments:
   #   model: a PLS model (object of class pls) 
   #   ncomp: number of components to use as optimal
   #
   # Returns:
   #   model: the same model but with selected components 
   
   if (ncomp > model$ncomp || ncomp < 0)
      stop('Wrong number of selected components!')
   
   model$ncomp.selected = ncomp      
   model$calres$ncomp.selected = ncomp
   
   if (!is.null(model$cvres)) 
      model$cvres$ncomp.selected = ncomp
   
   if (!is.null(model$testres)) 
      model$testres$ncomp.selected = ncomp
   
   model
}   

predict.pls = function(model, x, y = NULL, cv = F)
{   
   # Applies a PLS model to a new data.
   #
   # Arguments:
   #   model: a PLS model (object of class PLS)
   #   x: a matrix with predictor values
   #   y: a vector with response values (optional)
   #   cv: logical, is it prediction for cross-validation or not
   #
   # Returns:
   #   res: PLS results (object of class plsres)
   
   nresp = dim(model$coeffs$values)[3]
   nobj = nrow(x)
   
   # preprocess x and calculate scores, total and full variance
   x = as.matrix(x)
   x  = prep.autoscale(x, center = model$xcenter, scale = model$xscale)
   xscores = x %*% (model$weights %*% solve(t(model$xloadings) %*% model$weights))  
   xresiduals = x - xscores %*% t(model$xloadings)
   
   # make predictions
   yp = array(0, dim = c(nrow(x), model$ncomp, nresp))
   for (i in 1:nresp)
      yp[, , i] = x %*% model$coeffs$values[, , i]
   
   # if y is provided, calculate y residuals
   if (!is.null(y))
   {   
      yy = prep.autoscale(y, center = model$ycenter, scale = model$yscale)
      yscores = as.matrix(yy) %*% model$yloadings   

      ypp = yp[, ncol(yp), , drop = F]
      dim(ypp) = dim(yp)[c(1, 3)]
      yresiduals = ypp - yy
      
      dimnames(yp) = list(rownames(x), colnames(model$coeffs$values), colnames(y))   
   }
   else
   {
      dimnames(yp) = list(rownames(x), colnames(model$coeffs$values), colnames(model$calres$y.ref))         
   }   
   
   # unscale predicted y values
   if (model$yscale != F)
      for (i in 1:nresp)
         yp[, , i] = sweep(yp[, , i, drop = F], 2L, model$yscale[i], '*', check.margin = F)
   
   # uncenter predicted y values
   if (length(model$ycenter) > 1 || model$ycenter != F)
      for (i in 1:nresp)
         yp[, , i] = sweep(yp[, , i, drop = F], 2L, model$ycenter[i], '+', check.margin = F)
   
   if (cv == F)
   {   
      # normal predictions      
      # calculate PLS decomposition for x and y data
      # and return all results
      rownames(xscores) = rownames(x)
      colnames(xscores) = paste("Comp", 1:model$ncomp)
      xdecomp = ldecomp(xscores, model$xloadings, xresiduals, 
                        totvar = sum(x^2),
                        ncomp.selected = model$ncomp.selected)
      if (!is.null(y))
      {   
         dimnames(yscores) = dimnames(xscores)
         ydist = ldecomp.getDistances(xscores, model$yloadings, yresiduals)
         ydecomp = ldecomp(yscores, model$yloadings, yresiduals, sum(yy^2),
                           model$ytnorm, model$ncomp.selected,
                           ydist$T2, ydist$Q2)
      }
      else
      {
         ydecomp = NULL
      }      
      
      res = plsres(yp, y.ref = y, ncomp.selected = model$ncomp.selected, 
                   xdecomp = xdecomp, ydecomp = ydecomp)    
   }      
   else
   {
      # predictions for cross-validation      
      # just return predictions, scores and residuals
      # decomposition will be calculated in crossval() function
      res = list(
         yp = yp,
         xscores = xscores,
         yscores = yscores,
         xresiduals = xresiduals,
         yresiduals = yresiduals
      )
   }         
}  

plotRMSE.pls = function(model, ny = 1, type = 'b',
                        main = 'RMSE', xlab = 'Components', ylab = NULL, 
                        show.legend = T, show.labels = F, ...)
{
   # Makes RMSE plot.
   #
   # Arguments:
   #   model: a PLS model (object of class PLS)  
   #   ny: number of response variable to make the plot for
   #   type: type of the plot('b', 'l' or 'h')
   #   main: main plot title
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   show.legend: logical, show or not legend for the plot
   
   ncomp = model$ncomp
   
   if (is.null(ylab))
   {   
      if (nrow(model$calres$rmse) == 1)
         ylab = 'RMSE'
      else if (is.null(rownames(model$calres$rmse)))
         ylab = sprintf('RMSE (#%d)', ny)
      else         
         ylab = sprintf('RMSE (%s)', rownames(model$calres$rmse)[ny])
   }
   
   legend = c('cal')
   data = cbind(1:ncomp, model$calres$rmse[ny, ])   
   labels = mdaplot.formatValues(model$calres$rmse[ny, ])
   
   if (!is.null(model$cvres)) 
   { 
      data = cbind(data, model$cvres$rmse[ny, ])
      labels = cbind(labels, mdaplot.formatValues(model$cvres$rmse[ny, ]))
      legend = c(legend, 'cv')
   }   
   
   if (!is.null(model$testres)) 
   { 
      data = cbind(data, model$testres$rmse[ny, ])
      labels = cbind(labels, mdaplot.formatValues(model$testres$rmse[ny, ]))
      legend = c(legend, 'test')
   }     
   
   if (show.legend == F)
      legend = NULL
   
   if (show.labels == F)
      labels = NULL
   
   mdaplotg(data, type = type, legend = legend, labels = labels, 
            main = main, xlab = xlab, ylab = ylab, ...)
}

plotXVariance.pls = function(model, type = 'b',
                             main = 'X variance', xlab = 'Components', 
                             ylab = 'Explained variance, %', 
                             show.legend = T, ...)
{
   # Description.
   #
   # Arguments:
   #   :  
   #
   # Returns:
   #   : 
   
   plotVariance(model, decomp = 'xdecomp', type = type, main = main,
                xlab = xlab, ylab = ylab, show.legend = show.legend, ...)
}

plotYVariance.pls = function(model, type = 'b',
                             main = 'Y variance', xlab = 'Components', 
                             ylab = 'Explained variance, %', 
                             show.legend = T, ...)
{
   # Description.
   #
   # Arguments:
   #   :  
   #
   # Returns:
   #   : 
   
   plotVariance(model, decomp = 'ydecomp', type = type, main = main,
                xlab = xlab, ylab = ylab, show.legend = show.legend, ...)
}


plotXCumVariance.pls = function(model, type = 'b',
                                main = 'X cumulative variance', 
                                xlab = 'Components', ylab = 'Explained variance, %', 
                                show.legend = T, ...)
{
   # Description.
   #
   # Arguments:
   #   :  
   #
   # Returns:
   #   : 
   
   plotVariance(model, decomp = 'xdecomp', variance = 'cumexpvar', 
                type = type, main = main, xlab = xlab, ylab = ylab, 
                show.legend = show.legend, ...)   
}

plotYCumVariance.pls = function(model, type = 'b',
                                main = 'Y cumulative variance', 
                                xlab = 'Components', ylab = 'Explained variance, %', 
                                show.legend = T, ...)
{
   # Description.
   #
   # Arguments:
   #   :  
   #
   # Returns:
   #   : 
   
   plotVariance(model, decomp = 'ydecomp', variance = 'cumexpvar', 
                type = type, main = main, xlab = xlab, ylab = ylab, 
                show.legend = show.legend, ...)   
}

plotVariance.pls = function(model, decomp = 'xdecomp', variance = 'expvar',
                            type = 'b',
                            main = 'X variance', xlab = 'Components', 
                            ylab = 'Explained variance, %', 
                            show.labels = F,
                            show.legend = T, ...)
{
   # Description.
   #
   # Arguments:
   #   :  
   #
   # Returns:
   #   : 
   
   ncomp = model$ncomp
   
   legend = c('cal')
   data = cbind(1:ncomp, model$calres[[decomp]][[variance]])   
   labels = mdaplot.formatValues(model$calres[[decomp]][[variance]])
   
   if (!is.null(model$cvres)) 
   { 
      data = cbind(data, model$cvres[[decomp]][[variance]])
      labels = cbind(labels, mdaplot.formatValues(model$cvres[[decomp]][[variance]]))
      legend = c(legend, 'cv')
   }   
   
   if (!is.null(model$testres)) 
   { 
      data = cbind(data, model$testres[[decomp]][[variance]])
      labels = cbind(labels, mdaplot.formatValues(model$testres[[decomp]][[variance]]))
      legend = c(legend, 'test')
   }     
   
   if (show.legend == F)
      legend = NULL
   
   if (show.labels == F)
      labels = NULL
   
   mdaplotg(data, legend = legend, type = type, main = main, 
            xlab = xlab, ylab = ylab, labels = labels, ...)
}

plotXScores.pls = function(model, comp = c(1, 2), main = 'X scores',
                           xlab = NULL, ylab = NULL,
                           show.axes = F, 
                           show.labels = F, show.legend = T, ...)
{
   # Description.
   #
   # Arguments:
   #   :  
   #
   # Returns:
   #   : 
   
   ncomp = length(comp)
   
   if (ncomp < 1 || ncomp > 2)
      stop('The plot can be made for one or two components only!')
      
   if (is.null(main))
      main = sprintf('XY scores (ncomp = %d)', ncomp)
   
   if (is.null(ylab))
   {   
      if (ncomp == 2)
         ylab = colnames(model$calres$xdecomp$scores)[comp[2]]
      else
         ylab = colnames(model$calres$xdecomp$scores)[comp]
   }      
   
   if (is.null(xlab))
   {   
      if (ncomp == 2)
         xlab = colnames(model$calres$xdecomp$scores)[comp[1]]
      else
         xlab = 'Objects'
   }
   
   if (ncomp == 2)
   {   
      cdata = cbind(model$calres$xdecomp$scores[, comp[1], drop = F], 
                    model$calres$xdecomp$scores[, comp[2], drop = F])
   }
   else
   {
      ncrows = nrow(model$calres$xdecomp$scores);
      cdata = cbind(1:ncrows, 
                    model$calres$xdecomp$scores[, comp, drop = F])
   }   
   
   data = list(cdata = cdata)
   legend = c('cal')
   
   if (!is.null(model$testres)) 
   { 
      if (ncomp == 2)
      {   
         tdata = cbind(model$testres$xdecomp$scores[, comp[1], drop = F], 
                    model$testres$xdecomp$scores[, comp[2], drop = F])
      }
      else
      {
         tdata = cbind(ncrows + (1:nrow(model$testres$xdecomp$scores)), 
                       model$testres$xdecomp$scores[, comp, drop = F])         
      }   
      data$tdata = tdata
      legend = c(legend, 'test')      
   }   
   
   if (show.legend == F)
      legend = NULL
   
   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, legend = legend, show.labels = show.labels, show.lines = show.lines,
            main = main, xlab = xlab, ylab = ylab, ...)   
}  

plotXYScores.pls = function(model, comp = 1, main = NULL,
                            xlab = 'X scores', ylab = 'Y scores',
                            show.lines = T,
                            show.labels = F, show.legend = T, ...)
{
   # Description.
   #
   # Arguments:
   #   :  
   #
   # Returns:
   #   : 
   
   if (is.null(comp)) 
      comp = model$ncomp.selected   
   else if (comp <= 0 || comp > model$ncomp) 
      stop('Wrong component number!')
   
   if (is.null(main))
      main = sprintf('XY scores (comp %d)', comp)
   
   cdata = cbind(model$calres$xdecomp$scores[, comp, drop = F], 
                 model$calres$ydecomp$scores[, comp, drop = F])
   colnames(cdata) = c('X scores', 'Y scores')
   data = list(cdata = cdata)
   legend = c('cal')
   
   if (!is.null(model$testres)) 
   { 
      tdata = cbind(model$testres$xdecomp$scores[, comp, drop = F], 
                    model$testres$ydecomp$scores[, comp, drop = F])
      data$tdata = tdata
      legend = c(legend, 'test')      
   }   
   
   if (show.legend == F)
      legend = NULL
   
   if (show.lines == T)
      show.lines = c(0, 0)
   
   mdaplotg(data, legend = legend, show.labels = show.labels, main = main,
            xlab = xlab, ylab = ylab, show.lines = show.lines, ...)   
}  

plotPredictions.pls = function(model, ncomp = NULL, ny = 1, main = NULL, xlab = NULL,
                               ylab = NULL, colmap = 'default', col = NULL,
                               legend.position = 'topleft', show.lines = T, show.labels = F,
                               show.legend = T, ...)
{
   # Makes predictions (predicted vs measured values) plot.
   #
   # Arguments:
   #   model: a PLS model (object of class pls)  
   #   ncomp: number of components to make the plot for (default: ncomp.selected)
   #   ny: which y variables to make the plot for
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   colmap: color map for marker colors
   #   col: color for markers
   #   legend.position: position of the plot legend
   #   show.lines: logical, show target lines for points or not
   #   show.labels: logical, show labels for data objects or not
   #   show.legend: logical, show plot legend or not
   
   # set default values for main title, x and y axis labels
   if (is.null(main))
   {   
      if (is.null(ncomp))
         main = 'Predictions'
      else
         main = sprintf('Predictions (ncomp = %d)', ncomp)
   }
   
   if (is.null(xlab))
   {
      if (!is.null(colnames(model$calres$y.ref)))
         xlab = sprintf('%s, measured', colnames(model$calres$y.ref)[ny])
      else
         xlab = 'y, measured'
   }  
   
   if (is.null(ylab))
   {
      if (!is.null(colnames(model$calres$y.ref)))
         ylab = sprintf('%s, predicted', colnames(model$calres$y.ref)[ny])
      else
         ylab = 'y, predicted'
   }  
   
   # check number of components
   if (is.null(ncomp)) 
      ncomp = model$ncomp.selected
   else if (ncomp <= 0 || ncomp > model$ncomp) 
      stop('Wrong value for number of components!')
   
   # make list with groups of measured and predicted values for each type of result
   cdata = cbind(model$calres$y.ref[, ny], model$calres$y.pred[, ncomp, ny])
   colnames(cdata) = c('y, measured', 'y, predicted')
   rownames(cdata) = rownames(model$calres$y.ref)
   legend = c('cal')
   data = list(cdata = cdata)
   
   if (!is.null(model$cvres)) 
   { 
      cvdata = cbind(model$cvres$y.ref[, ny], model$cvres$y.pred[, ncomp, ny])
      colnames(cvdata) = c('y, measured', 'y, predicted')
      rownames(cvdata) = rownames(model$cvres$y.ref)
      legend = c(legend, 'cv')
      data$cvdata = cvdata
   }   
   
   if (!is.null(model$testres)) 
   { 
      testdata = cbind(model$testres$y.ref[, ny], model$testres$y.pred[, ncomp, ny])
      colnames(testdata) = c('y, measured', 'y, predicted')
      rownames(testdata) = rownames(model$testres$y.ref)
      legend = c(legend, 'test')
      data$testdata = testdata
   }   
   
   if (show.legend == F)
      legend = NULL
   
   mdaplotg(data, legend = legend, show.labels = show.labels, main = main, 
            colmap = colmap, xlab = xlab, ylab = ylab, col = col, 
            legend.position = legend.position, ...)  
   
   if (show.lines == T)
      mdaplot.showRegressionLine(data, colmap = colmap, col = col)
}

plotYResiduals.pls = function(model, ncomp = NULL, ny = 1, type = 'p', 
                              main = NULL, 
                              ylab = NULL,
                              xlab = NULL,
                              show.lines = T,
                              show.labels = F,
                              show.legend = T,
                              ...)
{
   # Makes Y residuals plot.
   #
   # Arguments:
   #   model: a PLS model (object of class pls) 
   #   ncomp: number of components to make the plot for
   #   ny: number of response variable to make the plot for
   #   type: type of the plot
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   show.lines: logical, show confidence limits as lines or not 
   #   show.labels: logical, show labels for data objects or not
   #   show.legend: logical, show plot legend or not
   
   if (is.null(main))
   {   
      if (is.null(ncomp))
         main = 'Y residuals'
      else
         main = sprintf('Y residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(xlab))
   {
      if (type == 'p')
      {   
         if (ncol(model$calres$y.ref) == 1)
            xlab = 'y, measured'
         else
            xlab = sprintf('%s, measured', colnames(model$calres$y.ref)[ny])
      }
      else
      {
         xlab = 'Objects'
      }   
   }  

   if (is.null(ylab))
   {
      if (ncol(model$calres$y.ref) == 1)
         ylab = 'y residuals'
      else
         ylab = sprintf('y residuals (%s)', colnames(model$calres$y.ref)[ny])
   }  
   
   if (is.null(ncomp)) 
      ncomp = model$ncomp.selected
   else if (ncomp <= 0 || ncomp > model$ncomp) 
      stop('Wrong value for number of components!')
   
   if (type == 'p')
   {   
      cdata = cbind(model$calres$y.ref[, ny], model$calres$y.pred[, ncomp, ny] - model$calres$y.ref[, ny])
      colnames(cdata) = c('y values', 'y residuals')
   }
   else
   {   
      nobj = nrow(model$calres$y.ref)
      cdata = cbind(1:nrow(model$calres$y.ref), model$calres$y.pred[, ncomp, ny] - model$calres$y.ref[, ny])
      colnames(cdata) = c('Objects', 'y residuals')
   }
   rownames(cdata) = rownames(model$calres$y.pred)
   legend = c('cal')
   data = list(cdata = cdata)
   
   if (!is.null(model$cvres)) 
   { 
      if (type == 'p')
      {   
         cvdata = cbind(model$cvres$y.ref[, ny], 
                        model$cvres$y.pred[, ncomp, ny] - model$cvres$y.ref[, ny])
         colnames(cvdata) = c('y values', 'y residuals')
      }
      else
      {   
         cvdata = cbind(nobj + (1:nrow(model$cvres$y.ref)), 
                        model$cvres$y.pred[, ncomp, ny] - model$cvres$y.ref[, ny])
         colnames(cvdata) = c('Objects', 'y residuals')
         nobj = nobj + nrow(model$cvres$y.ref)
      }      
      rownames(cvdata) = rownames(model$cvres$y.pred)
      legend = c(legend, 'cv')
      data$cvdata = cvdata
   }   
   
   if (!is.null(model$testres)) 
   { 
      if (type == 'p')
      {   
         tdata = cbind(model$testres$y.ref[, ny], 
                       model$testres$y.pred[, ncomp, ny] - model$testres$y.ref[, ny])
         colnames(tdata) = c('y values', 'y residuals')
      }
      else
      {   
         tdata = cbind(nobj + (1:nrow(model$testres$y.ref)), 
                       model$testres$y.pred[, ncomp, ny] - model$testres$y.ref[, ny])
         colnames(tdata) = c('Objects', 'y residuals')
      }
      rownames(tdata) = rownames(model$testres$y.pred)
      legend = c(legend, 'test')
      data$tdata = tdata
   }   
   
   if (show.legend == F)
      legend = NULL
   
   if (show.lines == T)
      show.lines = c(NA, 0)
   
   mdaplotg(data, legend = legend, type = type, show.labels = show.labels, xlab = xlab, ylab = ylab, 
            main = main, show.lines = show.lines, ...)      
}

plotRegcoeffs.pls = function(model, ncomp = NULL, ny = 1, main = NULL, ylab = NULL, ...)
{
   # Makes regression coefficients plot.
   #
   # Arguments:
   #   model: a PLS model (object of class pls)  
   #   ncomp: number of components to make the plot for
   #   ny: number of response variable to make the plot for
   #   main: main title for the plot
   
   if (is.null(main))
   {   
      if (is.null(ncomp))
         main = 'Regression coefficients'
      else
         main = sprintf('Regression coefficients (ncomp = %d)', ncomp)
   }
   
   if (is.null(ylab))
   {   
      if (ncol(model$calres$y.ref) == 1)
         ylab = 'Coefficients'
      else
         ylab = sprintf('Coefficients (%s)', colnames(model$calres$y.ref)[ny])
   }
   
   if (is.null(ncomp)) 
      ncomp = model$ncomp.selected
   else if (ncomp <= 0 || ncomp > model$ncomp) 
      stop('Wrong value for number of components!')
   
   plot(model$coeffs, ncomp = ncomp, ny = ny, main = main, ylab = ylab, ...)
}

plotXLoadings.pls = function(model, comp = c(1, 2), type = 'p', main = 'X loadings', 
                             ylab = NULL, xlab = NULL, show.axes = F, ...)
{
   # Makes X loadings plot.
   #
   # Arguments:
   #   model: a PLS model (object of class pls)  
   #   comp: which components to make the plot for
   #   type: type of the plot
   #   main: main title for the plot
   
   ncomp = length(comp)
   
   if (min(comp) < 1 || max(comp) > model$ncomp)
      stop('Wrong components number!')
   
   if (type == 'p' && ncomp != 2)
      stop('Scatter plot can be made only for two components!')

   if (type == 'p')
   {
      data = cbind(model$xloadings[, comp[1], drop = F],
                   model$xloadings[, comp[2], drop = F])      
      legend = NULL
   }   
   else
   {
      data = cbind(1:nrow(model$xloadings),
                   model$xloadings[, comp, drop = F])  
      legend = colnames(model$xloadings[, comp, drop = F])
   }   
   
   if (is.null(ylab))
   {   
      if (type == 'p')
         ylab = colnames(model$xloadings)[comp[2]]
      else   
         ylab = 'Loadings'
   }

   if (is.null(xlab))
   {   
      if (type == 'p')
         xlab = colnames(model$xloadings)[comp[1]]
      else   
         xlab = 'Variables'
   }

   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, main = main, type = type, ylab = ylab, xlab = xlab, legend = legend, 
            show.lines = show.lines, ...)
}

plotXYLoadings.pls = function(model, comp = c(1, 2), main = 'XY loadings', 
                             ylab = NULL, xlab = NULL, show.axes = F, ...)
{
   # Makes XY loadings plot.
   #
   # Arguments:
   #   model: a PLS model (object of class pls)  
   #   comp: which components to make the plot for
   #   main: main title for the plot
   
   if (length(comp) != 2)
      stop('This plot can be made for only two components!')
   
   xdata = cbind(model$xloadings[, comp[1], drop = F],
                 model$xloadings[, comp[2], drop = F])      
   ydata = cbind(model$yloadings[, comp[1], drop = F],
                 model$yloadings[, comp[2], drop = F])      
   data = list(xdata = xdata, ydata = ydata)
   
   legend = c('X', 'Y')
   
   if (is.null(ylab))
      ylab = colnames(model$xloadings)[comp[2]]
   
   if (is.null(xlab))
      xlab = colnames(model$xloadings)[comp[1]]
   
   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, main = main, type = 'p', ylab = ylab, xlab = xlab, legend = legend, 
            show.lines = show.lines, ...)
}

plotXResiduals.pls = function(model, ncomp = NULL, 
                              main = NULL, xlab = 'T2', ylab = 'Q2',
                              show.labels = F, show.legend = T, ...)
{
   # Makes X residuals plot (T2 vs Q2).
   #
   # Arguments:
   #   model: a PLS model (object of class PLS)  
   #   ncomp: number of components to make the plot for
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   show.labels: logical, show labels for data objects or not 
   #   show.legend: logical, show legend for the plot or not
   
   if (is.null(main))
   {   
      if (is.null(ncomp))
         main = 'X residuals'
      else
         main = sprintf('X residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp)) 
      ncomp = model$ncomp.selected
   else if (ncomp <= 0 || ncomp > model$ncomp) 
      stop('Wrong value for number of components!')
   
   cdata = cbind(model$calres$xdecomp$T2[, ncomp], model$calres$xdecomp$Q2[, ncomp])
   
   colnames(cdata) = c('T2', 'Q2')
   rownames(cdata) = rownames(model$calres$xdecomp$scores)
   
   data = list(cdata = cdata)
   legend = NULL
   
   if (!is.null(model$testres))
   {
      tdata = cbind(model$testres$xdecomp$T2[, ncomp], model$testres$xdecomp$Q2[, ncomp])      
      colnames(tdata) = c('T2', 'Q2')
      rownames(tdata) = rownames(model$testres$xdecomp$scores)
      
      data$tdata = tdata
      legend = c('cal', 'test')
   }      
   
   mdaplotg(data, main = main, xlab = xlab, ylab = ylab,
            show.labels = show.labels, legend = legend, ...)
} 

## makes a plot with regression results ##
plot.pls = function(model, ncomp = NULL, ny = 1, show.legend = T, show.labels = F)
{
   # Makes a plot for PLS model overview.
   #
   # Arguments:
   #   model: a PLS model (object of class pls)  
   #   ncomp: number of components to show the summary for (default: ncomp.selected)
   #   ny: number of response variable to show the summary for
   #   show.legend: logical, show the plot legend or not
   #   show.labels: logical, show data object labels or not
   
   if (!is.null(ncomp) && (ncomp <= 0 || ncomp > model$ncomp)) 
      stop('Wrong value for number of components!')
   
   par(mfrow = c(2, 2))      
   plotXResiduals(model, ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   
   plotRegcoeffs(model, ncomp = ncomp, ny = ny, show.labels = F)   
   plotRMSE(model, ny = ny, show.legend = show.legend)   
   plotPredictions(model, ncomp = ncomp, ny = ny, show.labels = show.labels, show.legend = show.legend)   
   par(mfrow = c(1, 1))
}

summary.pls = function(model, ncomp = NULL, ny = NULL)
{
   # Shows summary information for PLS model.
   #
   # Arguments:
   #   model: a PLS model (object of class pls)  
   #   ncomp: number of components to show the summary for (default: ncomp.selected)
   #   ny: number of response variable to show the summary for
   
   if (is.null(ncomp))
      ncomp = model$ncomp.selected
   else if (ncomp <= 0 || ncomp > model$ncomp)
      stop('Wrong value for number of components!')
   
   if (is.null(ny))
      ny = 1:ncol(model$calres$y.ref)
   
   cat('\nPLS model (class pls) summary\n')
   cat('\nPerformance and validation:\n')
   cat(sprintf('Number of selected components: %d\n', ncomp))
   
   for (y in ny)
   {   
      if (ncol(model$calres$y.ref) > 1)
         cat(sprintf('\nResponse variable #%d (%s)\n', y, colnames(model$calres$y.ref)[y]))
      
      data = as.matrix(model$calres, ncomp = ncomp, ny = y)
      rownames(data) = 'Cal'
      
      if (!is.null(model$cvres))
      {
         data = rbind(data, as.matrix(model$cvres, ncomp = ncomp, ny = y))      
         rownames(data)[2] = 'CV'
      }
      
      if (!is.null(model$testres))
      {
         data = rbind(data, as.matrix(model$testres, ncomp = ncomp, ny = y))
         rownames(data)[nrow(data)] = 'Test'
      }   
      
      data = data[, -c(1, 3), drop = F]
      
      data[, 1:2] = round(data[, 1:2], 2)      
      data[, 4:5] = round(data[, 4:5], 2)  
      data[, 3] = mdaplot.formatValues(data[, 3], round.only = T)
      data[, 6] = round(data[, 6], 4)      
      data[, 7] = round(data[, 7], 2)      
      
      print(data)
   }   
   cat('\n')
}

print.pls = function(model, ...)
{
   # Prints information about the PLS model.
   #
   # Arguments:
   #   model: a PLS model (object of class pls)  
   
   cat('\nPLS model (class pls)\n')
   cat('\nCall:\n')
   print(model$call)
   cat('\nMajor fields:\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$coeffs - vector with regression coefficients\n')
   cat('$xloadings - vector with x loadings\n')
   cat('$yloadings - vector with Y loadings\n')
   cat('$weights - vector with weights\n')
   cat('$calres - results for calibration set\n')
   if (!is.null(model$cvres))
   {
      cat('$cvres - results for cross-validation\n')      
   }   
   if (!is.null(model$testres))
   {
      cat('$testres - results for test set\n')      
   }   
   cat('\nTry summary(model) and plot(model) to see the model performance.\n')   
}
