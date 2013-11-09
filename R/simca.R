## class and methods for SIMCA classification ##

simca = function(data, classname, ncomp = 20, center = T, scale = F, cv = NULL, test.data = NULL, 
                 test.c = NULL, alpha = 0.05, method = 'svd', info = '', ...)
{
   # Calibrate and validate a SIMCA classification model for one class
   #
   # Arguments:
   #   data: a matrix with data values
   #   classname: a short text with name of the class
   #   ncomp: maximum number of components to calculate
   #   center: logical, mean center the data values or not 
   #   scale: logical, standardize the data values or not 
   #   cv: number of segments for random cross-validation (1 - for full CV)
   #   test.data: a matrix with data values for test set validation
   #   test.c: a matrix with class values for test set validation
   #   alpha: a significance level for Q2 residuals
   #   method: method to find principal component space (only SVD is supported so far)
   #   info: a text with information about the model
   #
   # Returns:
   #   model: a SIMCA model (object of simca class)
   
   data = as.matrix(data)
   
   # check if data has missing values
   if (sum(is.na(data)) > 0)
   {
      warning('Data has missing values, will try to fix using pca.mvreplace.')
      data = pca.mvreplace(data, center = center, scale = scale)
   }   
   
   if (!is.character(classname))
      stop('Argument "classname" must be a text!')
   
   if (length(classname) > 20)
      stop('Argument "classname" must have up to 20 symbols!')
   
   # correct maximum number of components
   ncomp = min(ncomp, ncol(data), nrow(data) - 1)
   
   # calibrate model  
   model = pca.cal(data, ncomp, center = center, scale = scale, method = method)
   model$ncomp = ncomp
   model$ncomp.selected = model$ncomp
   model$classname = classname
   model$info = info
   model$alpha = alpha
   
   # calculate and assign limit values for T2 and Q2 residuals
   lim = ldecomp.getResLimits(model$eigenvals, nrow(data), model$ncomp.selected, model$alpha)
   model$T2lim = lim$T2lim
   model$Q2lim = lim$Q2lim   

   model$call = match.call()   
   class(model) = c("simca", "pca")
   
   # apply model to calibration set
   model$calres = predict.simca(model, data, c.ref = rep(1, nrow(data)))
   model$modpower = model$calres$modpower
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = simca.crossval(model, data, cv, center = center, scale = scale)
   
   # apply model to test set if provided
   if (!is.null(test.data))
   {
      if (is.null(test.c))
         test.c = matrix(T, nrow(test.data), 1)
      
      model$testres = predict.simca(model, test.data, c.ref = test.c)
   }
   
   model
}

predict.simca = function(model, data, c.ref = NULL, cv = F)
{
   # Apply the SIMCA model to a new data set
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   data: a matrix with new data values
   #   c.ref: reference class values for the data set (optional)
   #   cv: logical, is it CV predictions or not 
   #
   # Returns:
   #   res: results of SIMCA classification (object of simcares class)
   
   data = as.matrix(data)
   pres = predict.pca(model, data, cv)     
   pres$Q2lim = model$Q2lim
   pres$T2lim = model$T2lim
   
   c.pred = simca.classify(model, pres)
   
   # check c.ref values and add dimnames
   if (!is.null(c.ref))
   {   
      if (is.logical(c.ref))
         c.ref = c.ref * 2 - 1
      c.ref = as.matrix(c.ref)
      rownames(c.ref) = rownames(data)
      colnames(c.ref) = model$classname
   } 
   
   cres = classres(c.pred, c.ref = c.ref)
   res = simcares(pres, cres)
   
   res
}

simca.classify = function(model, res)
{
   # Make classification based on calculated Q2 and T2 statistics
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   res: results of projection of new data to PC space 
   #
   # Returns:
   #   c.pred: predicted class values
   
   ncomp = model$ncomp
   c.pred = array(0, dim = c(nrow(res$Q2), ncomp, 1))
   dimnames(c.pred) = list(rownames(res$Q2), paste('Comp', 1:ncomp), model$classname)
   
   for (i in 1:ncomp)
   {
      c.pred[, i, 1] = 
         (res$T2[, i] - model$T2lim[1, i]) < 0.00000001 & 
         (res$Q2[, i] - model$Q2lim[1, i]) < 0.00000001
   }   
   c.pred = c.pred * 2 - 1
   
   c.pred
}  

simca.crossval = function(model, data, cv, center = T, scale = F)
{
   # Cross-validation of a SIMCA model
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)  
   #   data: a matrix with data values
   #   cv: number of segments for cross-validation (1 for full CV)
   #   center: logical, mean center data values or not
   #   scale: logical, standardize data values or not
   #
   # Returns:
   #   res: results of cross-validation (object of class simcares) 
   
   ncomp = model$ncomp   
   nobj = nrow(data)
   nvar = ncol(data)
   
   # get matrix with indices for cv segments
   idx = crossval(nobj, cv)
   
   seglen = ncol(idx);
   
   Q2 = matrix(0, ncol = ncomp, nrow = nobj)   
   T2 = matrix(0, ncol = ncomp, nrow = nobj)   
   c.pred = array(0, dim = c(nobj, ncomp, 1))
   c.ref = matrix(1, ncol = 1, nrow = nobj)
   
   # loop over segments
   for (i in 1:nrow(idx))
   {
      ind = na.exclude(idx[i,])
      
      if (length(ind) > 0)
      {   
         datac = data[-ind, , drop = F]
         datat = data[ind, , drop = F]
         
         m = pca.cal(datac, ncomp, center, scale)               
         lim = ldecomp.getResLimits(m$eigenvals, nrow(datac), ncomp, model$alpha)
         m$T2lim = lim$T2lim
         m$Q2lim = lim$Q2lim
         m$ncomp = ncomp         
         
         res = predict.pca(m, datat, cv = T)
         Q2[ind, ] = res$Q2
         T2[ind, ] = res$T2
         c.pred[ind, , ] = simca.classify(m, res)
      }
   }  
   
   rownames(Q2) = rownames(T2) = rownames(c.pred) = rownames(c.ref) = rownames(data)
   colnames(Q2) = colnames(T2) = colnames(c.pred) = colnames(model$scores)
   colnames(c.ref) = model$classname
   
   pres = pcares(NULL, NULL, NULL, model$calres$totvar, model$tnorm, model$ncomp.selected, T2, Q2)
   cres = classres(c.pred, c.ref = c.ref)   
   res = simcares(pres, cres)
   
   res
}  

plotPredictions.simca = function(model, which = 'crossval', ncomp = NULL, type = 'h',
                                 main = NULL, xlab = NULL, 
                                 ylab = NULL, legend = NULL, ...)
{
   # makes a plot with predictions
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   which: which results to show (not used so far)
   #   ncomp: number of components to show the predictions for (default - selected)
   #   type: plot type
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   legend: legend strings
   #
   
   nc = 1
   
   if (is.null(main))
   {
      if (is.null(ncomp))
         main = 'Predictions'
      else
         main = sprintf('Predictions (ncomp = %d)', ncomp)      
   }   
   
   ncomp = getSelectedComponents(model, ncomp)
   
   if (ncomp > model$ncomp)
      stop('Wrong value for argument "ncomp"!')
   
   if (is.null(ylab))
      ylab = 'Predicted values'
   
   if (is.null(xlab))
      xlab = 'Objects'   
   
   if (is.null(model$p.pred))
   {   
      y = model$c.pred[ , ncomp, nc, drop = F]
   }   
   else
   {   
      y = model$p.pred[ , ncomp, nc, drop = F]
   }
      
   if (is.null(model$calres$p.pred))
   {   
      if (!is.null(model$cvres$c.pred))
      {   
         y = cbind(model$calres$c.pred[ , ncomp, nc, drop = F], 
                   model$cvres$c.pred[ , ncomp, nc, drop = F])
         if (is.null(legend))
            legend = c('Cal', 'CV')          
      }
      else             
         y = model$calres$c.pred[ , ncomp, nc, drop = F]
   }   
   else
   {   
      if (!is.null(model$cvres$p.pred))
      {   
         y = cbind(model$calres$p.pred[ , ncomp, nc, drop = F], 
                   model$cvres$p.pred[ , ncomp, nc, drop = F])
         if (is.null(legend))
            legend = c('Cal', 'CV')          
      }
      else             
         y = model$calres$p.pred[ , ncomp, nc, drop = F]
   }
   
   data = cbind(1:nrow(y), y)
   rownames(data) = rownames(y)
   mdaplotg(data, type = type, main = main, xlab = xlab, ylab = ylab, legend = legend, ...)            
}

# plotSpecificity.simca = function(model, type = 'h', legend = NULL, main = 'Specificity', 
#                                  xlab = 'Components', ylab = '', ylim = c(0, 1.15), ...)
# {
#    # makes a plot with specificity values vs. number of PCs
#    #
#    # Arguments:
#    #   model: a SIMCA model (object of class simca)
#    #   type: plot type
#    #   legend: legend strings
#    #   main: main title for the plot
#    #   xlab: label for x axis
#    #   ylab: label for y axis
#    #   ylim: limits for y axis
#    #
#    
#    plotPerformance(model, which = 'specificity', type = type, legend = legend, main = main,
#                    xlab = xlab, ylab = ylab, ylim = ylim, ...)
# }

plotSensitivity.simca = function(model, type = 'h', legend = NULL, main = 'Sensitivity', 
                                 xlab = 'Components', ylab = '', ylim = c(0, 1.15), ...)
{
   # makes a plot with sensitivity values vs. number of PCs
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   type: plot type
   #   legend: legend strings
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   ylim: limits for y axis
   #

   plotPerformance(model, which = 'sensitivity', type = type, legend = legend, main = main,
                   xlab = xlab, ylab = ylab, ylim = ylim, ...)
}

plotPerformance.simca = function(model, which = 'specificity', type = 'h', legend = NULL, 
                                 main = 'Specificity', xlab = 'Components', ylab = '', 
                                 ylim = c(0, 1.15), ...)
{
   # makes a plot with specificity or sensitivity values vs. number of PCs
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   which: which parameter to make the plot for ('specificity', 'sensitivity')
   #   type: plot type
   #   legend: legend strings
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   ylim: limits for y axis
   #
   
   data = cbind(1:model$ncomp, model$calres[[which]][1, ])
   labels = mdaplot.formatValues(model$calres[[which]][1, ])
   legend_str = 'cal'
   
   if (!is.null(model$cvres))
   {
      data = cbind(data, model$cvres[[which]][1, ])   
      labels = cbind(labels, mdaplot.formatValues(model$cvres[[which]][1, ]))
      legend_str = c(legend_str, 'cv')
   }   
   
   if (!is.null(model$testres))
   {
      data = cbind(data, model$testres[[which]][1, ])   
      labels = cbind(labels, mdaplot.formatValues(model$testres[[which]][1, ]))
      legend_str = c(legend_str, 'test')
   }
   
   if (!is.null(legend))
      legend_str = legend
   
   mdaplotg(data, type = type, main = main, xlab = xlab, ylab = ylab, legend = legend_str,
            ylim = ylim, labels = labels, ...)
}

plotModellingPower.simca = function(model, ncomp = NULL, type = 'h', legend = NULL, 
                                    xlab = 'Variables', ylab = '', 
                                    main = NULL, ...)
{
   # makes a plot with modelling power of variables
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   ncomp: number of components to make the plot for
   #   type: plot type
   #   legend: legend strings
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #
   
   if (is.null(main))
   {
      if (is.null(ncomp))
         main = 'Modelling power'
      else
         main = sprintf('Modelling power (ncomp = %d)', ncomp)      
   }   

   ncomp = getSelectedComponents(model, ncomp)
   
   nvar = nrow(model$modpower)
   if (is.null(type))
   {   
      if (nvar < 20)
         type = 'h'
      else
         type = 'l'
   }
   
   data = cbind(1:nvar, model$modpower[, ncomp, drop = F])
   mdaplot(data, type = type, xlab  = xlab, ylab = ylab, main = main, ...)
}   

plot.simca = function(model, ncomp = NULL, ...)
{
   # makes a plot with overview of a SIMCA model
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   ncomp: number of components to make the plot for
   #
   
   par(mfrow = c(2, 2))
   plotScores(model, ...)
   plotModellingPower(model, ncomp = ncomp, main = 'Modelling power', 
                      show.labels = ncol(model$modpower) < 10, ...)
   plotResiduals(model, main = 'Residuals', ncomp = ncomp, ...)
   plotCumVariance(model, ...)
   par(mfrow = c(1, 1))
}  

print.simca = function(model, ...)
{
   cat('\nSIMCA one class model (class simca)\n')
   
   cat('\nCall:\n')
   print(model$call)
   
   cat('\nMajor fields:\n')   
   cat('$classname - name of the class\n')
   cat('$alpha - significance level\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$loadings - matrix with loadings\n')
   cat('$eigenvals - eigenvalues for components\n')
   cat('$center - values for centering data\n')
   cat('$scale - values for scaling data\n')
   cat('$info - information about the model\n')
   cat('$cv - number of segments for cross-validation\n')
   cat('$calres - results (scores, etc) for calibration set\n')
   
   if (!is.null(model$cvres))
   {
      cat('$cvres - results for cross-validation\n')      
   }   
   if (!is.null(model$testres))
   {
      cat('$testres - results for test set\n')      
   }       
}  

summary.simca = function(model, ...)
{
   ncomp = model$ncomp   
   cat(sprintf('\nSIMCA model for class "%s" summary\n\n', model$classname))
   cat(sprintf('Significance level (alpha): %.2f\n\n', model$alpha))
   data = cbind(round(model$calres$expvar, 2),
                round(model$calres$cumexpvar, 2),
                round(model$calres$sensitivity[1, ], 2)
   )   
   colnames(data) = c('Expvar', 'Cumexpvar', 'Sens (cal)')
   
   if (!is.null(model$cvres))
   {
      cnames = colnames(data)
      data = cbind(data,
                   round(model$cvres$sensitivity[1, ], 2)
      )
      colnames(data) = c(cnames, 'Sens (cv)')
   }   
   
   if (!is.null(model$testres))
   {
      cnames = colnames(data)
      if (is.null(model$testres$specificity[1, ]) || min(model$testres$specificity[1, ]) == 1)
      {
         data = cbind(data,
                      round(model$testres$sensitivity[1, ], 2)
         )
         colnames(data) = c(cnames, 'Sens (test)')         
      }
      else  
      {   
         data = cbind(data,
                      round(model$testres$specificity[1, ], 2),
                      round(model$testres$sensitivity[1, ], 2)
         )
         colnames(data) = c(cnames, 'Spec (test)', 'Sens (test)')
      }
   }   
   
   print(data)   
}  
