## class and methods for SIMCA classification ##

simca = function(x, classname, ncomp = 20, center = T, scale = F, cv = NULL, x.test = NULL, 
                 c.test = NULL, alpha = 0.05, method = 'svd', info = '', ...)
{
   # Calibrate and validate a SIMCA classification model for one class
   #
   # Arguments:
   #   x: a matrix with data values
   #   classname: a short text with name of the class
   #   ncomp: maximum number of components to calculate
   #   center: logical, mean center the data values or not 
   #   scale: logical, standardize the data values or not 
   #   cv: number of segments for random cross-validation (1 - for full CV)
   #   x.test: a matrix with data values for test set validation
   #   c.test: a matrix with class values for test set validation
   #   alpha: a significance level for Q2 residuals
   #   method: method to find principal component space (only SVD is supported so far)
   #   info: a text with information about the model
   #
   # Returns:
   #   model: a SIMCA model (object of simca class)
   
   x = as.matrix(x)
   
   # check if data has missing values
   if (sum(is.na(x)) > 0)
   {
      warning('Data has missing values, will try to fix using pca.mvreplace.')
      x = pca.mvreplace(x, center = center, scale = scale)
   }   
   
   if (!is.character(classname))
      stop('Argument "classname" must be a text!')
   
   if (length(classname) > 20)
      stop('Argument "classname" must have up to 20 symbols!')
   
   # correct maximum number of components
   ncomp = min(ncomp, ncol(x), nrow(x) - 1)
   
   # calibrate model  
   model = pca.cal(x, ncomp, center = center, scale = scale, method = method)
   model$ncomp = ncomp
   model$ncomp.selected = model$ncomp
   model$nclasses = 1
   model$classname = classname
   model$info = info
   model$alpha = alpha
   
   # calculate and assign limit values for T2 and Q2 residuals
   lim = ldecomp.getResLimits(model$eigenvals, nrow(x), model$ncomp.selected, model$alpha)
   model$T2lim = lim$T2lim
   model$Q2lim = lim$Q2lim   

   model$call = match.call()   
   class(model) = c("simca", "classmodel", "pca")
   
   # apply model to calibration set
   model$calres = predict.simca(model, x, c.ref = rep(1, nrow(x)))
   model$modpower = model$calres$modpower
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = simca.crossval(model, x, cv, center = center, scale = scale)
   
   # apply model to test set if provided
   if (!is.null(x.test))
   {
      if (is.null(c.test))
         c.test = matrix(T, nrow(x.test), 1)
      
      model$testres = predict.simca(model, x.test, c.ref = c.test)
   }
   
   model
}

predict.simca = function(model, x, c.ref = NULL, cv = F)
{
   # Apply the SIMCA model to a new data set
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   x: a matrix with new data values
   #   c.ref: reference class values for the data set (optional)
   #   cv: logical, is it CV predictions or not 
   #
   # Returns:
   #   res: results of SIMCA classification (object of simcares class)
   
   x = as.matrix(x)
   
   if (is.null(rownames(x)))
      rownames(x) = 1:nrow(x)
   
   if (is.null(colnames(x)))
      colnames(x) = paste('v', 1:ncol(x), sep = '')
   
   pres = predict.pca(model, x, cv)     
   pres$Q2lim = model$Q2lim
   pres$T2lim = model$T2lim
   
   c.pred = simca.classify(model, pres)
   
   # check c.ref values and add dimnames
   if (!is.null(c.ref))
   {   
      if (is.character(c.ref))
         c.ref = c.ref == model$classname
      
      if (is.logical(c.ref))
         c.ref = c.ref * 2 - 1
      
      c.ref = as.matrix(c.ref)
      rownames(c.ref) = rownames(x)
      colnames(c.ref) = model$classname
   } 
   
   cres = classres(c.pred, c.ref = c.ref, ncomp.selected = model$ncomp.selected)
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

simca.crossval = function(model, x, cv, center = T, scale = F)
{
   # Cross-validation of a SIMCA model
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)  
   #   x: a matrix with data values
   #   cv: number of segments for cross-validation (1 for full CV)
   #   center: logical, mean center data values or not
   #   scale: logical, standardize data values or not
   #
   # Returns:
   #   res: results of cross-validation (object of class simcares) 
   
   ncomp = model$ncomp   
   nobj = nrow(x)
   nvar = ncol(x)
   
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
         x.cal = x[-ind, , drop = F]
         x.val = x[ind, , drop = F]
         
         m = pca.cal(x.cal, ncomp, center, scale)               
         lim = ldecomp.getResLimits(m$eigenvals, nrow(x.cal), ncomp, model$alpha)
         m$T2lim = lim$T2lim
         m$Q2lim = lim$Q2lim
         m$ncomp = ncomp         
         m$classname = model$classname
         m$ncomp.selected = model$ncomp.selected

         res = predict.pca(m, x.val, cv = T)
         Q2[ind, ] = res$Q2
         T2[ind, ] = res$T2
         
         c.pred[ind, , ] = simca.classify(m, res)
      }
   }  
  
   dimnames(c.pred) = list(rownames(x), colnames(model$loadings), model$classname)
   rownames(Q2) = rownames(T2) = rownames(c.pred) = rownames(c.ref) = rownames(x)
   colnames(Q2) = colnames(T2) = colnames(c.pred) = colnames(model$loadings)
   pres = pcares(NULL, NULL, NULL, model$calres$totvar, model$tnorm, model$ncomp.selected, T2, Q2)
   cres = classres(c.pred, c.ref = c.ref)   
   res = simcares(pres, cres)
   
   res
}  

getSelectedComponents.simca = function(obj, ncomp = NULL)
{
   if (is.null(ncomp))
   {   
      if (is.null(obj$ncomp.selected))
         ncomp = 1
      else
         ncomp = obj$ncomp.selected
   }   

   ncomp
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
   
   if (!is.null(model$info))
      cat(sprintf('Info: %s\n', model$info))
   
   cat(sprintf('Significance level (alpha): %.2f\n', model$alpha))
   cat(sprintf('Selected number of components: %d\n\n', model$ncomp.selected))
   
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
