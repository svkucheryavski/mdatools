## class and methods for SIMCA classification ##

simca = function(x, classname, ncomp = 15, center = T, scale = F, cv = NULL, x.test = NULL, 
                 c.test = NULL, alpha = 0.05, method = 'svd', info = '')
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

#' SIMCA predictions
#' 
#' @description
#' Applies SIMCA model to a new data set
#' 
#' @param object
#' a SIMCA model (object of class \code{simca})
#' @param x
#' a matrix with x values (predictors)
#' @param c.ref
#' a vector with reference class values
#' @param cv
#' logical, are predictions for cross-validation or not
#' @param ...
#' other arguments
#' 
#' @return
#' SIMCA results (an object of class \code{simcares})
#'
#' @details
#' See examples in help for \code{\link{simca}} function.
#'  
predict.simca = function(object, x, c.ref = NULL, cv = F, ...)
{
   x = as.matrix(x)
   
   if (is.null(rownames(x)))
      rownames(x) = 1:nrow(x)
   
   if (is.null(colnames(x)))
      colnames(x) = paste('v', 1:ncol(x), sep = '')
   
   pres = predict.pca(object, x, cv)     
   pres$Q2lim = object$Q2lim
   pres$T2lim = object$T2lim
   
   c.pred = simca.classify(object, pres)
   
   # check c.ref values and add dimnames
   if (!is.null(c.ref))
   {   
      if (is.character(c.ref))
         c.ref = c.ref == object$classname
      
      if (is.logical(c.ref))
         c.ref = c.ref * 2 - 1
      
      c.ref = as.matrix(c.ref)
      rownames(c.ref) = rownames(x)
      colnames(c.ref) = object$classname
   } 
   
   cres = classres(c.pred, c.ref = c.ref, ncomp.selected = object$ncomp.selected)
   res = simcares(pres, cres)
   
   res
}

#' SIMCA classification
#' 
#' @description
#' Make classification based on calculated T2 and Q2 values and corresponding limits
#' 
#' @param model
#' a SIMCA model (object of class \code{simca})
#' @param res
#' results of projection data to PCA space
#' 
#' @return
#' vector with predicted class values (\code{c.pred})
#'
#' @details
#' This is a service function for SIMCA class, do not use it manually.
#'  
simca.classify = function(model, res)
{
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

#' Cross-validation of a SIMCA model
#' 
#' @description
#' Does the cross-validation of a SIMCA model
#' 
#' @param model
#' a SIMCA model (object of class \code{simca})
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param cv
#' number of segments (if cv = 1, full cross-validation will be used)
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#'
#' @return
#' object of class \code{simcares} with results of cross-validation
#'  
simca.crossval = function(model, x, cv, center = T, scale = F)
{
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

#' Modelling power plot for SIMCA model
#' 
#' @description
#' Shows a plot with modelling power values for each predictor
#' 
#' @param obj
#' a SIMCA model (object of class \code{simca})
#' @param ncomp
#' number of components to show the values for
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
plotModellingPower.simca = function(obj, ncomp = NULL, type = 'h', main = NULL, 
                                    xlab = 'Variables', ylab = '', ...)
{
   if (is.null(main))
   {
      if (is.null(ncomp))
         main = 'Modelling power'
      else
         main = sprintf('Modelling power (ncomp = %d)', ncomp)      
   }   

   ncomp = getSelectedComponents(obj, ncomp)
   
   nvar = nrow(obj$modpower)
   if (is.null(type))
   {   
      if (nvar < 20)
         type = 'h'
      else
         type = 'l'
   }
   
   data = cbind(1:nvar, obj$modpower[, ncomp, drop = F])
   mdaplot(data, type = type, xlab  = xlab, ylab = ylab, main = main, ...)
}   

#' Model overview plot for SIMCA
#' 
#' @description
#' Shows a set of plots for SIMCA model.
#' 
#' @param x
#' a SIMCA model (object of class \code{simca})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
plot.simca = function(x, ncomp = NULL, ...)
{
   obj = x
   
   par(mfrow = c(2, 2))
   plotScores(obj, ...)
   plotModellingPower(obj, ncomp = ncomp, main = 'Modelling power', 
                      show.labels = ncol(obj$modpower) < 10, ...)
   plotResiduals(obj, main = 'Residuals', ncomp = ncomp, ...)
   plotCumVariance(obj, ...)
   par(mfrow = c(1, 1))
}  


#' Summary method for SIMCA model object
#' 
#' @method summary simca
#' @S3method summary simca
#'
#' @description
#' Shows performance statistics for the model.
#' 
#' @param object
#' a SIMCA model (object of class \code{simca})
#' @param ...
#' other arguments
#' 
summary.simca = function(object, ...)
{
   obj = object
   
   ncomp = obj$ncomp   
   cat(sprintf('\nSIMCA model for class "%s" summary\n\n', obj$classname))
   
   if (!is.null(obj$info))
      cat(sprintf('Info: %s\n', obj$info))
   
   cat(sprintf('Significance level (alpha): %.2f\n', obj$alpha))
   cat(sprintf('Selected number of components: %d\n\n', obj$ncomp.selected))
   
   data = cbind(round(obj$calres$expvar, 2),
                round(obj$calres$cumexpvar, 2),
                round(obj$calres$sensitivity[1, ], 2)
   )   
   colnames(data) = c('Expvar', 'Cumexpvar', 'Sens (cal)')
   
   if (!is.null(obj$cvres))
   {
      cnames = colnames(data)
      data = cbind(data,
                   round(obj$cvres$sensitivity[1, ], 2)
      )
      colnames(data) = c(cnames, 'Sens (cv)')
   }   
   
   if (!is.null(obj$testres))
   {
      cnames = colnames(data)
      if (is.null(obj$testres$specificity[1, ]) || min(obj$testres$specificity[1, ]) == 1)
      {
         data = cbind(data,
                      round(obj$testres$sensitivity[1, ], 2)
         )
         colnames(data) = c(cnames, 'Sens (test)')         
      }
      else  
      {   
         data = cbind(data,
                      round(obj$testres$specificity[1, ], 2),
                      round(obj$testres$sensitivity[1, ], 2)
         )
         colnames(data) = c(cnames, 'Spec (test)', 'Sens (test)')
      }
   }   
   
   print(data)   
}  

#' Print method for SIMCA model object
#' 
#' @method print simca
#' @S3method print simca
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a SIMCA model (object of class \code{simca})
#' @param ...
#' other arguments
#' 
print.simca = function(x, ...)
{
   obj = x
   
   cat('\nSIMCA one class model (class simca)\n')
   
   cat('\nCall:\n')
   print(obj$call)
   
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
   
   if (!is.null(obj$cvres))
   {
      cat('$cvres - results for cross-validation\n')      
   }   
   if (!is.null(obj$testres))
   {
      cat('$testres - results for test set\n')      
   }       
}  

