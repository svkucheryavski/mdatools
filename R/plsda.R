# class and methods for Partial Least Squares Discriminant Analysis #

plsda = function(x, c, ncomp = 15, center = T, scale = F, cv = NULL, 
               x.test = NULL, c.test = NULL, method = 'simpls', alpha = 0.05, info = '')
{
   # Calibrate and validate a PLS-DA model.
   #
   # Arguments:
   #   x: a matrix with predictor values    
   #   c: a vector with class values 
   #   ncomp: maximum number of components to calculate
   #   center: logical, center or not x and y data
   #   scale: logical, standardize or not x data
   #   cv: number of segments for cross-validation (1 - full CV)
   #   x.test: a matrix with predictor values for test set validation
   #   c.test: a vector with class values for test set validation
   #   method: a method to calculate PLS model
   #   alpha: a sigificance limit for Q2 values
   #   info: a short string with information about the model
   #
   # Returns:
   #   model: a PLS-DA model (object of class pls) 
   
   x = as.matrix(x)
   c = as.matrix(c)
   y = plsda.getReferenceValues(c)

   # correct maximum number of components
   ncomp = min(ncol(x), nrow(x) - 1, ncomp)
   
   # build a model and apply to calibration set
   model = pls.cal(x, y, ncomp, center = center, scale = scale, method = method)
   model$ncomp.selected = ncomp
   model$alpha = alpha   
   model$classnames = unique(c)
   model$nclasses = length(model$classnames)
   model$calres = predict.plsda(model, x, c)
   model = selectCompNum.pls(model, ncomp)
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = plsda.crossval(model, x, c, cv, center = center, scale = scale)    
   
   # do test set validation if provided
   if (!is.null(x.test) && !is.null(c.test))
   {
      x.test = as.matrix(x.test)
      model$testres = predict.plsda(model, x.test, c.test)
   }
   
   model$call = match.call()
   model$info = info
   
   class(model) = c("plsda", "classmodel", "pls")
   
   model
}

#' Reference values for PLS-DA
#' 
#' @description
#' Generates matrix with reference y values (-1 and +1) for a 
#' vector with class values
#' 
#' @param c 
#' vector with class values (discrete)
#' @param classnames
#' vector with names for the classes
#' 
#' @return
#' the generated matrix with one column for each class
#' 
plsda.getReferenceValues = function(c, classnames = NULL)
{
   # generate matrix with y values
   
   if (is.null(classnames))
      classnames = unique(c)

   nclasses = length(classnames)
   y = matrix(-1, nrow = length(c), ncol = nclasses)
   
   for (i in 1:nclasses)
      y[c == classnames[i], i] = 1
   
   rownames(y) = rownames(c)
   colnames(y) = classnames

   y
}

#' PLS-DA predictions
#' 
#' @description
#' Applies PLS-DA model to a new data set
#' 
#' @param object
#' a PLS-DA model (object of class \code{plsda})
#' @param x
#' a matrix with x values (predictors)
#' @param c
#' a vector with reference class values
#' @param cv
#' logical, are predictions for cross-validation or not
#' @param ...
#' other arguments
#' 
#' @return
#' PLS-DA results (an object of class \code{plsdares})
#'
#' @details
#' See examples in help for \code{\link{plsda}} function.
#'  
predict.plsda = function(object, x, c = NULL, cv = F, ...)
{   
   y = plsda.getReferenceValues(c, object$classnames)
   plsres = predict.pls(object, x, y)
   cres = classify.plsda(object, plsres$y.pred, c)

   res = plsdares(plsres, cres)

   res
}  

#' PLS-DA classification
#' 
#' @description
#' Converts PLS predictions of y values to predictions of classes
#' 
#' @param model
#' a PLS-DA model (object of class \code{plsda})
#' @param y
#' a matrix with predicted y values
#' @param c.ref
#' a vector with reference class values
#' 
#' @return
#' Classification results (an object of class \code{classres})
#'
#' @details
#' This is a service function for PLS-DA class, do not use it manually.
#'  
classify.plsda = function(model, y, c.ref = NULL)
{
   c.pred = array(-1, dim(y))
   c.pred[y >= 0] = 1
   dimnames(c.pred) = list(rownames(y), paste('Comp', 1:model$ncomp), model$classnames)
   cres = classres(c.pred, c.ref = c.ref, p.pred = y, ncomp.selected = model$ncomp.selected)

   cres
}

#' Cross-validation of a PLS-DA model
#' 
#' @description
#' Does the cross-validation of a PLS-DA model
#' 
#' @param model
#' a PLS-DA model (object of class \code{plsda})
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param c
#' a vetor with c values (classes from calibration set)
#' @param cv
#' number of segments (if cv = 1, full cross-validation will be used)
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#'
#' @return
#' object of class \code{plsdares} with results of cross-validation
#'  
plsda.crossval = function(model, x, c, cv, center = T, scale = F)
{
   y = plsda.getReferenceValues(c, model$classnames)
   plsres = pls.crossval(model, x, y, cv, center, scale)
   cres = classify.plsda(model, plsres$y.pred, c)

   res = plsdares(plsres, cres)

   res
}

#' Model overview plot for PLS-DA
#' 
#' @description
#' Shows a set of plots (x residuals, regression coefficients, misclassification ratio and predictions) 
#' for PLS-DA model.
#' 
#' @param x
#' a PLS-DA model (object of class \code{plsda})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param nc
#' which class to show the plots
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{plsda}} function.
#' 
plot.plsda = function(x, ncomp = NULL, nc = 1, show.legend = T, show.labels = F, ...)
{
   model = x
   
   if (!is.null(ncomp) && (ncomp <= 0 || ncomp > model$ncomp)) 
      stop('Wrong value for number of components!')
   
   par(mfrow = c(2, 2))      
   plotXResiduals(model, ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   
   plotRegcoeffs(model, ncomp = ncomp, ny = nc, show.labels = show.labels)   
   plotMisclassified(model, nc = nc, show.legend = show.legend)   
   
   if (!is.null(model$cvres))
      plotPredictions(model, res = 'cvres', ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   
   else
      plotPredictions(model, ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   

   par(mfrow = c(1, 1))
}

#' Summary method for PLS-DA model object
#' 
#' @method summary plsda
#' @S3method summary plsda
#'
#' @description
#' Shows some statistics for the model.
#' 
#' @param object
#' a PLS-DA model (object of class \code{plsda})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param nc
#' which class to show the summary for (if NULL, will be shown for all)
#' @param ...
#' other arguments
#' 
summary.plsda = function(object, ncomp = NULL, nc = NULL, ...)
{
   obj = object
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   else if (ncomp <= 0 || ncomp > obj$ncomp)
      stop('Wrong value for number of components!')
   
   if (is.null(nc))
      nc = 1:obj$nclasses
   
   cat('\nPLS-DA model (class plsda) summary statistics\n\n')
   cat(sprintf('Number of selected components: %d\n', ncomp))
   
   if (!is.null(obj$info))
      cat(sprintf('Info: %s\n', obj$info))
      
   for (n in nc)
   {   
      cat(sprintf('\nClass #%d (%s)\n', n, obj$classnames[n]))
      
      data = as.matrix(obj$calres, ncomp = ncomp, nc = n)
      rownames(data) = 'Cal'
      
      if (!is.null(obj$cvres))
      {
         data = rbind(data, as.matrix(obj$cvres, ncomp = ncomp, nc = n))      
         rownames(data)[2] = 'CV'
      }
      
      if (!is.null(obj$testres))
      {
         data = rbind(data, as.matrix(obj$testres, ncomp = ncomp, nc = n))
         rownames(data)[nrow(data)] = 'Test'
      }   
      
      data[, 1:4] = round(data[, 1:4], 2)      
      print(data)
   }   
   cat('\n')
}

#' Print method for PLS-DA model object
#' 
#' @method print plsda
#' @S3method print plsda
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a PLS-DA model (object of class \code{plsda})
#' @param ...
#' other arguments
#' 
print.plsda = function(x, ...)
{
   obj = x
   
   cat('\nPLS-DA model (class plsda)\n')
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$coeffs - vector with regression coefficients\n')
   cat('$xloadings - vector with x loadings\n')
   cat('$yloadings - vector with Y loadings\n')
   cat('$weights - vector with weights\n')
   cat('$calres - results for calibration set\n')
   if (!is.null(obj$cvres))
   {
      cat('$cvres - results for cross-validation\n')      
   }   
   if (!is.null(obj$testres))
   {
      cat('$testres - results for test set\n')      
   }   
   cat('\nTry summary(model) and plot(model) to see the model performance.\n')   
}
