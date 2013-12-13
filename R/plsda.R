# class and methods for Partial Least Squares Discriminant Analysis #

plsda = function(x, c, ncomp = 20, center = T, scale = F, cv = NULL, 
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

predict.plsda = function(model, x, c = NULL, cv = F)
{   
   # Applies a PLS-DA model to a new data.
   #
   # Arguments:
   #   model: a PLS model (object of class PLS)
   #   x: a matrix with predictor values
   #   c: a vector with reference class values (optional)
   #   cv: logical, is it prediction for cross-validation or not
   #
   # Returns:
   #   res: PLS results (object of class plsres)

   y = plsda.getReferenceValues(c, model$classnames)
   plsres = predict.pls(model, x, y)
   cres = classify.plsda(model, plsres$y.pred, c)

   res = plsdares(plsres, cres)

   res
}  

classify.plsda = function(model, y, c.ref = NULL)
{
   c.pred = array(-1, dim(y))
   c.pred[y >= 0] = 1
   dimnames(c.pred) = list(rownames(y), paste('Comp', 1:model$ncomp), model$classnames)
   cres = classres(c.pred, c.ref = c.ref, p.pred = y, ncomp.selected = model$ncomp.selected)

   cres
}

plsda.crossval = function(model, x, c, cv, center = T, scale = F)
{
   y = plsda.getReferenceValues(c, model$classnames)
   plsres = pls.crossval(model, x, y, cv, center, scale)
   cres = classify.plsda(model, plsres$y.pred, c)

   res = plsdares(plsres, cres)

   res
}

plot.plsda = function(model, ncomp = NULL, nc = 1, show.legend = T, show.labels = F)
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
   plotRegcoeffs(model, ncomp = ncomp, ny = nc, show.labels = show.labels)   
   plotMisclassified(model, nc = nc, show.legend = show.legend)   
   
   if (!is.null(model$cvres))
      plotPredictions(model, res = 'cvres', ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   
   else
      plotPredictions(model, ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   

   par(mfrow = c(1, 1))
}

summary.plsda = function(model, ncomp = NULL, nc = NULL)
{
   # Shows summary information for PLSDA model.
   #
   # Arguments:
   #   model: a PLSDA model (object of class plsda)  
   #   ncomp: number of components to show the summary for (default: ncomp.selected)
   #   nc: number of class to show the summary for
   
   if (is.null(ncomp))
      ncomp = model$ncomp.selected
   else if (ncomp <= 0 || ncomp > model$ncomp)
      stop('Wrong value for number of components!')
   
   if (is.null(nc))
      nc = 1:model$nclasses
   
   cat('\nPLS-DA model (class plsda) summary statistics\n\n')
   cat(sprintf('Number of selected components: %d\n', ncomp))
   
   if (!is.null(model$info))
      cat(sprintf('Info: %s\n', model$info))
      
   for (n in nc)
   {   
      cat(sprintf('\nClass #%d (%s)\n', n, model$classnames[n]))
      
      data = as.matrix(model$calres, ncomp = ncomp, nc = n)
      rownames(data) = 'Cal'
      
      if (!is.null(model$cvres))
      {
         data = rbind(data, as.matrix(model$cvres, ncomp = ncomp, nc = n))      
         rownames(data)[2] = 'CV'
      }
      
      if (!is.null(model$testres))
      {
         data = rbind(data, as.matrix(model$testres, ncomp = ncomp, nc = n))
         rownames(data)[nrow(data)] = 'Test'
      }   
      
      data[, 1:4] = round(data[, 1:4], 2)      
      print(data)
   }   
   cat('\n')
}

print.plsda = function(model, ...)
{
   # Prints information about the PLS-DA model.
   #
   # Arguments:
   #   model: a PLS-DA model (object of class plsda)  
   
   cat('\nPLS-DA model (class plsda)\n')
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
