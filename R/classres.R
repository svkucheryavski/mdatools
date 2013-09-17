#  methods for classification results #
classres = function(c.pred, c.ref = NULL, ncomp.selected = NULL)
{
   # Class for storing and visualising of classification results 
   #
   # Arguments:
   #   c.pred: vector or matrix with predicted class
   #   c.ref: vector with reference (true) class
   #   ncomp.selected: if c.pred calculated for different components, which to use as default
   #
   # Returns:
   # a list (object of "classres" class) with following fields
   #   c.pred: a matrix with predicted class
   #   c.ref: a vector with reference (true) class
   #   rmse: root mean squared error for predicted vs measured values
   #   slope: slope for predicted vs measured values
   #   r2: coefficient of determination for predicted vs measured values
   #   bias: bias for predicted vs measured values
   #   rpd: RPD values
   
   
   if (!is.null(c.ref))
   {
      c.ref = as.matrix(c.ref)      
      obj = classres.getClassificationPerformance(c.ref, c.pred)
      obj$c.ref = c.ref
   }
   else
   {
      obj = list()
   }   
   
   obj$c.pred = c.pred
   obj$ncomp.selected = ncomp.selected
   
   obj$call = match.call()   
   class(obj) = "classres"
   
   obj
}


plotPredictions.classres = function(obj, nc = 1, ncomp = NULL, main = 'Predictions', 
                                  xlab = NULL, ylab = NULL, 
                                  show.line = T, ...)
{
   # Makes plot with predicted class values vs object number 
   # for selected class 
   #
   # Arguments:
   #   obj: object of "classres" class
   #   nc: number of class c to make the plot for
   #   ncomp: which column of c.pred to use
   #   main: main title of the plot
   #   xlab: text for x axis label
   #   ylab: text for y axis label
   #   show.line: logical, show or not separation line for classes
   #   col: color for plot and target line
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   
   if (is.null(ylab))
   {   
      if (~is.null(dim(obj$c.pred)) && dim(obj$c.pred)[3] > 1)
         ylab = sprintf('%s, predicted', dimnames(obj$c.pred)[[3]][nc])
      else
         ylab = 'class, predicted'
   }
   
   if (is.null(xlab))
      xlab = 'Objects'
      
   data = cbind(1:nrow(obj$c.pred[, , nc]), 
                obj$c.meas[, ncomp, nc, drop = F],
                obj$c.pred[, ncomp, nc, drop = F])     
   mdaplotg(data, type = 'p', main = main, xlab = xlab, ylab = ylab, ...)
}

classres.getClassificationPerformance = function(c.ref, c.pref)
{
   # Calculates and returns a list with classification performance parameters
   #
   # Arguments:
   #  c.ref - matrix with reference values (nobj x nclasses)
   #  c.pred - array with predicted values (nobj x ncomp x nclasses) 
   #
   
   nclasses = ncol(c.ref)
   ncomp = ncol(c.pred)
   
   tp = matrix(0, nrow = nclasses, ncol = ncomp)
   fp = matrix(0, nrow = nclasses, ncol = ncomp)
   fn = matrix(0, nrow = nclasses, ncol = ncomp)
   specificity = matrix(0, nrow = nclasses, ncol = ncomp)
   sensitivity = matrix(0, nrow = nclasses, ncol = ncomp)
   f1 = matrix(0, nrow = nclasses, ncol = ncomp)
   
   for (i in 1:nclasses)         
      fn[i, ] = colSums((c.ref[, i, drop = F] == 1) & (c.pred[, , i] == 0))
      fp[i, ] = colSums((c.ref[, i, drop = F] == 0) & (c.pred[, , i] == 1))
      tp[i, ] = colSums((c.ref[, i, drop = F] == 1) & (c.pred[, , i] == 1))
      sensitivity[i, ] = tp[i, ] / (tp[i, ] + fn[i, ])
      specificity[i, ] = tp[i, ] / (tp[i, ] + fp[i, ])
      f1[i, ] = 2 * sensitivity[i, ] * specificity[i, ] / (sensitivity[i, ] + specificity[i, ])
   end
   
   obj = list()
   obj$fn = fn
   obj$fp = fp
   obj$tp = tp
   obj$sensitivity = sensitivity
   obj$specificity = specificity
   obj$f1 = f1
   
   obj
}   

plot.classres = function(obj, ny = 1, ...)
{
   # Plot method for "regres" objects
   #
   # Arguments:
   #   obj: object of "regres" class
   #   ny: number of response variable y to make the plot for   
   
   plotPredictions.classres(obj, ny = ny, ...)
}   

as.matrix.classres = function(obj, ncomp = NULL, ny = 1)
{
   # as.matrix method for "classres" objects
   #
   # Arguments:
   #   obj: object of "regres" class
   #   ny: number of response variable y to make the plot for   
   
   if (!is.null(obj$c.ref))
   {  
      if (is.null(ncomp))
         res = cbind(obj$tp[nc, ], obj$fp[nc, ], obj$fn[nc, ], obj$specificity[nc, ], obj$sensitivity[nc, ])
      else
         res = cbind(obj$tp[nc, ncomp], obj$fp[nc, ncomp], obj$fn[nc, ncomp], 
                     obj$specificity[nc, ncomp], obj$sensitivity[nc, ncomp])
      
      colnames(res) = c('TP', 'FN', 'FP', 'Specificity', 'Sensitivity')
   }
   else
   {
      res = NULL
   }   
   
   res
}

print.classres = function(obj, ...)
{
   # Print method for "classres" object
   #
   # Arguments:
   #   obj: object of "classres" class
   
   cat('\nClassification results (class classres)\n')
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$c.pred - matrix or vector with predicted class\n')
   if (!is.null(obj$c.ref))
   {   
      cat('$c.ref - vector with reference (true) class\n')
      cat('$tp - number of true positives\n')
      cat('$fp - number of false positives\n')
      cat('$fn - number of false negatives\n')
      cat('$specificity - specificity\n')
      cat('$sensitivity - sensitivity\n')
   }
   
   if (ncol(obj$c.pred) > 1)   
      cat('$ncomp.selected - number of selected components for each class\n')
}   

summary.calres = function(obj, ncomp = NULL, nc = NULL, ...)
{
   # Summary method for "calres" object
   #
   # Arguments:
   #   obj: object of "calres" class
   #   ncomp: which column of c.pred to use
   
   cat('\nClassiciation results (class classres) summary\n')
   if (!is.null(obj$c.ref))
   {         
      if (is.null(ncomp))
         ncomp = obj$ncomp.selected
      
      if (is.null(nc))
         nc = 1:ncol(obj$c.ref)
      
      if (!is.null(ncomp))
         cat(sprintf('\nNumber of selected components: %d\n\n', ncomp))
      
      for (i in nc)
      {   
         cat(sprintf('\nClass "%s":\n', colnames(obj$c.ref)[i]))
         res = as.matrix.regres(obj, nc = i, ncomp = ncomp)
         rownames(res) = ncomp
         print(res)
      }
      
   }      
   else
   {
      cat('No reference data provided to calculate prediction performance.')
   }   
}   

