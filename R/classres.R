#  methods for classification results #
classres = function(c.pred, c.ref = NULL, p.pred = NULL, ncomp.selected = NULL)
{
   # Class for storing and visualising of classification results 
   #
   # Arguments:
   #   c.pred: vector or matrix with predicted class
   #   c.ref: vector with reference (true) class
   #   p.pred: vector with predicted probabilities for the classification decision
   #   ncomp.selected: if c.pred calculated for different components, which to use as default
   #
   # Returns:
   # a list (object of "classres" class) with following fields
   #   c.pred: a matrix with predicted class
   #   c.ref: a matrix with reference (true) class
   #   p.pred: a matrix with predicted probabilities
   #   tp: a matrix with true positives for each class and component
   #   fp: a matrix with false positives for each class and component
   #   fn: a matrix with false negatives for each class and component
   #   specificity: a matrix with specificity values
   #   sensitivity: a matrix with sensitivity values
   
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
   obj$p.pred = p.pred
   
   obj$ncomp = dim(c.pred)[2]
   if (is.null(ncomp.selected) && obj$ncomp > 1)
      obj$ncomp.selected = obj$ncomp
   else
      obj$ncomp.selected == NULL
   
   obj$call = match.call()   
   class(obj) = "classres"
   
   obj
}


classres.getClassificationPerformance = function(c.ref, c.pred)
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
   {   
      fn[i, ] = colSums((c.ref[, i] == 1) & (c.pred[, , i, drop = F] == -1))
      fp[i, ] = colSums((c.ref[, i] == -1) & (c.pred[, , i, drop = F] == 1))
      tp[i, ] = colSums((c.ref[, i] == 1) & (c.pred[, , i, drop = F] == 1))
      
      sensitivity[i, ] = tp[i, ] / (tp[i, ] + fn[i, ])
      specificity[i, ] = tp[i, ] / (tp[i, ] + fp[i, ])
      f1[i, ] = 2 * sensitivity[i, ] * specificity[i, ] / (sensitivity[i, ] + specificity[i, ])
   }

   rownames(fn) = rownames(fp) = rownames(tp) = rownames(f1) = 
      rownames(sensitivity) = rownames(specificity) = colnames(c.ref)
   colnames(fn) = colnames(fp) = colnames(tp) = colnames(f1) = 
      colnames(sensitivity) = colnames(specificity) = dimnames(c.pred)[[2]]
   
   obj = list()
   obj$fn = fn
   obj$fp = fp
   obj$tp = tp
   obj$sensitivity = sensitivity
   obj$specificity = specificity
   obj$f1 = f1
   
   obj
}   

showPredictions.classres = function(obj, ncomp = NULL)
{
   #  Shows table with predictions for selected or user specified number of components
   #
   # Arguments:
   #   obj: object of "classres" class
   #   ncomp: for which components to show the results for

   ncomp = getSelectedComponents(obj, ncomp)
   
   pred = obj$c.pred[ , ncomp, , drop = F]
   dim(pred) = c(dim(obj$c.pred)[1], dim(obj$c.pred)[3]);
   dimnames(pred) = list(dimnames(obj$c.pred)[[1]], dimnames(obj$c.pred)[[3]])
   
   print(pred)
}  

plotSensitivity.classres = function(obj, ...)
{
   plotPerformance(obj, which = 'sensitivity', ...)
}   

plotSpecificity.classres = function(obj, ...)
{
   plotPerformance(obj, which = 'specificity', ...)
}   

plotPerformance.classres = function(obj, nc = 1, which = 'all', main = NULL, xlab = 'Complexity', 
                                    type = 'h', legend = NULL, ylim = NULL, ...)
{
   # Makes plot with snesitivity and specificity vs number of components
   #
   # Arguments:
   #   obj: object of "classres" class
   #   nc: number of class to make the plot for
   #   which: which of the two performance parameters to show
   #   main: main title of the plot
   #   xlab: text for x axis label
   #   ylab: text for y axis label
   #   legend: legend for the plot items
   #   ylim: limits for y axis
   
   
   if (length(obj$sensitivity) < 2)
      warning('There are values only for one component!')
   else
   {   
      if (dim(obj$c.pred)[3] > 1)
         titlestr = sprintf('for class "%s"', dimnames(obj$c.pred)[[3]][nc]) 
      else
         titlestr = '';
      
      if (which == 'all')
      {   
         if (is.null(obj$sensitivity) || is.null(obj$specificity))
            stop('The sensitivity and specificity values are not available!')

         if (is.null(main))
            main = paste('Prediction performance', titlestr);
      
         if (is.null(legend))
            legend = c('sensitivity', 'specificity')
      
         if (is.null(ylim))
            ylim = c(0, 1.1)
      
         data = cbind(1:length(obj$sensitivity[nc, ]), obj$sensitivity[nc, ], obj$specificity[nc, ])
         labels = round(cbind(obj$sensitivity[nc, ], obj$specificity[nc, ]), 2)
      }
      else if (which == 'specificity')
      {
         if (is.null(obj$specificity))
            stop('The specificity values are not available!')
         
         if (is.null(main))
            main = paste('Specificity', titlestr)
                  
         if (is.null(ylim))
            ylim = c(0, 1.1)
         
         data = cbind(1:length(obj$specificity[nc, ]), obj$specificity[nc, ])         
         labels = matrix(round(obj$specificity[nc, ], 2), ncol = 1)
      }   
      else if (which == 'sensitivity')
      {
         if (is.null(obj$specificity))
            stop('The sensitivity values are not available!')
         
         if (is.null(main))
            main = paste('Sensitivity', titlestr)
         
         if (is.null(ylim))
            ylim = c(0, 1.1)
         
         data = cbind(1:length(obj$sensitivity[nc, ]), obj$sensitivity[nc, ])                  
         labels = matrix(round(obj$sensitivity[nc, ], 2), ncol = 1)
      }   
      
      mdaplotg(data, type = type, legend = legend, main = main, xlab = xlab, ylim = ylim,
               labels = labels, ...)
   }
}

plotPredictions.classres = function(obj, nc = 1, type = 'h', ncomp = NULL, main = NULL, 
                                    xlab = NULL, ylab = NULL, legend = NULL, ylim = c(-1.2, 1.2),
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
   #   legend: legend for the plot groups
   #   show.line: logical, show or not separation line for classes
   
   if (is.null(ylab))
      ylab = 'Predicted values'
   
   if (is.null(xlab))
      xlab = 'Objects'
   
   if (is.null(main))
   {   
      if (!is.null(dim(obj$c.pred)) && dim(obj$c.pred)[3] > 1)
         main = sprintf('Predictions for class "%s"', dimnames(obj$c.pred)[[3]][nc])
      else
         main = 'Predictions'
      
      if (!is.null(ncomp))
         main = sprintf('%s (ncomp = %d)', main, ncomp)
   }
   
   ncomp = getSelectedComponents(obj, ncomp)
      
   if (is.null(obj$p.pred))
   {   
      y = obj$c.pred[ , ncomp, nc, drop = F]
   }   
   else
   {   
      y = obj$p.pred[ , ncomp, nc, drop = F]
   }

   obj_idx = 1:length(y)
   
   if (!is.null(obj$c.ref))
   {   
      if (is.null(legend))
         legend = c(colnames(obj$c.ref)[nc], 'Others')
   
      members_idx = obj$c.ref[, nc] == 1;
      members_num = sum(members_idx == T);
      
      members = cbind(obj_idx[members_idx], y[members_idx])
      rownames(members) = rownames(y)[members_idx]
      
      if (sum(!members_idx) > 0)
      {         
         nonmembers = cbind(obj_idx[!members_idx], y[!members_idx])
         rownames(nonmembers) = rownames(y)[!members_idx]
         data = list(members, nonmembers)
         c1 = mdaplot.getColors(n = 1)
         c2 = mdaplot.getColors(n = 2, colmap = 'gray')[1]
         mdaplotg(data, col = c(c1, c2), type = type, main = main, xlab = xlab, ylab = ylab, 
                  legend = legend, ylim = ylim, ...)
      }    
      else
      {   
         mdaplot(members, type = type, main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
      }   
   }
   else
   {
      data = cbind(1:length(y), y)
      rownames(data) = rownames(y)
      mdaplot(data, type = type, main = main, xlab = xlab, ylab = ylab, ...)      
   }   
}

plot.classres = function(obj, nc = 1, ...)
{
   # Plot method for "regres" objects
   #
   # Arguments:
   #   obj: object of "regres" class
   #   nc: number of response variable y to make the plot for   
   
   plotPredictions.classres(obj, nc = nc, ...)
}   

as.matrix.classres = function(obj, ncomp = NULL, nc = 1)
{
   # as.matrix method for "classres" objects
   #
   # Arguments:
   #   obj: object of "regres" class
   #   nc: number of response variable y to make the plot for   
   
   if (!is.null(obj$c.ref))
   {  
      if (is.null(ncomp))
         res = cbind(obj$tp[nc, ], obj$fp[nc, ], obj$fn[nc, ], obj$specificity[nc, ], obj$sensitivity[nc, ])
      else
         res = cbind(obj$tp[nc, ncomp], obj$fp[nc, ncomp], obj$fn[nc, ncomp], 
                     round(obj$specificity[nc, ncomp], 3), round(obj$sensitivity[nc, ncomp], 3))
      
      res[, 4:5] = round(res[, 4:5], 3)
      colnames(res) = c('TP', 'FP', 'FN', 'Spec', 'Sens')
   }
   else
   {
      res = NULL
   }   
   
   res
}

print.classres = function(obj, str = NULL, ...)
{
   # Print method for "classres" object
   #
   # Arguments:
   #   obj: object of "classres" class
   
   if (is.null(str))
      str = 'Classification results (class classres)\nMajor fields:'
   
   if (nchar(str) > 0)
      cat(sprintf('\n%s\n', str))   
   
   cat('$c.pred - predicted class values\n')
   if (!is.null(obj$c.ref))
   {   
      cat('$c.ref - reference (true) class values\n')
      cat('$tp - number of true positives\n')
      cat('$fp - number of false positives\n')
      cat('$fn - number of false negatives\n')
      cat('$specificity - specificity of predictions\n')
      cat('$sensitivity - sensitivity of predictions\n')
   }
   
}   

summary.classres = function(obj, ncomp = NULL, nc = NULL, ...)
{
   # Summary method for "calres" object
   #
   # Arguments:
   #   obj: object of "calres" class
   #   ncomp: which column of c.pred to use
   
   cat('\nClassiciation results (class classres) summary\n')
   if (!is.null(obj$c.ref))
   {         
      
      if (is.null(nc))
         nc = 1:ncol(obj$c.ref)

      if (!is.null(ncomp))
      {   
         cat(sprintf('\nNumber of selected components: %d', ncomp))
      }   
      else
      {   
         if (is.null(obj$ncomp.selected))
         {   
            ncomp = 1
         }   
         else
         {   
            ncomp = obj$ncomp.selected
            cat(sprintf('\nNumber of selected components: %d', ncomp))
         }   
      }
   
      cat(sprintf('\nNumber of classes: %d\n', ncol(obj$c.ref)))
      
      for (i in nc)
      {   
         cat(sprintf('\nClass "%s":\n', colnames(obj$c.ref)[i]))
         res = as.matrix.classres(obj, nc = i)    
         print(res)
      }      
   }      
   else
   {
      cat('No reference data provided to calculate prediction performance.')
   }   
}   

