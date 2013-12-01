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
   #
   # In fact ncomp.selected is a vector with a value for each class
   #

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
   obj$ncomp.selected = ncomp.selected
   obj$nclasses = dim(c.pred)[3]
   obj$ncomp = dim(c.pred)[2]

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

   ncomp = dim(c.pred)[2]
   nclasses = dim(c.pred)[3]

   tp = matrix(0, nrow = nclasses, ncol = ncomp)
   fp = matrix(0, nrow = nclasses, ncol = ncomp)
   fn = matrix(0, nrow = nclasses, ncol = ncomp)
   specificity = matrix(0, nrow = nclasses, ncol = ncomp)
   sensitivity = matrix(0, nrow = nclasses, ncol = ncomp)
   f1 = matrix(0, nrow = nclasses, ncol = ncomp)

   for (i in 1:nclasses)         
   {   
      if (is.numeric(c.ref))
      {   
         fn[i, ] = colSums((c.ref[, 1] == i) & (c.pred[, , i, drop = F] == -1))
         fp[i, ] = colSums((c.ref[, 1] != i) & (c.pred[, , i, drop = F] == 1))
         tp[i, ] = colSums((c.ref[, 1] == i) & (c.pred[, , i, drop = F] == 1))
      }
      else
      {
         cname = dimnames(c.pred)[[3]][i]
         fn[i, ] = colSums((c.ref[, 1] == cname) & (c.pred[, , i, drop = F] == -1))
         fp[i, ] = colSums((c.ref[, 1] != cname) & (c.pred[, , i, drop = F] == 1))
         tp[i, ] = colSums((c.ref[, 1] == cname) & (c.pred[, , i, drop = F] == 1))         
      }   

      sensitivity[i, ] = tp[i, ] / (tp[i, ] + fn[i, ])
      specificity[i, ] = tp[i, ] / (tp[i, ] + fp[i, ])
      f1[i, ] = 2 * sensitivity[i, ] * specificity[i, ] / (sensitivity[i, ] + specificity[i, ])
   }

   rownames(fn) = rownames(fp) = rownames(tp) = rownames(f1) = 
   rownames(sensitivity) = rownames(specificity) = dimnames(c.pred)[[3]]
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

getSelectedComponents.classres = function(obj, ncomp = NULL)
{
   if (is.null(ncomp))
   {   
      if (is.null(obj$ncomp.selected) || dim(obj$c.pred)[2] == 1)
         ncomp = 1
      else
         ncomp = obj$ncomp.selected
   }   

   if (length(ncomp) == 1)
      ncomp = matrix(ncomp, nrow = 1, ncol = obj$nclasses)
   
   ncomp
}

showPredictions.classres = function(obj, ncomp = NULL)
{
   #  Shows table with predictions for selected or user specified number of components
   #
   # Arguments:
   #   obj: object of "classres" class
   #   ncomp: for which components to show the results for

   ncomp = getSelectedComponents.classres(obj, ncomp)

   if (obj$nclasses == 1)
      pred = obj$c.pred[ , ncomp[1], , drop = F]
   else
   {  
      if (length(ncomp) == 1)
         ncomp = matrix(ncomp, nrow = 1, ncol = obj$nclasses)

   pred = NULL
   for (i in 1:obj$nclasses)
      pred = cbind(pred, obj$c.pred[, ncomp[i], i, drop = F])
   }

   dim(pred) = c(nrow(pred), obj$nclasses)
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

plotPredictions.classres = function(obj, nc = NULL, ncomp = NULL, type = 'p', main = NULL, 
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

   if (is.null(xlab))
      xlab = 'Objects'

   if (is.null(nc))
      nc = 1:obj$nclasses

   if (is.null(main))
   {   
      main = 'Predictions'

      if (!is.null(ncomp) && length(ncomp) == 1)
         main = sprintf('%s (ncomp = %d)', main, ncomp)
   }

   ncomp = getSelectedComponents.classres(obj, ncomp)
   
   if (max(ncomp) > dim(obj$c.pred)[2])
      stop('Wrong value for ncomp parameter!')

   if (is.null(ylab))
      ylab = 'Classes'

   if (!is.null(obj$c.ref))
   {   
      cdata = NULL
      refdata = NULL
      nullidx = matrix(TRUE, nrow = nrow(obj$c.pred), ncol = 1)

      for (i in 1:length(nc))
      {
         nullidx = obj$c.pred[, ncomp[i], nc[i]] == -1 & nullidx

         objn = which(obj$c.pred[, ncomp[i], nc[i]] == 1)
         vals = matrix(nc[i], nrow = length(objn), ncol = 1)

         data = cbind(objn, vals)
         rownames(data) = rownames(obj$c.pred)[objn]

         refdata = c(refdata, obj$c.ref[objn])
         cdata = rbind(cdata, data)
      }

      data = cbind(which(nullidx == T), 0)
      rownames(data) = rownames(obj$c.pred)[nullidx]
      cdata = rbind(cdata, data)
      refdata = c(refdata, obj$c.ref[nullidx == T])

      pdata = list()
      refc = unique(obj$c.ref)
      legend_str = NULL;

      for (i in 1:length(refc))
      {
         idx = refdata == refc[i]         
         data = cdata[idx, 1:2]

         pdata[[i]] = data
         legend_str[i] = toString(refc[i])
      }  

      if (is.null(legend) && length(legend_str) > 1)
         legend = legend_str

      mdaplotg(pdata, type = 'p', main = main, xlab = xlab, ylab = ylab, legend = legend,
               yticks = c(0, nc), yticklabels = c('None', dimnames(obj$c.pred)[[3]][nc]), ...)
   }
   else
   {
      cdata = NULL
      nullidx = matrix(TRUE, nrow = nrow(obj$c.pred), ncol = 1)

      for (i in 1:length(nc))
      {
         nullidx = obj$c.pred[, ncomp[i], nc[i]] == -1 & nullidx

         objn = which(obj$c.pred[, ncomp[i], nc[i]] == 1)
         vals = matrix(nc[i], nrow = length(objn), ncol = 1)

         data = cbind(objn, vals)            
         rownames(data) = rownames(obj$c.pred)[objn]

         cdata = rbind(cdata, data)
      }

      data = cbind(which(nullidx == T), 0)
      rownames(data) = rownames(obj$c.pred)[nullidx]
      cdata = rbind(cdata, data)

      mdaplot(cdata, type = 'p', main = main, xlab = xlab, ylab = ylab,
              yticks = c(0, nc), yticklabels = c('None', dimnames(obj$c.pred)[[3]][nc]), ...)         
   }   
}

plot.classres = function(obj, nc = NULL, ...)
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
         nc = 1:dim(obj$c.pred)[3]
      show(nc)
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

