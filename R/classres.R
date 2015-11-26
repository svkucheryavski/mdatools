#' Results of classification
#' @description 
#' \code{classres} is used to store results classification for one or multiple classes.
#'
#'      
#' @param c.pred
#' matrix with predicted values (+1 or -1) for each class.
#' @param c.ref
#' matrix with reference values for each class.
#' @param p.pred
#' matrix with probability values for each class.
#' @param ncomp.selected
#' vector with selected number of components for each class.
#'
#' @details 
#' There is no need to create a \code{classres} object manually, it is created automatically when 
#' build a classification model (e.g. using \code{\link{simca}} or \code{\link{plsda}}) or apply the 
#' model to new data. For any classification method from \code{mdatools}, a class using to represent 
#' results of classification (e.g. \code{\link{simcares}}) inherits fields and methods of 
#' \code{classres}.
#'
#' @return
#' \item{c.pred}{predicted class values (+1 or -1).}
#' \item{c.ref}{reference (true) class values if provided.}
#' 
#' The following fields are available only if reference values were provided.
#' \item{tp}{number of true positives.}
#' \item{fp}{nmber of false positives.}
#' \item{fn}{number of false negatives.}
#' \item{specificity}{specificity of predictions.}
#' \item{sensitivity}{sensitivity of predictions.}
#'
#' @seealso 
#' Methods \code{classres} class:
#' \tabular{ll}{
#'  \code{\link{showPredictions.classres}} \tab shows table with predicted values.\cr
#'  \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
#'  \code{\link{plotSensitivity.classres}} \tab makes plot with sensitivity vs. components values.\cr
#'  \code{\link{plotSpecificity.classres}} \tab makes plot with specificity vs. components values.\cr
#'  \code{\link{plotMisclassified.classres}} \tab makes plot with misclassified ratio values.\cr
#'  \code{\link{plotPerformance.classres}} \tab makes plot with misclassified ration, specificity and sensitivity values.\cr
#' }
classres = function(c.pred, c.ref = NULL, p.pred = NULL, ncomp.selected = NULL)
{
   if (!is.null(c.ref))
   {
      c.ref = as.matrix(c.ref)
      obj = getClassificationPerformance(c.ref, c.pred)
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
   obj$classnames = dimnames(c.pred)[[3]]

   obj$call = match.call()   
   class(obj) = "classres"

   obj
}

#' Calculation of  classification performance parameters
#'
#' @description
#' Calculates and returns performance parameters for classification result (e.g. number of false
#' negatives, false positives, sensitivity, specificity, etc.).
#' 
#' @param c.ref
#' reference class values for objects (vector with numeric or text values)
#' @param c.pred
#' predicted class values for objects (array nobj x ncomponents x nclasses) 
#'
#' @return
#' Returns a list with following fields:
#' \tabular{ll}{
#'    \code{$fn} \tab number of false negatives (nclasses x ncomponents) \cr
#'    \code{$fp} \tab number of false positives (nclasses x ncomponents) \cr
#'    \code{$tp} \tab number of true positives (nclasses x ncomponents) \cr
#'    \code{$sensitivity} \tab sensitivity values (nclasses x ncomponents) \cr
#'    \code{$specificity} \tab specificity values (nclasses x ncomponents) \cr
#'    \code{$sensitivity} \tab misclassified ratio values (nclasses x ncomponents) \cr
#' }
#'
#' @details
#' The function is called automatically when a classification result with reference values is 
#' created, for example when applying a \code{plsda} or \code{simca} models.
#' 
getClassificationPerformance = function(c.ref, c.pred)
{
   ncomp = dim(c.pred)[2]
   nobj = dim(c.pred)[1]
   nclasses = dim(c.pred)[3]

   tp = matrix(0, nrow = nclasses, ncol = ncomp)
   fp = matrix(0, nrow = nclasses, ncol = ncomp)
   fn = matrix(0, nrow = nclasses, ncol = ncomp)
   tn = matrix(0, nrow = nclasses, ncol = ncomp)
   
   specificity = matrix(0, nrow = nclasses + 1, ncol = ncomp)
   sensitivity = matrix(0, nrow = nclasses + 1, ncol = ncomp)
   misclassified = matrix(0, nrow = nclasses + 1, ncol = ncomp)

   classnames = dimnames(c.pred)[[3]]
   
   for (i in 1:nclasses)         
   {   
      fn[i, ] = colSums((c.ref[, 1] == classnames[i]) & (c.pred[, , i, drop = F] == -1))
      fp[i, ] = colSums((c.ref[, 1] != classnames[i]) & (c.pred[, , i, drop = F] == 1))
      tp[i, ] = colSums((c.ref[, 1] == classnames[i]) & (c.pred[, , i, drop = F] == 1))
      tn[i, ] = colSums((c.ref[, 1] != classnames[i]) & (c.pred[, , i, drop = F] == -1))
      
      sensitivity[i, ] = tp[i, ] / (tp[i, ] + fn[i, ])
      specificity[i, ] = tn[i, ] / (tn[i, ] + fp[i, ])
      misclassified[i, ] = (fp[i, ] + fn[i, ]) / nobj
   }

   sensitivity[nclasses + 1, ] = colSums(tp) / (colSums(tp) + colSums(fn))
   specificity[nclasses + 1, ] = colSums(tn) / (colSums(tn) + colSums(fp))
   misclassified[nclasses + 1, ] = (colSums(fp) + colSums(fn)) / (nclasses * nobj)

   rownames(fn) = rownames(fp) = rownames(tp) = dimnames(c.pred)[[3]]
   rownames(sensitivity) = rownames(specificity) = rownames(misclassified) = c(dimnames(c.pred)[[3]], 'Total')
   colnames(fn) = colnames(fp) = colnames(tp) =  colnames(sensitivity) = colnames(specificity) = dimnames(c.pred)[[2]]

   obj = list()
   obj$fn = fn
   obj$fp = fp
   obj$tp = tp
   obj$tn = tn
   obj$sensitivity = sensitivity
   obj$specificity = specificity
   obj$misclassified = misclassified

   obj
}   

#' Get selected components
#' 
#' @description
#' Returns number of components depending on user selection and object properites
#' 
#' @param obj
#' object with classification results (e.g. \code{plsdares} or \code{simcamres}).
#' @param ncomp
#' number of components specified by user.
#' 
#' @details
#' This is a technical function used for selection proper value for number of components in 
#' plotting functions.
#' 
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

#' Show predicted class values
#' 
#' @description
#' Shows a table with predicted class values for classification result.
#' 
#' @param obj
#' object with classification results (e.g. \code{plsdares} or \code{simcamres}).
#' @param ncomp
#' number of components to show the predictions for (NULL - use selected for a model).
#' @param ...
#' other parameters
#' 
#' @details
#' The function prints a matrix where every column is a class and every row is an data object.
#' The matrix has either -1 (does not belong to the class) or +1 (belongs to the class) values.
#' 
#' @export
showPredictions.classres = function(obj, ncomp = NULL, ...)
{
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


#' Sensitivity plot for classification results
#' 
#' @description
#' Makes a plot with sensitivity values vs. model complexity (e.g. number of components) for 
#' classification results.
#' 
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#' 
#' @export
plotSensitivity.classres = function(obj, nc = NULL, ...)
{
   plotPerformance(obj, nc = nc, param = 'sensitivity', ...)
}   

#' Specificity plot for classification results
#' 
#' @description
#' Makes a plot with specificity values vs. model complexity (e.g. number of components) for
#' classification results.
#' 
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#' 
#' @export
plotSpecificity.classres = function(obj, nc = NULL, ...)
{
   plotPerformance(obj, nc = nc, param = 'specificity', ...)
}   

#' Misclassified ratio plot for classification results
#' 
#' @description
#' Makes a plot with misclassified ratio values vs. model complexity (e.g. number of components) for
#' classification results.
#' 
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#' 
#' @export
plotMisclassified.classres = function(obj, nc = NULL, ...)
{
   plotPerformance(obj, nc = nc, param = 'misclassified', ...)
}   


#' Performance plot for classification results
#' 
#' @description
#' Makes a plot with classification performance parameters vs. model complexity (e.g. number of 
#' components) for classification results.
#' 
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param param
#' which performance parameter to make the plot for (\code{'sensitivity'}, \code{'specificity'}, 
#' \code{'misclassified'}, \code{'all'}).
#' @param type
#' type of the plot
#' @param legend
#' vector with legend items
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ylim
#' vector with two values - limits for y axis
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#' 
#' @export
plotPerformance.classres = function(obj, nc = NULL, param = 'all', type = 'h', legend = NULL, 
                                    main = NULL, xlab = 'Components', 
                                    ylab = '', ylim = c(0, 1.1), ...)
{
   if (is.null(nc))
   {
      nc =  obj$nclasses + 1
      
      if (obj$nclasses > 1)
         classname = '(all classes)'
      else
         classname = ''
   }
   else
   {
      if (nc > obj$nclasses || nc < 1)
         stop('Wrong value for argument "nc"!')

      classname = sprintf('(%s)', dimnames(obj$c.pred)[[3]][nc])
   }

   if (param == 'all')
   {   
      if (is.null(main))
      {  
         main = sprintf('Prediction performance %s', classname);
      }

      if (is.null(legend))
         legend = c('sensitivity', 'specificity', 'misclassified')

      data = cbind(1:length(obj$sensitivity[nc, ]), obj$sensitivity[nc, ], obj$specificity[nc, ],
                   obj$misclassified[nc, ])
      labels = round(cbind(obj$sensitivity[nc, ], obj$specificity[nc, ], obj$misclassified[nc, ]), 2)
      
      mdaplotg(data, type = type, legend = legend, main = main, xlab = xlab, ylim = ylim,
               ylab = ylab, labels = labels, ...)
   }
   else
   {
      if (is.null(main))
         main = sprintf('%s%s %s', toupper(substring(param, 1, 1)), substring(param, 2, ), classname)
      
      data = cbind(1:length(obj[[param]][nc, ]), obj[[param]][nc, ])
      labels = round(obj[[param]][nc, ], 2)
      mdaplot(data, type = type, main = main, xlab = xlab, ylab = ylab, ylim = ylim, labels = labels, ...)
   }
}

#' Prediction plot for classification results
#' 
#' @description
#' Makes a plot with predicted class values for classification results.
#' 
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ncomp
#' which number of components to make the plot for (can be one value for all classes or vector with
#' separate values for each, if NULL - model selected number will be used).
#' @param type
#' type of the plot
#' @param legend
#' vector with legend items
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ylim
#' vector with two values - limits for y axis
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} or \code{\link{mdaplot}} function 
#' can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#' 
#' @export
plotPredictions.classres = function(obj, nc = NULL, ncomp = NULL, type = 'p', legend = NULL, 
                                    main = NULL, xlab = 'Objects', ylab = NULL, ylim = c(-1.2, 1.2),
                                    ...)
{
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

         if (length(objn) > 0)
         {   
            vals = matrix(nc[i], nrow = length(objn), ncol = 1)

            data = cbind(objn, vals)
            rownames(data) = rownames(obj$c.pred)[objn]

            refdata = c(refdata, obj$c.ref[objn])
            cdata = rbind(cdata, data)
         }   
      }

      if (length(which(nullidx == T)) > 0)
      {  
         data = cbind(which(nullidx == T), 0)
         rownames(data) = rownames(obj$c.pred)[nullidx]
         cdata = rbind(cdata, data)
         refdata = c(refdata, obj$c.ref[nullidx == T])
      }

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

#' Plot function for classification results
#' 
#' @description
#' Generic plot function for classification results. Shows predicted class values.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' other arguments
#' 
#' @export
plot.classres = function(x, nc = NULL, ...)
{
   plotPredictions.classres(x, nc = nc, ...)
}   

#' as.matrix method for classification results
#' 
#' @description
#' Generic \code{as.matrix} function for classification results. Returns matrix with performance 
#' values for specific class.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param ncomp
#' model complexity (number of components) to show the parameters for.
#' @param nc
#' if there are several classes, which class to show the parameters for.
#' @param ...
#' other arguments
#' 
#' @export
as.matrix.classres = function(x, ncomp = NULL, nc = 1, ...)
{
   obj = x
   
   if (!is.null(obj$c.ref))
   {  
      if (is.null(ncomp))
         res = cbind(obj$tp[nc, ], obj$fp[nc, ], obj$tn[nc, ], obj$fn[nc, ], 
                     obj$specificity[nc, ], obj$sensitivity[nc, ])
      else
         res = cbind(obj$tp[nc, ncomp], obj$fp[nc, ncomp], obj$tn[nc, ncomp], obj$fn[nc, ncomp], 
                     round(obj$specificity[nc, ncomp], 3), round(obj$sensitivity[nc, ncomp], 3))

   res[, 4:5] = round(res[, 4:5], 3)
   colnames(res) = c('TP', 'FP', 'TN', 'FN', 'Spec', 'Sens')
   }
   else
   {
      res = NULL
   }   

   res
}

#' Print information about classification result object
#' 
#' @description
#' Generic \code{print} function for classification results. Prints information about major fields
#' of the object.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param str
#' User specified text (e.g. to be used for particular method, like PLS-DA, etc).
#' @param ...
#' other arguments
#' 
#' @export
print.classres = function(x, str = NULL, ...)
{
   if (is.null(str))
      str = 'Classification results (class classres)\nMajor fields:'

   if (nchar(str) > 0)
      cat(sprintf('\n%s\n', str))   

   cat('$c.pred - predicted class values\n')
   if (!is.null(x$c.ref))
   {   
      cat('$c.ref - reference (true) class values\n')
      cat('$tp - number of true positives\n')
      cat('$fp - number of false positives\n')
      cat('$fn - number of false negatives\n')
      cat('$specificity - specificity of predictions\n')
      cat('$sensitivity - sensitivity of predictions\n')
      cat('$misclassified - misclassification ratio for predictions\n')
   }
}   

#' Summary statistics about classification result object
#' 
#' @description
#' Generic \code{summary} function for classification results. Prints performance values for the 
#' results.
#'
#' @param object
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param ncomp
#' which number of components to make the plot for (can be one value for all classes or vector with
#' separate values for each, if NULL - model selected number will be used).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' other arguments
#' 
#' @export
summary.classres = function(object, ncomp = NULL, nc = NULL, ...)
{
   cat('\nClassiciation results (class classres) summary\n')
   if (!is.null(object$c.ref))
   {         

      if (is.null(nc))
         nc = 1:dim(object$c.pred)[3]
      show(nc)
      if (!is.null(ncomp))
      {   
         cat(sprintf('\nNumber of selected components: %d', ncomp))
      }   
      else
      {   
         if (is.null(object$ncomp.selected))
         {   
            ncomp = 1
         }   
         else
         {   
            ncomp = object$ncomp.selected
            cat(sprintf('\nNumber of selected components: %d', ncomp))
         }   
      }

      cat(sprintf('\nNumber of classes: %d\n', ncol(object$c.ref)))

      for (i in nc)
      {   
         cat(sprintf('\nClass "%s":\n', colnames(object$c.ref)[i]))
         res = as.matrix.classres(object, nc = i)    
         print(res)
      }      
   }      
   else
   {
      cat('No reference data provided to calculate prediction performance.')
   }   
}   