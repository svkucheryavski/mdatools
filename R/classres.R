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
classres = function(c.pred, c.ref = NULL, p.pred = NULL, ncomp.selected = NULL) {
   if (!is.null(c.ref)) {
      attrs = mda.getattr(c.ref)
      c.ref = as.matrix(c.ref)
      c.ref = mda.setattr(c.ref, attrs)
      obj = getClassificationPerformance(c.ref, c.pred)
      obj$c.ref = c.ref
   } else {
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
   # remove excluded rows for correct calculation of performance
   attrs = mda.getattr(c.pred)
   if (length(attrs$exclrows) > 0) {
      c.pred = c.pred[-attrs$exclrows, , , drop = F]
      if (nrow(c.ref) > nrow(c.pred))
         c.ref = c.ref[-attrs$exclrows, , drop = F]
   }
 
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
getSelectedComponents.classres = function(obj, ncomp = NULL) {
   if (is.null(ncomp)) {   
      if (is.null(obj$ncomp.selected) || dim(obj$c.pred)[2] == 1)
         ncomp = 1
      else
         ncomp = obj$ncomp.selected
   }   
   
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
   if (is.null(nc)) {
      nc =  obj$nclasses + 1
      
      if (obj$nclasses > 1)
         classname = '(all classes)'
      else
         classname = ''
   } else {
      if (nc > obj$nclasses || nc < 1)
         stop('Wrong value for argument "nc"!')

      classname = sprintf('(%s)', dimnames(obj$c.pred)[[3]][nc])
   }

   if (param == 'all') {   
      if (is.null(main))
         main = sprintf('Prediction performance %s', classname);

      data = list()
      if (!any(is.na(obj$sensitivity[nc, ])))
          data$sensitivity = obj$sensitivity[nc, ]
      if (!any(is.na(obj$specificity[nc, ])))
         data$specificity = obj$specificity[nc, ]
      if (!any(is.na(obj$misclassified[nc, ])))
         data$misclassified = obj$misclassified[nc, ]

      mdaplotg(data, type = type, legend = legend, main = main, xticks = 1:obj$ncomp, 
               xlab = xlab, ylim = ylim, ylab = ylab, ...)
   } else {
      if (is.null(main))
         main = sprintf('%s%s %s', toupper(substring(param, 1, 1)), substring(param, 2), classname)
      
      data = obj[[param]][nc, , drop = F]
      if (any(is.na(data)))
         stop('This performance parameter has NA values!')
      
      mdaplot(data, type = type, main = main, xticks = 1:obj$ncomp, xlab = xlab, ylab = ylab, ylim = ylim, ...)
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
#' which number of components to make the plot for (one value, if NULL - model selected number will be used).
#' This parameter shal not be used for multiclass models or results as predictions in this case are made only
#' for optimal number of components
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param ylab
#' label for y axis
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} or \code{\link{mdaplot}} function 
#' can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#' 
#' @export
plotPredictions.classres = function(obj, nc = NULL, ncomp = NULL, type = 'p', main = NULL, ylab = '', ...) {
  
   if (is.null(obj))
      stop('No classification results were provided!')
   
   # get classnames
   classnames = dimnames(obj$c.pred)[[3]]

   # if vector with class number was not provided show plot for all of them
   if (is.null(nc))
      nc = 1:obj$nclasses

   # if vector with names were provided instead of numbers - convert to numbers
   if (is.character(nc))
      nc = which(classnames %in% nc)
   
   # take unique values and check correctness
   nc = unique(nc)
   if (length(nc) == 0 || min(nc)< 1 || max(nc) > dim(obj$c.pred)[3])
      stop('Incorrect class numbers or names!')
   nclasses = length(nc)   
   
   # set main title
   if (is.null(main)) {   
      main = 'Predictions'
      if (!is.null(ncomp) && length(ncomp) == 1)
         main = sprintf('%s (ncomp = %d)', main, ncomp)
   }

   ncomp = getSelectedComponents.classres(obj, ncomp)
  
   if (max(ncomp) > dim(obj$c.pred)[2])
      stop('Wrong value for ncomp parameter!')

   # extract predicted values for particular component
   attrs = mda.getattr(obj$c.pred)
   c.pred = as.matrix(obj$c.pred[, ncomp, nc])
   
   if (nclasses > 1) {
      # prepare matrix for the results
      pdata = matrix(0, nrow = nrow(c.pred), ncol = nclasses + 1)
      # find objects that were not classified for any of the class and set rows = 1 in first column
      pdata[, 1] = apply(c.pred, 1, function(x)(all(x == -1)))
      # set values for the other 
      for (i in 1:nclasses)
         pdata[, i + 1] = (c.pred[, i] == 1) * (i + 1)
      # unfold the matrix
      dim(pdata) = NULL
      # add evector with indices of objects
      pdata = cbind(rep(1:nrow(c.pred), nclasses + 1), pdata)
      # remove all rows with zeros in the second column
      pdata = pdata[pdata[, 2] != 0, ]
      # assign proper rownames
      rownames(pdata) = rownames(c.pred)[pdata[, 1]]
      # assign other attributes
      pdata = mda.setattr(pdata, attrs)
      # exclude rows
      attr(pdata, 'exclrows') = NULL
      pdata = mda.exclrows(pdata, pdata[, 1] %in% attrs$exclrows)
      row.ind = pdata[, 1]
   } else {
      pdata = cbind((1:nrow(c.pred)), (c.pred == 1) + 1)
      pdata = mda.setattr(pdata, attrs)
   }
   
   colnames(pdata) = c(ifelse(is.null(attrs$yaxis.name), 'Objects', attrs$yaxis.name), 'Classes')
   
   if (is.null(obj$c.ref)) {   
      mdaplot(pdata, type = 'p', main = main, ylab = ylab, yticks = c(1:(nclasses + 1)), 
              yticklabels = c('None', classnames[nc]), ylim = c(0.8, nclasses + 1.2), ...)
   } else {
      pdata.g = as.factor(obj$c.ref[pdata[, 1]])
      mdaplotg(pdata, type = 'p', main = main, ylab = ylab, yticks = c(1:(nclasses + 1)), groupby = pdata.g,
              yticklabels = c('None', classnames[nc]), ylim = c(0.8, nclasses + 1.2), ...)
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
plot.classres = function(x, nc = NULL, ...){
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

   if (is.null(obj$c.ref))
      return()
   
   if (is.null(ncomp))
      res = cbind(obj$tp[nc, ], obj$fp[nc, ], obj$tn[nc, ], obj$fn[nc, ], 
               round(obj$specificity[nc, ], 3), round(obj$sensitivity[nc, ], 3))
   else
      res = cbind(obj$tp[nc, ncomp], obj$fp[nc, ncomp], obj$tn[nc, ncomp], obj$fn[nc, ncomp], 
               round(obj$specificity[nc, ncomp], 3), round(obj$sensitivity[nc, ncomp], 3))
   
   colnames(res) = c('TP', 'FP', 'TN', 'FN', 'Spec', 'Sens')
   
   if (any(is.na(obj$specificity)))
      res = res[, -5]

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