## class and methods for classification models ##

checkReferenceValues.classmodel = function(model, c.ref, x) {
   if (is.null(c.ref))
      return()
   
   attrs = mda.getattr(c.ref)
   
   if (is.logical(c.ref))
      c.ref = ifelse(c.ref, model$classname, 'None')
   if (is.data.frame(c.ref))
      c.ref = as.matrix(c.ref)
   if (!is.matrix(c.ref))
      c.ref = matrix(c.ref, ncol = 1)
   if (ncol(c.ref) != 1)
      stop('Class reference values should be a vector or matrix/data frame with one column')
   if(nrow(c.ref) != nrow(x))
      stop('Matrix with predictors and classes should have the same number of rows!')
   if (!is.character(c.ref))
      stop('Matrix/vector with reference class values should be either logical or character!')

   c.ref = mda.setattr(c.ref, attrs)
   c.ref
}

#' Predictions plot for classification model
#' 
#' @description
#' Makes a plot with class predictions for a classification model.
#'
#' @param obj 
#' a classification model (object of class \code{simca}, \code{plsda}, etc.). if \code{NULL} value 
#' is specified, the result will be selected automatically by checking the nearest available from 
#' test, cv and calibration results.
#' @param res 
#' which result to make the plot for (\code{'calres'}, \code{'cvres'} or \code{'testres'}).
#' @param nc 
#' if there are several classes, which class to make the plot for (NULL - for all).
#' @param ncomp 
#' what number of components to make the plot for (NULL - for selected in the model).
#' @param main
#' main title for the plot
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#' 
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#'   
#' @export  
plotPredictions.classmodel = function(obj, res = NULL, nc = NULL, ncomp = NULL, main = NULL, ...) {   
   if (is.null(res))
   {
      if (!is.null(obj$testres))
         res = 'testres'
      else if (!is.null(obj$cvres))
         res = 'cvres'
      else
         res = 'calres'
   }
   resnames = c('calres', 'cvres', 'testres')
   resstrings = c('calbiration', 'cross-validation', 'test set')
   resobj = obj[[res]]

   if (!(res %in% resnames))
      stop('Wrong value for argument "res" (use "calres", "cvres" or "testres")!')
   else if (is.null(resobj))
      stop('The result does not exist!')
   else
   {
      if (is.null(main))
      {
         if (is.null(ncomp))
            main = sprintf('Predictions for %s results', resstrings[resnames == res])
         else
            main = sprintf('Predictions for %s results (ncomp = %d)', resstrings[resnames == res], ncomp)
      }
      plotPredictions.classres(resobj, nc = nc, ncomp = ncomp, main = main, ...)
   }
}

#' Specificity plot for classification model
#' 
#' @description
#' Makes a plot with specificity values vs. model complexity (e.g. number of components)
#' 
#' @param obj
#' classification model (object of class \code{plsda}, \code{simca}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#' 
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#' 
#' @export
plotSpecificity.classmodel = function(obj, nc = NULL, ...)
{   
   plotPerformance(obj, nc = nc, param = 'specificity', ...)
}

#' Sensitivity plot for classification model
#' 
#' @description
#' Makes a plot with sensitivity values vs. model complexity (e.g. number of components)
#' 
#' @param obj
#' classification model (object of class \code{plsda}, \code{simca}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#' 
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#' 
#' @export
plotSensitivity.classmodel = function(obj, nc = NULL, ...)
{
   plotPerformance(obj, nc = nc, param = 'sensitivity', ...)
}


#' Misclassified ratio plot for classification model
#' 
#' @description
#' Makes a plot with misclassified ratio values vs. model complexity (e.g. number of components)
#' 
#' @param obj
#' classification model (object of class \code{plsda}, \code{simca}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#' 
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#' 
#' @export
plotMisclassified.classmodel = function(obj, nc = NULL, ...)
{
   plotPerformance(obj, nc = nc, param = 'misclassified', ...)
}


#' Performance plot for classification model
#' 
#' @description
#' Makes a plot with sensitivity values vs. model complexity (e.g. number of components)
#' 
#' @param obj
#' classification model (object of class \code{plsda}, \code{simca}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param param
#' which parameter to make the plot for (\code{'specificity'}, \code{'sensitivity'}, 
#' or \code{'misclassified'})
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ylim
#' vector with two values - limits for y axis
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#' 
#' @export
plotPerformance.classmodel = function(obj, nc = NULL, param = 'specificity', type = 'h', 
                                 main = NULL, xlab = 'Components', ylab = '', 
                                 ylim = c(0, 1.15), ...)
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

      classname = sprintf('(%s)', obj$classnames[nc])
   }
   
   data = list(cal = obj$calres[[param]][nc, ])

   if (any(is.nan(data$cal))) {
      warning('Specificity values are not available!')
   } else {
      if (!is.null(obj$cvres))
         data$cv = obj$cvres[[param]][nc, ]   
      
      if (!is.null(obj$testres))
         data$test = obj$testres[[param]][nc, ]   
      
      if (is.null(main))
         main = sprintf('%s%s %s', toupper(substring(param, 1, 1)), substring(param, 2), 
                        toString(classname))
      mdaplotg(data, type = type, main = main, xticks = 1:obj$ncomp, xlab = xlab, ylab = ylab, ylim = ylim, ...)
   }
   
}

