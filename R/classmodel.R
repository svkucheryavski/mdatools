## class and methods for classification models ##

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
plotPredictions.classmodel = function(obj, res = NULL, nc = NULL, ncomp = NULL, main = NULL, ...)
{   
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
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#' 
plotPerformance.classmodel = function(obj, nc = NULL, param = 'specificity', type = 'h', legend = NULL, 
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
   
   data = cbind(1:obj$ncomp, obj$calres[[param]][nc, ])
   labels = matrix(mdaplot.formatValues(obj$calres[[param]][nc, ]), ncol = 1)
   legend_str = 'cal'
   
   if (!is.null(obj$cvres))
   {
      data = cbind(data, obj$cvres[[param]][nc, ])   
      labels = cbind(labels, mdaplot.formatValues(obj$cvres[[param]][nc, ]))
      legend_str = c(legend_str, 'cv')
   }   
   
   if (!is.null(obj$testres))
   {
      data = cbind(data, obj$testres[[param]][nc, ])   
      labels = cbind(labels, mdaplot.formatValues(obj$testres[[param]][nc, ]))
      legend_str = c(legend_str, 'test')
   }
  
   if (is.null(main))
      main = sprintf('%s%s %s', toupper(substring(param, 1, 1)), substring(param, 2, ), toString(classname))

   if (!is.null(legend))
      legend_str = legend
   
   mdaplotg(data, type = type, main = main, xlab = xlab, ylab = ylab, legend = legend_str,
            ylim = ylim, labels = labels, ...)
}

