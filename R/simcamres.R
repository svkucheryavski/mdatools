## class and methods for SIMCA multi class classification results ##

simcamres = function(cres, T2, Q2, T2lim, Q2lim)
{
   # Creates an object of simcamres class. 
   #
   # Arguments:
   #  cres: an object of classres class (results for classification)
   #  T2: T2 values for the objects and classes (selected component only)
   #  Q2: Q2 values for the objects and classes (selected component only)
   #  T2lim: T2 limits for the classes (selected component only)
   #  Q2lim: Q2 limits for the classes (selected component only)

   res = cres
   res$T2 = T2
   res$Q2 = Q2
   res$T2lim = T2lim
   res$Q2lim = Q2lim
   res$classnames = dimnames(cres$c.pred)[[3]]
   class(res) = c('simcamres', 'classres')   
   
   res
}   

#' Residuals plot for SIMCAM results
#' 
#' @description
#' Shows a plot with Q2 vs. T2 residuals for SIMCAM results
#' 
#' @param obj
#' SIMCAM results (object of class \code{simcamres})
#' @param nc
#' which class (SIMCA model) to show the plot for
#' @param show.limits
#' logical, show or not lines with statistical limits for the residuals
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param legend
#' vector with legend items
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcamres}} function.
#' 
plotResiduals.simcamres = function(obj, nc = 1, show.limits = T, type = 'p', main = NULL, 
                                  xlab = 'T2', ylab = 'Q2', legend = NULL, ...)
{
   # set main title
   if (is.null(main))
      main = sprintf('Residuals (%s)', obj$classnames[nc])

   if (!is.null(obj$c.ref))
   {   
      classes = unique(obj$c.ref)
      data = list()
      for (i in 1:length(classes))
         data[[i]] = cbind(obj$T2[obj$c.ref == classes[i], nc], obj$Q2[obj$c.ref == classes[i], nc])
      
      if (is.null(legend))
         legend = classes
      
      if (show.limits == T)
         show.lines = c(obj$T2lim[nc], obj$Q2lim[nc])
      else
         show.lines = F
      
      mdaplotg(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines,
               legend = legend, ...)
   }
   else
   {

      data = cbind(obj$T2[, nc], obj$Q2[, nc])
      if (show.limits == T)
         show.lines = c(obj$T2lim[nc], obj$Q2lim[nc])
      else
         show.lines = F
      
      mdaplot(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines, ...)
   }
}

#' Cooman's plot for SIMCAM results
#' 
#' @description
#' Shows a Cooman's plot for a pair of SIMCA models
#' 
#' @param obj
#' SIMCAM results (object of class \code{simcamres})
#' @param nc
#' vector with two values - classes (SIMCA models) to show the plot for
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.limits
#' logical, show or not lines with statistical limits for the residuals
#' @param legend
#' vector with legend items 
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcamres}} function.
#' 
plotCooman.simcamres = function(obj, nc = c(1, 2), type = 'p', main = "Cooman's plot", xlab = NULL, 
                                ylab = NULL, show.limits = T, legend = NULL, ...)
{
   # set labels for axes
   if (is.null(xlab))
      xlab = sprintf('Distance to class %s', obj$classnames[nc[1]])

   if (is.null(ylab))
      ylab = sprintf('Distance to class %s', obj$classnames[nc[2]])
   
   if (!is.null(obj$c.ref))
   {   
      classes = unique(obj$c.ref)
      data = list()
      for (i in 1:length(classes))
         data[[i]] = cbind(sqrt(obj$Q2[obj$c.ref == classes[i], nc[1]]), 
                           sqrt(obj$Q2[obj$c.ref == classes[i], nc[2]])
                           )
      
      if (is.null(legend))
         legend = classes
      
      if (show.limits == T)
         show.lines = c(obj$Q2lim[nc[1]], obj$Q2lim[nc[2]])
      else
         show.lines = F
      
      mdaplotg(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines,
               legend = legend, ...)
   }
   else
   {
      data = cbind(sqrt(obj$Q2[, nc[1]]), sqrt(obj$Q2[, nc[2]]))
      
      if (show.limits == T)
         show.lines = c(obj$Q2lim[nc[1]], obj$Q2lim[nc[2]])
      else
         show.lines = F
      
      mdaplot(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines, ...)
   }
}

#' Model overview plot for SIMCAM results
#' 
#' @description
#' Just shows a prediction plot for SIMCAM results.
#' 
#' @param x
#' SIMCAM results (object of class \code{simcamres})
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{simcamres}} function.
#' 
plot.simcamres = function(x, ...)
{
   plotPredictions(x)
}

#' Summary method for SIMCAM results object
#' 
#' @method summary simcamres
#' @S3method summary simcamres
#' 
#' @description
#' Shows performance statistics for the results.
#' 
#' @param object
#' SIMCAM results (object of class \code{simcamres})
#' @param ...
#' other arguments
summary.simcamres = function(object, ...)
{
   obj = object
   
   cat('\nSummary for SIMCA multiple classes classification result\n')
   if (!is.null(obj$c.ref))
   {   
      classres = NULL
      for (i in 1:obj$nclasses)
         classres = rbind(classres, as.matrix.classres(obj, nc = i))
      rownames(classres) = paste(obj$classnames, ' (', obj$ncomp.selected, ' comp)', sep = '')   
      print(classres)
   }
   else
   {
      cat('\nReference values are not provided.\n')
   }
}  

#' Print method for SIMCAM results object
#' 
#' @method print simcamres
#' @S3method print simcamres
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' SIMCAM results (object of class \code{simcamres})
#' @param ...
#' other arguments
#' 
print.simcamres = function(x, ...)
{   
   obj = x
   
   cat('Result for SIMCA multiple classes classification (class simcamres)\n')
   print.classres(obj, '')
   cat('\n')
}
