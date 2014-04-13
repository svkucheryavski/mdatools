plsdares = function(plsres, cres)
{
   # Creates an object of plsda class. 
   #
   # Arguments:
   #  plsres: an object of plsres class (results for PLS regression)
   #  cres: an object of classres class (results for classification)
   
   obj = c(plsres, cres)
   class(obj) = c('plsdares', 'classres', 'plsres')   
   
   obj$call = match.call()   
   
   obj
}   

#' Overview plot for PLS-DA results
#' 
#' @description
#' Shows a set of plots (x residuals, y variance, classification performance and predictions) 
#' for PLS-DA results.
#' 
#' @param x
#' PLS-DA results (object of class \code{plsdares})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param nc
#' which class to show the summary for (if NULL, will be shown for all)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plot.plsdares = function(x, nc = NULL, ncomp = NULL, show.labels = F, ...)
{
   obj = x
   
   par(mfrow = c(2, 2))
   plotXResiduals.plsres(obj, ncomp = ncomp, show.labels = show.labels)
   plotYVariance.plsres(obj, ncomp = ncomp, show.labels = show.labels)
   plotPerformance(obj, nc = nc, ncomp = ncomp, show.labels = show.labels)
   plotPredictions(obj, nc = nc, ncomp = ncomp, show.labels = show.labels)
   par(mfrow = c(1, 1))
}

#' as.matrix method for PLS-DA results
#' 
#' @method as.matrix plsdares
#' @S3method as.matrix plsdares
#'
#' @description
#' Returns a matrix with model performance statistics for PLS-DA results
#' 
#' @param x
#' PLS-DA results (object of class \code{plsdares})
#' @param ncomp
#' number of components to calculate the statistics for
#' @param nc
#' for which class to calculate the statistics for
#' @param ...
#' other arguments
#' 
as.matrix.plsdares = function(x, ncomp = NULL, nc = NULL, ...)
{
   obj = x
   
   plsmat = as.matrix.plsres(obj, ncomp = ncomp, ny = nc)
   classmat = as.matrix.classres(obj, ncomp = ncomp, nc = nc)
   mat = cbind(plsmat[, 1:4, drop = F], classmat)
   
   mat
}

#' Summary method for PLS-DA results object
#' 
#' @method summary plsdares
#' @S3method summary plsdares
#' 
#' @description
#' Shows performance statistics for the results.
#' 
#' @param object
#' PLS-DA results (object of class \code{plsdares})
#' @param nc
#' which class to show the summary for (if NULL, will be shown for all)
#' @param ...
#' other arguments
#' 
summary.plsdares = function(object, nc = NULL, ...)
{
   obj = object
   
   if (is.null(nc))
      nc = 1:obj$nclasses
   cat('\nPLS-DA results (class plsdares) summary:\n');
   cat(sprintf('Number of selected components: %d\n', obj$ncomp.selected))
   for (n in nc)
   {
      cat(sprintf('\nClass #%d (%s)\n', n, obj$classnames[n]))
      
      mat = as.matrix(obj, nc = n)
      mat[, 1:4] = round(mat[, 1:4], 2)
      print(mat)
      cat('\n')
   }
}

#' Print method for PLS-DA results object
#' 
#' @method print plsdares
#' @S3method print plsdares
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' PLS-DA results (object of class \code{plsdares})
#' @param ...
#' other arguments
#' 
print.plsdares = function(x, ...)
{
   obj = x
   
   cat('\nPLS-DA results (class plsdares)\n')
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$ncomp.selected - number of selected components\n')
   cat('$c.pred - array with predicted class values\n')
   if (!is.null(obj$c.ref))
   {   
      cat('$c.ref - vector with reference class values\n')
      cat('$tp - number of true positives\n')
      cat('$fp - number of false positives\n')
      cat('$fn - number of false negatives\n')
      cat('$specificity - specificity of predictions\n')
      cat('$sensitivity - sensitivity of predictions\n')
      cat('$misclassified - misclassification ratio for predictions\n')
      cat('$ydecomp - decomposition of y values (ldecomp object)\n')
   }
   cat('$xdecomp - decomposition of x values (ldecomp object)\n')
}
