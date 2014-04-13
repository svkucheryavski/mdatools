## class and methods for PCA results ##

pcares = function(scores, loadings, residuals, totvar, tnorm = NULL, ncomp.selected = NULL, ...)
{
   # Creates an object of pcares class. In fact the class is a wrapper for ldecomp and
   # uses its methods and attributes.
   
   res = ldecomp(scores, loadings, residuals, totvar, tnorm = tnorm, ncomp.selected = ncomp.selected, ...)
   class(res) = c('pcares', 'ldecomp')   
   
   res
}   


#' Pplot method for PCA results object
#' 
#' @description
#' Show several plots to give an overview about the PCA results
#' 
#' @param x
#' PCA results (object of class \code{pcares})
#' @param comp
#' which components to show the scores plot for (can be one value or vector with two values).
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#' 
plot.pcares = function(x, comp = c(1, 2), show.labels = T, ...)
{   
   par(mfrow = c(2, 2))
   plotScores(x, comp = comp, show.labels = show.labels, ...)
   plotResiduals(x, show.labels = show.labels, ...)
   plotVariance(x, show.labels = show.labels, ...)
   plotCumVariance(x, show.labels = show.labels, ...)
   par(mfrow = c(1, 1))
}


#' Summary method for PCA results object
#' 
#' @method summary pcares
#' @S3method summary pcares
#'
#' @description
#' Shows some statistics (explained variance, eigenvalues) about the results.
#' 
#' @param object
#' PCA results (object of class \code{pcares})
#' @param ...
#' other arguments
#' 
summary.pcares = function(object, ...)
{
   summary.ldecomp(object, 'Summary for PCA results', ...)
}


#' Print method for PCA results object
#' 
#' @method print pcares
#' @S3method print pcares
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' PCA results (object of class \code{pcares})
#' @param ...
#' other arguments
#' 
print.pcares = function(x, ...)
{   
   print.ldecomp(x, 'Results for PCA decomposition (class pcares)', ...)
   cat('\n')
}