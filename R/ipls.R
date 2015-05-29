## class and methods for iPLS ##

#' Interval iPLS variable selection
#' 
#' @description 
#' Applies iPLS alrogithm to find variable intervals most important for
#' prediction
#' 
#'  @param x
#'  a matrix with predictor values
#'  @param y
#'  a vector with response values
#'  @param ncomp
#'  maximum number of components to use in PLS model
#'  @param nint
#'  number of intervals
#'  @param center
#'  logical, center or not the data values
#'  @param scale
#'  logical, standardize or not the data values
#'  @param cv
#'  which cross-validation to use
#'  @param method
#'  iPLS method (\code{'forward'} or \code{'backward'})
#'  
#'  @return 
#'  
#'  @details 
#'  
ipls = function(x, y, ncomp = 10, nint = NULL, center = T, scale = F, cv = 10, method = 'forward')
{

   model$call = match.call()   
   class(model) = "ipls"
   
   model
}

#' Runs the iPLS algorithm
#' 
#' @param x
#'  
ipls.fit = function(x, y, ncomp, nint, center, scale, method = 'forward', cv = 10)
{
   
}
      
            
#' Overview plot for iPLS model
#' 
#' @description
#' Shows a set of plots () for iPLS model.
#' 
#' @param x
#' a  (object of class \code{pca})
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{ipls}} function.
#' 
plot.ipls = function(x, ...)
{   
}

#' Print method for iPLS
#' 
#' @method print ipls
#' @S3method print ipls
#'
#' @description
#' Prints information about the iPLS object structure
#' 
#' @param x
#' a iPLS (object of class \code{ipls})
#' @param ...
#' other arguments
#' 
print.ipls = function(x, ...)
{
}

#' Summary method for iPLS model object
#' 
#' @method summary ipls
#' @S3method summary ipls
#'
#' @description
#' Shows some statistics () for the .
#' 
#' @param object
#' a iPLS (object of class \code{ipls})
#' @param ...
#' other arguments
#' 
summary.ipls = function(object, ...)
{
}
