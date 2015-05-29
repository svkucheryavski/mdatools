## class and methods for Spectral Unmixing with Purity based method (SUPURE) ##

#' Makes unmixing of set of spectra to pure components.
#'
#' @prop data
#' a matrix with spectral values
#' @prop ncomp
#' number of pure components to identify
#' @prop suggest
#' a sequence with variable numbers suggested as pure variables for each component
#' @prop savgol.dorder
#' a derivative order 
#' @prop savgol.width
#' a width of filter used to smooth signal prior to take a derivative
#' @prop savgol.porder
#' a polynomial order used for smoothing
#' @prop wavelength
#' a vector with wavelength if necessary
#' 
#' @returns
#'  A \code{supure} object - list with several fields, including:
#'   
supure = function(data, ncomp, suggest = NULL, savgol.deriv = 2, 
                  savgol.width = NULL, savgol.porder = 2, wavelength = NULL)
{
   
   res$call = match.call()   
   class(res) = "supure"
   
   res
}

#' Overview plot for supure object
#' 
#' @description
#' Shows a set of plots () for .
#' 
#' @param x
#' a  (object of class \code{pca})
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{supure}} function.
#' 
plot.supure = function(x, ...)
{   
}

#' Print method for XXX
#' 
#' @method print supure
#' @S3method print supure
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a XXX (object of class \code{supure})
#' @param ...
#' other arguments
#' 
print.supure = function(x, ...)
{
}

#' Summary method for XXX model object
#' 
#' @method summary supure
#' @S3method summary supure
#'
#' @description
#' Shows some statistics () for the .
#' 
#' @param object
#' a XXX (object of class \code{supure})
#' @param ...
#' other arguments
#' 
summary.supure = function(object, ...)
{
}
