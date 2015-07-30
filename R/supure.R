## class and methods for Spectral Unmixing with Purity based method (SUPURE) ##

#' Makes unmixing of set of spectra to pure components.
#'
#' @param data
#' a matrix with spectral values
#' @param  ncomp
#' number of pure components to identify
#' @param  suggest
#' a sequence with variable numbers suggested as pure variables for each component
#' @param  savgol.deriv
#' a derivative order 
#' @param  savgol.width
#' a width of filter used to smooth signal prior to take a derivative
#' @param  savgol.porder
#' a polynomial order used for smoothing
#' @param  wavelength
#' a vector with wavelength if necessary
#' 
#' @return
#'  A \code{supure} object - list with several fields, including:
#'  
#' @export 
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
#' @export
plot.supure = function(x, ...)
{   
}

#' Print method for SUPURE object
#' 
#' @method print supure
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' unmixing results (object of class \code{supure})
#' @param ...
#' other arguments
#'
#' @export 
print.supure = function(x, ...)
{
}

#' Summary method for SUPURE object
#' 
#' @description
#' Shows some statistics () for the .
#' 
#' @param object
#' unmixing results (object of class \code{supure})
#' @param ...
#' other arguments
#' 
#' @export 
summary.supure = function(object, ...)
{
}
