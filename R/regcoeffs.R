#' Regression coefficients
#' 
#' @description
#' A class for storing and visualising of regression coefficients for any regression model.
#' 
#' @param coeffs
#' vector or matrix with regression coefficients
#' 
regcoeffs = function(coeffs)
{ 
   # A class for storing and visualising of regression coefficients 
   # Can be used for any regression model
   #
   # Arguments:
   #   coeffs: vector or matrix with regression coefficients
   #
   # Returns:
   # a list, object of the class "regcoeffs", with following fields
   #   values: matrix with regression coefficient values    
  

   regcoeffs = list(values = coeffs)
   regcoeffs$call = match.call()
   
   class(regcoeffs) = "regcoeffs"
   
   regcoeffs
}

#' as.matrix method for regression coefficients class
#' 
#' @description
#' returns matrix with regression coeffocoents for given response number and amount of components
#' 
#' @param x
#' regression coefficients object (class \code{regcoeffs})
#' @param ncomp
#' number of components to return the coefficients for
#' @param ny
#' number of response variable to return the coefficients for
#' @param ...
#' other arguments
#' 
as.matrix.regcoeffs = function(x, ncomp = 1, ny = 1, ...)
{
   return (x$values[, ncomp, ny, drop = F])
}

#' print method for regression coefficients class
#' 
#' @description
#' prints regression coeffocoent values for given response number and amount of components
#' 
#' @param x
#' regression coefficients object (class \code{regcoeffs})
#' @param ncomp
#' number of components to return the coefficients for
#' @param ny
#' number of response variable to return the coefficients for
#' @param digits
#' decimal digits round the coefficients to
#' @param ...
#' other arguments
#' 
print.regcoeffs = function(x, ncomp = 1, ny = 1, digits = 3, ...)
{
   obj = x
   
   cat('\nRegression coefficients (class plsregcoeffs)\n')
   print(round(obj$values[, ncomp, ny, drop = F], digits))
}   

#' Regression coefficients plot
#' 
#' @description
#' Shows plot with regression coefficient values for every predictor variable (x)
#' 
#' @param x
#' regression coefficients object (class \code{regcoeffs})
#' @param ncomp
#' number of components to return the coefficients for
#' @param ny
#' number of response variable to return the coefficients for
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.line
#' logical, show or not line for 0 value
#' @param ...
#' other arguments
#' 
plot.regcoeffs = function(x, ncomp = 1, ny = 1, type = NULL, main = 'Regression coefficients',
                          xlab = 'Variables', ylab = 'Coefficients', show.line = T, ...)
{
   obj = x
   
   coeffs = obj$values[, ncomp, ny, drop = F]
   ncoeff = length(coeffs)
   
   if (show.line == T)
      show.line = c(NA, 0)
   
   if (is.null(type))
   {   
      if (ncoeff < 30)
         type = 'b'
      else
         type = 'l'
   }
  
   data = cbind(1:ncoeff, coeffs)
   rownames(data) = rownames(obj$values)

   mdaplot(data, type = type, main = main, xlab = xlab, ylab = ylab, 
           show.grid = T, show.lines = show.line, ...)
}   
