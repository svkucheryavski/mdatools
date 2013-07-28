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

as.matrix.regcoeffs = function(obj, ncomp = 1, ny = 1)
{
   # Returns values of regression coefficients as a matrix 
   #
   # Arguments:
   #   obj: object of class "regcoeffs" 
   #   ncomp.selected: which column of the matrix to show (number of components for PLS/PCR) 
   #
   # Returns:
   #   res: matrix with regression coefficients
   
   return (obj$values[, ncomp, ny, drop = F])
}

print.regcoeffs = function(obj, ncomp.selected = 1, ny = 1, digits = 3)
{
   # Shows regression coefficient values
   #
   # Arguments:
   #   obj: object of class "regcoeffs" 
   #   ncomp.selected: which column of the matrix to show (number of components for PLS/PCR) 
   #   digits: how many decimal numbers to show
   
   cat('\nRegression coefficients (class plsregcoeffs)\n')
   print(round(obj$values[, ncomp.selected, ny, drop = F], digits))
}   

plot.regcoeffs = function(obj, ncomp = 1, ny = 1, main = 'Regression coefficients', type = NULL,
                          xlab = 'Variables', ylab = 'Coefficients', show.lines = T, ...)
{
   # Makes regression coefficients plot
   #
   # Arguments:
   #   obj: object of class "regcoeffs" 
   #   ncomp.selected: which column of the matrix to show (number of components for PLS/PCR) 
   #   main: main title for the plot
   #   xlab: label for X axis
   #   ylab: label for Y axus
   #   show.line: logical, show or not zero line
   
   coeffs = obj$values[, ncomp, ny, drop = F]
   ncoeff = length(coeffs)
   
   if (show.lines == T)
      show.lines = c(NA, 0)
   
   if (is.null(type))
   {   
      if (ncoeff < 30)
         type = 'b'
      else
         type = 'l'
   }
   
   mdaplot(cbind(1:ncoeff, coeffs), type = type, main = main, xlab = xlab, ylab = ylab, 
           show.grid = T, show.line = c(NA, 0), ...)
}   
