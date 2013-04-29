# class and methods for regression coefficients #
regcoeffs = function(coeffs, ...) UseMethod("regcoeffs")

## default method ##
regcoeffs.default = function(coeffs)
{ 
   regcoeffs = list(values = coeffs)
   regcoeffs$call = match.call()
   
   class(regcoeffs) = "regcoeffs"
   
   regcoeffs
}

as.matrix.regcoeffs = function(regcoeffs, ...)
{
   return (regcoeffs$values)
}

print.regcoeffs = function(regcoeffs, digits = 3, ...)
{
   cat('\nRegression coefficients (class regcoeffs)\n')
   print(round(regcoeffs$values, digits))
}   

plot.regcoeffs = function(regcoeffs, main = 'Regression coefficients',
                          xlab = 'Variables', ylab = 'Coefficients',
                          pch = 16, col = 'blue', ...)
{
   
   
   # remove Intercept
   ncoeff = length(regcoeffs$values)
   coeffs = regcoeffs$values[2:ncoeff]
   ncoeff = ncoeff - 1
   
   # select limits for y axis
   ylim = max(abs(coeffs))
   
   # chose plot type depending on number of coefficients
   if (ncoeff < 30) { type = 'b' }else{ type = 'l' } 
   
   # show plot
   plot(coeffs, type = type, col = col, pch = pch,
        main = main,
        xlab = xlab,
        ylab = ylab,
        ylim = c(-ylim, ylim),
        axes = F
   )
   abline(c(0, 0), c(0, ncoeff))
   
   if (is.null(dim(coeffs))) { names = names(coeffs)} 
   else {names = rownames(coeffs)}
   
   # show axes and labels if needed
   if (ncoeff > 20)
   {   
      atx = seq(1, ncoeff, ncoeff/10)
   }
   else
   {
      atx = 1:ncoeff
   }   
   axis(1, at = atx, labels = names[atx], cex.axis = 0.85)      
   axis(2, cex.axis = 0.85)      
   if (ncoeff < 30)
   {
      text(1:ncoeff, coeffs, names, cex = 0.6, pos = 3, col = 'gray')
   }
   grid()
   box()
}   