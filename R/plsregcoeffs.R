# class and methods for PLS regression coefficients #
plsregcoeffs = function(coeffs, ...) UseMethod("plsregcoeffs")

## default method ##
plsregcoeffs.default = function(coeffs)
{ 
   plsregcoeffs = list(values = as.matrix(coeffs))
   plsregcoeffs$call = match.call()
   
   class(plsregcoeffs) = "plsregcoeffs"
   
   plsregcoeffs
}

as.matrix.plsregcoeffs = function(plsregcoeffs, nlv = 0, ...)
{
   if (nlv == 0) 
   { 
      return (plsregcoeffs$values)
   }
   else
   {   
      return (plsregcoeffs$values[, nlv])
   }   
}

print.plsregcoeffs = function(coeffs, nlv = 0, digits = 3, ...)
{
   show(coeffs)
   cat('\nRegression coefficients (class plsregcoeffs)\n')
   if (nlv == 0) { nlv = ncol(coeffs.values)}
   print(round(coeffs$values[, nlv], digits))
}   

plot.plsregcoeffs = function(plsregcoeffs, nlv = 0, main = 'Regression coefficients',
                          xlab = 'Variables', ylab = 'Coefficients',
                          pch = 16, col = 'blue', ...)
{
   
   if (nlv == 0) { nlv = ncol(plsregcoeffs$values)}
   
   coeffs = plsregcoeffs$values[, nlv]
   ncoeff = length(coeffs)
   
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