#' Regression coefficients

#' @description
#' class for storing and visualisation of regression coefficients
#' for regression models
#' 
#' @param coeffs
#' vector or matrix with regression coefficients
#' @param ci.coeffs
#' array (nobj x ncomp x ny x cv) with regression coefficients for 
#' computing confidence intervals (e.g. from jack-knifing)
#' @param ci.alpha
#' significance level for computing of the confidence intervals 
#' 
#' @return
#' a list (object of \code{regcoeffs} class) with fields, including:
#' \tabular{ll}{
#'    \code{values} \tab an array (nvar x ncomp x ny) with regression coefficients \cr
#'    \code{ci} \tab an array (nvar x ncomp x ny) with confidence intervals for coefficients\cr
#'    \code{p.values} \tab an array (nvar x ncomp x ny) with p-values for coefficients \cr
#' }
#' last two fields are available if proper values for calculation of the statistics were provided.
#' 
regcoeffs = function(coeffs, ci.coeffs = NULL, ci.alpha = 0.1)
{   

   regcoeffs = list(values = coeffs)
   if (!is.null(ci.coeffs))
   {
      stat = regcoeffs.getStat(coeffs, ci.coeffs, ci.alpha)
      regcoeffs$ci = stat$ci
      regcoeffs$t.values = stat$t.values
      regcoeffs$p.values = stat$p.values
      regcoeffs$alpha = ci.alpha
   }   
   regcoeffs$call = match.call()
   
   class(regcoeffs) = "regcoeffs"
   
   regcoeffs
}

#' Confidence intervals and p-values for regression coeffificents
#' 
#' @description
#' calculates confidence intervals and t-test based p-values for 
#' regression coefficients based on jack-knifing procedure
#' 
#' @param obj
#' regression coefficients array for a model
#' @param ci.coeffs
#' array with regression coefficients for calculation of condifence intervals
#' @param ci.alpha
#' significance level to calculate the confidence intervals
#' 
#' @return
#' a list with statistics (\code{$ci} - array with confidence intervals, 
#' \code{$p.values} - array with p-values, \code{$t.values} - array with t-values)
#' 
regcoeffs.getStat = function(obj, ci.coeffs, ci.alpha = 0.1)
{
   s = dim(ci.coeffs)
   nvar = s[1]
   ncomp = s[2]
   ny = s[3]
   nobj = s[4]

   t = qt(1 - ci.alpha/2, nobj - 1)

   ci = array(0, dim = c(nvar, ncomp, ny, 2))
   t.values = array(0, dim = c(nvar, ncomp, ny))
   p.values = array(0, dim = c(nvar, ncomp, ny))
   for (y in 1:ny)
   {
      for (comp in 1:ncomp)
      {
         coeffs = ci.coeffs[, comp, y, ]
         m = apply(coeffs, 1, mean)
         ssq = apply(t(scale(t(coeffs), center = m, scale = FALSE))^2, 1, sum)
         se = sqrt( (nobj - 1)/nobj * ssq )
         ci[, comp, y, ] = cbind(m - t * se, m + t * se)
         tvals = m/se
         tmin = apply(cbind(tvals, -tvals), 1, min)
         t.values[, comp, y] = tvals
         p.values[, comp, y] = 2 * pt(tmin, nobj - 1)
      }   
   }   
   
   dimnames(t.values) = dimnames(obj)
   dimnames(p.values) = dimnames(obj)
   
   stat = list(
      ci = ci,
      t.values = t.values,
      p.values = p.values
      )
}

#' as.matrix method for regression coefficients class
#' 
#' @method as.matrix regcoeffs
#' @S3method as.matrix regcoeffs
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
#' @method print regcoeffs
#' @S3method print regcoeffs
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
#' @param col
#' vector with colors for the plot (vector or one value)
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.line
#' logical, show or not line for 0 value
#' @param show.ci
#' logical, show or not confidence intervals if they are available
#' @param ...
#' other arguments
#' 
plot.regcoeffs = function(x, ncomp = 1, ny = 1, type = NULL, col = NULL, 
                          main = 'Regression coefficients',
                          xlab = 'Variables', ylab = 'Coefficients', show.line = T, 
                          show.ci = F, ...)
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
   if (show.ci == F || is.null(obj$ci))
   {   
      rownames(data) = rownames(obj$values)
      mdaplot(data, type = type, main = main, xlab = xlab, ylab = ylab, 
              show.grid = T, show.lines = show.line, ...)
   }
   else
   {
      ci = obj$ci[, ncomp, ny, ]
      rownames(data) = rownames(obj$values)
      
      if (type == 'l')
      {   
         if (is.null(col))
         {
            cc = mdaplot.getColors(ngroups = 2)
            cg = mdaplot.getColors(ngroups = 4, colmap = 'gray')
            col = c(cg[1], cc[1], cc[1])            
         }
         
         type = c('l', 'l', 'l')
         mdata = cbind(1:ncoeff, coeffs, ci[, 1], ci[, 2])
      }   
      else
      {
         mdata = list()
         mdata[[1]] = data 
         mdata[[2]] = cbind(1:ncoeff, coeffs, ci[, 1], ci[, 2])         
         type = c(type, 'e')

         if (is.null(col))
         {
            cc = mdaplot.getColors(ngroups = 2)
            cg = mdaplot.getColors(ngroups = 4, colmap = 'gray')
            col = c(cg[1], cc[1])            
         }   
      }   
      mdaplotg(mdata, type = type, main = main, xlab = xlab, ylab = ylab, 
              show.grid = T, show.lines = show.line, colmap = col, ...)
   }   
}   
