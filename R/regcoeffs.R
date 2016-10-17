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
regcoeffs = function(coeffs, ci.coeffs = NULL, ci.alpha = 0.1) {   
   regcoeffs = list()
   regcoeffs$values = coeffs
   if (!is.null(ci.coeffs)) {
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
#' @param coeffs.values
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
regcoeffs.getStat = function(coeffs.values, ci.coeffs, ci.alpha = 0.1) {

   # get attributes   
   attrs = mda.getattr(coeffs.values)
   exclvars = attrs$exclrows
   nexclvar = length(exclvars)

   # get dimensions and set t-value
   nvar = dim(ci.coeffs)[1]
   ncomp = dim(ci.coeffs)[2]
   ny = dim(ci.coeffs)[3]
   nobj = dim(ci.coeffs)[4]
   t = qt(1 - ci.alpha/2, nobj - 1)
 
   # set up matrices and calculate statistics 
   ci = array(0, dim = c(nvar, ncomp, ny, 2))
   t.values = array(0, dim = c(nvar, ncomp, ny))
   p.values = array(0, dim = c(nvar, ncomp, ny))
   for (y in 1:ny) {
      for (comp in 1:ncomp) {
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

   if (nexclvar > 0) {
      ci.out = array(0, dim = c(nvar + nexclvar, ncomp, ny, 2))
      t.values.out = array(0, dim = c(nvar + nexclvar, ncomp, ny))
      p.values.out = array(0, dim = c(nvar + nexclvar, ncomp, ny))
      ci.out[-exclvars, , ,] = ci
      t.values.out[-exclvars, , ] = t.values
      p.values.out[-exclvars, , ] = p.values
   } else {
      ci.out = ci
      t.values.out = t.values
      p.values.out = p.values
   }  
   
   dimnames(t.values.out) = dimnames(p.values.out) = dimnames(coeffs.values)
   dimnames(ci.out) = c(dimnames(coeffs.values), list(c('Lo', 'Up')))
   t.values.out = mda.setattr(t.values.out, attrs)
   p.values.out = mda.setattr(p.values.out, attrs)
   ci.out = mda.setattr(ci.out, attrs)
   attr(t.values.out, 'name') = 't-values (Jack-knife)'
   attr(p.values.out, 'name') = 'p-values (Jack-knife)'
   attr(t.values.out, 'name') = sprintf('%d%% confidence interval', round((1 - ci.alpha)*100))
   
   stat = list(
      ci = ci.out,
      t.values = t.values.out,
      p.values = p.values.out
   )
   
   stat
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
#' @export
as.matrix.regcoeffs = function(x, ncomp = 1, ny = 1, ...) {
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
#' @export
print.regcoeffs = function(x, ncomp = 1, ny = 1, digits = 3, ...) {
   obj = x
   
   cat('\nRegression coefficients (class regcoeffs)\n')
   cat('\nCall:\n')
   print(obj$call)
   cat('\nMajor fields:\n')
   cat('$values - array with regression coefficients\n')
   cat('$ci - matrix with confidence intervals\n')
   cat('$t - vector with t-values\n')
   cat('$p - vector with p-values\n')
   cat('$alpha - significance level used to calculate the coefficients\n')
   cat('\nThe last four fields available only if Jack-Knife was used.\n')
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
#' @param ylab
#' label for y axis
#' @param show.line
#' logical, show or not line for 0 value
#' @param show.ci
#' logical, show or not confidence intervals if they are available
#' @param ...
#' other arguments
#' 
#' @export
plot.regcoeffs = function(x, ncomp = 1, ny = 1, type = NULL, col = NULL, main = NULL, ylab = NULL, 
                          show.line = T, show.ci = T, ...) {

   obj = x
   attrs = mda.getattr(obj$values)
   
   if (is.null(main)) {   
      if (is.null(ncomp))
         main = 'Regression coefficients'
      else
         main = sprintf('Regression coefficients (ncomp = %d)', ncomp)
   }
   
   if (is.null(ylab)) {   
      if (dim(obj$values)[3] == 1 || is.null(dimnames(obj$values)[[3]]))
         ylab = 'Coefficients'
      else
         ylab = sprintf('Coefficients (%s)', dimnames(obj$values)[[3]][ny])
   }
   
   if (show.line == T)
      show.line = c(NA, 0)
   ncoeff = nrow(obj$values)
   
   if (is.null(type)) {   
      if (ncoeff < 30)
         type = 'b'
      else
         type = 'l'
   }
  
   data = matrix(obj$values[, ncomp, ny, drop = F], ncol = 1)
   data = mda.setattr(data, attrs)
   rownames(data) = rownames(obj$values)
   if (show.ci == TRUE && !is.null(obj$ci)) {   
      ci.col = mdaplot.getColors(1)
      main.col   = 'lightgray'
      err.margin = matrix(obj$ci[, ncomp, ny, 2], ncol = 1) - data
      err.maring = mda.setattr(err.margin, attrs)
      attr(data, 'name') = 'Regression coefficients'      
      
      if (type == 'l')
         mdaplotg(list(data, data + err.margin, data - err.margin), type = c('l', 'l', 'l'), 
                  main = main, ylab = ylab, show.legend = F, colmap = c(main.col, ci.col, ci.col), 
                  show.grid = T, show.lines = show.line, ...)
     else
        mdaplotg(list(data, mda.t(mda.cbind(data, err.margin))), type = c(type, 'e'), 
                 main = main, ylab = ylab, show.legend = F, 
                 colmap = c(main.col, ci.col), show.grid = T, show.lines = show.line, ...)
      
   } else {
      main.col = ifelse(is.null(col), mdaplot.getColors(1), col)
      mdaplot(data, type = type, show.grid = T, main = main, ylab = ylab, show.lines = show.line, 
              col = main.col, ...)
   }
}   
