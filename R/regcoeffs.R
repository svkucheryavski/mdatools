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
regcoeffs = function(coeffs, ci.coeffs = NULL) {   
   regcoeffs = list()
   regcoeffs$values = coeffs
   if (!is.null(ci.coeffs)) {
      stat = regcoeffs.getStat(coeffs, ci.coeffs)
      regcoeffs$se = stat$se
      regcoeffs$t.values = stat$t.values
      regcoeffs$p.values = stat$p.values
      regcoeffs$niter = stat$niter
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
#' 
#' @return
#' a list with statistics (\code{$ci} - array with confidence intervals, 
#' \code{$p.values} - array with p-values, \code{$t.values} - array with t-values)
#' 
regcoeffs.getStat = function(coeffs.values, ci.coeffs) {

   # get attributes   
   attrs = mda.getattr(coeffs.values)
   exclvars = attrs$exclrows
   nexclvar = length(exclvars)

   # get dimensions and set t-value
   nvar = dim(ci.coeffs)[1]
   ncomp = dim(ci.coeffs)[2]
   ny = dim(ci.coeffs)[3]
   nobj = dim(ci.coeffs)[4]
 
   # set up matrices and calculate statistics 
   t.values = array(0, dim = c(nvar, ncomp, ny))
   p.values = array(0, dim = c(nvar, ncomp, ny))
   se = array(0, dim = c(nvar, ncomp, ny))
   for (y in 1:ny) {
      for (comp in 1:ncomp) {
         coeffs = ci.coeffs[, comp, y, ]
         m = apply(coeffs, 1, mean)
         ssq = apply(t(scale(t(coeffs), center = m, scale = FALSE))^2, 1, sum)
         se[, comp, y] = sqrt( (nobj - 1)/nobj * ssq )
         tvals = m/se[, comp, y]
         tmin = apply(cbind(tvals, -tvals), 1, min)
         t.values[, comp, y] = tvals
         p.values[, comp, y] = 2 * pt(tmin, nobj - 1)
      }   
   }   

   if (nexclvar > 0) {
      se.out = array(0, dim = c(nvar + nexclvar, ncomp, ny))
      t.values.out = array(0, dim = c(nvar + nexclvar, ncomp, ny))
      p.values.out = array(1, dim = c(nvar + nexclvar, ncomp, ny))
      se.out[-exclvars, , ] = se
      t.values.out[-exclvars, , ] = t.values
      p.values.out[-exclvars, , ] = p.values
   } else {
      se.out = se
      t.values.out = t.values
      p.values.out = p.values
   }  
   
   dimnames(t.values.out) = dimnames(p.values.out) = dimnames(coeffs.values)
   dimnames(se.out) = c(dimnames(coeffs.values))
   t.values.out = mda.setattr(t.values.out, attrs)
   p.values.out = mda.setattr(p.values.out, attrs)
   se.out = mda.setattr(se.out, attrs)
   attr(se.out, 'name') = 'Standard error (Jack-knife)'
   attr(t.values.out, 'name') = 't-values (Jack-knife)'
   attr(p.values.out, 'name') = 'p-values (Jack-knife)'

   stat = list(
      se = se.out,
      t.values = t.values.out,
      p.values = p.values.out,
      niter = nobj
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
   cat('$se - matrix with standard errors\n')
   cat('$t - matrix with t-values\n')
   cat('$p - matrix with p-values\n')
   cat('\nThe last three fields available only if Jack-Knife was used.\n')
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
#' @param alpha
#' significance level for confidence intervals (a number between 0 and 1, e.g. for 95\% alpha = 0.05)
#' @param ...
#' other arguments
#' 
#' @export
plot.regcoeffs = function(x, ncomp = 1, ny = 1, type = NULL, col = NULL, main = NULL, ylab = NULL, 
                          show.line = T, show.ci = T, alpha = 0.05, ...) {

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
   if (show.ci == TRUE && !is.null(obj$se)) {   
      ci.col = mdaplot.getColors(1)
      main.col   = 'lightgray'
      err.margin = matrix(obj$se[, ncomp, ny], ncol = 1) * qt(1 - alpha/2, obj$niter)
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

#' Summary method for regcoeffs object
#' 
#' @description
#' Shows estimated coefficients and statistics (if available).
#' 
#' @param object
#' object of class \code{regcoeffs}
#' @param ncomp
#' how many components to use 
#' @param ny
#' which y variable to show the summary for 
#' @param alpha
#' significance level for confidence interval (if statistics available) 
#' @param ...
#' other arguments
#' 
#' @details 
#' Statistcs are shown if Jack-Knifing was used when model is calibrated.
#' 
#' @export
summary.regcoeffs = function(object, ncomp = 1, ny = 1, alpha = 0.05, ...) {
   obj = object
   nxvar = dim(obj$values)[1]
   nyvar = dim(obj$values)[3]
   
   if (!is.numeric(ncomp) || length(ncomp) != 1 || ncomp <= 0 || ncomp > dim(obj$values)[2]) 
      stop('Wrong value for number of components!')
   
   if (!is.numeric(ny) || length(ny) != 1 || ny <= 0 || ny > nyvar) 
      stop('Wrong value for number of components!')

   if (alpha <= 0 || alpha >= 1)
      stop('Wrong value for alpha (must be between 0 and 1)')
   
   attrs = mda.getattr(obj$values)
   coeffs = obj$values[, ncomp, ny, drop = F]
   dim(coeffs) = c(nxvar, 1)
   
   if (!is.null(obj$p.values)) {
      
      ci.CL = (1 - alpha) * 100
      t = qt(1 - alpha/2, obj$niter - 1)
      
      se = matrix(obj$se[, ncomp, ny], ncol = 1)
      ci.lo = coeffs - t * se
      ci.up = coeffs + t * se
      
      coeffs = data.frame(coeffs, 
                     round(obj$t.values[, ncomp, ny], 1),
                     se,
                     round(obj$p.values[, ncomp, ny], 3),
                     ci.lo,
                     ci.up
      )
      colnames(coeffs) = c('Estimated coefficients', 't-value', 'SE', 'p-value', 
                           sprintf('%d%% CI (lo)', ci.CL),
                           sprintf('%d%% CI (up)', ci.CL)
      )
   } else {
      colnames(coeffs) = 'Estimated coefficients'
   }
   
   rownames(coeffs) = dimnames(obj$values)[[1]]
   attr(coeffs, 'exclrows') = attrs$exclrows
   attr(coeffs, 'name') = paste('Regression coefficients for', dimnames(obj$values)[[3]][ny])
   
   mda.show(coeffs)
   if (ncol(coeffs) > 1)
      cat(sprintf('\nDegrees of freedom: %d, t-value: %.2f\n', obj$niter - 1, t))
   cat('\n\n')
}
