## class and methods for linear decomposition X = TP' + E ##

#' Linear decomposition of data
#' 
#' @description
#' Creates an object of ldecomp class.
#'
#' @param scores
#' matrix with score values (nobj x ncomp).
#' @param loadings
#' matrix with loading values (nvar x ncomp).
#' @param residuals
#' matrix with data residuals 
#' @param totvar
#' full variance of original data, preprocessed and centered
#' @param tnorm
#' singular values for score normalization
#' @param ncomp.selected
#' number of selected components
#' @param T2
#' matrix with calculated T2 values (e.g. for CV)
#' @param Q2
#' matrix with calculated Q2 values (e.g. for CV)
#' @param cal
#' logical, true if data is for calibration of a LDECOMP based model
#'
#' @return
#' Returns an object (list) of \code{ldecomp} class with following fields:
#' \item{scores }{matrix with score values (nobj x ncomp).}
#' \item{residuals }{matrix with data residuals (nobj x nvar).}
#' \item{T2 }{matrix with T2 distances (nobj x ncomp).}
#' \item{Q2 }{matrix with Q2 distances (nobj x ncomp).}
#' \item{tnorm }{vector with singular values used for scores normalization.}
#' \item{ncomp.selected }{selected number of components.}
#' \item{expvar }{explained variance for each component.}
#' \item{cumexpvar }{cumulative explained variance.}
#' \item{modpower}{modelling power of variables.}
#'
#' @details
#' \code{ldecomp} is a general class for decomposition X = TP' + E. Here, X is a data matrix, 
#' T - matrix with scores, P - matrix with loadings and E - matrix with residuals. It is used, 
#' for example, for PCA results (\code{\link{pcares}}), in PLS and other methods. The class also 
#' includes methods for calculation and plotting residuals, variances, and so on.
#'
#' There is no need to use the \code{ldecomp} manually. For example, when build PCA model 
#' with \code{\link{pca}} or apply it to a new data, the results will automatically inherit 
#' all methods of \code{ldecomp}.
#'
ldecomp = function(scores = NULL, loadings = NULL, residuals = NULL, 
                   totvar, tnorm = NULL, ncomp.selected = NULL,
                   T2 = NULL, Q2 = NULL, cal = TRUE)
{
   if (!is.null(scores))
   {   
      scores = as.matrix(scores)
      rownames(scores) = rownames(residuals)
      colnames(scores) = paste('Comp', 1:ncol(scores))
   }
   
   obj = list(
      scores = scores,
      residuals = residuals, 
      totvar = totvar
   )
   
   if (is.null(ncomp.selected))
      obj$ncomp.selected = ncol(scores)
   else
      obj$ncomp.selected = ncomp.selected
   
   # calculate residual distances and explained variance
   if (is.null(Q2) && is.null(T2) && !is.null(scores) && !is.null(loadings) && !is.null(residuals))
   {   
      res = ldecomp.getDistances(scores, loadings, residuals, tnorm, cal)
      
      if (is.null(Q2))
         obj$Q2 = res$Q2
      
      if (is.null(T2))
         obj$T2 = res$T2
      
      if (is.null(tnorm))
         obj$tnorm = res$tnorm
      
      obj$modpower = res$modpower
   }
   else
   {
      obj$Q2 = Q2
      obj$T2 = T2
      obj$tnorm = tnorm
   }   
   
   var = ldecomp.getVariances(obj$Q2, totvar)   
   obj$expvar = var$expvar
   obj$cumexpvar = var$cumexpvar
   
   obj$call = match.call()
   
   class(obj) = "ldecomp"
   
   obj
}

#' Residuals distances for linear decomposition
#'
#' @description
#' Computes residual distances (Q2 and T2) and modelling power for a data decomposition X = TP' + E.
#' 
#' @param scores
#' matrix with scores (T).
#' @param loadings
#' matrix with loadings (P).
#' @param residuals
#' matrix with residuals (E).
#' @param tnorm
#' vector with singular values for scores normalisation (if NULL will be calculated from \code{scores}).
#' @param cal
#' logical, are these results for calibration set or not
#' 
#' @details
#' The distances are calculated for every 1:n components, where n goes from 1 to ncomp 
#' (number of columns in scores and loadings).
#' 
#' @return
#' Returns a list with Q2, Q2var, T2 and modelling power values for each component.
#'  
ldecomp.getDistances = function(scores, loadings, residuals, tnorm = NULL, cal = TRUE)
{
   ncomp = ncol(scores)
   nobj = nrow(scores)
   nvar = nrow(loadings)
   
   T2 = matrix(0, nrow = nobj, ncol = ncomp)
   Q2 = matrix(0, nrow = nobj, ncol = ncomp)
   modpower = matrix(0, nrow = nvar, ncol = ncomp)
      
   # calculate normalized scores
   if (is.null(tnorm))
      tnorm = sqrt(colSums(scores ^ 2)/(nrow(scores) - 1));
   
   scoresn = sweep(scores, 2L, tnorm, '/', check.margin = F);  

   # calculate variance for data columns
   data = scores %*% t(loadings) + residuals;
   
   if (nobj > 1 && cal == TRUE)
      datasd = sqrt(colSums(data^2)/(nobj - 1))
   
   # calculate distances for each set of components
   for (i in 1:ncomp)
   {
      
      exp = scores[, 1:i, drop = F] %*% t(loadings[, 1:i, drop = F]);
      res = data - exp;
      
      Q2[, i] = rowSums(res^2)
      T2[, i] = rowSums(scoresn[, 1:i, drop = F]^2)
      
      if (nobj > i && cal == TRUE)
         modpower[, i] = 1 - sqrt(colSums(res^2)/(nobj - i - 1))/datasd
   }   
   
   # set dimnames and return results
   colnames(Q2) = colnames(T2) = colnames(modpower) = colnames(scores)
   rownames(Q2) = rownames(T2) = rownames(scores)
   rownames(modpower) = rownames(loadings)
      
   res = list(
      Q2 = Q2,
      T2 = T2,
      modpower = modpower,
      tnorm = tnorm
   )
}


#' Explained variance for linear decomposition
#' 
#' @description
#' Computes explained variance and cumulative explained variance for a data decomposition X = TP' + E.
#'
#' @param Q2
#' Q2 values (squared residuals distance from object to component space).
#' @param totvar
#' Total variance of the original data (after preprocessing).
#' 
#' @return
#' Returns a list with two vectors.
#' 
ldecomp.getVariances = function(Q2, totvar)
{   
   cumresvar = colSums(Q2) / totvar * 100
   cumexpvar = 100 - cumresvar
   expvar = c(cumexpvar[1], diff(cumexpvar))
   
   res = list(
      expvar = expvar,
      cumexpvar = cumexpvar
   )
}

#' Statistical limits for Q2 and T2 residuals
#' 
#' @description
#' Computes statisticsl limits for Q2 and T2 residuals
#' 
#' @param eigenvals
#' vector with eigenvalues 
#' @param nobj
#' number of objects in data
#' @param ncomp
#' number of calculated components
#' @param alpha
#' significance level
#'
#' @details
#' T2 limits are calculated using Hotelling statistics. 
#' 
#' @return
#' Returns a list with two vectors:  \code{T2lim} and \code{Q2lim}.
#' 
ldecomp.getResLimits = function(eigenvals, nobj, ncomp, alpha = 0.05)
{   
   T2lim = matrix(0, nrow = 1, ncol = ncomp)
   for (i in 1:ncomp)
   {
      if (nobj == ncomp)
         T2lim[1, i] = 0
      else
         T2lim[1, i] = (i * (nobj - 1) / (nobj - i)) * qf(1 - alpha, i, nobj - i);  
   }
   
   # calculate Q2 limit using F statistics
   Q2lim = matrix(0, nrow = 1, ncol = ncomp)
   conflim = 100 - alpha * 100;   
   nvar = length(eigenvals)
   
   for (i in 1:ncomp)
   {   
      if (i < nvar)
      {         
         evals = eigenvals[(i + 1):nvar]         
         
         cl = 2 * conflim - 100
         t1 = sum(evals)
         t2 = sum(evals^2)
         t3 = sum(evals^3)
         h0 = 1 - 2 * t1 * t3/3/(t2^2);
         
         if (h0 < 0.001)
            h0 = 0.001
         
         ca = sqrt(2) * erfinv(cl/100)
         h1 = ca * sqrt(2 * t2 * h0^2)/t1
         h2 = t2 * h0 * (h0 - 1)/(t1^2)
         Q2lim[1, i] = t1 * (1 + h1 + h2)^(1/h0)
      }
      else
         Q2lim[1, i] = 0
   }
   
   colnames(T2lim) = colnames(Q2lim) = paste('Comp', 1:ncomp)
   res = list(
      T2lim = T2lim,
      Q2lim = Q2lim
   )   
}   

#' Cumulative explained variance plot for linear decomposition
#' 
#' @description
#' Shows a plot with cumulative explained variance values vs. number of components.
#' 
#' @param obj
#' object of \code{ldecomp} class.
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
plotCumVariance.ldecomp = function(obj, type = 'b', main = 'Cumulative variance',
                                   xlab = 'Components', ylab = 'Explained variance, %',
                                   show.labels = F, ...)
{
   data = cbind(1:length(obj$cumexpvar), obj$cumexpvar)
   if (type != 'h')
   {   
      data = rbind(c(0, 0), data)
      rownames(data) = round(c(0, obj$cumexpvar), 2)
   }
   else
   {
      rownames(data) = round(obj$cumexpvar, 2)      
   }  
   
   colnames(data) = c(xlab, ylab)
   mdaplot(data, main = main, xlab = xlab, ylab = ylab, type = type, show.labels = show.labels, ...)
}

#' Explained variance plot for linear decomposition
#' 
#' @description
#' Shows a plot with explained variance values vs. number of components.
#' 
#' @param obj
#' object of \code{ldecomp} class.
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
plotVariance.ldecomp = function(obj, type = 'b', main = 'Variance',
                                xlab = 'Components', ylab = 'Explained variance, %',
                                show.labels = F, ...)
{
   data = cbind(1:length(obj$expvar), obj$expvar)
   colnames(data) = c(xlab, ylab)
   rownames(data) = round(obj$expvar, 2)
   mdaplot(data, main = main, xlab = xlab, ylab = ylab, show.labels = show.labels, type = type, ...)
}


#' Scores plot for linear decomposition
#' 
#' @description
#' Shows a plot with scores values for data objects.
#' 
#' @param obj
#' object of \code{ldecomp} class.
#' @param comp
#' which components to show the plot for (can be one value or vector with two values).
#' @param main
#' main title for the plot
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
plotScores.ldecomp = function(obj, comp = c(1, 2), main = 'Scores', 
                              show.labels = F, show.axes = F, ...)
{
   if (is.null(obj$scores))
   {
      warning('Scores values are not specified!')
   }   
   else
   {   
      if (length(comp) == 1)
      {   
         # scores vs objects
         data = cbind(1:nrow(obj$scores), obj$scores[, comp])      
         colnames(data) = c('Objects', colnames(obj$scores)[comp])
         rownames(data) = rownames(obj$scores)
         
         mdaplot(data, main = main, show.labels = show.labels, ...)
      }
      else if (length(comp) == 2)
      {
         # scores vs scores
         data = obj$scores[, c(comp[1], comp[2])]   
         
         if (show.axes == T)
            show.lines = c(0, 0)      
         else
            show.lines = F
         
         mdaplot(data, main = main, show.labels = show.labels, show.lines = show.lines, ...)
      }
      else
      {
         stop('Wrong number of components!')
      }   
   }
}  

#' Residuals plot for linear decomposition
#' 
#' @description
#' Shows a plot with T2 vs Q2 values for data objects.
#' 
#' @param obj
#' object of \code{ldecomp} class.
#' @param ncomp
#' what number of components to show the plot for (if NULL, model selected value will be used).
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.limits
#' logical, show or not lines for statistical limits of the residuals
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
plotResiduals.ldecomp = function(obj, ncomp = NULL, main = NULL, xlab = 'T2', ylab = 'Q2', 
                                 show.labels = F, show.limits = T, ...)
{
   if (is.null(main))
   {
      if (is.null(ncomp))
         main = 'Residuals'
      else
         main = sprintf('Residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
      
   if (show.limits == T)
      show.lines = c(obj$T2lim[1, ncomp], obj$Q2lim[1, ncomp])
   else
      show.lines = F

   data = cbind(obj$T2[, ncomp], obj$Q2[, ncomp])
   colnames(data) = c(xlab, ylab)
   mdaplot(data, main = main, xlab = xlab, ylab = ylab, show.labels = show.labels, 
           show.lines = show.lines, ...)
}  

#' Print method for linear decomposition
#'
#' @method print ldecomp
#' @S3method print ldecomp
#' 
#' @description
#' Generic \code{print} function for linear decomposition. Prints information about the \code{ldecomp}
#' object.
#' 
#' @param x
#' object of class \code{ldecomp}
#' @param str
#' user specified text to show as a description of the object
#' @param ...
#' other arguments
#' 
print.ldecomp = function(x, str = NULL, ...)
{   
   if (is.null(str))
      str ='Results of data decomposition (class ldecomp)'
   
   if (nchar(str) > 0)   
      cat(sprintf('\n%s\n', str))
   
   cat('\nMajor fields:\n')   
   cat('$scores - matrix with score values\n')
   cat('$T2 - matrix with T2 distances\n')
   cat('$Q2 - matrix with Q2 residuals\n')
   cat('$ncomp.selected - selected number of components\n')
   cat('$expvar - explained variance for each component\n')
   cat('$cumexpvar - cumulative explained variance\n')
}

#' as.matrix method for ldecomp object
#' 
#' @method as.matrix ldecomp
#' @S3method as.matrix ldecomp
#'
#' @description
#' Generic \code{as.matrix} function for linear decomposition. Returns a matrix with information 
#' about the decomposition.
#' 
#' @param x
#' object of class \code{ldecomp}
#' @param ...
#' other arguments
#' 
as.matrix.ldecomp = function(x, ...)
{
   data = cbind(x$expvar, x$cumexpvar)   
   colnames(data) = c('Expvar', 'Cumexpvar')   
   data
}  

#' Summary statistics for linear decomposition
#'
#' @method summary ldecomp
#' @S3method summary ldecomp
#' 
#' @description
#' Generic \code{summary} function for linear decomposition. Prints statistic about the decomposition.
#' 
#' @param object
#' object of class \code{ldecomp}
#' @param str
#' user specified text to show as a description of the object
#' @param ...
#' other arguments
#' 
summary.ldecomp = function(object, str = NULL, ...)
{
   if (is.null(str))
      str ='Summary for data decomposition (class ldecomp)'
   
   cat('\n')
   cat(str, '\n')
   cat(sprintf('\nSelected components: %d\n\n', object$ncomp.selected))      
   
   data = as.matrix(object)
   print(round(data, 2))   
}

#' Inverse error function 
#' 
#' @param x
#' a matrix or vector with data values
erfinv = function (x) qnorm((1 + x)/2)/sqrt(2)


