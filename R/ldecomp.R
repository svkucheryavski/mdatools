## class and methods for linear decomposition X = TP' + E ##

ldecomp = function(scores = NULL, loadings = NULL, residuals = NULL, 
                   totvar, tnorm = NULL, ncomp.selected = NULL,
                   T2 = NULL, Q2 = NULL)
{
   # Creates an object of ldecomp class.
   #
   # The object is needed to store and visualise results for decomposition X = TP' + E 
   # In case of cross-validated results, only distances and variances are stored
   #
   # Arguments:
   #   scores: matrix with score values (nobj x ncomp).
   #   loadings: matrix with loading values (nvar x ncomp).
   #   residuals: matrix with data residuals 
   #   totvar: full variance of original data, preprocessed and centered
   #   tnorm: singular values for score normalization
   #   ncomp.selected: number of selected components
   #
   # Returns:
   #  object (list) of class ldecomp with following fields:   
   #   obj$ncomp.selected: selected number of components
   #   obj$scores: matrix with score values (nobj x ncomp).
   #   obj$Q2: matrix with Q2 residuals (nobj x ncomp).
   #   obj$T2: matrix with T2 distances (nobj x ncomp)  
   #   obj$totvar: total variance of the data
   
   if (!is.null(scores))
   {   
      scores = as.matrix(scores)
      rownames(scores) = rownames(residuals)
      colnames(scores) = paste('Comp', 1:ncol(scores))
   }
   
   obj = list(
      scores = scores,
      totvar = totvar
   )
   
   if (is.null(ncomp.selected))
      obj$ncomp.selected = ncol(scores)
   else
      obj$ncomp.selected = ncomp.selected
   
   # calculate residual distances and explained variance
   if (is.null(Q2) && is.null(T2) && !is.null(scores) && !is.null(loadings) && !is.null(residuals))
   {   
      res = ldecomp.getDistances(scores, loadings, residuals, tnorm)
      
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

ldecomp.getDistances = function(scores, loadings, residuals, tnorm = NULL)
{
   # Computes residual distances (Q2 and T2) for a decomposition.
   # The distances are calculated for every 1:n components, where n
   # goes from 1 to ncomp (number of columns in scores and loadings)
   #
   # Arguments:
   #   scores: matrix with scores (nobj x ncomp).
   #   loadings: matrix with loadings (nvar x ncomp)
   #   residuals: matrix with data residuals 
   #   tnorm: singular values for normalizing scores
   #
   # Returns:
   #   res$Q2: matrix with Q2 residuals (nobj x ncomp).
   #   res$T2: matrix with T2 distances (nobj x ncomp)   
   #   res$mpower: modelling power (nvar x ncomp)
   
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
   datasd = apply(data, 2, sd)
   
   # calculate distances for each set of components
   for (i in 1:ncomp)
   {
      
      exp = scores[, 1:i, drop = F] %*% t(loadings[, 1:i, drop = F]);
      res = data - exp;
      
      Q2[, i] = rowSums(res^2)
      T2[, i] = rowSums(scoresn[, 1:i, drop = F]^2)
      modpower[, i] = 1 - apply(res, 2, sd)/datasd
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

ldecomp.getVariances = function(Q2, totvar)
{   
   # Computes explained variance and cumulative explained variance 
   # for every component of a decomposition.
   #
   # Arguments:
   #   Q2: matrix with Q2 values
   #   totvar: total variance of the data
   #
   # Returns:
   #   res$expvar: vector with explained variance for every component
   #   res$cumexpvar: vector with cumulative explained variance
   
   cumresvar = colSums(Q2) / totvar * 100
   cumexpvar = 100 - cumresvar
   expvar = c(cumexpvar[1], diff(cumexpvar))
   
   res = list(
      expvar = expvar,
      cumexpvar = cumexpvar
   )
}

ldecomp.getResLimits = function(eigenvals, nobj, ncomp, alpha = 0.05)
{   
   # Computes statistical limits for Q2 residuals and T2 distances.
   #
   # Arguments:
   #   eigenvals: vector with eigenvalues for a model.
   #   nobj: number of objects in calibration data
   #   ncomp: number of selected components 
   #   alpha: significance level for Q2 limits
   #
   # Returns:
   #   res$Q2lim: limit for Q2 residuals
   #   res$T2lim: limit for T2 distances
   
   # calculate T2 limit using Hotelling statistics
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

plotCumVariance.ldecomp = function(obj, type = 'b', main = 'Cumulative variance',
                                   xlab = 'Components', ylab = 'Explained variance, %',
                                   show.labels = F, ...)
{
   # Shows cumulative explained variance plot.
   #
   # Arguments:
   #   obj: object of ldecomp class.
   #   type: type of the plot
   #   main: main title for the plot
   #   xlab: text for x axis label
   #   ylab: text for y axis label
   #   show.labels: show or not labels for plot points
   
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

plotVariance.ldecomp = function(obj, type = 'b', main = 'Variance',
                                xlab = 'Components', ylab = 'Explained variance, %',
                                show.labels = F, ...)
{
   # Shows explained variance plot.
   #
   # Arguments:
   #   obj: object of ldecomp class
   #   type: type of the plot
   #   main: main title for the plot
   #   xlab: text for x axis label
   #   ylab: text for y axis label
   #   show.labels: logical, show or not labels for plot points
   
   
   data = cbind(1:length(obj$expvar), obj$expvar)
   colnames(data) = c(xlab, ylab)
   rownames(data) = round(obj$expvar, 2)
   mdaplot(data, main = main, xlab = xlab, ylab = ylab, show.labels = show.labels, type = type, ...)
}

plotScores.ldecomp = function(obj, comp = c(1, 2), main = 'Scores', 
                              show.labels = F, show.axes = F, ...)
{
   # Shows scores plot.
   #
   # Arguments:
   #   obj: object of ldecomp class.
   #   comp: which components to show on x and y axis. 
   #   cgroup: variable for color grouping of plot points.
   #   main: main title for the plot
   #   show.labels: show or not labels for plot points.
   #   show.colorbar: show or not a colorbar legend if cgroup is provided.
   #   show.axes: show or not axes crossing (0, 0) point.
   
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

plotResiduals.ldecomp = function(obj, ncomp = NULL, main = NULL, xlab = 'T2', ylab = 'Q2', 
                                 show.labels = F, show.limits = T, ...)
{
   # Shows T2 vs Q2 residuals plot.
   #
   # Arguments:
   #   obj: object of ldecomp class.
   #   ncomp: number of components for the residuals 
   #   main: main title for the plot
   #   xlab: text for x axis label
   #   ylab: text for y axis label
   #   show.labels: show or not labels for plot points.
   
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

print.ldecomp = function(obj, str = NULL, ...)
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

as.matrix.ldecomp = function(obj)
{
   data = cbind(obj$expvar, obj$cumexpvar)   
   colnames(data) = c('Expvar', 'Cumexpvar')   
   data
}  

summary.ldecomp = function(obj, str = NULL)
{
   if (is.null(str))
      str ='Summary for data decomposition (class ldecomp)'
   
   cat('\n')
   cat(str, '\n')
   cat(sprintf('\nSelected components: %d\n\n', obj$ncomp.selected))      
   
   print(round(data, 2))   
}

erfinv = function (x) qnorm((1 + x)/2)/sqrt(2)


