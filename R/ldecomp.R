# class and methods for linear decomposition X = TP' + E #
ldecomp = function(scores, loadings, residuals, fullvar, ...) UseMethod("ldecomp")

ldecomp.default = function(scores, loadings, residuals, fullvar, tnorm = NULL, ncomp.selected = NULL)
{
   # Creates an object of ldecomp class. 
   #
   # Arguments:
   #   scores: matrix with score values (nobj x ncomp).
   #   loadings: matrix with loading values (nvar x ncomp).
   #   residuals: matrix with data residuals 
   #   fullvar: full variance of original data, preprocessed and centered
   #   tnorm: singular values for score normalization
   #   ncomp.selected: number of selected components
   #
   # Returns:
   #  object (list) of class ldecomp with following fields:   
   #   obj$scores: matrix with score values (nobj x ncomp).
   #   obj$residuals: matrix with residuals (nobj x nvar).
   #   obj$fullvar: full variance of original data
   #   obj$Q2: matrix with Q2 residuals (nobj x ncomp).
   #   obj$T2: matrix with T2 distances (nobj x ncomp)   
   #   obj$ncomp.selected: selected number of components
   #   obj$expvar: explained variance for each component
   #   obj$cumexpvar: cumulative explained variance
   
   scores = as.matrix(scores)
   loadings = as.matrix(loadings)
   residuals = as.matrix(residuals)
   
   # check dimension   
   if (ncol(scores) != ncol(loadings) || 
          nrow(scores) != nrow(residuals) || nrow(loadings) != ncol(residuals))
      stop('Dimensions of scores, loadings and data do not correspond to each other!')
   
   # set names for scores and loadings
   rownames(scores) = rownames(residuals)
   colnames(scores) = paste('Comp', 1:ncol(scores))
   rownames(loadings) = colnames(residuals)
   colnames(loadings) = paste('Comp', 1:ncol(loadings))
   
   if (is.null(ncomp.selected))
      ncomp.selected = ncol(scores)
   
   # calculate residual distances and explained variance
   obj = ldecomp.getDistances(scores, loadings, residuals, tnorm)
   var = ldecomp.getVariances(obj$Q2, fullvar)
   
   obj$expvar = var$expvar
   obj$cumexpvar = var$cumexpvar
   obj$scores = scores
   obj$residuals = residuals
   obj$fullvar = fullvar
   obj$ncomp.selected = ncomp.selected
   obj$call = match.call()
   
   class(obj) = "ldecomp"
   
   return (obj)
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
   #
   # Returns:
   #   res$Q2: matrix with Q2 residuals (nobj x ncomp).
   #   res$T2: matrix with T2 distances (nobj x ncomp)   

   ncomp = ncol(scores)
   nobj = nrow(scores)
   T2 = matrix(0, nrow = nobj, ncol = ncomp)
   Q2 = matrix(0, nrow = nobj, ncol = ncomp)

   # calculate normalized scores
   if (is.null(tnorm))
      tnorm = sqrt(colSums(scores ^ 2)/(nrow(scores) - 1));
   scoresn = sweep(scores, 2L, tnorm, '/', check.margin = F);  

   # calculate distances for each set of components
   for (i in 1:ncomp)
   {
      if (i < ncomp)
         res = residuals + scores[, (i + 1):ncomp, drop = F] %*% t(loadings[, (i + 1):ncomp, drop = F])
      else
         res = residuals
      
      Q2[, i] = rowSums(res^2)
      T2[, i] = rowSums(scoresn[, 1:i, drop = F]^2)
   }   
   
   # set dimnames and return results
   colnames(Q2) = colnames(T2) = colnames(scores)
   rownames(Q2) = rownames(T2) = rownames(scores)
   
   res = list(
      Q2 = Q2,
      T2 = T2,
      tnorm = tnorm
   )
}

ldecomp.getVariances = function(Q2, fullvar)
{   
   # Computes explained variance and cumulative explained variance 
   # for every component of a decomposition.
   #
   # Arguments:
   #   scores: matrix with scores.
   #   loadings: matrix with loadings.
   #   residuals: matrix with residuals.
   #   Q2: matrix with Q2 values
   # Returns:
   #   res$expvar: vector with explained variance for every component
   #   res$cumexpvar: vector with cumulative explained variance

   cumresvar = colSums(Q2) / fullvar * 100
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
   if (nobj == ncomp)
      T2lim = 0
   else
      T2lim = (ncomp * (nobj - 1) / (nobj - ncomp)) * qf(1 - alpha, ncomp, nobj - ncomp);  
   
   # calculate Q2 limit using F statistics
   conflim = 100 - alpha * 100;
   
   nvar = length(eigenvals)
   
   if (ncomp < nvar)
   {         
      eigenvals = eigenvals[(ncomp + 1):nvar]         
      
      cl = 2 * conflim - 100
      t1 = sum(eigenvals)
      t2 = sum(eigenvals^2)
      t3 = sum(eigenvals^3)
      h0 = 1 - 2 * t1 * t3/3/(t2^2);
      
      if (h0 < 0.001)
         h0 = 0.001
      
      ca = sqrt(2) * erfinv(cl/100)
      h1 = ca * sqrt(2 * t2 * h0^2)/t1
      h2 = t2 * h0 * (h0 - 1)/(t1^2)
      Q2lim = t1 * (1 + h1 + h2)^(1/h0)
   }
   else
      Q2lim = 0
   
   res = list(
      T2lim = T2lim,
      Q2lim = Q2lim
   )   
}   

plotCumVariance.ldecomp = function(obj, show.labels = F)
{
   # Shows cumulative explained variance plot.
   #
   # Arguments:
   #   obj: object of ldecomp class.
   #   show.labels: show or not labels for plot points
   
   data = cbind(1:length(obj$cumexpvar), obj$cumexpvar)
   data = rbind(c(0, 0), data)
   colnames(data) = c('Components', 'Explained variance, %')
   rownames(data) = round(c(0, obj$cumexpvar), 1)
   mdaplots.linescatter(data, main = 'Cumulative variance', show.labels = show.labels)
}

plotVariance.ldecomp = function(obj, show.labels = F)
{
   # Shows explained variance plot.
   #
   # Arguments:
   #   obj: object of ldecomp class.
   #   show.labels: show or not labels for plot points
   #
   
   data = cbind(1:length(obj$expvar), obj$expvar)
   colnames(data) = c('Components', 'Explained variance, %')
   rownames(data) = round(obj$expvar, 1)
   mdaplots.linescatter(data, main = 'Variance', show.labels = show.labels)
}

plotScores.ldecomp = function(obj, comp = c(1, 2), cgroup = NULL, 
                              show.labels = F, show.colorbar = T,
                              show.axes = T)
{
   # Shows scores plot.
   #
   # Arguments:
   #   obj: object of ldecomp class.
   #   comp: which components to show on x and y axis. 
   #   cgroup: variable for color grouping of plot points.
   #   show.labels: show or not labels for plot points.
   #   show.colorbar: show or not a colorbar legend if cgroup is provided.
   #   show.axes: show or not axes crossing (0, 0) point.
   
   if (length(comp) == 1)
   {   
      # scores vs objects
      data = cbind(1:nrow(obj$scores), obj$scores[, comp])
      colnames(data) = c('Objects', colnames(obj$scores)[comp])
      rownames(data) = rownames(obj$scores)
      mdaplots.scatter(data, main = 'Scores', cgroup = cgroup, 
                       show.labels = show.labels, 
                       show.colorbar = show.colorbar
                       )
   }
   else if (length(comp) == 2)
   {
      # scores vs scores
      data = obj$scores[, c(comp[1], comp[2])]   

      if (show.axes == T)
         show.lines = c(0, 0)      
      else
         show.lines = F
      
      mdaplots.scatter(data, main = 'Scores', cgroup = cgroup, 
                       show.labels = show.labels, 
                       show.colorbar = show.colorbar,
                       show.lines = show.lines)
   }
   else
   {
      stop('Wrong number of components!')
   }   
}  

plotResiduals.ldecomp = function(obj, ncomp = NULL, cgroup = NULL, 
                                 show.labels = F, show.colorbar = T)
{
   # Shows T2 vs Q2 residuals plot.
   #
   # Arguments:
   #   obj: object of ldecomp class.
   #   ncomp: number of components for the residuals 
   #   cgroup: variable for color grouping of plot points.
   #   show.labels: show or not labels for plot points.
   #   show.colorbar: show or not a colorbar legend if cgroup is provided.
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   
   data = cbind(obj$T2[, ncomp], obj$Q2[, ncomp])
   colnames(data) = c('T2', 'Q2')
   mdaplots.scatter(data, main = sprintf('Residuals (ncomp = %d)', ncomp), 
                    cgroup = cgroup, show.labels = show.labels,
                    show.colorbar = show.colorbar)
}  

print.ldecomp = function(obj, str = NULL)
{   
   if (is.null(str))
      str ='Results of data decomposition (class ldecomp)'
   
   cat('\n')
   cat(str, '\n')
   cat('\nMajor fields:\n')   
   cat('$scores - matrix with score values (nobj x ncomp)\n')
   cat('$T2 - matrix with T2 distances (nobj x ncomp)\n')
   cat('$Q2 - matrix with Q2 residuals (nobj x ncomp)\n')
   cat('$ncomp.selected - selected number of components\n')
   cat('$expvar - explained variance for each component\n')
   cat('$cumexpvar - cumulative explained variance\n\n')
}

summary.ldecomp = function(obj, str = NULL)
{
   if (is.null(str))
      str ='Summary for data decomposition (class ldecomp)'

   cat('\n')
   cat(str, '\n')
   cat(sprintf('\nSelected components: %d\n\n', obj$ncomp.selected))      
   
   data = cbind(round(obj$expvar, 2),
                round(obj$cumexpvar, 2))
   
   colnames(data) = c('Exp. var', 'Cum. exp. var')
   show(data)   
}

erfinv = function (x) qnorm((1 + x)/2)/sqrt(2)


