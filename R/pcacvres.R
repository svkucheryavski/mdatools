# class and methods for PCA cross-validates results #
pcacvres = function(T2, Q2, fullvar, ...) UseMethod("pcacvres")

pcacvres.default = function(T2, Q2, fullvar, ncomp.selected = NULL, ...)
{
   # Creates an object of pcacvres class. The returned object also inherits class
   # ldecomp and some of its methods.
   #
   # Arguments:
   #   T2: matrix with T2 distances (nobj x ncomp) for cross-validation.
   #   Q2: matrix with Q2 residuals (nobj x ncomp) for cross-validation.
   #   fullvar: full variance of data 
   #   ncomp.selected: number of selected components
   #
   # Returns:
   #  list with cross-validation results (object of class pcacvres)   
   #   obj$Q2: matrix with Q2 residuals (nobj x ncomp).
   #   obj$T2: matrix with T2 distances (nobj x ncomp)   
   #   obj$ncomp.selected: selected number of components
   #   obj$expvar: explained variance for each component
   #   obj$cumexpvar: cumulative explained variance
   
   if (is.null(ncomp.selected))
      ncomp.selected = ncol(T2)
   
   obj = list(
      T2 = T2,
      Q2 = Q2
      )
   
   var = ldecomp.getVariances(obj$Q2, fullvar)
   
   obj$expvar = var$expvar
   obj$cumexpvar = var$cumexpvar
   obj$ncomp.selected = ncomp.selected
   
   obj$call = match.call()   
   class(obj) = c('pcacvres', 'ldecomp')   
   
   return (obj)
}   

plotScores.pcacvres = function(obj, ...)
{
   # stub functon to show that scores plot is not available for cv results
   
   stop('Scores plot is not available for cross-validated results.')
}   

print.pcacvres = function(obj)
{   
   cat('\nCross-validation results for PCA model (class pcacvres) \n')

   cat('\nMajor fields:\n')   
   cat('$T2 - matrix with T2 distances (nobj x ncomp)\n')
   cat('$Q2 - matrix with Q2 residuals (nobj x ncomp)\n')
   cat('$ncomp.selected - selected number of components\n')
   cat('$expvar - explained variance for each component\n')
   cat('$cumexpvar - cumulative explained variance\n\n')
}

summary.pcacvres = function(obj)
{
   cat('\nSummary for cross-validation of PCA model\n')
   cat(sprintf('Selected components: %d\n', obj$ncomp.selected))      

   data = cbind(round(obj$expvar, 2),
                round(obj$cumexpvar, 2))
   
   colnames(data) = c('Exp. var', 'Cum. exp. var')
   show(data)   
}   