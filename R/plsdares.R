plsdares = function(plsres, cres)
{
   # Creates an object of plsda class. 
   #
   # Arguments:
   #  plsres: an object of plsres class (results for PLS regression)
   #  cres: an object of classres class (results for classification)
   
   obj = c(plsres, cres)
   class(obj) = c('plsdares', 'classres', 'plsres')   
   
   obj$call = match.call()   
   
   obj
}   

as.matrix.plsdares = function(obj, ncomp = NULL, nc = NULL)
{
   plsmat = as.matrix.plsres(obj, ncomp = ncomp, ny = nc)
   classmat = as.matrix.classres(obj, ncomp = ncomp, nc = nc)
   mat = cbind(plsmat[, 1:4, drop = F], classmat)

   mat
}

plot.plsdares = function(obj, nc = NULL, ncomp = NULL, ...)
{
   par(mfrow = c(2, 2))
   plotXResiduals.plsres(obj, ncomp = ncomp, ...)
   plotYVariance.plsres(obj, ncomp = ncomp, ...)
   plotPerformance(obj, nc = nc, ncomp = ncomp, ...)
   plotPredictions(obj, nc = nc, ncomp = ncomp, ...)
   par(mfrow = c(1, 1))
}

summary.plsdares = function(obj, nc = NULL)
{
   if (is.null(nc))
      nc = 1:obj$nclasses
   cat('\nPLSDA results (class plsdares) summary:\n');
   cat(sprintf('Number of selected components: %d\n', obj$ncomp.selected))
   for (n in nc)
   {
      cat(sprintf('\nClass #%d (%s)\n', n, obj$classnames[n]))
      
      mat = as.matrix(obj, nc = n)
      mat[, 1:4] = round(mat[, 1:4], 2)
      print(mat)
      cat('\n')
   }
}

print.plsdares = function(obj)
{

   cat('\nPLSDA results (class plsdares)\n')
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$ncomp.selected - number of selected components\n')
   cat('$c.pred - array with predicted class values\n')
   if (!is.null(obj$c.ref))
   {   
      cat('$c.ref - vector with reference class values\n')
      cat('$tp - number of true positives\n')
      cat('$fp - number of false positives\n')
      cat('$fn - number of false negatives\n')
      cat('$specificity - specificity of predictions\n')
      cat('$sensitivity - sensitivity of predictions\n')
      cat('$misclassified - misclassification ratio for predictions\n')
      cat('$ydecomp - decomposition of y values (ldecomp object)\n')
   }
   cat('$xdecomp - decomposition of x values (ldecomp object)\n')
}
