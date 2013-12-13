## class and methods for PCA results ##

pcares = function(scores, loadings, residuals, totvar, tnorm = NULL, ncomp.selected = NULL, ...)
{
   # Creates an object of pcares class. In fact the class is a wrapper for ldecomp and
   # uses its methods and attributes.
   
   res = ldecomp(scores, loadings, residuals, totvar, tnorm = tnorm, ncomp.selected = ncomp.selected, ...)
   class(res) = c('pcares', 'ldecomp')   
   
   res
}   


print.pcares = function(obj, ...)
{   
   print.ldecomp(obj, 'Results for PCA decomposition (class pcares)', ...)
   cat('\n')
}

summary.pcares = function(obj, ...)
{
   summary.ldecomp(obj, 'Summary for PCA results', ...)
}