#  methods for classification results #
classres = function(res, ...) UseMethod("classres")

classres.getPredictionPerformance = function(res)
{
   fn = NULL
   fp = NULL
   tp = NULL
   sensitivity = NULL
   specificity = NULL
   f1 = NULL
   
   if (!is.null(res.creference))
   {
      fn = sum((res.creference == 1) & (res.cpredictions == 0))
      fp = sum((res.creference == 0) & (res.cpredictions == 1))
      tp = sum((res.creference == 1) & (res.cpredictions == 1))
      
      sensitivity = tp / (tp + fn)
      specificity = tp / (tp + fp)
      f1 = 2 * sensitivity * specificity / (sensitivity + specificity)
   }
   
   res$fn = fn
   res$fp = fp
   res$tp = tp
   res$sensitivity = sensitivity
   res$specificity = specificity
   res$f1 = f1
   
   return (res)
}   