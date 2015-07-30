#' Partial Least Squares Discriminant Analysis
#'
#' @description 
#' \code{plsda} is used to calibrate, validate and use of partial least squares discrimination 
#' analysis (PLS-DA) model.
#'
#' @param x
#' matrix with predictors.
#' @param c
#' vector with class values (either class number or class name as text for each object).
#' @param ncomp 
#' maximum number of components to calculate.
#' @param center 
#' logical, center or not predictors and response values.
#' @param scale 
#' logical, scale (standardize) or not predictors and response values.
#' @param cv
#' number of segments for cross-validation (if cv = 1, full cross-validation will be used).
#' @param x.test
#' matrix with predictors for test set.
#' @param c.test
#' vector with reference class values for test set (same format as calibration values).
#' @param method
#' method for calculating PLS model.
#' @param alpha
#' significance level for calculating statistical limits for residuals.
#' @param coeffs.ci
#' method to calculate p-values and confidence intervals for regression coefficients (so far only 
#' jack-knifing is availavle: \code{='jk'}).
#' @param coeffs.alpha
#' significance level for calculating confidence intervals for regression coefficients.
#' @param info
#' short text with information about the model.
#' @param light
#' run normal or light (faster) version of PLS without calculationg some performance statistics.
#' @param ncomp.selcrit
#' criterion for selecting optimal number of components (\code{'min'} for first local minimum of 
#' RMSECV and \code{'wold'} for Wold's rule.)
#' 
#' @return
#' Returns an object of \code{plsda} class with following fields (most inherited from class 
#' \code{pls}):
#' \item{ncomp }{number of components included to the model.} 
#' \item{ncomp.selected }{selected (optimal) number of components.} 
#' \item{xloadings }{matrix with loading values for x decomposition.} 
#' \item{yloadings }{matrix with loading values for y (c)  decomposition.} 
#' \item{weights }{matrix with PLS weights.} 
#' \item{coeffs }{matrix with regression coefficients calculated for each component.}   
#' \item{info }{information about the model, provided by user when build the model.} 
#' \item{calres }{an object of class \code{\link{plsdares}} with PLS-DA results for a calibration 
#' data.} 
#' \item{testres }{an object of class \code{\link{plsdares}} with PLS-DA results for a test data, 
#' if it was provided.} 
#' \item{cvres }{an object of class \code{\link{plsdares}} with PLS-DA results for cross-validation,
#' if this option was chosen.} 
#'
#' @details 
#' The \code{plsda} class is based on \code{pls} with extra functions and plots covering 
#' classification functionality. All plots for \code{pls} can be used. E.g. of you want to see the 
#' real predicted values (y in PLS) instead of classes use \code{plotPredictions.pls(model)} instead
#' of \code{plotPredictions(model)}.
#' 
#' Calculation of confidence intervals and p-values for regression coefficients are available
#' only by jack-knifing so far. See help for \code{\link{regcoeffs}} objects for details.
#'
#' @seealso 
#' Specific methods for \code{plsda} class:
#' \tabular{ll}{
#'  \code{print.plsda} \tab prints information about a \code{pls} object.\cr
#'  \code{summary.plsda} \tab shows performance statistics for the model.\cr
#'  \code{plot.plsda} \tab shows plot overview of the model.\cr
#'  \code{\link{predict.plsda}} \tab applies PLS-DA model to a new data.\cr
#' }
#' 
#' Methods, inherited from \code{classmodel} class:
#' \tabular{ll}{
#'  \code{\link{plotPredictions.classmodel}} \tab shows plot with predicted values.\cr
#'  \code{\link{plotSensitivity.classmodel}} \tab shows sensitivity plot.\cr
#'  \code{\link{plotSpecificity.classmodel}} \tab shows specificity plot.\cr
#'  \code{\link{plotMisclassified.classmodel}} \tab shows misclassified ratio plot.\cr
#' }
#' 
#' See also methods for class \code{\link{pls}}.
#' 
#' @author 
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#' 
#' @examples
#' ### Examples for PLS-DA model class
#' 
#' library(mdatools)
#' 
#' ## 1. Make a PLS-DA model with full cross-validation and show model overview
#' 
#' # make a calibration set from iris data (3 classes)
#' # use names of classes as class vector
#' x.cal = iris[seq(1, nrow(iris), 2), 1:4] 
#' c.cal = iris[seq(1, nrow(iris), 2), 5]
#' 
#' model = plsda(x.cal, c.cal, ncomp = 3, cv = 1, info = 'IRIS data example')
#' model = selectCompNum(model, 1)
#' 
#' # show summary and basic model plots
#' # misclassification will be shown only for first class
#' summary(model)
#' plot(model)
#' 
#' # summary and model plots for second class
#' summary(model, nc = 2)
#' plot(model, nc = 2)
#' 
#' # summary and model plot for specific class and number of components
#' summary(model, nc = 3, ncomp = 3)
#' plot(model, nc = 3, ncomp = 3)
#' 
#' ## 2. Show performance plots for a model
#' par(mfrow = c(2, 2))
#' plotSpecificity(model)
#' plotSensitivity(model)
#' plotMisclassified(model)
#' plotMisclassified(model, nc = 2)
#' par(mfrow = c(1, 1))
#' 
#' ## 3. Show both class and y values predictions
#' par(mfrow = c(2, 2))
#' plotPredictions(model)
#' plotPredictions(model, res = 'calres', ncomp = 2, nc = 2)
#' plotPredictions(structure(model, class = "pls"))
#' plotPredictions(structure(model, class = "pls"), ncomp = 2, ny = 2)
#' par(mfrow = c(1, 1))
#' 
#' ## 4. All plots from ordinary PLS can be used, e.g.:
#' par(mfrow = c(2, 2))
#' plotXYScores(model)
#' plotYVariance(model)
#' plotXResiduals(model)
#' plotRegcoeffs(model, ny = 2)
#' par(mfrow = c(1, 1))
#' 
#' @export
plsda = function(x, c, ncomp = 15, center = T, scale = F, cv = NULL, 
               x.test = NULL, c.test = NULL, method = 'simpls', alpha = 0.05, 
               coeffs.ci = NULL, coeffs.alpha = 0.1, info = '', light = F,
               ncomp.selcrit = 'min')
{
   x = as.matrix(x)
   c = as.matrix(c)
   y = plsda.getReferenceValues(c)

   # set LOO cross-validation if jack.knife is selected
   jack.knife = F
   if (!is.null(coeffs.ci) && coeffs.ci == 'jk')
   {
      jack.knife = T
      if (is.null(cv))
      {   
         cv = 1      
      }   
      else
      {
         if (is.list(cv))
         {   
            if (length(cv) == 1 && cv[[1]] == 'loo')
               cvnseg = nrow(y)
            else   
               cvnseg = cv[[2]]
         }   
         else
         {
            cvnseg = cv
         }   
         
         if (cvnseg > 1 && cvnseg < 10)
            warning('Number of segments in cross-validation is too small for jack-knifing!')   
      }   
   }   
   
   # correct maximum number of components
   if (!is.null(cv))
   {
      if (!is.numeric(cv))
         nseg = cv[[2]]
      else
         nseg = cv
      
      if (nseg == 1)
         nobj.cv = 1
      else
         nobj.cv = ceiling(nrow(x)/nseg)  
   }
   else
   {   
      nobj.cv = 0
   }
   
   ncomp = min(ncol(x), nrow(x) - 1 - nobj.cv, ncomp)
   
   # build a model and apply to calibration set
   model = pls.cal(x, y, ncomp, center = center, scale = scale, method = method, light = light)
   model$light = light
   model$ncomp.selcrit = ncomp.selcrit
   model$alpha = alpha   
   model$classnames = unique(c)
   model$nclasses = length(model$classnames)
   model$calres = predict.plsda(model, x, c)
   
   # do cross-validation if needed
   if (!is.null(cv)) 
   {
      res = plsda.crossval(model, x, c, cv, center = center, scale = scale, jack.knife = jack.knife)    
      if (jack.knife == TRUE)
      {   
         model$coeffs = regcoeffs(model$coeffs$values, res$jkcoeffs, coeffs.alpha)
         res[['jkcoeffs']] = NULL
         model$cvres = res
      }
      else
      {
         model$cvres = res;
      }   
      
   }
   # do test set validation if provided
   if (!is.null(x.test) && !is.null(c.test))
   {
      x.test = as.matrix(x.test)
      model$testres = predict.plsda(model, x.test, c.test)
   }
   
   model = selectCompNum.pls(model, ncomp)
   model$call = match.call()
   model$info = info
   
   class(model) = c("plsda", "classmodel", "pls")
   
   model
}

#' Reference values for PLS-DA
#' 
#' @description
#' Generates matrix with reference y values (-1 and +1) for a 
#' vector with class values
#' 
#' @param c 
#' vector with class values (discrete)
#' @param classnames
#' vector with names for the classes
#' 
#' @return
#' the generated matrix with one column for each class
#' 
#' @export
plsda.getReferenceValues = function(c, classnames = NULL)
{
   # generate matrix with y values
   
   if (is.null(classnames))
      classnames = unique(c)

   nclasses = length(classnames)
   y = matrix(-1, nrow = length(c), ncol = nclasses)
   
   for (i in 1:nclasses)
      y[c == classnames[i], i] = 1
   
   rownames(y) = rownames(c)
   colnames(y) = classnames

   y
}

#' PLS-DA predictions
#' 
#' @description
#' Applies PLS-DA model to a new data set
#' 
#' @param object
#' a PLS-DA model (object of class \code{plsda})
#' @param x
#' a matrix with x values (predictors)
#' @param c
#' a vector with reference class values
#' @param cv
#' logical, are predictions for cross-validation or not
#' @param ...
#' other arguments
#' 
#' @return
#' PLS-DA results (an object of class \code{plsdares})
#'
#' @details
#' See examples in help for \code{\link{plsda}} function.
#'  
#' @export 
predict.plsda = function(object, x, c = NULL, cv = F, ...)
{   
   y = plsda.getReferenceValues(c, object$classnames)
   plsres = predict.pls(object, x, y)
   cres = classify.plsda(object, plsres$y.pred, c)

   res = plsdares(plsres, cres)

   res
}  

#' PLS-DA classification
#' 
#' @description
#' Converts PLS predictions of y values to predictions of classes
#' 
#' @param model
#' a PLS-DA model (object of class \code{plsda})
#' @param y
#' a matrix with predicted y values
#' @param c.ref
#' a vector with reference class values
#' 
#' @return
#' Classification results (an object of class \code{classres})
#'
#' @details
#' This is a service function for PLS-DA class, do not use it manually.
#'  
classify.plsda = function(model, y, c.ref = NULL)
{
   c.pred = array(-1, dim(y))
   c.pred[y >= 0] = 1
   dimnames(c.pred) = list(rownames(y), paste('Comp', 1:model$ncomp), model$classnames)
   cres = classres(c.pred, c.ref = c.ref, p.pred = y, ncomp.selected = model$ncomp.selected)

   cres
}

#' Cross-validation of a PLS-DA model
#' 
#' @description
#' Does the cross-validation of a PLS-DA model
#' 
#' @param model
#' a PLS-DA model (object of class \code{plsda})
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param c
#' a vetor with c values (classes from calibration set)
#' @param cv
#' number of segments (if cv = 1, full cross-validation will be used)
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#' @param jack.knife
#' logical, do jack-knifing or not
#'
#' @return
#' object of class \code{plsdares} with results of cross-validation
#'  
plsda.crossval = function(model, x, c, cv, center = T, scale = F, jack.knife = T)
{
   y = plsda.getReferenceValues(c, model$classnames)
   plsres = pls.crossval(model, x, y, cv, center, scale, jack.knife)
   cres = classify.plsda(model, plsres$y.pred, c)

   res = plsdares(plsres, cres)

   res
}

#' Model overview plot for PLS-DA
#' 
#' @description
#' Shows a set of plots (x residuals, regression coefficients, misclassification ratio and 
#' predictions) for PLS-DA model.
#' 
#' @param x
#' a PLS-DA model (object of class \code{plsda})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param nc
#' which class to show the plots
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{plsda}} function.
#' 
#' @export
plot.plsda = function(x, ncomp = NULL, nc = 1, show.legend = T, show.labels = F, ...)
{
   model = x
   
   if (!is.null(ncomp) && (ncomp <= 0 || ncomp > model$ncomp)) 
      stop('Wrong value for number of components!')
   
   par(mfrow = c(2, 2))      
   plotXResiduals(model, ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   
   plotRegcoeffs(model, ncomp = ncomp, ny = nc, show.labels = show.labels)   
   plotMisclassified(model, nc = nc, show.legend = show.legend)   
   
   if (!is.null(model$cvres))
      plotPredictions(model, res = 'cvres', ncomp = ncomp, show.labels = show.labels, 
                      show.legend = show.legend)   
   else
      plotPredictions(model, ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   

   par(mfrow = c(1, 1))
}

#' Summary method for PLS-DA model object
#' 
#' @description
#' Shows some statistics for the model.
#' 
#' @param object
#' a PLS-DA model (object of class \code{plsda})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param nc
#' which class to show the summary for (if NULL, will be shown for all)
#' @param ...
#' other arguments
#' 
#' @export
summary.plsda = function(object, ncomp = NULL, nc = NULL, ...)
{
   obj = object
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   else if (ncomp <= 0 || ncomp > obj$ncomp)
      stop('Wrong value for number of components!')
   
   if (is.null(nc))
      nc = 1:obj$nclasses
   
   cat('\nPLS-DA model (class plsda) summary statistics\n\n')
   cat(sprintf('Number of selected components: %d\n', ncomp))
   
   if (!is.null(obj$info))
      cat(sprintf('Info: %s\n', obj$info))
      
   for (n in nc)
   {   
      cat(sprintf('\nClass #%d (%s)\n', n, obj$classnames[n]))
      
      data = as.matrix(obj$calres, ncomp = ncomp, nc = n)
      rownames(data) = 'Cal'
      
      if (!is.null(obj$cvres))
      {
         data = rbind(data, as.matrix(obj$cvres, ncomp = ncomp, nc = n))      
         rownames(data)[2] = 'CV'
      }
      
      if (!is.null(obj$testres))
      {
         data = rbind(data, as.matrix(obj$testres, ncomp = ncomp, nc = n))
         rownames(data)[nrow(data)] = 'Test'
      }   
      
      data[, 1:4] = round(data[, 1:4], 2)      
      print(data)
   }   
   cat('\n')
}

#' Print method for PLS-DA model object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a PLS-DA model (object of class \code{plsda})
#' @param ...
#' other arguments
#' 
#' @export
print.plsda = function(x, ...)
{
   obj = x
   
   cat('\nPLS-DA model (class plsda)\n')
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$coeffs - vector with regression coefficients\n')
   cat('$xloadings - vector with x loadings\n')
   cat('$yloadings - vector with Y loadings\n')
   cat('$weights - vector with weights\n')
   cat('$calres - results for calibration set\n')
   if (!is.null(obj$cvres))
   {
      cat('$cvres - results for cross-validation\n')      
   }   
   if (!is.null(obj$testres))
   {
      cat('$testres - results for test set\n')      
   }   
   cat('\nTry summary(model) and plot(model) to see the model performance.\n')   
}
