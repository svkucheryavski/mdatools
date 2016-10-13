#' Partial Least Squares Discriminant Analysis
#'
#' @description 
#' \code{plsda} is used to calibrate, validate and use of partial least squares discrimination 
#' analysis (PLS-DA) model.
#'
#' @param x
#' matrix with predictors.
#' @param c
#' vector with class values (should be a factor with either class number or class name for each object).
#' @param ncomp 
#' maximum number of components to calculate.
#' @param center 
#' logical, center or not predictors and response values.
#' @param scale 
#' logical, scale (standardize) or not predictors and response values.
#' @param cv
#' number of segments for cross-validation (if cv = 1, full cross-validation will be used).
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
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
plsda = function(x, c, ncomp = 15, center = T, scale = F, cv = NULL, exclcols = NULL, 
                 exclrows = NULL, x.test = NULL, c.test = NULL, method = 'simpls', alpha = 0.05, 
                 coeffs.ci = NULL, coeffs.alpha = 0.1, info = '', light = F, 
                 ncomp.selcrit = 'min') {
   
   # build a model and apply to calibration set
   model = plsda.cal(x, c, ncomp, center = center, scale = scale, method = method, light = light,
                     alpha = alpha, coeffs.ci = coeffs.ci, coeffs.alpha = coeffs.alpha, info = info,
                     ncomp.selcrit = ncomp.selcrit, cv = cv, exclcols = exclcols, exclrows = exclrows)
   
   # do test set validation if provided
   if (!is.null(x.test) && !is.null(c.test))
      model$testres = predict.plsda(model, x.test, c.test)

   # select optimal number of components and return object with model
   model = selectCompNum.pls(model, ncomp)
   model
}

#' Calibrate PLS-DA model
#' 
#' @param x
#' matrix with predictors.
#' @param c
#' vector with reference class values.
#' @param ncomp 
#' maximum number of components to calculate.
#' @param center 
#' logical, center or not predictors and response values.
#' @param scale 
#' logical, scale (standardize) or not predictors and response values.
#' @param cv
#' number of segments for cross-validation (if cv = 1, full cross-validation will be used).
#' @param method
#' method for calculating PLS model.
#' @param light
#' run normal or light (faster) version of PLS without calculationg some performance statistics.
#' @param alpha
#' significance level for calculating statistical limits for residuals.
#' @param coeffs.ci
#' method to calculate p-values and confidence intervals for regression coefficients (so far only 
#' jack-knifing is availavle: \code{='jk'}).
#' @param coeffs.alpha
#' significance level for calculating confidence intervals for regression coefficients.
#' @param info
#' short text with information about the model.
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param ncomp.selcrit
#' criterion for selecting optimal number of components (\code{'min'} for first local minimum of 
#' 
#' @export
plsda.cal = function(x, c, ncomp, center, scale, cv, method, light, alpha, coeffs.ci, coeffs.alpha,
                     info, exclcols = NULL, exclrows = NULL, ncomp.selcrit) {
   
   c = checkReferenceValues.classmodel(model, c, x)
   y = mda.df2mat(as.factor(c), full = TRUE)
   y[y == 0] = -1      
   
   if (length(exclcols) > 0)
      x = mda.exclcols(x, exclcols)
   
   if (length(exclrows) > 0) {
      x = mda.exclrows(x, exclrows)
      y = mda.exclrows(y, exclrows)
      c = mda.exclrows(c, exclrows)
   }
   
   # build a model and apply to calibration set
   model = pls.cal(x, y, ncomp, center = center, scale = scale, method = method, coeffs.ci = coeffs.ci,
                    coeffs.alpha = coeffs.alpha, info = info, light = light, alpha = alpha, cv = cv, 
                    ncomp.selcrit = ncomp.selcrit)
   model$classnames = unique(c)
   model$nclasses = length(model$classnames)
   model$calres = predict.plsda(model, x, c)
   
   # do cross-validation if needed
   if (!is.null(cv)) 
      model$cvres = plsda.crossval(model, x, c, center, scale, method)    
   
   # combine everything to model object
   class(model) = c("plsda", "classmodel", "pls")
   model
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
#' @param c.ref
#' a vector with reference class values (should be a factor)
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
predict.plsda = function(object, x, c.ref = NULL, ...) {   
   
   y.ref = NULL
   if (!is.null(c.ref)) {
      c.ref = checkReferenceValues.classmodel(object, c.ref, x)
      y.ref = mda.df2mat(as.factor(c.ref), full = TRUE)
      y.ref[y.ref == 0] = -1      
   }

   # do PLS predictions
   plsres = predict.pls(object, x, y.ref)

   # classify objects and set attributes   
   c.pred = classify.plsda(object, plsres$y.pred)   

   # combine everything to plsdares object
   cres = classres(c.pred, c.ref = c.ref, p.pred = plsres$y.pred, ncomp.selected = object$ncomp.selected)
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
#' 
#' @return
#' Classification results (an object of class \code{classres})
#'
#' @details
#' This is a service function for PLS-DA class, do not use it manually.
#'  
classify.plsda = function(model, y) {
   c.pred = array(-1, dim(y))
   c.pred[y >= 0] = 1
   
   c.pred = mda.setattr(c.pred, mda.getattr(y))
   attr(c.pred, 'name') = 'Class, predicted values'
   dimnames(c.pred) = dimnames(y)
   
   c.pred
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
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#' @param method
#' method for calculating PLS model.
#'
#' @return
#' object of class \code{plsdares} with results of cross-validation
#'  
plsda.crossval = function(model, x, c, center, scale, method) {
  
   c = checkReferenceValues.classmodel(model, c, x)
   attrs = mda.getattr(x)
   if (length(attrs$exclrows) > 0) {
      c = c[-attrs$exclrows, , drop = F]
      x = x[-attrs$exclrows, , drop = F]
      attr(c, 'exclrows') = NULL
      attr(x, 'exclrows') = NULL
   }
   
   y = mda.df2mat(as.factor(c), full = TRUE)
   y[y == 0] = -1      
   
   if (!is.null(model$coeffs$p.values))
      jack.knife = TRUE
   else
      jack.knife = FALSE
   plsres = pls.crossval(model, x, y, model$cv, center, scale, method, jack.knife)
   c.pred = classify.plsda(model, plsres$y.pred)
   cres = classres(c.pred, c.ref = c, p.pred = plsres$y.pred, ncomp.selected = model$ncomp.selected)

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
plot.plsda = function(x, ncomp = NULL, nc = 1, show.legend = T, show.labels = F, ...) {
   model = x
   
   if (!is.null(ncomp) && (ncomp <= 0 || ncomp > model$ncomp)) 
      stop('Wrong value for number of components!')
   
   par(mfrow = c(2, 2))      
   plotXResiduals(model, ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   
   plotRegcoeffs(model, ncomp = ncomp, ny = nc, show.labels = show.labels)   
   plotMisclassified(model, nc = nc, show.legend = show.legend)   
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
summary.plsda = function(object, ncomp = NULL, nc = NULL, ...) {
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
      
   for (n in nc) {   
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
print.plsda = function(x, ...) {
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
