#' SIMCA one-class classification
#' 
#' @description
#' \code{simca} is used to make SIMCA (Soft Independent Modelling of Class Analogies) model for 
#' one-class classification.
#' 
#' @param x
#' a numerical matrix with data values.
#' @param classname
#' short text (up to 20 symbols) with class name.
#' @param ncomp
#' maximum number of components to calculate.
#' @param center
#' logical, do mean centering of data or not.
#' @param scale
#' logical, do sdandardization of data or not.
#' @param cv
#' number of segments for random cross-validation (1 for full cross-validation).
#' @param x.test
#' a numerical matrix with test data.
#' @param c.test
#' a vector with text values (names of classes) of test data objects.
#' @param alpha
#' significance level for calculating limit for T2 and Q residuals.
#' @param method
#' method to compute principal components.
#' @param info
#' text with information about the model.
#' 
#' @details 
#' SIMCA is in fact PCA model with additional functionality, so \code{simca} class inherits most 
#' of the functionality of \code{\link{pca}} class. 
#'
#' @return 
#' Returns an object of \code{simca} class with following fields:
#' \item{classname }{a short text with class name.} 
#' \item{modpower }{a matrix with modelling power of variables.} 
#' \item{calres }{an object of class \code{\link{simcares}} with classification results for a 
#' calibration data.} 
#' \item{testres }{an object of class \code{\link{simcares}} with classification results for a test 
#' data, if it was provided.} 
#' \item{cvres }{an object of class \code{\link{simcares}} with classification results for 
#' cross-validation, if this option was chosen.} 
#' 
#' Fields, inherited from \code{\link{pca}} class:
#' \item{ncomp }{number of components included to the model.} 
#' \item{ncomp.selected }{selected (optimal) number of components.} 
#' \item{loadings }{matrix with loading values (nvar x ncomp).} 
#' \item{eigenvals }{vector with eigenvalues for all existent components.} 
#' \item{expvar }{vector with explained variance for each component (in percent).} 
#' \item{cumexpvar }{vector with cumulative explained variance for each component (in percent).} 
#' \item{T2lim }{statistical limit for T2 distance.} 
#' \item{Qlim }{statistical limit for Q residuals.} 
#' \item{info }{information about the model, provided by user when build the model.} 
#'
#' @references 
#' S. Wold, M. Sjostrom. "SIMCA: A method for analyzing chemical data in terms of similarity and 
#' analogy" in B.R. Kowalski (ed.), Chemometrics Theory and Application, American Chemical Society 
#' Symposium Series 52, Wash., D.C., American Chemical Society, p. 243-282.
#' 
#' @seealso 
#' Methods for \code{simca} objects:
#' \tabular{ll}{
#'  \code{print.simca} \tab shows information about the object.\cr
#'  \code{summary.simca} \tab shows summary statistics for the model.\cr
#'  \code{plot.simca} \tab makes an overview of SIMCA model with four plots.\cr
#'  \code{\link{predict.simca}} \tab applies SIMCA model to a new data.\cr
#'  \code{\link{plotModellingPower.simca}} \tab shows plot with modelling power of variables.\cr
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
#' Methods, inherited from \code{\link{pca}} class:
#' \tabular{ll}{
#'  \code{\link{selectCompNum.pca}} \tab set number of optimal components in the model\cr
#'  \code{\link{plotScores.pca}} \tab shows scores plot.\cr
#'  \code{\link{plotLoadings.pca}} \tab shows loadings plot.\cr
#'  \code{\link{plotVariance.pca}} \tab shows explained variance plot.\cr
#'  \code{\link{plotCumVariance.pca}} \tab shows cumulative explained variance plot.\cr
#'  \code{\link{plotResiduals.pca}} \tab shows Q vs. T2 residuals plot.\cr
#' }
#' 
#' @examples
#' ## make a SIMCA model for Iris setosa class with full cross-validation
#' library(mdatools)
#' 
#' data = iris[, 1:4]
#' class = iris[, 5]
#' 
#' # take first 20 objects of setosa as calibration set 
#' se = data[1:20, ]
#' 
#' # make SIMCA model and apply to test set
#' model = simca(se, 'setosa', cv = 1)
#' model = selectCompNum(model, 1)
#' 
#' # show infromation, summary and plot overview
#' print(model)
#' summary(model)
#' plot(model)
#' 
#' # show predictions 
#' par(mfrow = c(2, 1))
#' plotPredictions(model, show.labels = TRUE)
#' plotPredictions(model, res = 'calres', ncomp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' # show performance, modelling power and residuals for ncomp = 2
#' par(mfrow = c(2, 2))
#' plotSensitivity(model)
#' plotMisclassified(model)
#' plotModellingPower(model, ncomp = 2, show.labels = TRUE)
#' plotResiduals(model, ncomp = 2)
#' par(mfrow = c(1, 1))
#'
#' @export
simca = function(x, classname, ncomp = 15, center = T, scale = F, cv = NULL, x.test = NULL, 
                 c.test = NULL, alpha = 0.05, method = 'svd', info = '')
{
   x = as.matrix(x)
   
   # check if data has missing values
   if (sum(is.na(x)) > 0)
   {
      warning('Data has missing values, will try to fix using pca.mvreplace.')
      x = pca.mvreplace(x, center = center, scale = scale)
   }   
   
   if (!is.character(classname))
      stop('Argument "classname" must be a text!')
   
   if (length(classname) > 20)
      stop('Argument "classname" must have up to 20 symbols!')
   
   # correct maximum number of components
   ncomp = min(ncomp, ncol(x), nrow(x) - 1)
   
   # calibrate model  
   model = pca.cal(x, ncomp, center = center, scale = scale, method = method)
   model$ncomp = ncomp
   model$ncomp.selected = model$ncomp
   model$nclasses = 1
   model$classname = classname
   model$info = info
   model$alpha = alpha
   
   # calculate and assign limit values for T2 and Q residuals
   lim = ldecomp.getResLimits(model$eigenvals, nrow(x), model$ncomp.selected, model$alpha)
   model$T2lim = lim$T2lim
   model$Qlim = lim$Qlim   

   model$call = match.call()   
   class(model) = c("simca", "classmodel", "pca")
   
   # apply model to calibration set
   model$calres = predict.simca(model, x, c.ref = rep(classname, nrow(x)))
   model$modpower = model$calres$modpower
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = simca.crossval(model, x, cv, center = center, scale = scale)
   
   # apply model to test set if provided
   if (!is.null(x.test))
   {
      if (is.null(c.test))
         c.test = rep(classname, nrow(x.test))
      model$testres = predict.simca(model, x.test, c.ref = c.test)
   }
   
   model
}

#' SIMCA predictions
#' 
#' @description
#' Applies SIMCA model to a new data set
#' 
#' @param object
#' a SIMCA model (object of class \code{simca})
#' @param x
#' a matrix with x values (predictors)
#' @param c.ref
#' a vector with reference class names (same as class names for models)
#' @param cv
#' logical, are predictions for cross-validation or not
#' @param ...
#' other arguments
#' 
#' @return
#' SIMCA results (an object of class \code{simcares})
#'
#' @details
#' See examples in help for \code{\link{simca}} function.
#' 
#' @export
predict.simca = function(object, x, c.ref = NULL, cv = F, ...)
{
   x = as.matrix(x)
   
   if (is.null(rownames(x)))
      rownames(x) = 1:nrow(x)
   
   if (is.null(colnames(x)))
      colnames(x) = paste('v', 1:ncol(x), sep = '')
   
   pres = predict.pca(object, x, cv)     
   pres$Qlim = object$Qlim
   pres$T2lim = object$T2lim
   
   c.pred = simca.classify(object, pres)
   
   # check c.ref values and add dimnames
   if (!is.null(c.ref))
   {   
      c.ref = as.matrix(c.ref)
      rownames(c.ref) = rownames(x)
      colnames(c.ref) = object$classname
   } 
   
   cres = classres(c.pred, c.ref = c.ref, ncomp.selected = object$ncomp.selected)
   res = simcares(pres, cres)
   
   res
}

#' SIMCA classification
#' 
#' @description
#' Make classification based on calculated T2 and Q values and corresponding limits
#' 
#' @param model
#' a SIMCA model (object of class \code{simca})
#' @param res
#' results of projection data to PCA space
#' 
#' @return
#' vector with predicted class values (\code{c.pred})
#'
#' @details
#' This is a service function for SIMCA class, do not use it manually.
#'  
simca.classify = function(model, res)
{
   ncomp = model$ncomp
   c.pred = array(0, dim = c(nrow(res$Q), ncomp, 1))
   dimnames(c.pred) = list(rownames(res$Q), paste('Comp', 1:ncomp), model$classname)
   
   for (i in 1:ncomp)
   {
      c.pred[, i, 1] = 
         res$T2[, i] <= model$T2lim[1, i] & res$Q[, i] <= model$Qlim[1, i]
   }   
   c.pred = c.pred * 2 - 1
  
   c.pred
}  

#' Cross-validation of a SIMCA model
#' 
#' @description
#' Does the cross-validation of a SIMCA model
#' 
#' @param model
#' a SIMCA model (object of class \code{simca})
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param cv
#' number of segments (if cv = 1, full cross-validation will be used)
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#'
#' @return
#' object of class \code{simcares} with results of cross-validation
#'  
simca.crossval = function(model, x, cv, center = T, scale = F)
{
   ncomp = model$ncomp   
   nobj = nrow(x)
   
   # get matrix with indices for cv segments
   idx = crossval(nobj, cv)
   nseg = nrow(idx);
   nrep = dim(idx)[3]
   
   Q = matrix(0, ncol = ncomp, nrow = nobj)   
   T2 = matrix(0, ncol = ncomp, nrow = nobj)   
   Qlim = matrix(0, ncol = ncomp, nrow = 1)   
   T2lim = matrix(0, ncol = ncomp, nrow = 1)   
   
   c.pred = array(0, dim = c(nobj, ncomp, 1))
   c.ref = matrix(model$classname, ncol = 1, nrow = nobj)
   
   # loop over segments
   for (iRep in 1:nrep)
   {
      
      for (iSeg in 1:nseg)
      {
         ind = na.exclude(idx[iSeg, ,iRep])
      
         if (length(ind) > 0)
         {   
            x.cal = x[-ind, , drop = F]
            x.val = x[ind, , drop = F]
         
            m = pca.cal(x.cal, ncomp, center, scale)               
            res = predict.pca(m, x.val, cv = T)
            Q[ind, ] = Q[ind, ] + res$Q
            T2[ind, ] = T2[ind, ] + res$T2
                        
            lim = ldecomp.getResLimits(m$eigenvals, nrow(x.cal), ncomp, model$alpha)
            T2lim = T2lim + lim$T2lim
            Qlim = Qlim + lim$Qlim         
         }
      }  
   }
   
   Q = Q / nrep;
   T2 = T2 / nrep;
   Qlim = Qlim / nrep;
   T2lim = T2lim / nrep;
   m = list(Qlim = Qlim, T2lim = T2lim, classname = model$classname, ncomp = model$ncomp)
   r = list(Q = Q, T2 = T2, classname = model$classname)
   
   c.pred = simca.classify(m, r)
   
   
   dimnames(c.pred) = list(rownames(x), colnames(model$loadings), model$classname)
   rownames(Q) = rownames(T2) = rownames(c.pred) = rownames(c.ref) = rownames(x)
   colnames(Q) = colnames(T2) = colnames(c.pred) = colnames(model$loadings)
   pres = pcares(NULL, NULL, NULL, model$calres$totvar, model$tnorm, model$ncomp.selected, T2, Q)
   cres = classres(c.pred, c.ref = c.ref)   
   res = simcares(pres, cres)
   
   res
}  

#' Modelling power plot for SIMCA model
#' 
#' @description
#' Shows a plot with modelling power values for each predictor
#' 
#' @param obj
#' a SIMCA model (object of class \code{simca})
#' @param ncomp
#' number of components to show the values for
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @export
plotModellingPower.simca = function(obj, ncomp = NULL, type = 'h', main = NULL, 
                                    xlab = 'Variables', ylab = '', ...)
{
   if (is.null(main))
   {
      if (is.null(ncomp))
         main = 'Modelling power'
      else
         main = sprintf('Modelling power (ncomp = %d)', ncomp)      
   }   

   ncomp = getSelectedComponents(obj, ncomp)
   
   nvar = nrow(obj$modpower)
   if (is.null(type))
   {   
      if (nvar < 20)
         type = 'h'
      else
         type = 'l'
   }
   
   data = cbind(1:nvar, obj$modpower[, ncomp, drop = F])
   mdaplot(data, type = type, xlab  = xlab, ylab = ylab, main = main, ...)
}   

#' Model overview plot for SIMCA
#' 
#' @description
#' Shows a set of plots for SIMCA model.
#' 
#' @param x
#' a SIMCA model (object of class \code{simca})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
#' @export
plot.simca = function(x, ncomp = NULL, ...)
{
   obj = x
   
   par(mfrow = c(2, 2))
   plotScores(obj, ...)
   plotModellingPower(obj, ncomp = ncomp, main = 'Modelling power', 
                      show.labels = ncol(obj$modpower) < 10, ...)
   plotResiduals(obj, main = 'Residuals', ncomp = ncomp, ...)
   plotCumVariance(obj, ...)
   par(mfrow = c(1, 1))
}  


#' Summary method for SIMCA model object
#' 
#' @description
#' Shows performance statistics for the model.
#' 
#' @param object
#' a SIMCA model (object of class \code{simca})
#' @param ...
#' other arguments
#' 
#' @export
summary.simca = function(object, ...)
{
   obj = object
   
   cat(sprintf('\nSIMCA model for class "%s" summary\n\n', obj$classname))
   
   if (!is.null(obj$info))
      cat(sprintf('Info: %s\n', obj$info))
   
   cat(sprintf('Significance level (alpha): %.2f\n', obj$alpha))
   cat(sprintf('Selected number of components: %d\n\n', obj$ncomp.selected))
   
   data = cbind(round(obj$calres$expvar, 2),
                round(obj$calres$cumexpvar, 2),
                round(obj$calres$sensitivity[1, ], 2)
   )   
   colnames(data) = c('Expvar', 'Cumexpvar', 'Sens (cal)')
   
   if (!is.null(obj$cvres))
   {
      cnames = colnames(data)
      data = cbind(data,
                   round(obj$cvres$sensitivity[1, ], 2)
      )
      colnames(data) = c(cnames, 'Sens (cv)')
   }   
   
   if (!is.null(obj$testres))
   {
      cnames = colnames(data)
      if (is.null(obj$testres$specificity[1, ]) || min(obj$testres$specificity[1, ]) == 1)
      {
         data = cbind(data,
                      round(obj$testres$sensitivity[1, ], 2)
         )
         colnames(data) = c(cnames, 'Sens (test)')         
      }
      else  
      {   
         data = cbind(data,
                      round(obj$testres$specificity[1, ], 2),
                      round(obj$testres$sensitivity[1, ], 2)
         )
         colnames(data) = c(cnames, 'Spec (test)', 'Sens (test)')
      }
   }   
   
   print(data)   
}  

#' Print method for SIMCA model object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a SIMCA model (object of class \code{simca})
#' @param ...
#' other arguments
#' 
#' @export
print.simca = function(x, ...)
{
   obj = x
   
   cat('\nSIMCA one class model (class simca)\n')
   
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$classname - name of the class\n')
   cat('$alpha - significance level\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$loadings - matrix with loadings\n')
   cat('$eigenvals - eigenvalues for components\n')
   cat('$center - values for centering data\n')
   cat('$scale - values for scaling data\n')
   cat('$info - information about the model\n')
   cat('$cv - number of segments for cross-validation\n')
   cat('$calres - results (scores, etc) for calibration set\n')
   
   if (!is.null(obj$cvres))
   {
      cat('$cvres - results for cross-validation\n')      
   }   
   if (!is.null(obj$testres))
   {
      cat('$testres - results for test set\n')      
   }       
}  

