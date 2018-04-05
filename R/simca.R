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
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values)
#' @param x.test
#' a numerical matrix with test data.
#' @param c.test
#' a vector with classes of test data objects (can be text with names of classes or logical).
#' @param method
#' method to compute principal components.
#' @param rand
#' vector with parameters for randomized PCA methods (if NULL, conventional PCA is used instead)
#' @param lim.type
#' which method to use for calculation of critical limits for residuals (see details)
#' @param alpha
#' significance level for calculating critical limits for T2 and Q residuals.
#' @param gamma
#' significance level for calculating outlier limits for T2 and Q residuals.
#' @param info
#' text with information about the model.
#' 
#' @details 
#' SIMCA is in fact PCA model with additional functionality, so \code{simca} class inherits most 
#' of the functionality of \code{\link{pca}} class. It uses critical limits calculated for Q and T2 
#' residuals calculated for PCA model for making classification decistion.
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
simca = function(x, classname, ncomp = 15, center = T, scale = F, cv = NULL, exclcols = NULL,
                 exclrows = NULL, x.test = NULL,  c.test = NULL, method = 'svd', rand = NULL, 
                 lim.type = 'jm', alpha = 0.05, gamma = 0.01, info = '') {
   
   if (!is.character(classname))
      stop('Argument "classname" must be a text!')
   
   if (length(classname) > 20)
      stop('Argument "classname" must have up to 20 symbols!')

   # add proper attributes if some rows or columns must be excluded
   if (length(exclcols) > 0)
      x = mda.exclcols(x, exclcols)
   
   if (length(exclrows) > 0)
      x = mda.exclrows(x, exclrows)
   
   # calibrate model  
   model = pca.cal(x, ncomp, center = center, scale = scale, method = method, rand = rand,
                   lim.type = lim.type, alpha = alpha, gamma = gamma, info = info, cv = NULL)
   model$nclasses = 1
   model$classname = classname
   model$call = match.call()   
   class(model) = c("simca", "classmodel", "pca")
   
   # apply model to calibration set
   model$calres = predict.simca(model, x, c.ref = rep(classname, nrow(x)), cal = TRUE)
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = simca.crossval(model, x, cv, center = center, scale = scale)
   
   # apply model to test set if provided
   if (!is.null(x.test)){
      # if classes are not provided we assume the object are from the same class
      if (is.null(c.test)) 
         c.test = matrix(classname, nrow = nrow(x.test), ncol = 1) 
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
#' @param cal
#' logical, are predictions for calibration set or not
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
predict.simca = function(object, x, c.ref = NULL, cal = FALSE, ...) {
   if (cal == FALSE) {
      pres = predict.pca(object, x, cal = cal)     
   } else {
      pres = object$calres
   }
  
   # check reference values
   c.ref = checkReferenceValues.classmodel(object, c.ref, x)
  
   # do predictions and set attributes
   cres = simca.classify(object, pres)
   dimnames(cres$c.pred) = dimnames(cres$p.pred) = list(
      rownames(x), 
      colnames(object$loadings), 
      object$classname
   )
   cres$c.pred = mda.setattr(cres$c.pred, mda.getattr(x), 'row')
   attr(cres$c.pred, 'name') = 'Class, predicted values'
   cres$p.pred = mda.setattr(cres$p.pred, mda.getattr(x), 'row')
   attr(cres$p.pred, 'name') = 'Class, probabilities'
   
   # combine everything to simcares object
   cres = classres(
      c.pred = cres$c.pred, 
      p.pred = cres$p.pred,
      c.ref = c.ref, 
      ncomp.selected = object$ncomp.selected
   )
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
simca.classify = function(model, res) {
   ncomp = model$ncomp
   nobj = nrow(res$Q)
   c.pred = array(0, dim = c(nobj, ncomp, 1))
   p.pred = array(0, dim = c(nobj, ncomp, 1))
   
   for (i in 1:ncomp) {
      p.pred[, i, 1] = getProbabilities(model, i, res$Q[, i], res$T2[, i])
      c.pred[, i, 1] = p.pred[, i, 1] >= 0.5
   }   
   c.pred = c.pred * 2 - 1
   res = list(c.pred = c.pred, p.pred = p.pred)
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
simca.crossval = function(model, x, cv, center = T, scale = F) {
   ncomp = model$ncomp   
   
   # convert data to a matrix 
   attrs = mda.getattr(x)
   x = mda.df2mat(x)
   x.nrows = nrow(x)
   x.ncols = ncol(x)
   
   # remove excluded rows 
   if (length(attrs$exclrows) > 0)
      x = x[-attrs$exclrows, , drop = F]
   
   # remove excluded columns
   if (length(attrs$exclcols) > 0)
      x = x[, -attrs$exclcols, drop = F]
   
   # get matrix with indices for cv segments
   nobj = nrow(x)
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
   for (iRep in 1:nrep) {
      for (iSeg in 1:nseg) {
         ind = na.exclude(idx[iSeg, ,iRep])
      
         if (length(ind) > 0) {   
            x.cal = x[-ind, , drop = F]
            x.val = x[ind, , drop = F]
            
            # autoscale calibration set
            x.cal = prep.autoscale(x.cal, center = center, scale = scale)
            
            # get loadings
            m = pca.run(x.cal, ncomp, model$method)               
            
            # apply autoscaling to the validation set
            x.val = prep.autoscale(x.val, center = attr(x.cal, 'prep:center'), 
                                   scale = attr(x.cal, 'prep:scale'))
            
            # get scores
            scores = x.val %*% m$loadings
            residuals = x.val - tcrossprod(scores, m$loadings)
            
            # compute distances
            res = ldecomp.getDistances(scores, m$loadings, residuals, model$tnorm)
            Q[ind, ] = Q[ind, ] + res$Q
            T2[ind, ] = T2[ind, ] + res$T2
         }
      }  
   }
   
   Q = Q / nrep;
   T2 = T2 / nrep;
   r = list(Q = Q, T2 = T2, classname = model$classname)
   
   # classify data
   cres = simca.classify(model, r)
   dimnames(cres$c.pred) = dimnames(cres$p.pred) = list(
      rownames(x), 
      colnames(model$loadings), 
      model$classname
   )
   
   # add names
   rownames(Q) = rownames(T2) = rownames(c.pred) = rownames(c.ref) = rownames(x)
   colnames(Q) = colnames(T2) = colnames(c.pred) = colnames(model$loadings)
   
   # add attributes
   attr(Q, 'name') = 'Squared residual distance (Q)'
   attr(Q, 'xaxis.name') = attr(model$calres$scores, 'xaxis.name')
   attr(Q, 'yaxis.name') = attr(model$calres$scores, 'yaxis.name')
   attr(Q, 'yaxis.values') = attr(model$calres$scores, 'yaxis.values')
   
   attr(T2, 'name') = 'T2 residuals'
   attr(T2, 'xaxis.name') = attr(model$calres$scores, 'xaxis.name')
   attr(T2, 'yaxis.name') = attr(model$calres$scores, 'yaxis.name')
   attr(T2, 'yaxis.values') = attr(model$calres$scores, 'yaxis.values')
   
   attr(c.ref, 'name') = 'Class, reference values'
   attr(c.ref, 'xaxis.name') = attr(model$calres$scores, 'xaxis.name')
   attr(c.ref, 'yaxis.name') = attr(model$calres$scores, 'yaxis.name')
   attr(c.ref, 'yaxis.values') = attr(model$calres$scores, 'yaxis.values')
   
   attr(cres$c.pred, 'name') = 'Class, predicted values'
   attr(cres$c.pred, 'xaxis.name') = attr(model$calres$scores, 'xaxis.name')
   attr(cres$c.pred, 'yaxis.name') = attr(model$calres$scores, 'yaxis.name')
   attr(cres$c.pred, 'yaxis.values') = attr(model$calres$scores, 'yaxis.values')
   
   attr(cres$c.pred, 'name') = 'Class, probabilities'
   attr(cres$c.pred, 'xaxis.name') = attr(model$calres$scores, 'xaxis.name')
   attr(cres$c.pred, 'yaxis.name') = attr(model$calres$scores, 'yaxis.name')
   attr(cres$c.pred, 'yaxis.values') = attr(model$calres$scores, 'yaxis.values')
   
   # compute variance
   var = ldecomp.getVariances(Q, model$calres$totvar)
   
   # in CV results there are no scores nor residuals, only residual distances and variances
   pres = pcares(ncomp.selected = model$ncomp.selected, dist = list(T2 = T2, Q = Q), var = var)
   pres$Qlim = model$Qlim
   pres$T2lim = model$T2lim
   pres$lim.type = model$lim.type
   pres$alpha = model$alpha
   pres$gamma = model$gamma

   # combine results together   
   cres = classres(c.pred = cres$c.pred, p.pred = cres$p.pred, c.ref = c.ref)   
   res = simcares(pres, cres)
   
   res
}

#' Probability of class belonging for PCA/SIMCA results
#' 
#' @description 
#' Computes probability of class belonging for each object based on Q and T2 residuls
#' 
#' @param obj
#' object with SIMCA model
#' @param ncomp
#' number of components to compute the probabilities for.
#' @param Q
#' vector with Q values for selected component
#' @param T2
#' vector with T2 values for selected component
#' @param ...
#' other arguments
#'  
#' @export
getProbabilities.simca = function(obj, ncomp, Q, T2, ...) {
   p = NULL
   Qlim = obj$Qlim[, ncomp]
   T2lim = obj$T2lim[, ncomp]
   lim.type = obj$lim.type
   alpha = obj$alpha
   
   # compute class probability based statistic probability
   p = function(x, alpha){ ifelse((p = (1 - x) / alpha * 0.5) > 1, 1, p) }
   
   # probability for T2 by Hotelling
   T2.p = p(reslim.hotelling(ncomp, T2 = T2, T2lim = T2lim, return = 'prob'), alpha)

   # probability for Q by selected method
   if (obj$lim.type == 'chisq') {
      Q.p = p(reslim.chisq(Q = Q, Qlim = Qlim, return = 'prob'), alpha)
      p = pmin(Q.p, T2.p)
   } else if (obj$lim.type == 'jm') {
      Q.p = p(reslim.jm(obj$eigenvals, Q = Q, ncomp = ncomp, return = 'prob'), alpha)
      p = pmin(Q.p, T2.p)
   } else if (substr(obj$lim.type, 1, 2) == 'dd') {
      p = p(reslim.dd(Q, T2, Qlim = Qlim, T2lim = T2lim, type = lim.type, return = 'prob'), alpha)
   } else {
      stop('Wrong value for "type" parameter!')
   }
   
   p
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
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @export
plotModellingPower.simca = function(obj, ncomp = NULL, type = 'h', main = NULL, ylab = '', ...){
   if (is.null(main)) {
      if (is.null(ncomp))
         main = 'Modelling power'
      else
         main = sprintf('Modelling power (ncomp = %d)', ncomp)      
   }   

   ncomp = getSelectedComponents(obj, ncomp)
   
   nvar = nrow(obj$modpower)
   if (is.null(type)) {   
      if (nvar < 20)
         type = 'h'
      else
         type = 'l'
   }
   
   data = mda.subset(obj$modpower, select = ncomp)
   mdaplot(mda.t(data), type = type, ylab = ylab, main = main, ...)
}   

#' Shows extreme plot for SIMCA model
#' 
#' @param obj
#' SIMCA model
#' @param ncomp
#' Number of components to show the plot for
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other arguments
#'  
#' @description 
#' The plot shows the number of extreme objects rejected by the model vs. the expected number, 
#' which depends on significance level and total number of objects. The light blue area shows 95% 
#' tolerance limits. See more details in [1].
#'  
#' @references 
#' 1. A.L. Pomerantsev, O.Ye. Rodionova, Concept and role of extreme objects in PCA/SIMCA,
#' Journal of Chemometrics, 28 (2014) pp. 429-438.
#' 
#' @export   
plotExtreme.simca = function(obj, ncomp = NULL, main = NULL, xlab = 'Expected', 
                             ylab = 'Observed', ...) {
   
   ncomp = getSelectedComponents(obj, ncomp)
   
   if (is.null(main)) {
      if (is.null(ncomp))
         main = 'Extreme plot'
      else
         main = sprintf('Extreme plot (ncomp = %d)', ncomp)      
   }   
   
   Q = obj$calres$Q[, ncomp]
   T2 = obj$calres$T2[, ncomp]
   attrs = mda.getattr(obj$calres$Q)
   if (length(attrs$exclrows) > 0) {
      Q = Q[-attrs$exclrows]
      T2 = T2[-attrs$exclrows]
   }
   
   T2.p = reslim.hotelling(ncomp, T2 = T2, T2lim = obj$T2lim[, ncomp], return = 'prob')
   if (obj$lim.type == 'chisq') {
      Q.p = reslim.chisq(Q, Qlim = obj$Qlim[, ncomp], return = 'prob')
      p = pmax(Q.p, T2.p)
   } else if (obj$lim.type == 'jm') {
      Q.p = reslim.jm(obj$eigenvals, Q, ncomp, return = 'prob')
      p = pmax(Q.p, T2.p)
   } else if (substr(obj$lim.type, 1, 2) == 'dd') {
      p = reslim.dd(Q, T2, type = obj$lim.type, Qlim = obj$Qlim[, ncomp], 
                    T2lim = obj$T2lim[, ncomp], return = 'prob')
   } else {
      stop('Wrong value for "type" parameter!')
   }
   
   p = sort(p)
   nobj = length(p)
   N = Np = Nm = rep(0, nobj)
   for (i in 1:nobj) {
      alpha = i / nobj
      N[i] = sum((1 - p) < alpha)
      D = 2 * sqrt(i * (1 - alpha))
      Np[i] = i + D;
      Nm[i] = i - D;
   }
   
   expected = 1:nobj
   plot(expected, expected, type = 'l', col = 'blue', main = main, xlab = xlab, ylab = ylab)
   for (i in 1:nobj){
      lines(c(expected[i],expected[i]),c(Nm[i],Np[i]), col = 'lightblue')
   }
   lines(expected, Nm, col = 'lightblue')
   lines(expected, Np, col = 'lightblue')
   points(expected, N, type = 'p', pch = 16, col = 'red')
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
plot.simca = function(x, ncomp = NULL, ...) {
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
summary.simca = function(object, ...) {
   obj = object
   
   cat(sprintf('\nSIMCA model for class "%s" summary\n\n', obj$classname))
   
   if (!is.null(obj$info))
      cat(sprintf('Info: %s\n', obj$info))
   
   
   cat(sprintf('Method for critical limits: %s\n', obj$lim.type))
   cat(sprintf('Significance level (alpha): %.2f\n', obj$alpha))
   cat(sprintf('Selected number of components: %d\n\n', obj$ncomp.selected))
   
   data = cbind(obj$calres$expvar, obj$calres$cumexpvar, obj$calres$sensitivity[1, ])
   colnames(data) = c('Expvar', 'Cumexpvar', 'Sens (cal)')
   
   if (!is.null(obj$cvres)) {
      cnames = colnames(data)
      data = cbind(data, obj$cvres$expvar, obj$cvres$sensitivity[1, ])
      colnames(data) = c(cnames, 'Expvar (cv)', 'Sens (cv)')
   }   
   
   if (!is.null(obj$testres)) {
      cnames = colnames(data)
      if (!is.null(obj$testres$c.ref) && !any(is.nan(obj$testres$specificity))) {
         data = cbind(
            data, 
            obj$testres$expvar, 
            obj$testres$specificity[1, ], 
            obj$testres$sensitivity[1, ]
         )
         colnames(data) = c(cnames, 'Expvar (test)', 'Spec (test)', 'Sens (test)')
      } else {
         data = cbind(data, obj$testres$expvar, obj$testres$sensitivity[1, ])
         colnames(data) = c(cnames, 'Expvar (test)', 'Sens (test)')
      }
   }   
   
   print(round(data, 2))   
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
print.simca = function(x, ...) {
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

