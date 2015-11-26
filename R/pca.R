#' Principal Component Analysis
#'
#' @description 
#' \code{pca} is used to build and explore a principal component analysis (PCA) model.
#'
#' @param x
#' a numerical matrix with calibration data.
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
#' @param alpha
#' significance level for calculating limit for Q residuals.
#' @param method
#' method to compute principal components.
#' @param info
#' a short text line with model description.
#'
#' @details 
#' So far only SVD (Singular Value Decompisition) method is available, more coming soon. 
#'   
#' By default \code{pca} uses number of components (\code{ncomp}) as a minimum of number of 
#' objects - 1, number of variables and default or provided value. Besides that, there is also 
#' a parameter for selecting an optimal number of components (\code{ncomp.selected}). The optimal 
#' number of components is used to build a residuals plot (with Q residuals vs. Hotelling T2 
#' values), calculate confidence limits for Q residuals, as well as for SIMCA classification. 
#'   
#' If data contains missing values (NA) the \code{pca} will use an iterative algorithm to fit the 
#' values with most probable ones. The algorithm is implemented in a function 
#' \code{\link{pca.mvreplace}}. The same center and scale options will be used. You can also
#' do this step manually before calling \code{pca} and play with extra options.
#' 
#' @return 
#' Returns an object of \code{pca} class with following fields:
#' \item{ncomp }{number of components included to the model.} 
#' \item{ncomp.selected }{selected (optimal) number of components.} 
#' \item{loadings }{matrix with loading values (nvar x ncomp).} 
#' \item{eigenvals }{vector with eigenvalues for all existent components.} 
#' \item{expvar }{vector with explained variance for each component (in percent).} 
#' \item{cumexpvar }{vector with cumulative explained variance for each component (in percent).} 
#' \item{T2lim }{statistical limit for T2 distance.} 
#' \item{Qlim }{statistical limit for Q residuals.} 
#' \item{info }{information about the model, provided by user when build the model.} 
#' \item{calres }{an object of class \code{\link{pcares}} with PCA results for a calibration data.} 
#' \item{testres }{an object of class \code{\link{pcares}} with PCA results for a test data, if it 
#' was provided.} 
#' \item{cvres }{an object of class \code{\link{pcares}} with PCA results for cross-validation, 
#' if this option was chosen.} 
#' 
#' @author 
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @seealso 
#' Methods for \code{pca} objects:
#' \tabular{ll}{
#'    \code{plot.pca} \tab makes an overview of PCA model with four plots.\cr
#'    \code{summary.pca} \tab shows some statistics for the model.\cr
#'    \code{\link{selectCompNum.pca}} \tab set number of optimal components in the model\cr
#'    \code{\link{predict.pca}} \tab applies PCA model to a new data.\cr
#'    \code{\link{plotScores.pca}} \tab shows scores plot.\cr
#'    \code{\link{plotLoadings.pca}} \tab shows loadings plot.\cr
#'    \code{\link{plotVariance.pca}} \tab shows explained variance plot.\cr
#'    \code{\link{plotCumVariance.pca}} \tab shows cumulative explained variance plot.\cr
#'    \code{\link{plotResiduals.pca}} \tab shows Q vs. T2 residuals plot.\cr
#' }
#'  Most of the methods for plotting data are also available for PCA results (\code{\link{pcares}})
#'  objects. Also check \code{\link{pca.mvreplace}}, which replaces missing values in a data matrix 
#'  with approximated using iterative PCA decomposition.
#'  
#' @examples 
#' library(mdatools)
#' ### Examples for PCA class
#' 
#' ## 1. Make PCA model for People data with autoscaling
#' ## and full cross-validation
#' 
#' data(people)
#' model = pca(people, scale = TRUE, cv = 1, info = 'Simple PCA model')
#' model = selectCompNum(model, 4)
#' summary(model)
#' plot(model, show.labels = TRUE)
#' 
#' ## 2. Add missing values, make a new model and show plots
#' peoplemv = people
#' peoplemv[2, 7] = NA
#' peoplemv[6, 2] = NA
#' peoplemv[10, 4] = NA
#' peoplemv[22, 12] = NA
#' 
#' modelmv = pca(peoplemv, scale = TRUE, info = 'Model with missing values')
#' modelmv = selectCompNum(modelmv, 4)
#' summary(modelmv)
#' plot(modelmv, show.labels = TRUE)
#' 
#' ## 3. Show scores and loadings plots for the model
#' par(mfrow = c(2, 2))
#' plotScores(model, comp = c(1, 3), show.labels = TRUE)
#' plotScores(model, comp = 2, type = 'h', show.labels = TRUE)
#' plotLoadings(model, comp = c(1, 3), show.labels = TRUE)
#' plotLoadings(model, comp = c(1, 2), type = 'h', show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' ## 4. Show residuals and variance plots for the model
#' par(mfrow = c(2, 2))
#' plotVariance(model, type = 'h')
#' plotCumVariance(model, show.labels = TRUE, legend.position = 'bottomright')
#' plotResiduals(model, show.labels = TRUE)
#' plotResiduals(model, ncomp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' @export   
pca = function(x, ncomp = 15, center = T, scale = F, cv = NULL, x.test = NULL, 
               alpha = 0.05, method = 'svd', info = '')
{
   x = as.matrix(x)

   # check if data has missing values
   if (sum(is.na(x)) > 0)
   {
      warning('Data has missing values, will try to fix using pca.mvreplace.')
      x = pca.mvreplace(x, center = center, scale = scale)
   }   
   
   # correct maximum number of components
   ncomp = min(ncomp, ncol(x), nrow(x) - 1)

   # calibrate model  
   model = pca.cal(x, ncomp, center = center, scale = scale, method = method)
   model$ncomp = ncomp
   model$ncomp.selected = model$ncomp
   model$info = info
   model$alpha = alpha
   
   # apply model to calibration set
   model$calres = predict.pca(model, x)
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = pca.crossval(model, x, cv, center = center, scale = scale)
   
   # apply model to test set if provided
   if (!is.null(x.test))
      model$testres = predict.pca(model, x.test)
   
   # calculate and assign limit values for T2 and Q residuals
   lim = ldecomp.getResLimits(model$eigenvals, nrow(x), model$ncomp, model$alpha)
   model$T2lim = lim$T2lim
   model$Qlim = lim$Qlim
   
   model$call = match.call()   
   class(model) = "pca"
   
   model
}

#' Get calibration data
#' 
#' @description
#' Get data, used for calibration of the PCA model
#' 
#' @param obj
#' PCA model (object of class \code{pca})
#' @param ...
#' other parameters
#' 
getCalibrationData.pca = function(obj, ...)
{
   x = obj$calres$scores %*% t(obj$loadings) + obj$calres$residuals
   
   if (is.numeric(attr(x, 'prep:scale')))
      x = sweep(x, 2L, attr(x, 'prep:scale'), '*', check.margin = F)
   
   if (is.numeric(attr(x, 'prep:center')))
      x = sweep(x, 2L, attr(x, 'prep:center'), '+', check.margin = F)
   
   x
}

#' Select optimal number of components for PCA model
#' 
#' @description
#' Allows user to select optimal number of components for PCA model
#' 
#' @param model
#' PCA model (object of class \code{pca})
#' @param ncomp
#' number of components to select
#' 
#' @return
#' the same model with selected number of components
#' 
#' @export
selectCompNum.pca = function(model, ncomp)
{
   if (ncomp < 1 || ncomp > model$ncomp)
      stop('Wrong number of selected components!')
   
   model$ncomp.selected = ncomp   
   
   model$calres$ncomp.selected = ncomp
   
   if (!is.null(model$testres))
      model$testres$ncomp.selected = ncomp

   if (!is.null(model$cvres))
      model$cvres$ncomp.selected = ncomp

   model
}

#' Replace missing values in data
#' 
#' \code{pca.mvreplace} is used to replace missing values in a data matrix with 
#' approximated by iterative PCA decomposition.
#'
#' @param x
#' a matrix with data, containing missing values.
#' @param center
#' logical, do centering of data values or not.
#' @param scale
#' logical, do standardization of data values or not.
#' @param maxncomp
#' maximum number of components in PCA model.
#' @param expvarlim
#' minimum amount of variance, explained by chosen components (used for selection of optimal number 
#' of components in PCA models).
#' @param covlim
#' convergence criterion.
#' @param maxiter
#' maximum number of iterations if convergence criterion is not met.
#'
#' @details 
#' The function uses iterative PCA modeling of the data to approximate and impute missing values.  
#' The result is most optimal for data sets with low or moderate level of noise and with number of
#' missing values less than 10\% for small dataset and up to 20\% for large data.
#'
#' @return 
#' Returns the same matrix \code{x} where missing values are replaced with approximated.
#' 
#' @references 
#' Philip R.C. Nelson, Paul A. Taylor, John F. MacGregor. Missing data methods in PCA and PLS: 
#' Score calculations with incomplete observations. Chemometrics and Intelligent Laboratory 
#' Systems, 35 (1), 1996.
#'
#' @author 
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @examples
#' library(mdatools)
#' 
#' ## A very simple example of imputing missing values in a data with no noise
#' 
#' # generate a matrix with values
#' s = 1:6
#' odata = cbind(s, 2*s, 4*s)
#' 
#' # make a matrix with missing values
#' mdata = odata
#' mdata[5, 2] = mdata[2, 3] = NA
#' 
#' # replace missing values with approximated
#' rdata = pca.mvreplace(mdata, scale = TRUE)
#' 
#' # show all matrices together
#' show(cbind(odata, mdata, round(rdata, 2)))
#' 
#' @export
pca.mvreplace = function(x, center = T, scale = F, maxncomp = 7,
                         expvarlim = 0.95, covlim = 10^-6, maxiter = 100)
{
   x.rep = x
   mvidx = is.na(x.rep)

   # calculate number of missing values for every variable
   # and make initial estimates with mean values
   for (i in 1:ncol(x))
   {
      mv = is.na(x[, i])
      
      if (sum(mv)/length(x[, i]) > 0.2)
         stop(sprintf('To many missing values in column #%d', i))
      
      x.rep[mv, i] = mean(x[, i], na.rm = T)        
   }  
   
   # autoscale 
   x.rep = scale(x.rep, center = center, scale = scale)
   
   if (scale == T)
      gsd = attr(x.rep, 'scaled:scale')
   
   if (center == T)
      gmean = attr(x.rep, 'scaled:center');         
   
   x = x.rep
   
   n = 1
   scoresp = 0
   scores = 1
   cond = 1
   while (cond > covlim && n < maxiter)
   {    
      n = n + 1
      
      # rescale data on every iteration
      x.rep = scale(x.rep, center = T, scale = F)
      lmean = attr(x.rep, 'scaled:center')
      
      res = pca.svd(x.rep, maxncomp)
      
      expvar = cumsum(res$eigenvals/sum(res$eigenvals))
      ncomp = min(which(expvar >= expvarlim), maxncomp)
            
      if (ncomp == 0)
         ncomp = 1
      if (ncomp == length(expvar))
         ncomp = ncomp - 1
      
      # get and trancate scores and loadings and reestimate the values
      scoresp = scores
      loadings = res$loadings[, 1:ncomp]      
      scores = x.rep %*% loadings
      x.new = scores %*% t(loadings)   
      
      # remove centering
      x.new = sweep(x.new, 2L, lmean, '+', check.margin = F)

      x.rep = x
      x.rep[mvidx] = x.new[mvidx]
      
      if (n > 2)
      {
         # calculate difference between scores for convergence 
         ncompcond = min(ncol(scores), ncol(scoresp))
         cond = sum((scores[, 1:ncompcond] - scoresp[, 1:ncompcond])^2)
      }      
   }   
   
   # rescale the data back and return
   if (scale == T)
      x.rep = sweep(x.rep, 2L, gsd, '*', check.margin = F)
   
   if (center == T)
      x.rep = sweep(x.rep, 2L, gmean, '+', check.margin = F)

   x.rep
}

#' PCA model calibration
#' 
#' @description
#' Calibrates (builds) a PCA model for given data and parameters
#' 
#' @param x
#' matrix with data values
#' @param ncomp
#' number of principal components to calculate
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#' @param method
#' algorithm for compiting PC space (only 'svd' is supported so far)
#' 
#' @return
#' an object with calibrated PCA model
#' 
pca.cal = function(x, ncomp, center = T, scale = F, method = 'svd')
{
   x = prep.autoscale(x, center = center, scale = scale)
   model = pca.svd(x, ncomp)
   
   model$tnorm = sqrt(colSums(model$scores ^ 2)/(nrow(model$scores) - 1));   
   
   rownames(model$loadings) = colnames(x)
   colnames(model$loadings) = paste('Comp', 1:ncol(model$loadings))
   model$center = attr(x, 'prep:center')
   model$scale = attr(x, 'prep:scale')
      
   model
}  

#' Singular Values Decomposition based PCA algorithm
#' 
#' @description
#' Computes principal component space using Singular Values Decomposition
#' 
#' @param x
#' a matrix with data values (preprocessed)
#' @param ncomp
#' number of components to calculate
#' 
#' @return
#' a list with scores, loadings and eigencalues for the components
#' 
pca.svd = function(x, ncomp = NULL)
{
   if (is.null(ncomp)) 
      ncomp = min(ncol(x), nrow(x) - 1)
   else
      ncomp = min(ncomp, ncol(x), nrow(x) - 1)
   
   s = svd(x)
   loadings = s$v[, 1:ncomp]
      
   res = list(
      loadings = loadings,
      scores = x %*% loadings,
      eigenvals = (s$d^2)/(nrow(x) - 1)
   )
   
   res
}

#' NIPALS based PCA algorithm
#' 
#' @description
#' Calculates principal component space using non-linear iterative partial least squares algorithm 
#' (NIPALS)
#' 
#' @param x
#' a matrix with data values (preprocessed)
#' @param ncomp
#' number of components to calculate
#' 
#' @return
#' a list with scores, loadings and eigencalues for the components
#' 
#' @references
#' Geladi, Paul; Kowalski, Bruce (1986), "Partial Least Squares 
#' Regression:A Tutorial", Analytica Chimica Acta 185: 1-17 
#'    
pca.nipals = function(x, ncomp)
{
   nobj = nrow(x)
   nvar = ncol(x)   
   ncomp = min(ncomp, nobj - 1, nvar)
   
   scores = matrix(0, nrow = nobj, ncol = ncomp)
   loadings = matrix(0, nrow = nvar, ncol = ncomp)
   eigenvals = rep(0, ncomp)
   
   E = x
   for (i in 1:ncomp)
   {      
      ind = which.max(apply(E, 2, sd))
      t = E[, ind, drop = F]
      tau = 99999
      th = 9999

      while (th > 0.000001)
      {      
         p = (t(E) %*% t) / as.vector((t(t) %*% t))
         p = p / as.vector(t(p) %*% p) ^ 0.5
         t = (E %*% p)/as.vector(t(p) %*% p)
         th = abs(tau - as.vector(t(t) %*% t))
         tau = as.vector(t(t) %*% t)
      }
            
      E = E - t %*% t(p)
      scores[, i] = t
      loadings[, i] = p
      eigenvals[i] = tau / (nobj - 1)
   }

   res = list(
      loadings = loadings,
      scores = scores,
      eigenvals = eigenvals
   )   
}

#' Cross-validation of a PCA model
#' 
#' @description
#' Does the cross-validation of a PCA model
#' 
#' @param model
#' a PCA model (object of class \code{pca})
#' @param x
#' a matrix with data values (calibration set)
#' @param cv
#' number of segments (if cv = 1, full cross-validation will be used)
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#'
#' @return
#' object of class \code{pcares} with results of cross-validation
#'  
pca.crossval = function(model, x, cv, center = T, scale = F)
{      
   ncomp = model$ncomp   
   nobj = nrow(x)
   
   # get matrix with indices for cv segments
   idx = crossval(nobj, cv)
   nrep = dim(idx)[3]
      
   Q = matrix(0, ncol = ncomp, nrow = nobj)   
   T2 = matrix(0, ncol = ncomp, nrow = nobj)   
   
   # loop over repetitions and segments
   
   for (iRep in 1:nrep)
   {   
      for (iSeg in 1:nrow(idx))
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
         }
      }  
   }
   
   Q = Q / nrep
   T2 = T2 / nrep
   rownames(Q) = rownames(T2) = rownames(x)
   colnames(Q) = colnames(T2) = colnames(model$scores)

   # in CV results there are no scores only residuals and variances
   res = pcares(NULL, NULL, NULL, model$calres$totvar, model$tnorm, model$ncomp.selected,
                T2, Q)
   res$Qlim = model$Qlim
   res$T2lim = model$T2lim
   
   res
}  

#' PCA predictions
#' 
#' @description
#' Applies PCA model to a new data
#' 
#' @param object
#' a PCA model (object of class \code{pca})
#' @param x
#' a matrix with data values
#' @param cv
#' logical, are predictions for cross-validation or not
#' @param ...
#' other arguments
#' 
#' @return
#' PCA results (an object of class \code{pcares})
#'  
#' @export
predict.pca = function(object, x, cv = F, ...)
{
   x = prep.autoscale(x, object$center, object$scale)
   scores = x %*% object$loadings
   residuals = x - scores %*% t(object$loadings)
   
   if (cv == F)
   {   
      totvar = sum(x^2)
      res = pcares(scores, object$loadings, residuals, totvar, object$tnorm, object$ncomp.selected)
      res$Qlim = object$Qlim
      res$T2lim = object$T2lim
   }   
   else
   {
      res = ldecomp.getDistances(scores, object$loadings, residuals, object$tnorm)   
   }

   res
}  

#' Explained variance plot for PCA
#' 
#' @description
#' Shows a plot with explained variance or cumulative explained variance for components.
#' 
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param type
#' type of the plot ('b', 'l', 'h')
#' @param variance
#' which variance to use ('expvar', 'cumexpvar')
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{pca}} function.
#' 
#' @export
plotVariance.pca = function(obj, type = 'b', variance = 'expvar', 
                            main = 'Variance', xlab = 'Components', 
                            ylab = 'Explained variance, %',
                            show.legend = T, show.labels = F, ...)
{
   data = cbind(1:length(obj$calres[[variance]]), obj$calres[[variance]])
   labels = mdaplot.formatValues(obj$calres[[variance]])
   legend  = 'cal'
   
   if (!is.null(obj$cvres))
   {
      data = cbind(data, obj$cvres[[variance]])
      labels = cbind(labels, mdaplot.formatValues(obj$cvres[[variance]]))
      legend = c(legend, 'cv')
   }      

   if (!is.null(obj$testres))
   {
      data = cbind(data, obj$testres[[variance]])
      labels = cbind(labels, mdaplot.formatValues(obj$testres[[variance]]))
      legend = c(legend, 'test')
   }      
      
   if (show.legend == F)
      legend = NULL
   
   if (show.labels == F)
      labels = NULL
   
   mdaplotg(data, main = main, xlab = xlab, ylab = ylab, labels = labels, legend = legend, 
            type = type, show.labels = show.labels, ...)   
}

#' Cumulative explained variance plot for PCA
#' 
#' @description
#' Shows a plot with cumulative explained variance for components.
#' 
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{pca}} function.
#' 
#' @export
plotCumVariance.pca = function(obj, xlab = 'Components', ylab = 'Explained variance, %', 
                               main = 'Cumulative variance', ...)
{
   plotVariance.pca(obj, variance = 'cumexpvar', xlab = xlab, ylab = ylab, main = main, ...)   
}

#' Scores plot for PCA
#' 
#' @description
#' Shows a scores plot for selected components.
#' 
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param comp
#' a value or vector with several values - number of components to show the plot for
#' @param type
#' type of the plot ('b', 'l', 'h')
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{pca}} function.
#' 
#' @export
plotScores.pca = function(obj, comp = c(1, 2), type = 'p', main = 'Scores', xlab = NULL, 
                          ylab = NULL, show.labels = F, show.legend = T,
                          show.axes = T, ...)
{
   ncomp = length(comp)
   legend = NULL
   
   if (ncomp == 1)
   {   
      # scores vs objects plot
      
      if (comp > obj$ncomp || comp < 1)
         stop('Wrong number of components!')
      
      if (is.null(xlab))
         xlab = 'Objects'
      
      if (is.null(ylab))
         ylab = colnames(obj$calres$scores)[comp]
      
      nobj.cal = nrow(obj$calres$scores)      
      cdata = cbind(1:nobj.cal, obj$calres$scores[, comp])      
      colnames(cdata) = c(xlab, ylab)
      rownames(cdata) = rownames(obj$calres$scores)
      
      data = list(cdata = cdata)
      
      if (!is.null(obj$testres))
      {
         nobj.test = nrow(obj$testres$scores)
         tdata = cbind((nobj.cal + 1):(nobj.cal + nobj.test), obj$testres$scores[, comp])      
         rownames(tdata) = rownames(obj$testres$scores)
         data$tdata = tdata
         if (show.legend == T)
            legend = c('cal', 'test')
      }   
      
      mdaplotg(data, type = type, main = main, show.labels = show.labels, legend = legend, 
               xlab = xlab, ylab = ylab, ...)
   }
   else if (ncomp == 2)
   {
      # scores vs scores plot
      
      if (comp[1] > obj$ncomp || comp[1] < 1 || comp[2] > obj$ncomp || comp[2] < 1)
         stop('Wrong component numbers!')
      
      if (is.null(xlab))
         xlab = colnames(obj$calres$scores)[comp[1]]
      
      if (is.null(ylab))
         ylab = colnames(obj$calres$scores)[comp[2]]

      cdata = cbind(obj$calres$scores[, comp[1]], obj$calres$scores[, comp[2]])      
      rownames(cdata) = rownames(obj$calres$scores)
      
      data = list(cdata = cdata)
      
      if (!is.null(obj$testres))
      {
         tdata = cbind(obj$testres$scores[, comp[1]], obj$testres$scores[, comp[2]])      
         colnames(tdata) = c(xlab, ylab)
         rownames(tdata) = rownames(obj$testres$scores)
         data$tdata = tdata
         if (show.legend == T)
            legend = c('cal', 'test')
      }   
      
      if (show.axes == T)
         show.lines = c(0, 0)      
      else
         show.lines = F

      mdaplotg(data, type = type, main = main, show.labels = show.labels, legend = legend, 
               show.lines = show.lines, xlab = xlab, ylab = ylab, ...)
      }
   else
   {
      stop('Wrong number of components!')
   }   
}  


#' Residuals plot for PCA
#' 
#' @description
#' Shows a plot with Q residuals vs. Hotelling T2 values for selected number of components.
#' 
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param show.limits
#' logical, show or not lines with statistical limits for the residuals
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{pca}} function.
#' 
#' @export
plotResiduals.pca = function(obj, ncomp = NULL, main = NULL, xlab = 'T2',
                             ylab = 'Squared residual distance (Q)', show.labels = F, 
                             show.legend = T, show.limits = T, ...)
{
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
      show.lines = c(obj$T2lim[1, ncomp], obj$Qlim[1, ncomp])
   else
      show.lines = F
   
   if (ncomp > obj$ncomp || ncomp < 1)
      stop('Wrong number of components!')

   cdata = cbind(obj$calres$T2[, ncomp], obj$calres$Q[, ncomp])
   rownames(cdata) = rownames(obj$calres$scores)
   legend = 'cal'   
   data = list(cdata = cdata)

   if (!is.null(obj$cvres))
   {
      cvdata = cbind(obj$cvres$T2[, ncomp], obj$cvres$Q[, ncomp])      
      rownames(cvdata) = rownames(obj$cvres$T2)      
      data$cvdata = cvdata
      legend = c(legend, 'cv')
   }      
   
   if (!is.null(obj$testres))
   {
      tdata = cbind(obj$testres$T2[, ncomp], obj$testres$Q[, ncomp])      
      rownames(tdata) = rownames(obj$testres$scores)      
      data$tdata = tdata
      legend = c(legend, 'test')
   }      

   if (show.legend == F)
      legend = NULL
   
   mdaplotg(data, main = main, xlab = xlab, ylab = ylab,
            show.labels = show.labels, legend = legend, show.lines = show.lines, ...)
}  

#' Loadings plot for PCA
#' 
#' @description
#' Shows a loadings plot for selected components.
#' 
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param comp
#' a value or vector with several values - number of components to show the plot for
#' @param type
#' type of the plot ('b', 'l', 'h')
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{pca}} function.
#' 
#' @export
plotLoadings.pca = function(obj, comp = c(1, 2), type = NULL, main = 'Loadings', xlab = NULL, 
                            ylab = NULL, show.labels = T, show.legend = T,  show.axes = T, ...)
{   
   ncomp = length(comp)
   
   if (ncomp == 2 && (type == 'p' || is.null(type)))
   {
      # scatter plot
      
      data = obj$loadings[, c(comp[1], comp[2])]      
      mdaplot(data, show.labels = show.labels, main = main, xlab = xlab, ylab = ylab, 
              show.lines = c(0, 0), ...)      
   }  
   else if (ncomp < 1 | ncomp > 8 )
   {
      stop ('Number of components must be between 1 and 8!')
   }  
   else
   {
      # loadings vs objects
      
      if (is.null(type))
         type = 'b'
      
      if (is.null(xlab))
         xlab = 'Variables'
      
      if (is.null(ylab))
         ylab = 'Loadings'
      
      data = cbind(1:nrow(obj$loadings), obj$loadings[, comp, drop = F])            
      rownames(data) = rownames(obj$loadings)
      
      if (show.legend == T)
         legend = colnames(data)[-1];

      mdaplotg(data, legend = legend, type = type, show.labels = show.labels, 
               main = main, ylab = ylab, xlab = xlab, ...)
   }   
}

#' Model overview plot for PCA
#' 
#' @description
#' Shows a set of plots (scores, loadings, residuals and explained variance) for PCA model.
#' 
#' @param x
#' a PCA model (object of class \code{pca})
#' @param comp
#' vector with two values - number of components to show the scores and loadings plots for
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{pca}} function.
#' 
#' @export
plot.pca = function(x, comp = c(1, 2), show.labels = F, show.legend = T, ...)
{   
   obj = x
   
   par(mfrow = c(2, 2))
   plotScores(obj, comp = comp, show.labels = show.labels, show.legend = show.legend)
   plotLoadings(obj, comp = comp, show.labels = show.labels, show.legend = show.legend)
   plotResiduals(obj, ncomp = obj$ncomp.selected,  show.labels = show.labels, 
                 show.legend = show.legend, show.limits = T)
   plotCumVariance(obj, show.legend = show.legend)
   par(mfrow = c(1, 1))
}

#' Print method for PCA model object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a PCA model (object of class \code{pca})
#' @param ...
#' other arguments
#' 
#' @export
print.pca = function(x, ...)
{
   obj = x
   
   cat('\nPCA model (class pca)\n')
   
   if (length(obj$info) > 1)
   {
      cat('\nInfo:\n')
      cat(obj$info)      
   }   
   
   cat('\n\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$loadings - matrix with loadings\n')
   cat('$eigenvals - eigenvalues for components\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$center - values for centering data\n')
   cat('$scale - values for scaling data\n')
   cat('$cv - number of segments for cross-validation\n')
   cat('$alpha - significance level for Q residuals\n')
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

#' Summary method for PCA model object
#' 
#' @description
#' Shows some statistics (explained variance, eigenvalues) for the model.
#' 
#' @param object
#' a PCA model (object of class \code{pca})
#' @param ...
#' other arguments
#' 
#' @export
summary.pca = function(object, ...)
{
   obj = object
   
   cat('\nPCA model (class pca) summary\n')

   if (length(obj$info) > 0)
      cat(sprintf('\nInfo:\n%s\n\n', obj$info))
   
   data = cbind(round(obj$eigenvals[1:obj$ncomp], 3), 
                round(obj$calres$expvar, 2),
                round(obj$calres$cumexpvar, 2))
   
   colnames(data) = c('Eigvals', 'Expvar', 'Cumexpvar')
   show(data)
}
