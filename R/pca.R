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
#' method to compute principal components ('svd', 'nipals').
#' @param info
#' a short text line with model description.
#'
#' @details 
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
   # calibrate and cross-validate model  
   model = pca.cal(x, ncomp, center = center, scale = scale, method = method, cv = cv, alpha = alpha, info = info)
   model$call = match.call()   
   
   # apply model to test set if provided
   if (!is.null(x.test))
      model$testres = predict.pca(model, x.test)
   
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
   x = tcrossprod(obj$calres$scores, obj$loadings) + obj$calres$residuals
   
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
      x.new = tcrossprod(scores, loadings)   
      
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

#' Runs one of the selected PCA methods
#' 
#' @param x
#' data matrix 
#' @param ncomp
#' number of components 
#' @param method
#' name of PCA methods ('svd', 'nipals')
#' 
#' @export
pca.run = function(x, ncomp, method) {
   # compute loadings, scores and eigenvalues for data without excluded elements
   if (method == 'svd')
      res = pca.svd(x, ncomp)
   else if(method == 'nipals')
      res = pca.nipals(x, ncomp)
   else
      stop('Wrong value for PCA method!')
   
   res
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
#' algorithm for compiting PC space (only 'svd' and 'nipals' are supported so far)
#' 
#' @return
#' an object with calibrated PCA model
#' 
pca.cal = function(x, ncomp, center, scale, method, cv, alpha, info)
{
   # prepare empty list for model object
   model = list()
   
   # convert data to a matrix 
   x = mda.df2mat(x)
   x.nrows = nrow(x)
   x.ncols = ncol(x)
   
   # check if data has missing values
   if (sum(is.na(x)) > 0)
      stop('Data has missing values, try to fix this using pca.mvreplace.')
   
   # get attributes
   attrs = mda.getattr(x)
   
   # correct maximum number of components
   ncols = x.ncols - length(attrs$exclcols) 
   nrows = x.nrows - length(attrs$exclrows) 
   ncomp = min(ncomp, ncols, nrows - 1)
  
   # prepare data for model calibration and cross-validation
   x.cal = x
   
   # remove excluded rows 
   if (length(attrs$exclrows) > 0)
      x.cal = x.cal[-attrs$exclrows, , drop = F]
   
   # autoscale and save the mean and std values for predictions 
   x.cal = prep.autoscale(x.cal, center = center, scale = scale)
   model$center = attr(x.cal, 'prep:center')
   model$scale = attr(x.cal, 'prep:scale')
   
   # remove excluded columns
   if (length(attrs$exclcols) > 0)
      x.cal = x.cal[, -attrs$exclcols, drop = F]
   
   # compute loadings, scores and eigenvalues for data without excluded elements
   res = pca.run(x.cal, ncomp, method)
   
   # correct loadings for missing columns in x 
   # corresponding rows in loadings will be set to 0 and excluded
   loadings = matrix(0, nrow = x.ncols, ncol = ncomp)
   
   if (length(attrs$exclcols) > 0) {
      loadings[-attrs$exclcols, ] = res$loadings
      loadings = mda.exclrows(loadings, attrs$exclcols)      
   } else {
      loadings = res$loadings
   }

   # set names and attributes for the loadings
   rownames(loadings) = colnames(x)
   colnames(loadings) = paste('Comp', 1:ncol(loadings))
   attr(loadings, 'name') = 'Loadings'
   attr(loadings, 'xaxis.name') = 'Components'
   attr(loadings, 'yaxis.name') = attrs$xaxis.name
   attr(loadings, 'yaxis.values') = attrs$xaxis.values
  
   # add calculated loadings and eigenvalues to the model object
   model$loadings = loadings
   model$eigenvals = res$eigenvals
   model$method = method
   
   # calculate tnorm using data without excluded values
   model$tnorm = sqrt(colSums(res$scores ^ 2)/(nrow(res$scores) - 1));   

   # setup other fields and return the model   
   model$ncomp = ncomp
   model$ncomp.selected = model$ncomp
   model$info = info
   model$alpha = alpha

   # compute statistical limits for the distances
   lim = ldecomp.getResLimits(res$eigenvals, nobj = nrows, ncomp = ncomp, alpha = alpha)
   model$Qlim = lim$Qlim
   model$T2lim = lim$T2lim

   class(model) = "pca"
   
   # get calibration results
   model$calres = predict(model, x, cal = TRUE)
   model$calres$Qlim = model$Qlim
   model$calres$T2lim = model$T2lim
   model$modpower = model$calres$modpower
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = pca.crossval(model, x.cal, cv, center = center, scale = scale)

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
         p = crossprod(E, t) / as.vector(crossprod(t))
         p = p / as.vector(crossprod(p)) ^ 0.5
         t = (E %*% p)/as.vector(crossprod(p))
         th = abs(tau - as.vector(crossprod(t)))
         tau = as.vector(crossprod(t))
      }
      
      E = E - tcrossprod(t, p)
      scores[, i] = t
      loadings[, i] = p
      eigenvals[i] = tau / (nobj - 1)
   }
   
   s = svd(E)
   
   res = list(
      loadings = loadings,
      scores = scores,
      eigenvals = c(eigenvals, (s$d[1:(nvar - ncomp + 1)]^2)/(nrow(x) - 1))
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
         
            # autoscale calibration set
            x.cal = prep.autoscale(x.cal, center = center, scale = scale)
            
            # get loadings
            m = pca.run(x.cal, ncomp, model$method)               
            
            # get scores
            scores = x.val %*% m$loadings
            residuals = x.val - tcrossprod(scores, m$loadings)
             
            # comput distances
            res = ldecomp.getDistances(scores, m$loadings, residuals, model$tnorm, cv = TRUE)
            Q[ind, ] = Q[ind, ] + res$Q
            T2[ind, ] = T2[ind, ] + res$T2
         }
      }  
   }
   
   # prepare results
   Q = Q / nrep
   T2 = T2 / nrep
   rownames(Q) = rownames(T2) = rownames(x)
   colnames(Q) = colnames(T2) = colnames(model$calres$scores)

   # add attributes
   attr(Q, 'name') = 'Squared residual distance (Q)'
   attr(Q, 'xaxis.name') = attr(model$calres$scores, 'xaxis.name')
   attr(Q, 'yaxis.name') = attr(model$calres$scores, 'yaxis.name')
   attr(Q, 'yaxis.values') = attr(model$calres$scores, 'yaxis.values')
   
   attr(T2, 'name') = 'T2 residuals'
   attr(T2, 'xaxis.name') = attr(model$calres$scores, 'xaxis.name')
   attr(T2, 'yaxis.name') = attr(model$calres$scores, 'yaxis.name')
   attr(T2, 'yaxis.values') = attr(model$calres$scores, 'yaxis.values')
   
   # compute variance
   var = ldecomp.getVariances(Q, model$calres$totvar)
   
   # in CV results there are no scores nor residuals, only residual distances and variances
   res = pcares(NULL, NULL, model$ncomp.selected, T2, Q, var$expvar, var$cumexpvar)
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
#' @param ...
#' other arguments
#' 
#' @return
#' PCA results (an object of class \code{pcares})
#'  
#' @export
predict.pca = function(object, x, cal = FALSE, ...)
{
   # get attributes
   attrs = mda.getattr(x)
 
   # compute scores
   x = prep.autoscale(x, center = object$center, scale = object$scale)
   scores = x %*% object$loadings
   
   # set names and attributes
   rownames(scores) = rownames(x)
   colnames(scores) = colnames(object$loadings)
   scores = mda.setattr(scores, attrs, type = 'row')
   attr(scores, 'name') = 'Scores'
   attr(scores, 'xaxis.name') = 'Components'
   
   # calculate residuals and set all attributes from x
   residuals = x - tcrossprod(scores, object$loadings)
   residuals = mda.setattr(residuals, attrs)
   attr(scores, 'name') = 'Residuals'
   
   # calculate residual distances
   dist = ldecomp.getDistances(scores, object$loadings, residuals, object$tnorm, cal = cal) 
   
   # compute total variance 
   if (length(attrs$exclrows) > 0)
      x = x[-attrs$exclrows, , drop = F]
   
   if (length(attrs$exclcols) > 0)
      x = x[, -attrs$exclcols, drop = F]
   
   totvar = sum(x^2)
   
   # calculate explained variance
   var = ldecomp.getVariances(dist$Q, totvar) 
   
   # create and return the results object
   res = pcares(scores, residuals, object$ncomp.selected, dist$T2, dist$Q, var$expvar, var$cumexpvar)
   res$modpower = dist$modpower
   res$Qlim = object$Qlim
   res$T2lim = object$T2lim
   res$totvar = totvar
   
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
                            show.legend = T, ...)
{

   data = list()
   data$cal = obj$calres[[variance]]
   
   if (!is.null(obj$cvres))
      data$cv = obj$cvres[[variance]]
   
   if (!is.null(obj$testres))
      data$test = obj$testres[[variance]]
   
   mdaplotg(data, main = main, xlab = xlab, ylab = ylab, show.legend = show.legend, type = type, ...)   
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
                          ylab = NULL, show.labels = F, show.legend = NULL, cgroup = NULL,
                          show.axes = TRUE, ...)
{
   ncomp = length(comp)

   if (type != 'p') {
      plotScores(obj$calres, comp = comp, type = type, main = main, xlab = xlab, ylab = ylab, 
                 show.labels = show.labels, show.legend = show.legend, show.axes = show.axes, ...)
   } else {
      data = list() 
      data$cal = mda.subset(obj$calres$scores, select = comp)
      colnames(data$cal) = paste('Comp ', comp, ' (', round(obj$calres$expvar[comp], 2) , '%)', sep = '')
      
      if (!is.null(obj$testres)) {
         data$test = mda.subset(obj$testres$scores, select = comp)
         colnames(data$test) = paste('Comp ', comp, ' (', round(obj$calres$expvar[comp], 2) , '%)', sep = '')
      }
      
      if (show.axes == TRUE) {
         if (ncomp == 1)
            show.lines = c(NA, 0)
         else
            show.lines = c(0, 0)
      } else {
         show.lines = FALSE
      }
      
      if (is.null(show.legend)) {
         if (length(data) > 0)
            show.legend = TRUE
         else
            show.legend = FALSE 
      }
      
      if (length(data) == 1)
         mdaplot(data[[1]], type = type, main = main, show.labels = show.labels, show.lines = show.lines, 
                 xlab = xlab, ylab = ylab, cgroup = cgroup, ...)
      else
         mdaplotg(data, type = type, main = main, show.labels = show.labels, show.lines = show.lines, 
                  xlab = xlab, ylab = ylab, show.legend = show.legend, ...)
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
                             show.legend = T, show.limits = T, cgroup = NULL, ...)
{
   if (is.null(main)) {
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

   data = list()
   
   data$cal = mda.cbind(mda.subset(obj$calres$T2, select = ncomp), mda.subset(obj$calres$Q, select = ncomp))
   colnames(data$cal) = c(xlab, ylab)
   
   
   if (!is.null(obj$cvres)) {
      data$cv = mda.cbind(mda.subset(obj$cvres$T2, select = ncomp), mda.subset(obj$cvres$Q, select = ncomp))
      colnames(data$cv) = c(xlab, ylab)
   }      
   
   if (!is.null(obj$testres)) {
      data$test = mda.cbind(mda.subset(obj$testres$T2, select = ncomp), mda.subset(obj$testres$Q, select = ncomp))
      colnames(data$test) = c(xlab, ylab)
   }      

   if (length(data) == 1) {
      mdaplot(data[[1]], main = main, xlab = xlab, ylab = ylab, cgroup = cgroup,
               show.labels = show.labels, show.lines = show.lines, ...)
   } else {
      mdaplotg(data, main = main, xlab = xlab, ylab = ylab,
              show.labels = show.labels, show.legend = show.legend, show.lines = show.lines, ...)
   }
      
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
                            ylab = NULL, show.labels = NULL, show.legend = TRUE,  show.axes = TRUE, ...)
{   
   if (max(comp) > obj$ncomp || min(comp) < 1)
      stop('Wrong number of components!')
   
   if (is.null(type)) {
      if (length(comp) == 2)
         type = 'p'
      else
         type = 'l'
   }
   
   data = mda.subset(obj$loadings, select = comp)      
   
   if (type == 'p') {
      if (show.axes == TRUE) {
         if (length(comp) > 1)
            show.lines = c(0, 0)
         else
            show.lines = c(NA, 0)
      } else {
         show.lines = FALSE
      }
      
      if (is.null(show.labels))
         show.labels = TRUE
      
      mdaplot(data, type = type, show.labels = show.labels, show.lines = show.lines, 
               main = main, ylab = ylab, xlab = xlab, ...)
   } else {
      if (is.null(show.legend))
         show.legend = TRUE 
      
      if (is.null(show.labels))
         show.labels = FALSE
      
      if (show.axes == TRUE && type != 'h')
         show.lines = c(NA, 0)
      else
         show.lines = FALSE
      
      mdaplotg(mda.t(data), show.legend = show.legend, type = type, show.labels = show.labels, 
               show.lines = show.lines, main = main, ylab = ylab, xlab = xlab, ...)
   }  
}


#' PCA biplot
#' 
#' @export
plotBiplot.pca = function(obj, comp = c(1, 2), pch = c(16, NA), col = mdaplot.getColors(2), main = 'Biplot', 
                          lty = 1, lwd = 1, show.labels = FALSE, show.axes = TRUE, show.excluded = FALSE,
                          lab.col = c('#90A0D0', '#D09090'), ...) {
   
   if (length(comp) != 2)
      stop('Biplot can be made only for two principal components!')
   
   if (show.axes == TRUE)
      show.lines = c(0, 0)
   else
      show.lines = FALSE
   
   loadings = mda.subset(obj$loadings, select = comp)
   scores = mda.subset(obj$calres$scores, select = comp)
   attrs = mda.getattr(scores)
   
   loadsScaleFactor = sqrt(max(rowSums(loadings^2)))
   
   scoresScaleFactor = max(abs(scores))
   scores = (scores / scoresScaleFactor) * loadsScaleFactor
   scores = mda.setattr(scores, attrs)
   
   
   colnames(scores) = paste('Comp ', comp, ' (', round(obj$calres$expvar[comp], 2) , '%)', sep = '')
   
   mdaplotg(list(scores = scores, loadings = loadings), type = 'p', pch = pch, 
            show.legend = FALSE, show.labels = show.labels, lab.col = lab.col,
            main = main, col = col, show.lines = show.lines, show.excluded = show.excluded, ...)   
   
   if (show.excluded == TRUE && length(attr(loadings, 'exclrows')) > 0)
      loadings = loadings[-attr(loadings, 'exclrows'), , drop = F]
   
   segments(0, 0, loadings[, 1], loadings[, 2], col = col[2], lty = lty, lwd = lwd)
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
plot.pca = function(x, comp = c(1, 2), show.labels = FALSE, show.legend = TRUE, ...)
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
