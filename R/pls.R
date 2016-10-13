#' Partial Least Squares regression
#'
#' @description 
#' \code{pls} is used to calibrate, validate and use of partial least squares (PLS) 
#' regression model.
#' 
#' @param x
#' matrix with predictors.
#' @param y  
#' matrix with responses.
#' @param ncomp
#' maximum number of components to calculate.
#' @param center   
#' logical, center or not predictors and response values.
#' @param scale   
#' logical, scale (standardize) or not predictors and response values.
#' @param cv  
#' number of segments for cross-validation (if cv = 1, full cross-validation will be used).
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param x.test   
#' matrix with predictors for test set.
#' @param y.test  
#' matrix with responses for test set.
#' @param method
#' algorithm for computing PLS model (only 'simpls' is supported so far)
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
#' criterion for selecting optimal number of components (\code{'min'} for 
#' first local minimum of RMSECV and \code{'wold'} for Wold's rule.)
#' 
#' @return 
#' Returns an object of \code{pls} class with following fields:
#' \item{ncomp }{number of components included to the model.} 
#' \item{ncomp.selected }{selected (optimal) number of components.} 
#' \item{xloadings }{matrix with loading values for x decomposition.} 
#' \item{yloadings }{matrix with loading values for y decomposition.} 
#' \item{weights }{matrix with PLS weights.} 
#' \item{selratio }{array with selectivity ratio values.} 
#' \item{vipscores }{matrix with VIP scores values.} 
#' \item{coeffs }{object of class \code{\link{regcoeffs}} with regression coefficients calculated for each component.}   
#' \item{info }{information about the model, provided by user when build the model.} 
#' \item{calres }{an object of class \code{\link{plsres}} with PLS results for a calibration data.} 
#' \item{testres }{an object of class \code{\link{plsres}} with PLS results for a test data, if it was provided.} 
#' \item{cvres }{an object of class \code{\link{plsres}} with PLS results for cross-validation, if this option was chosen.} 
#'
#' @details 
#' So far only SIMPLS method [1] is available, more coming soon. Implementation works both with one
#' and multiple response variables.
#'
#' Like in \code{\link{pca}}, \code{pls} uses number of components (\code{ncomp}) as a minimum of 
#' number of objects - 1, number of x variables and the default or provided value. Regression 
#' coefficients, predictions and other results are calculated for each set of components from 1
#' to \code{ncomp}: 1, 1:2, 1:3, etc. The optimal number of components, (\code{ncomp.selected}), 
#' is found using Wold's R criterion, but can be adjusted by user using function
#' (\code{\link{selectCompNum.pls}}). The selected optimal number of components is used for all 
#' default operations - predictions, plots, etc. 
#'
#' Selectivity ratio [2] and VIP scores [3] are calculated for any PLS model authomatically, however
#' while selectivity ratio values are calculated for all computed components, the VIP scores are 
#' computed only for selected components (to save calculation time) and recalculated every time when 
#' \code{selectCompNum()} is called for the model. 
#' 
#' Calculation of confidence intervals and p-values for regression coefficients are available
#' only by jack-knifing so far. See help for \code{\link{regcoeffs}} objects for details. 
#'
#' @references     
#' 1. S. de Jong, Chemometrics and Intelligent Laboratory Systems 18 (1993) 251-263.
#' 2. Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), 35-48.
#' 3. Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), 103-112.
#'
#' @seealso    
#' Methods for \code{pls} objects:
#' \tabular{ll}{
#'  \code{print} \tab prints information about a \code{pls} object.\cr
#'  \code{\link{summary.pls}} \tab shows performance statistics for the model.\cr
#'  \code{\link{plot.pls}} \tab shows plot overview of the model.\cr
#'  \code{\link{pls.simpls}} \tab implementation of SIMPLS algorithm.\cr
#'  \code{\link{predict.pls}} \tab applies PLS model to a new data.\cr
#'  \code{\link{selectCompNum.pls}} \tab set number of optimal components in the model.\cr
#'  \code{\link{plotPredictions.pls}} \tab shows predicted vs. measured plot.\cr
#'  \code{\link{plotRegcoeffs.pls}} \tab shows regression coefficients plot.\cr      
#'  \code{\link{plotXScores.pls}} \tab shows scores plot for x decomposition.\cr
#'  \code{\link{plotXYScores.pls}} \tab shows scores plot for x and y decomposition.\cr
#'  \code{\link{plotXLoadings.pls}} \tab shows loadings plot for x decomposition.\cr
#'  \code{\link{plotXYLoadings.pls}} \tab shows loadings plot for x and y decomposition.\cr
#'  \code{\link{plotRMSE.pls}} \tab shows RMSE plot.\cr
#'  \code{\link{plotXVariance.pls}} \tab shows explained variance plot for x decomposition.\cr
#'  \code{\link{plotYVariance.pls}} \tab shows explained variance plot for y decomposition.\cr
#'  \code{\link{plotXCumVariance.pls}} \tab shows cumulative explained variance plot for y 
#'  decomposition.\cr
#'  \code{\link{plotYCumVariance.pls}} \tab shows cumulative explained variance plot for y 
#'  decomposition.\cr
#'  \code{\link{plotXResiduals.pls}} \tab shows T2 vs. Q plot for x decomposition.\cr
#'  \code{\link{plotYResiduals.pls}} \tab shows residuals plot for y values.\cr
#'  \code{\link{plotSelectivityRatio.pls}} \tab shows plot with selectivity ratio values.\cr
#'  \code{\link{plotVIPScores.pls}} \tab shows plot with VIP scores values.\cr
#'  \code{\link{getSelectivityRatio.pls}} \tab returns vector with selectivity ratio values.\cr
#'  \code{\link{getVIPScores.pls}} \tab returns vector with VIP scores values.\cr
#'  \code{\link{getRegcoeffs.pls}} \tab returns matrix with regression coefficients.\cr
#' }
#'      
#' Most of the methods for plotting data (except loadings and regression coefficients) are also 
#' available for PLS results 
#' (\code{\link{plsres}}) objects. There is also a randomization test for PLS-regression 
#' (\code{\link{randtest}}).   
#'
#' @author 
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @examples
#' ### Examples of using PLS model class
#' library(mdatools)   
#'   
#' ## 1. Make a PLS model for concentration of first component 
#' ## using full-cross validation and automatic detection of 
#' ## optimal number of components and show an overview
#' 
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#' 
#' model = pls(x, y, ncomp = 8, cv = 1)
#' summary(model)
#' plot(model)
#' 
#' ## 2. Make a PLS model for concentration of first component 
#' ## using test set and 10 segment cross-validation and show overview
#' 
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#' x.t = simdata$spectra.t
#' y.t = simdata$conc.t[, 1]
#' 
#' model = pls(x, y, ncomp = 8, cv = 10, x.test = x.t, y.test = y.t)
#' model = selectCompNum(model, 2)
#' summary(model)
#' plot(model)
#' 
#' ## 3. Make a PLS model for concentration of first component 
#' ## using only test set validation and show overview
#' 
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#' x.t = simdata$spectra.t
#' y.t = simdata$conc.t[, 1]
#' 
#' model = pls(x, y, ncomp = 6, x.test = x.t, y.test = y.t)
#' model = selectCompNum(model, 2)
#' summary(model)
#' plot(model)
#' 
#' ## 4. Show variance and error plots for a PLS model
#' par(mfrow = c(2, 2))
#' plotXCumVariance(model, type = 'h')
#' plotYCumVariance(model, type = 'b', show.labels = TRUE, legend.position = 'bottomright')
#' plotRMSE(model)
#' plotRMSE(model, type = 'h', show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' ## 5. Show scores plots for a PLS model
#' par(mfrow = c(2, 2))
#' plotXScores(model)
#' plotXScores(model, comp = c(1, 3), show.labels = TRUE)
#' plotXYScores(model)
#' plotXYScores(model, comp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' ## 6. Show loadings and coefficients plots for a PLS model
#' par(mfrow = c(2, 2))
#' plotXLoadings(model)
#' plotXLoadings(model, comp = c(1, 2), type = 'l')
#' plotXYLoadings(model, comp = c(1, 2), legend.position = 'topleft')
#' plotRegcoeffs(model)
#' par(mfrow = c(1, 1))
#' 
#' ## 7. Show predictions and residuals plots for a PLS model
#' par(mfrow = c(2, 2))
#' plotXResiduals(model, show.label = TRUE)
#' plotYResiduals(model, show.label = TRUE)
#' plotPredictions(model)
#' plotPredictions(model, ncomp = 4, xlab = 'C, reference', ylab = 'C, predictions')
#' par(mfrow = c(1, 1))
#' 
#' ## 8. Selectivity ratio and VIP scores plots
#' par(mfrow = c(2, 2))
#' plotSelectivityRatio(model)
#' plotSelectivityRatio(model, ncomp = 1)
#' par(mfrow = c(1, 1))
#' 
#' ## 9. Variable selection with selectivity ratio
#' selratio = getSelectivityRatio(model)
#' selvar = !(selratio < 8)
#' 
#' xsel = x[, selvar]
#' modelsel = pls(xsel, y, ncomp = 6, cv = 1)
#' modelsel = selectCompNum(modelsel, 3)
#' 
#' summary(model)
#' summary(modelsel)
#' 
#' ## 10. Calculate average spectrum and show the selected variables
#' i = 1:ncol(x)
#' ms = apply(x, 2, mean)
#' 
#' par(mfrow = c(2, 2))
#' 
#' plot(i, ms, type = 'p', pch = 16, col = 'red', main = 'Original variables')
#' plotPredictions(model)
#' 
#' plot(i, ms, type = 'p', pch = 16, col = 'lightgray', main = 'Selected variables')
#' points(i[selvar], ms[selvar], col = 'red', pch = 16)
#' plotPredictions(modelsel)
#' 
#' par(mfrow = c(1, 1))
#' 
#' @export   
pls = function(x, y, ncomp = 15, center = T, scale = F, cv = NULL, exclcols = NULL, exclrows = NULL,
               x.test = NULL, y.test = NULL, method = 'simpls', alpha = 0.05, coeffs.ci = NULL, 
               coeffs.alpha = 0.1, info = '', light = F, 
               ncomp.selcrit = 'min') {
   
   # build a model and apply to calibration set
   model = pls.cal(x, y, ncomp, center = center, scale = scale, method = method, coeffs.ci = coeffs.ci,
                   coeffs.alpha = coeffs.alpha, info = info, light = light, alpha = alpha, cv = cv, 
                   exclcols = exclcols, exclrows = exclrows, ncomp.selcrit = ncomp.selcrit)
   
   # do test set validation if provided
   if (!is.null(x.test) && !is.null(y.test)){
      model$testres = predict.pls(model, x.test, y.test)
   }

   # select optimal number of components
   model = selectCompNum(model)
   
   model
}

#' PLS model calibration
#' 
#' @description
#' Calibrates (builds) a PLS model for given data and parameters
#' 
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with y values (responses)
#' @param ncomp
#' number of components to calculate
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#' @param method
#' algorithm for computing PLS model (only 'simpls' is supported so far)
#' @param cv
#' logical, does calibration for cross-validation or not
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
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param ncomp.selcrit  
#' criterion for selecting optimal number of components (\code{'min'} for 
#' first local minimum of RMSECV and \code{'wold'} for Wold's rule.)
#' 
#' @return model
#' an object with calibrated PLS model
#' 
pls.cal = function(x, y, ncomp, center, scale, method, cv, alpha, coeffs.ci, coeffs.alpha, info, light,
                   exclcols = NULL, exclrows = NULL, ncomp.selcrit) {   
   # prepare empty list for model object
   model = list()
  
   # get attributes
   x.attrs = mda.getattr(x)
   y.attrs = mda.getattr(y)
   
   # if y is a vector convert it to a matrix
   if (is.null(dim(y)))
      y = matrix(y, ncol = 1)
   
   # check dimensions
   if (nrow(x) != nrow(y))
      stop('Number of rows for predictos and responses should be the same!')
   
   # exclude rows and columns if necessary
   if (length(exclcols) > 0) {
      x = mda.exclcols(x, exclcols)
      x.attrs$exclcols = attr(x, 'exclcols')
   }
   
   if (length(exclrows) > 0) {
      x = mda.exclrows(x, exclrows)
      y = mda.exclrows(y, exclrows)
      x.attrs$exclrows = attr(x, 'exclrows')
      y.attrs$exclrows = attr(y, 'exclrows')
   }
   
   # convert data to a matrix 
   x = mda.df2mat(x)
   y = mda.df2mat(y)
   x.nrows = nrow(x)
   x.ncols = ncol(x)
   y.nrows = nrow(y)
   y.ncols = ncol(y)
   
   # check if data has missing values
   if (any(is.na(x)))
      stop('Predictors have missing values, try to fix this using pca.mvreplace.')
   if (any(is.na(y)))
      stop('Responses have missing values, try to fix this using pca.mvreplace.')

   # set column names for predictors if missing   
   if (is.null(colnames(y)))
      colnames(y) = paste('y', 1:ncol(y), sep = '')
   
   # correct x-axis name
   if (is.null(x.attrs$xaxis.name))
      x.attrs$xaxis.name = 'Variables'
   if (is.null(x.attrs$yaxis.name))
      x.attrs$yaxis.name = 'Objects'
   
   # correct maximum number of components
   ncols = x.ncols - length(x.attrs$exclcols) 
   nrows = x.nrows - length(x.attrs$exclrows) 
   ncomp = min(ncomp, ncols, nrows - 1)
   
   # prepare data for model calibration and cross-validation
   x.cal = x
   y.cal = y
  
   # check excluded rows
   #if (length(x.attrs$exclrows) != length(y.attrs$exclrows) || any(x.attrs$exclrows != y.attrs$exclrows))
   #   stop('Excluded rows in response and predictors matrices are not the same!')
   
   # remove excluded rows 
   if (length(x.attrs$exclrows) > 0) {
      x.cal = x.cal[-x.attrs$exclrows, , drop = F]
      y.cal = y.cal[-x.attrs$exclrows, , drop = F]
   }
   
   # autoscale and save the mean and std values for predictions 
   x.cal = prep.autoscale(x.cal, center = center, scale = scale)
   model$xcenter = attr(x.cal, 'prep:center')
   model$xscale = attr(x.cal, 'prep:scale')
   
   y.cal = prep.autoscale(y.cal, center = center, scale = scale)
   model$ycenter = attr(y.cal, 'prep:center')
   model$yscale = attr(y.cal, 'prep:scale')
   
   # remove excluded columns
   if (length(x.attrs$exclcols) > 0)
      x.cal = x.cal[, -x.attrs$exclcols, drop = F]
   if (length(y.attrs$exclcols) > 0)
      y.cal = y.cal[, -y.attrs$exclcols, drop = F]
   
   # set LOO cross-validation if jack.knife is selected
   jack.knife = F
   if (!is.null(coeffs.ci) && coeffs.ci == 'jk') {
      jack.knife = T
      if (is.null(cv))    
         cv = 1      
   }   
   
   # correct maximum number of components
   if (!is.null(cv)) {
      if (!is.numeric(cv))
         nseg = cv[[2]]
      else
         nseg = cv
      
      if (nseg == 1)
         nobj.cv = 1
      else
         nobj.cv = ceiling(nrow(x)/nseg)  
   } else {   
      nobj.cv = 0
   }
   ncomp = min(ncol(x), nrow(x) - 1 - nobj.cv, ncomp)

   
   # compute model and redefine ncomp
   res = pls.run(x.cal, y.cal, method = method, ncomp = ncomp, cv = FALSE)
   ncomp = res$ncomp
   
   # correct results related to predictors for missing columns in x 
   # corresponding rows will be set to 0 and excluded
   xloadings = matrix(0, nrow = x.ncols, ncol = ncomp)
   weights = matrix(0, nrow = x.ncols, ncol = ncomp)
   coeffs = array(0, dim = c(x.ncols, ncomp, ncol(y.cal)))
   if (length(x.attrs$exclcols) > 0) {
      xloadings[-x.attrs$exclcols, ] = res$xloadings
      xloadings = mda.exclrows(xloadings, x.attrs$exclcols)      
      weights[-x.attrs$exclcols, ] = res$weights
      weights = mda.exclrows(weights, x.attrs$exclcols)      
      coeffs[-x.attrs$exclcols, , ] = res$coeffs
      coeffs = mda.exclrows(coeffs, x.attrs$exclcols)      
   } else {
      xloadings = res$xloadings
      weights = res$weights
      coeffs = res$coeffs
   }
   
   # set names and attributes 
   rownames(xloadings) = rownames(weights) = colnames(x)
   colnames(xloadings) = colnames(weights) = paste('Comp', 1:ncomp)
   attr(xloadings, 'name') = 'X loadings'
   attr(xloadings, 'xaxis.name') = attr(weights, 'xaxis.name') = attr(coeffs, 'xaxis.name') = 'Components'
   attr(xloadings, 'yaxis.name') = attr(weights, 'yaxis.name') = attr(coeffs, 'yaxis.name') = x.attrs$xaxis.name
   attr(xloadings, 'yaxis.values') = attr(weights, 'yaxis.values') = attr(coeffs, 'yaxis.values') = x.attrs$xaxis.values

   # do the same for response related results
   yloadings = matrix(0, nrow = y.ncols, ncol = ncomp)
   if (length(y.attrs$exclcols) > 0) {
      yloadings[-y.attrs$exclcols, ] = res$xloadings
      yloadings = mda.exclrows(xloadings, y.attrs$exclcols)      
   } else {
      yloadings = res$yloadings
   }
   
   # set names and attributes 
   dimnames(coeffs) = list(colnames(x), colnames(xloadings), colnames(y))
   rownames(yloadings) = colnames(y)
   colnames(yloadings) = colnames(xloadings)
   attr(yloadings, 'name') = 'Y loadings'
   attr(yloadings, 'xaxis.name') = 'Components'
   attr(yloadings, 'yaxis.name') = y.attrs$xaxis.name
   attr(yloadings, 'yaxis.values') = y.attrs$xaxis.values
  
   # set up model parameters
   model$xloadings = xloadings
   model$yloadings = yloadings
   model$weights = weights
   model$coeffs = regcoeffs(coeffs)
   model$method = method
   model$xtnorm = sqrt(colSums(res$xscores ^ 2)/(nrow(res$xscores) - 1))   
   model$ytnorm = sqrt(colSums(res$yscores ^ 2)/(nrow(res$yscores) - 1))   
   model$ncomp = ncomp
   model$alpha = alpha   
   model$light = light
   model$info = info
   model$exclrows = x.attrs$exclrows
   model$exclcols = x.attrs$exclcols
   model$cv = cv
   model$ncomp.selcrit = ncomp.selcrit
   
   model$call = match.call()
   class(model) = "pls"
   
   # do predictions for calibration set
   model$calres = predict.pls(model, x, y)
   
   
   # do cross-validation if needed
   if (!is.null(cv)) {   
      res = pls.crossval(model, x, y, cv, center = center, scale = scale, method = method, jack.knife = jack.knife)    
      if (jack.knife == T) {   
         model$coeffs = regcoeffs(model$coeffs$values, res$jkcoeffs, coeffs.alpha)
         res[['jkcoeffs']] = NULL
         model$cvres = res
      } else {
         model$cvres = res;
      }   
   }
   
   if (!light) {  
      # we calculate both for data without excluded rows and columns
      model$selratio = pls.calculateSelectivityRatio(model, x.cal)
   }
   
   model
}

#' Runs selected PLS algorithm
#' 
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param y
#' a matrix with y values (responses from calibration set)
#' @param ncomp
#' how many components to compute
#' @param method
#' algorithm for computing PLS model
#' @param cv
#' logical, is this for CV or not
#' 
#' @export
pls.run = function(x, y, ncomp, method, cv) {
   if (method == 'simpls')
      res = pls.simpls(x, y, ncomp, cv = cv)
   else
      stop('Method with this name is not supported!')
}

#' SIMPLS algorithm
#' 
#' @description
#' SIMPLS algorithm for calibration of PLS model
#' 
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with y values (responses)
#' @param ncomp
#' number of components to calculate
#' @param cv
#' logical, is model calibrated during cross-validation or not
#' 
#' @return 
#' a list with computed regression coefficients, loadings and scores for x and y matrices,
#' and weights.
#' 
#' @references
#' [1]. S. de Jong. SIMPLS: An Alternative approach to partial least squares regression. 
#' Chemometrics and Intelligent Laboratory Systems, 18, 1993 (251-263).
#' 
pls.simpls = function(x, y, ncomp, cv = FALSE) {
   x = as.matrix(x)
   y = as.matrix(y)
   
   # get names for objects, variables and components
   objnames = rownames(x);
   prednames = colnames(x);
   respnames = colnames(y);
   
   nobj = nrow(x)
   npred = ncol(x)
   nresp = ncol(y)
   
   # initial estimation
   A = t(x) %*% y
   M = t(x) %*% x
   C = diag(npred)
   
   # prepare space for results
   B = array(0, dim = c(npred, ncomp, nresp))
   W = matrix(0, nrow = npred, ncol = ncomp)
   P = matrix(0, nrow = npred, ncol = ncomp)
   Q = matrix(0, nrow = nresp, ncol = ncomp)
   
   # loop for each components
   for (n in 1:ncomp) {
      # get the dominate eigenvector of A'A
      e = eigen(crossprod(A))
      q = e$vectors[1:nresp]
      
      # calculate and store weights
      w = A %*% q
      c = crossprod(w, (M %*% w))
      w = w/sqrt(as.numeric(c))
      W[, n] = w
      
      # calculate and store x loadings
      p = M %*% w
      P[, n] = p
      
      # calculate and store y loadings
      q = crossprod(A, w)
      Q[, n] = q
      
      v = C %*% p
      v = v/sqrt(as.numeric(crossprod(v)))
      
      # calculate and store regression coefficients
      B[, n, ] = tcrossprod(W[, 1:n, drop = FALSE], Q[, 1:n, drop = FALSE])
      
      # recalculate matrices for the next compnonent
      C = C - tcrossprod(v)
      M = M - tcrossprod(p)
      A = C %*% A      
      
      if (cv == FALSE && e$value < 10^-12) {
         # stop cycle is egienvalue is almost zero
         break
      }
   }
   
   # truncate results if n is smaller than ncomp
   B = B[, 1:n, , drop = F]
   W = W[, 1:n, drop = F]
   P = P[, 1:n, drop = F]
   Q = Q[, 1:n, drop = F]
   
   # calculate x and y scores
   U = y %*% Q 
   TT = x %*% (W %*% solve(t(P) %*% W))  
   
   res = list(
      coeffs = B,
      weights = W,
      xloadings = P,
      xscores = TT,
      yloadings = Q,
      yscores = U,
      ncomp = n
   )  
   
   res
}  

#' Cross-validation of a PLS model
#' 
#' @description
#' Does the cross-validation of a PLS model
#' 
#' @param model
#' a PLS model (object of class \code{pls})
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param y
#' a matrix with y values (responses from calibration set)
#' @param cv
#' number of segments (if cv = 1, full cross-validation will be used)
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#' @param method
#' algorithm for computing PLS model
#' @param jack.knife
#' logical, do jack-knifing or not
#' 
#' @return
#' object of class \code{plsres} with results of cross-validation
#'  
pls.crossval = function(model, x, y, cv, center, scale, method, jack.knife = T) {
   # get attributes
   x.attrs = mda.getattr(x)
   y.attrs = mda.getattr(y)

   # remove excluded rows 
   if (length(x.attrs$exclrows) > 0){
      x = x[-x.attrs$exclrows, , drop = F]
      y = y[-x.attrs$exclrows, , drop = F]
   }
   
   # remove excluded columns 
   if (length(x.attrs$exclcols) > 0)
      x = x[, -x.attrs$exclcols, drop = F]
   if (length(y.attrs$exclcols) > 0)
      y = y[, -y.attrs$exclcols, drop = F]
   
   ncomp = model$ncomp
   nobj = nrow(x)
   nvar = ncol(x)
   nresp = ncol(y)
   
   # get matrix with indices for cv segments
   idx = crossval(nobj, cv)
   seglen = ncol(idx);
   nseg = nrow(idx);
   nrep = dim(idx)[3]

   yp.cv = array(0, dim = c(nobj, ncomp, nresp))
   Qx = matrix(0, ncol = ncomp, nrow = nobj)   
   T2x = matrix(0, ncol = ncomp, nrow = nobj)   
   Qy = matrix(0, ncol = ncomp, nrow = nobj)   
   T2y = matrix(0, ncol = ncomp, nrow = nobj)   
   jkcoeffs = array(0, dim = c(nvar, ncomp, nresp, nseg))
   
   # loop over segments and repetitions
   for (iRep in 1:nrep) {   
      for (iSeg in 1:nseg) {         
         ind = na.exclude(idx[iSeg, , iRep])
         
         if (length(ind) > 0) {   
            xc = x[-ind, , drop = F]
            yc = y[-ind, , drop = F]
            xt = x[ind, , drop = F]
            yt = y[ind, , drop = F]
           
            # autoscale calibration set
            xc = prep.autoscale(xc, center = center, scale = scale)
            yc = prep.autoscale(yc, center = center, scale = scale)
            
            # create a model            
            m = pls.run(xc, yc, ncomp, method = method, cv = TRUE)            
            m$xcenter = attr(xc, 'prep:center')
            m$ycenter = attr(yc, 'prep:center')
            m$xscale = attr(xc, 'prep:scale')
            m$yscale = attr(yc, 'prep:scale')
            
            # autoscale test set
            xt = prep.autoscale(xt, center = m$xcenter, scale = m$xscale)
            yt = prep.autoscale(yt, center = m$ycenter, scale = m$yscale)

            # get scores
            xscores = xt %*% (m$weights %*% solve(crossprod(m$xloadings, m$weights)))  
            yscores = as.matrix(yt) %*% m$yloadings   
            
            # make predictions 
            yp = apply(m$coeffs, 3, function(x, y)(y %*% x), xt)
            dim(yp) = c(nrow(xt), ncomp, ncol(yt))
            
            # get residuals
            xresiduals = xt - tcrossprod(xscores, m$xloadings)
            yresiduals = yt - yp[, ncol(yp), ]

            # unscale predicted y values
            if (scale == TRUE)
               yp = sweep(yp, 3, m$yscale, '*')

            # uncenter predicted y values
            if (center == TRUE)
               yp = sweep(yp, 3, m$ycenter, '+')
            
            # get distances
            xdist = ldecomp.getDistances(scores = xscores, loadings = m$xloadings, residuals = xresiduals, 
                                         tnorm = model$xtnorm)
            ydist = ldecomp.getDistances(scores = xscores, loadings = m$yloadings, residuals = yresiduals, 
                                         tnorm = model$ytnorm)
            
            # correct dimenstion for reg coeffs for JK
            dim(m$coeffs) = c(dim(m$coeffs), 1)
          
            # save results
            yp.cv[ind, , ] = yp.cv[ind, , , drop = F]  + yp
            Qx[ind, ]  = Qx[ind, , drop = F] + xdist$Q        
            T2x[ind, ]  = T2x[ind, , drop = F] + xdist$T2        
            Qy[ind, ]  = Qy[ind, , drop = F] + ydist$Q
            T2y[ind, ]  = T2y[ind, , drop = F] + ydist$T2       
            jkcoeffs[, , , iSeg] = jkcoeffs[, , , iSeg, drop = F] + m$coeffs
         }   
      }      
   }  
   
   # average results over repetitions
   yp.cv= yp.cv / nrep
   Qx = Qx / nrep
   T2x = T2x / nrep
   Qy = Qy / nrep
   T2y = T2y / nrep
   jkcoeffs = jkcoeffs / nrep
  
   # set up names
   dimnames(jkcoeffs) = list(colnames(x), colnames(model$coeffs$values), colnames(model$calres$y.ref), 1:nseg)
   dimnames(yp.cv) = list(rownames(x), colnames(model$coeffs$values), colnames(model$calres$y.ref))         
   
   # compute variance
   varx = ldecomp.getVariances(Qx, model$calres$xdecomp$totvar)
   vary = ldecomp.getVariances(Qy, model$calres$ydecomp$totvar)

   # get rid of some of the attributed
   x.attrs$exclrows = NULL
   x.attrs$exclcols = NULL
   y.attrs$exclrows = NULL
   y.attrs$exclcols = NULL
   
   # make pls results and return
   res = plsres(yp.cv, y.ref = y, ncomp.selected = model$ncomp.selected,
                xdecomp = ldecomp(ncomp.selected = model$ncomp.selected, dist = list(Q = Qx, T2 = T2x), 
                                  var = varx, loadings = model$xloadings, attrs = x.attrs),
                ydecomp = ldecomp(ncomp.selected = model$ncomp.selected, dist = list(Q = Qy, T2 = T2y), 
                                  var = vary, loadings = model$yloadings, attrs = y.attrs)
   )
   if (jack.knife == T)
      res$jkcoeffs = jkcoeffs
   
   res
}

#' Select optimal number of components for PLS model
#' 
#' @description
#' Allows user to select optimal number of components for PLS model
#' 
#' @param model
#' PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to select
#' 
#' @return
#' the same model with selected number of components
#'
#' @details
#' If number of components is not specified, the Wold's R criterion is used.
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
selectCompNum.pls = function(model, ncomp = NULL)
{
   if (!is.null(ncomp)) {
      # user defined number of components
      if (ncomp > model$ncomp || ncomp < 0)
         stop('Wrong number of selected components!')
   } else {
      # automatic estimation of the optimal number of components
      n = dim(model$coeffs$values)[3]
      
      if (!is.null(model$cvres))
         press = (model$cvres$rmse * n)^2
      else if (!is.null(model$testres))
         press = (model$testres$rmse * n)^2
      else {   
         press = (model$calres$rmse * n)^2
         warning('No validation results were found!')
      }
      
      if (model$ncomp.selcrit == 'wold') {   
         # using Wold's criterion
         ncomp = ncol(press)
         if (ncomp > 2)
         {
            r = press[, 2:ncomp] / press[, 1:(ncomp - 1)] 
            ind = which(r > 0.95, arr.ind = TRUE)
            if (length(ind) == 0)
               ncomp = 1
            else if (is.null(dim(ind)))
               ncomp = min(ind)
            else
               ncomp = min(ind[, 2])
         } else {
            ncomp = which.min(press)
         }   
      } else if (model$ncomp.selcrit == 'min') {
         # using first local minimum
         df = diff(as.vector(press)) > 0
         if (any(df))
            ncomp = which(df)[1]
         else
            ncomp = length(press)
      } else {
         stop('Wrong value for "ncomp.selcrit" argument!')
      }   
   }   
      
   model$ncomp.selected = ncomp      
   model$calres$ncomp.selected = ncomp

   if (!is.null(model$cvres)) 
      model$cvres$ncomp.selected = ncomp
   
   if (!is.null(model$testres)) 
      model$testres$ncomp.selected = ncomp

   if (!model$light)
      model$vipscores = pls.calculateVIPScores(model)
   
   model
}   

#' PLS predictions
#' 
#' @description
#' Applies PLS model to a new data set
#' 
#' @param object
#' a PLS model (object of class \code{pls})
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with reference y values (responses)
#' @param ...
#' other arguments
#' 
#' @return
#' PLS results (an object of class \code{plsres})
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'  
#' @export
predict.pls = function(object, x, y = NULL, ...) {   
   
   # preprocess x and calculate scores, total and full variance
   x.attrs = attributes(x)
   y.attrs = attributes(y)

   # correct dimension for y if it is a vector
   if (!is.null(y) && is.null(dim(y))) {
      y.rownames = names(y)
      y = matrix(y, ncol = 1)
      rownames(y) = y.rownames
   }
      
   # convert to matrices   
   x = mda.df2mat(x)
   y = mda.df2mat(y)
   
   # get dimensions
   nresp = dim(object$coeffs$values)[3]
   nobj = nrow(x)
   ncomp = object$ncomp
   
   # check dimensions
   if (ncol(x) != dim(object$coeffs$values)[1])
      stop('Wrong number of columns in matrix with predictors (x)!')
   
   if (!is.null(y)) {
      if (nrow(x) != nrow(y))
         stop('Matrices with predictors (x) and response values (y) should have the same number of rows!')
      if (ncol(y) != nresp)
         stop('Wrong number of columns in matrix with response values (y)!')
   }
   
   # autoscale x
   x = prep.autoscale(x, center = object$xcenter, scale = object$xscale)
 
   # compute x scores and residuals
   xscores = x %*% (object$weights %*% solve(crossprod(object$xloadings, object$weights)))  
   xresiduals = x - tcrossprod(xscores, object$xloadings)
   xdist = ldecomp.getDistances(xscores, object$xloadings, xresiduals, object$xtnorm) 
   
   
   # make predictions
   yp = apply(object$coeffs$values, 3, function(x, y)(y %*% x), x)
   dim(yp) = c(nrow(x), ncomp, dim(object$coeffs$values)[3])
   
   # if reference values are provided calculate and set up names for ydecomp
   y.ref = NULL
   ydecomp = NULL
   if (!is.null(y)) {   
      y.ref = y      
      
      # compute everything for yldecomp
      y = prep.autoscale(y, center = object$ycenter, scale = object$yscale)
      yscores = as.matrix(y) %*% object$yloadings 
      
      yresiduals = y - as.matrix(yp[, ncol(yp), ])
      ydist = ldecomp.getDistances(xscores, object$yloadings, yresiduals, object$ytnorm) 
      
      
      # compute total variance 
      if (length(y.attrs$exclrows) > 0)
         y = y[-y.attrs$exclrows, , drop = F]
      
      if (length(y.attrs$exclcols) > 0)
         y = y[, -y.attrs$exclcols, drop = F]
      
      ytotvar = sum(y^2)
      
      # create ydecomp object
      ydecomp = ldecomp(scores = yscores, residuals = yresiduals, loadings = object$yloadings, 
                        attrs = y.attrs, ncomp.selected = object$ncomp.selected, dist = ydist, 
                        totvar = ytotvar)
      ydecomp$totvar = ytotvar
   }
   
   # unscale predicted y values
   if (is.numeric(object$yscale))
      yp = sweep(yp, 3, object$yscale, '*')
   
   # uncenter predicted y values
   if (is.numeric(object$ycenter))
      yp = sweep(yp, 3, object$ycenter, '+')
   
   # set up all attributes and names
   dimnames(yp) = list(rownames(x), dimnames(object$coeffs$values)[[2]], dimnames(object$coeffs$values)[[3]])
   yp = mda.setattr(yp, y.attrs)
   attr(yp, 'name') = 'Response values, predicted'
   
   # we exclude rows always based on information from x (in case if y is NULL)
   attr(yp, 'exclrows') = NULL
   yp = mda.exclrows(yp, attr(x, 'exclrows'))
   
   # compute total variance 
   if (length(x.attrs$exclrows) > 0)
      x = x[-x.attrs$exclrows, , drop = F]
   
   if (length(object$exclcols) > 0)
      x = x[, -object$exclcols, drop = F]
   
   xtotvar = sum(x^2)
   
   # create xdecomp object
   xdecomp = ldecomp(scores = xscores, residuals = xresiduals, loadings = object$xloadings, attrs = x.attrs,
                     ncomp.selected = object$ncomp.selected, tnorm = object$xtnorm, totvar = xtotvar)
   xdecomp$totvar = xtotvar
  
   res = plsres(yp, y.ref = y.ref, ncomp.selected = object$ncomp.selected, 
                xdecomp = xdecomp, ydecomp = ydecomp)
   res
}  

#' Selectivity ratio calculation
#' 
#' @description
#' Calculates selectivity ration for each component and response variable in
#' the PLS model
#' 
#' @param model
#' a PLS model (object of class \code{pls})
#' @param x
#' predictor values from calibration set, preprocessed, centered and scaled
#' 
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#' 
#' @return
#' array \code{nvar x ncomp x ny} with selectivity ratio values
#' 
pls.calculateSelectivityRatio = function(model, x) {
   ny = dim(model$coeffs$values)[3]
   ncomp = dim(model$coeffs$values)[2]
   nvar = ncol(x)
   
   # get coefficients 
   coeffs = model$coeffs$values
   attrs = mda.getattr(coeffs)
   exclvars = attrs$exclrows
   nexclvar = length(exclvars)
   if (nexclvar > 0) {
      coeffs = coeffs[-exclvars, , , drop = FALSE]
      attrs$exclrow = NULL
   }
   
   selratio = array(0, dim = c(nvar, ncomp, ny))
   for (y in 1:ny) {   
      for (comp in 1:model$ncomp) {   
         b = coeffs[, comp, y]
         bnorm = sqrt(as.numeric(crossprod(b)))
         w = b/bnorm

         ttp = x %*% w
         ptp = crossprod(ttp, x) / as.numeric(crossprod(ttp))
         
         expvar = ttp %*% ptp
         resvar = apply(x - expvar, 2, var)
         expvar = apply(expvar, 2, var)
         
         selratio[, comp, y] = expvar / resvar        
      }
   }   

   if (nexclvar > 0) {
      selratio.out = array(0, dim = c(nvar + nexclvar, ncomp, ny))
      selratio.out[-exclvars, , ] = selratio
   } else {
      selratio.out = selratio
   }
   
   dimnames(selratio.out) = dimnames(model$coeffs$values)
   selratio.out = mda.setattr(selratio.out, attrs)   
   attr(selratio.out, 'name') = 'Selectivity ratio'
   selratio.out
}   

#' VIP scores calculation for PLS model
#'
#' @description
#' Calculates VIP (Variable Importance in Projection) scores for each component and 
#' response variable in the PLS model
#' 
#' @param object
#' a PLS model (object of class \code{pls})
#' 
#' @return
#' matrix \code{nvar x ny} with VIP score values for selected number of components
#' 
pls.calculateVIPScores = function(object) {
   
   # get values
   comp = object$ncomp.selected 
   coeffs = object$coeffs$values
   w = object$weights[, 1:comp, drop = F]
   xloads = object$xloadings[, 1:comp, drop = F];
   xscores = object$calres$xdecomp$scores[, 1:comp, drop = F];
   
   # remove hidden variables
   attrs = mda.getattr(coeffs)
   exclvars = attrs$exclrows
   nexclvar = length(exclvars)
   if (nexclvar > 0) {
      coeffs = coeffs[-exclvars, , , drop = FALSE]
      w = w[-exclvars, , drop = FALSE]
      xloads = xloads[-exclvars, , drop = FALSE]
   }

   # remove hidden objects
   if (length(attr(xscores, 'exclrows')) > 0)
      xscores = xscores[-attr(xscores, 'exclrows'), , drop = FALSE]

   ny = dim(coeffs)[3]
   nvar = dim(coeffs)[1]
   vipscores = matrix(0, nrow = nvar, ncol = ny)

   # regression coefficients for working with scores instead of x
   # T = X * WPW 
   # T * WPW' = X * WPW * WPW'
   # T * WPW' * (WPW * WPW')^-1 = X
   # YP = X * b 
   # YP = T * WPW' * (WPW * WPW')^-1 * b
   # YP = T * bT, where bT = WPW' * (WPW * WPW)^-1 * b
   wpw = (w %*% solve(t(xloads) %*% w))
   
   # normalise weights
   n = 1/sqrt(colSums(w^2))
   if (comp == 1)
      dim(n) = c(1, 1)
   else
      n = diag(n)                     
   wnorm = w %*% n
   
   for (y in 1:ny) {   
      b = coeffs[, comp, y, drop = F]
      dim(b) = c(dim(b)[1], 1) 
      
      bscores = ( t(wpw) %*% pinv(wpw %*% t(wpw)) ) %*% b
         
      TT = colSums(xscores^2)
      dim(TT) = c(1, comp)
      SS = bscores^2 * t(TT)
      vipscores[, y] = nvar * wnorm^2 %*% as.matrix(SS) / sum(SS)
   }   

   if (nexclvar > 0) {
      vipscores.out = matrix(0, nrow = nvar + nexclvar, ncol = ny)
      vipscores.out[-exclvars, ] = vipscores
   } else {
      vipscores.out = vipscores
   }
   
   rownames(vipscores.out) = dimnames(object$coeffs$values)[[1]]
   colnames(vipscores.out) = dimnames(object$coeffs$values)[[3]]
   vipscores.out = mda.setattr(vipscores.out, attrs)
   attr(vipscores.out, 'name') = 'VIP scores'
   vipscores.out   
}   

#' Selectivity ratio for PLS model
#' 
#' @description
#' Returns vector with selectivity ratio values for given number of components
#' and response variable
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to get the values for (if NULL user selected as optimal will be 
#' used)
#' @param ny
#' which response to get the values for (if y is multivariate)
#' @param ...
#' other parameters
#' 
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#' 
#' @return
#' vector with selectivity ratio values
#' 
#' @export
getSelectivityRatio.pls = function(obj, ncomp = NULL, ny = 1, ...) {
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   
   attrs = mda.getattr(obj$selratio)   
   selratio = obj$selratio[, ncomp, ny]
   dim(selratio) = c(dim(obj$selratio)[1], 1)
   selratio = mda.setattr(selratio, attrs)   
   rownames(selratio) = dimnames(obj$selratio)[[1]]
   colnames(selratio) = dimnames(obj$selratio)[[3]][[ny]]
   
   selratio
}  
 
#' VIP scores for PLS model
#' 
#' @description
#' Returns vector with VIP scores values for given number of components
#' and response variable
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ny
#' which response to get the values for (if y is multivariate)
#' @param ...
#' other parameters
#' 
#' @references
#' [1] Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), pp. 103-112.
#' 
#' @return
#' vector with VIP scores values
#' 
#' @export
getVIPScores.pls = function(obj, ny = 1, ...) {
   vipscores = mda.subset(obj$vipscores, select = ny)
}


#' Regression coefficients for PLS model'
#'
#' @description 
#' Returns a matrix with regression coefficients for
#' the PLS model which can be applied to a data directly
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to return the coefficients for
#' @param ny
#' if y is multivariate which variables you want to see the coefficients for
#' @param full
#' if TRUE the method also shows p-values for the coefficients (if available)
#' @param ...
#' other parameters
#'
#' @details 
#' The method recalculates the regression coefficients found by the PLS algorithm
#' taking into account centering and scaling of predictors and responses, so the 
#' matrix with coefficients can be applied directly to original data (yp = Xb).
#' 
#' If number of components is not specified, the optimal number, selected by user
#' or identified by a model will be used.
#'  
#' @return 
#' A matrix (n of predictors x n of responses) with regression coefficients.
#'  
#' @export
getRegcoeffs.pls = function(obj, ncomp = NULL, ny = NULL, full = FALSE, ...) {
   if (is.null(ncomp)) 
      ncomp = obj$ncomp.selected
   else if (length(ncomp) != 1 || ncomp <= 0 || ncomp > obj$ncomp) 
      stop('Wrong value for number of components!')
   
   if (is.null(ny))
      ny = 1:dim(obj$coeffs$values)[3]
   
   coeffs = obj$coeffs$values[, ncomp, ny, drop = F]
   coeffs = matrix(coeffs, nrow = dim(coeffs)[1], ncol = dim(coeffs)[3])
   xscale = obj$xscale
   if (is.logical(xscale))
      xscale = matrix(1, nrow = nrow(coeffs))
  
   xcenter = obj$xcenter 
   if (is.logical(xcenter))
      xcenter = matrix(0, nrow = nrow(coeffs))
      
   yscale = obj$yscale
   if (is.logical(yscale))
      yscale = matrix(1, nrow = ncol(coeffs))
  
   ycenter = obj$ycenter 
   if (is.logical(ycenter))
      ycenter = matrix(0, nrow = ncol(coeffs))
      
   # calculate intercept
   b0 = sweep(coeffs,  1, xcenter, '*')
   b0 = sweep(b0, 1, xscale, '/')
   b0 = apply(-b0, 2, sum)
   b0 = matrix(b0, nrow = 1)
   b0 = sweep(b0, 2, yscale, '*')
   b0 = sweep(b0, 2, ycenter, '+')
   
   # rescale coefficients
   coeffs = sweep(coeffs, 1, xscale, '/');
   coeffs = sweep(coeffs, 2, yscale, '*');
   coeffs = rbind(b0, coeffs)
   
   if (full == TRUE && !is.null(obj$coeffs$p.values)) {
      if (length(ny) != 1)
         stop('Full table can be shown only for selected response variable (ny)!')
      
      coeffs = cbind(coeffs, c(NA, obj$coeffs$p.values[, ncomp, ny]))
      colnames(coeffs) = c(dimnames(obj$coeffs$values)[[3]][ny], 'p-value')
   } else {
      colnames(coeffs) = dimnames(obj$coeffs$values)[[3]]
   }
   
   rownames(coeffs) = c('Intercept', dimnames(obj$coeffs$values)[[1]])
   coeffs
}   

#' VIP scores plot for PLS model
#' 
#' @description
#' Shows a plot with VIP scores values for given number of components
#' and response variable
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ny
#' which response to get the values for (if y is multivariate)
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @references
#' [1] Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), pp. 103-112.
#' 
#' @export
plotVIPScores.pls = function(obj, ny = 1, type = 'l', main = NULL, ylab = '', ...) {   
   
   main = getMainTitle(main, NULL, 'VIP scores')
   vipscores = getVIPScores.pls(obj, ny)
   mdaplot(mda.t(vipscores), type = type, main = main, ylab = ylab, ...)
}

#' Selectivity ratio plot for PLS model
#' 
#' @description
#' Shows a plot with selectivity ratio values for given number of components
#' and response variable
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to get the values for (if NULL user selected as optimal will be used)
#' @param ny
#' which response to get the values for (if y is multivariate)
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#' 
#' @export
plotSelectivityRatio.pls = function(obj, ncomp = NULL, ny = 1, type = 'l', 
                                    main = NULL, ylab = '', ...) {
   main = getMainTitle(main, ncomp, 'Selectivity ratio')
   ncomp = getSelectedComponents(obj, ncomp)
   
   selratio = getSelectivityRatio(obj, ncomp, ny)
   mdaplot(mda.t(selratio), type = type, main = main, ylab = ylab, ...)
} 


#' RMSE plot for PLS
#' 
#' @description
#' Shows plot with root mean squared error values vs. number of components for PLS model.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ny
#' number of response variable to make the plot for (if y is multivariate)
#' @param type
#' type of the plot('b', 'l' or 'h')
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param labels
#' what to show as labels (if this option is on)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotRMSE.pls = function(obj, ny = 1, type = 'b', main = 'RMSE', xlab = 'Components', ylab = NULL, 
                        labels = 'values', ...) {
   ncomp = obj$ncomp
   
   if (is.null(ylab)) {   
      if (nrow(obj$calres$rmse) == 1)
         ylab = 'RMSE'
      else if (is.null(rownames(obj$calres$rmse)))
         ylab = sprintf('RMSE (#%d)', ny)
      else         
         ylab = sprintf('RMSE (%s)', rownames(obj$calres$rmse)[ny])
   }
   
   data = list()
   data$cal = obj$calres$rmse[ny, , drop = F]   

   if (!is.null(obj$cvres))  
      data$cv = obj$cvres$rmse[ny, ]
   
   if (!is.null(obj$testres)) 
      data$test = obj$testres$rmse[ny, ]
  
   mdaplotg(data, type = type, labels = labels, main = main, xlab = xlab, ylab = ylab, ...)
}

#' Explained X variance plot for PLS
#' 
#' @description
#' Shows plot with explained X variance vs. number of components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot('b', 'l' or 'h')
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotXVariance.pls = function(obj, type = 'b',
                             main = 'X variance', xlab = 'Components', 
                             ylab = 'Explained variance, %', ...) {
   
   plotVariance(obj, decomp = 'xdecomp', type = type, main = main, xlab = xlab, ylab = ylab, ...)
}

#' Explained Y variance plot for PLS
#' 
#' @description
#' Shows plot with explained Y variance vs. number of components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot('b', 'l' or 'h')
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotYVariance.pls = function(obj, type = 'b',
                             main = 'Y variance', xlab = 'Components', 
                             ylab = 'Explained variance, %', ...) {
   
   plotVariance(obj, decomp = 'ydecomp', type = type, main = main, xlab = xlab, ylab = ylab, ...)
}

#' Cumulative explained X variance plot for PLS
#' 
#' @description
#' Shows plot with cumulative explained X variance vs. number of components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot('b', 'l' or 'h')
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotXCumVariance.pls = function(obj, type = 'b',
                                main = 'X cumulative variance', 
                                xlab = 'Components', ylab = 'Explained variance, %', ...) {
   
   plotVariance(obj, decomp = 'xdecomp', variance = 'cumexpvar', 
                type = type, main = main, xlab = xlab, ylab = ylab)   
}

#' Cumulative explained Y variance plot for PLS
#' 
#' @description
#' Shows plot with cumulative explained Y variance vs. number of components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot('b', 'l' or 'h')
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotYCumVariance.pls = function(obj, type = 'b',
                                main = 'Y cumulative variance', 
                                xlab = 'Components', ylab = 'Explained variance, %', ...) {
   
   plotVariance(obj, decomp = 'ydecomp', variance = 'cumexpvar', 
                type = type, main = main, xlab = xlab, ylab = ylab, ...)   
}

#' Variance plot for PLS
#' 
#' @description
#' Shows plot with variance values vs. number of components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param decomp
#' which decomposition to use ('xdecomp' for x or 'ydecomp' for y)
#' @param variance
#' which variance to use ('expvar', 'cumexpvar)
#' @param type
#' type of the plot('b', 'l' or 'h')
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param labels
#' what to show as labels for plot objects.
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotVariance.pls = function(obj, decomp = 'xdecomp', variance = 'expvar',
                            type = 'b', main = 'X variance', xlab = 'Components', 
                            ylab = 'Explained variance, %', labels = 'values', ...) {
   ncomp = obj$ncomp
   
   data = list()
   data$cal = obj$calres[[decomp]][[variance]]   

   if (!is.null(obj$cvres)) 
      data$cv = obj$cvres[[decomp]][[variance]]

   if (!is.null(obj$testres)) 
      data$test = obj$testres[[decomp]][[variance]]

   mdaplotg(data, type = type, main = main, xlab = xlab, ylab = ylab, labels = labels, ...)
}

#' X scores plot for PLS
#' 
#' @description
#' Shows plot with X scores values for selected components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param main
#' main plot title
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotXScores.pls = function(obj, comp = c(1, 2), main = 'X scores', show.axes = T, ...) {
   ncomp = length(comp)
   
   if (ncomp < 1 || ncomp > 2)
      stop('The plot can be made for one or two components only!')
      
   data = list()
   data$cal = mda.subset(obj$calres$xdecomp$scores, select = comp)
   colnames(data$cal) = paste('Comp ', comp, ' (', round(obj$calres$xdecomp$expvar[comp], 2), '%)');
   
   if (!is.null(obj$testres)) 
      data$test = mda.subset(obj$testres$xdecomp$scores, select = comp)
   
   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, show.lines = show.lines, main = main, ...)   
}  

#' XY scores plot for PLS
#' 
#' @description
#' Shows plot with X vs. Y scores values for selected component.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which component to show the plot for
#' @param main
#' main plot title
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXYScores.pls = function(obj, comp = 1, main = 'XY scores', show.axes = T, ...) {
   if (is.null(comp)) 
      comp = obj$ncomp.selected   
   else if (length(comp) != 1 || comp <= 0 || comp > obj$ncomp) 
      stop('Wrong component number!')
   
   data = list() 
   data$cal = cbind(
      obj$calres$xdecomp$scores[, comp, drop = F], 
      obj$calres$ydecomp$scores[, comp, drop = F]
   )
   data$cal = mda.setattr(data$cal, mda.getattr(obj$calres$xdecomp$scores))
   attr(data$cal, 'name') = sprintf('XY scores')
   rownames(data$cal) = rownames(obj$calres$xdecomp$scores)
   colnames(data$cal) = c(
      sprintf('X scores (Comp %d, %.2f%%)', comp, obj$calres$xdecomp$expvar[comp]), 
      sprintf('Y scores (Comp %d, %.2f%%)', comp, obj$calres$ydecomp$expvar[comp])
   )
   

   if (!is.null(obj$testres)){
      data$test = cbind(
         obj$testres$xdecomp$scores[, comp, drop = F], 
         obj$testres$ydecomp$scores[, comp, drop = F]
      )
      data$test = mda.setattr(data$test, mda.getattr(obj$testres$xdecomp$scores))
      rownames(data$test) = rownames(obj$testres$xdecomp$scores)
   }

   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, main = main, show.lines = show.lines, ...)   
}  

#' Predictions plot for PLS
#' 
#' @description
#' Shows plot with predicted vs. reference (measured) y values for selected components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' number of response variable to make the plot for (if y is multivariate)
#' @param main
#' main plot title
#' @param legend.position
#' position of legend on the plot (if shown)
#' @param show.line
#' logical, show or not line fit for the plot points
#' @param colmap
#' a colormap to use for coloring the plot items
#' @param col
#' a vector with color values for target lines fitted the points
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotPredictions.pls = function(obj, ncomp = NULL, ny = 1, main = NULL, legend.position = 'topleft', 
                               show.line = T, colmap = 'default', col = NULL, ...) {
   
   if (is.null(main)) {   
      if (is.null(ncomp))
         main = 'Predictions'
      else
         main = sprintf('Predictions (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   else if (length(ncomp) != 1 || ncomp < 1 || ncomp > ncol(obj$calres$y.pred))
      stop('Wrong number of components!')
   
   if (is.null(dimnames(obj$calres$y.pred)) || is.null(dimnames(obj$calres$y.pred)[[3]])) 
      yaxis.name = 'y, predicted'
   else
      yaxis.name = sprintf('%s, predicted', dimnames(obj$calres$y.pred)[[3]][ny])
   
   if (!is.null(dimnames(obj$calres$y.pred)) && !is.null(dimnames(obj$calres$y.pred)[[3]])) 
      xaxis.name = sprintf('%s, reference', dimnames(obj$calres$y.pred)[[3]][ny])
   else
      xaxis.name = 'y, reference'
   
   data = list()   
   attrs = mda.getattr(obj$calres$y.pred)
   data$cal = cbind(obj$calres$y.ref[, ny], obj$calres$y.pred[, ncomp, ny])
   data$cal = mda.setattr(data$cal, attrs)
   colnames(data$cal) = c(xaxis.name, yaxis.name)
   rownames(data$cal) = rownames(obj$calres$y.pred)

   if (!is.null(obj$cvres)) { 
      data$cv = cbind(obj$cvres$y.ref[, ny], obj$cvres$y.pred[, ncomp, ny])
      colnames(data$cv) = c(xaxis.name, yaxis.name)
      rownames(data$cv) = rownames(obj$cvres$y.pred)
   }   
   
   if (!is.null(obj$testres)) { 
      attrs = mda.getattr(obj$testres$y.pred)
      data$test = cbind(obj$testres$y.ref[, ny], obj$testres$y.pred[, ncomp, ny])
      data$test = mda.setattr(data$test, attrs)
      colnames(data$test) = c(xaxis.name, yaxis.name)
      rownames(data$test) = rownames(obj$testres$y.pred)
   }   
   
   mdaplotg(data, main = main, colmap = colmap, legend.position = legend.position, ...)  
   
   if (show.line == T)
      mdaplot.showRegressionLine(data, colmap = colmap, col = col)
}

#' Y residuals plot for PLS
#' 
#' @description
#' Shows plot with y residuals for selected components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' number of response variable to make the plot for (if y is multivariate)
#' @param main
#' main plot title
#' @param show.line
#' logical, show or not line for 0 value
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotYResiduals.pls = function(obj, ncomp = NULL, ny = 1, main = NULL, show.line = T, ...) {
   
   if (is.null(main)) {   
      if (is.null(ncomp))
         main = 'Y residuals'
      else
         main = sprintf('Y residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp)) 
      ncomp = obj$ncomp.selected
   else if (length(ncomp) != 1 || ncomp <= 0 || ncomp > obj$ncomp) 
      stop('Wrong value for number of components!')

   if (length(ny) != 1)
      stop('You can show prediction plot only for one selected response variable!')
   
   if (show.line == T)
      show.line = c(NA, 0)
   
   if (!is.null(dimnames(obj$calres$y.pred)) && !is.null(dimnames(obj$calres$y.pred)[[3]])) 
      xaxis.name = sprintf('%s, reference', dimnames(obj$calres$y.pred)[[3]][ny])
   else
      xaxis.name = 'y, reference'

   data = list()   
   attr = mda.getattr(obj$calres$y.pred)
   data$cal = cbind(obj$calres$y.ref[, ny], obj$calres$y.ref[, ny] - obj$calres$y.pred[, ncomp, ny])
   data$cal = mda.setattr(data$cal, attr)
   colnames(data$cal) = c(xaxis.name, 'Residuals')
   rownames(data$cal) = rownames(obj$calres$y.pred)
   
   if (!is.null(obj$cvres)) { 
      data$cv = cbind(obj$cvres$y.ref[, ny], obj$cvres$y.ref[, ny] - obj$cvres$y.pred[, ncomp, ny])
      colnames(data$cv) = c(xaxis.name, 'Residuals')
      rownames(data$cv) = rownames(obj$cvres$y.pred)
   }   
   
   if (!is.null(obj$testres)) { 
      attr = mda.getattr(obj$testres$y.pred)
      data$test = cbind(obj$testres$y.ref[, ny], obj$testres$y.ref[, ny] - obj$testres$y.pred[, ncomp, ny])
      data$test = mda.setattr(data$test, attr)
      colnames(data$test) = c(xaxis.name, 'Residuals')
      rownames(data$test) = rownames(obj$testres$y.pred)
   }   
   
   mdaplotg(data, type = 'p', main = main, show.lines = show.line, ...)      
}

#' Regression coefficient plot for PLS
#' 
#' @description
#' Shows plot with regression coefficient values for selected components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} and \code{plot.regcoeffs} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotRegcoeffs.pls = function(obj, ncomp = NULL, ...) {
   
   if (is.null(ncomp)) 
      ncomp = obj$ncomp.selected
   else if (ncomp <= 0 || ncomp > obj$ncomp) 
      stop('Wrong value for number of components!')
   
   plot(obj$coeffs, ncomp = ncomp, ...)
}

#' X loadings plot for PLS
#' 
#' @description
#' Shows plot with X loading values for selected components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotXLoadings.pls = function(obj, comp = c(1, 2), type = 'p', main = 'X loadings', show.axes = T, ...) {
   ncomp = length(comp)
   
   if (min(comp) < 1 || max(comp) > obj$ncomp)
      stop('Wrong components number!')
   
   if (type == 'p' && ncomp != 2)
      stop('Scatter plot can be made only for two components!')

   if (show.axes == T) {
      if (type == 'p')
         show.lines = c(0, 0)
      else
         show.lines = c(NA, 0)
   } else {
      show.lines = F
   }
   
   data = mda.subset(obj$xloadings, select = comp)
   if (type == 'p') {
      colnames(data) = paste('Comp ', comp, ' (', round(obj$calres$xdecomp$expvar[comp], 2) , '%)', sep = '')
      mdaplot(data, main = main, type = type, show.lines = show.lines, ...)
   } else {
      data = mda.t(data)
      mdaplotg(data, main = main, type = type, show.lines = show.lines, ...)
   }
}


#' XY loadings plot for PLS
#' 
#' @description
#' Shows plot with X and Y loading values for selected components.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param main
#' main plot title
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotXYLoadings.pls = function(obj, comp = c(1, 2), main = 'XY loadings', show.axes = F, ...) {
   
   if (length(comp) != 2)
      stop('This plot can be made for only two components!')
  
   data = list() 
   data$X = mda.subset(obj$xloadings, select = comp)
   data$Y = mda.subset(obj$yloadings, select = comp)

   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, main = main, type = 'p', show.lines = show.lines, ...)
}


#' X residuals plot for PLS
#' 
#' @description
#' Shows a plot with Q residuals vs. Hotelling T2 values for PLS decomposition of x data.
#' 
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
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
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plotXResiduals.pls = function(obj, ncomp = NULL, main = NULL, xlab = 'T2', 
                              ylab = 'Squared residual distance (Q)', ...) {
   
   if (is.null(main)) {   
      if (is.null(ncomp))
         main = 'X residuals'
      else
         main = sprintf('X residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp)) 
      ncomp = obj$ncomp.selected
   else if (length(ncomp) != 1 || ncomp <= 0 || ncomp > obj$ncomp) 
      stop('Wrong value for number of components!')
   
   data = list()
   attr = mda.getattr(obj$calres$xdecomp$Q)
   data$cal = cbind(obj$calres$xdecomp$T2[, ncomp, drop = F], obj$calres$xdecomp$Q[, ncomp, drop = F])
   data$cal = mda.setattr(data$cal, attr)
   rownames(data$cal) = rownames(obj$calres$xdecomp$Q)
   
   if (!is.null(obj$cvres)) { 
      data$cv = cbind(obj$cvres$xdecomp$T2[, ncomp, drop = F], obj$cvres$xdecomp$Q[, ncomp, drop = F])
      rownames(data$cv) = rownames(obj$cvres$xdecomp$Q)
   }   
   
   if (!is.null(obj$testres)) { 
      attr = mda.getattr(obj$testres$xdecomp$Q)
      data$test = cbind(obj$testres$xdecomp$T2[, ncomp, drop = F], obj$testres$xdecomp$Q[, ncomp, drop = F])
      data$test = mda.setattr(data$test, attr)
      rownames(data$test) = rownames(obj$testres$xdecomp$Q)
   }   
   
   mdaplotg(data, main = main, xlab = xlab, ylab = ylab, ...)
} 

#' Model overview plot for PLS
#' 
#' @description
#' Shows a set of plots (x residuals, regression coefficients, RMSE and predictions) for PLS model.
#' 
#' @param x
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' which y variable to show the summary for (if NULL, will be shown for all)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
#' @export
plot.pls = function(x, ncomp = NULL, ny = 1, show.legend = T, show.labels = F, ...) {
   obj = x
   
   if (!is.null(ncomp) && (ncomp <= 0 || ncomp > obj$ncomp)) 
      stop('Wrong value for number of components!')
   
   par(mfrow = c(2, 2))      
   plotXResiduals(obj, ncomp = ncomp, show.labels = show.labels, show.legend = show.legend)   
   plotRegcoeffs(obj, ncomp = ncomp, ny = ny, show.labels = F)   
   plotRMSE(obj, ny = ny, show.legend = show.legend)   
   plotPredictions(obj, ncomp = ncomp, ny = ny, show.labels = show.labels, show.legend = show.legend)   
   par(mfrow = c(1, 1))
}

#' Summary method for PLS model object
#' 
#' @description
#' Shows performance statistics for the model.
#' 
#' @param object
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' which y variable to show the summary for (if NULL, will be shown for all)
#' @param ...
#' other arguments
#' 
#' @export
summary.pls = function(object, ncomp = NULL, ny = NULL, ...) {
   obj = object
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   else if (ncomp <= 0 || ncomp > obj$ncomp)
      stop('Wrong value for number of components!')
   
   if (is.null(ny))
      ny = 1:ncol(obj$calres$y.ref)
   
   cat('\nPLS model (class pls) summary\n')
   cat('\nPerformance and validation:\n')
   cat(sprintf('Number of selected components: %d\n', ncomp))
   
   for (y in ny) {   
      if (ncol(obj$calres$y.ref) > 1)
         cat(sprintf('\nResponse variable #%d (%s)\n', y, colnames(obj$calres$y.ref)[y]))
      
      data = as.matrix(obj$calres, ncomp = ncomp, ny = y)
      rownames(data) = 'Cal'
      
      if (!is.null(obj$cvres))
      {
         data = rbind(data, as.matrix(obj$cvres, ncomp = ncomp, ny = y))      
         rownames(data)[nrow(data)] = 'CV'
      }
      
      if (!is.null(obj$testres))
      {
         data = rbind(data, as.matrix(obj$testres, ncomp = ncomp, ny = y))
         rownames(data)[nrow(data)] = 'Test'
      }   
      
      data = data[, -c(1, 3), drop = F]
      data[, 1:2] = round(data[, 1:2], 2)      
      data[, 3] = mdaplot.formatValues(data[, 3], round.only = T)
      data[, 4] = round(data[, 4], 2)  
      data[, 5] = round(data[, 5], 4)      
      data[, 6] = round(data[, 6], 2)      
      
      print(data)
   }   
   cat('\n')
}

#' Print method for PLS model object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a PLS model (object of class \code{pls})
#' @param ...
#' other arguments
#'
#' @export 
print.pls = function(x, ...) {
   obj = x
   
   cat('\nPLS model (class pls)\n')
   cat('\nCall:\n')
   print(obj$call)
   cat('\nMajor fields:\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$coeffs - object (regcoeffs) with regression coefficients\n')
   cat('$xloadings - vector with x loadings\n')
   cat('$yloadings - vector with y loadings\n')
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
