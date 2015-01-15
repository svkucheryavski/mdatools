# class and methods for Partial Least Squares regression #

pls = function(x, y, ncomp = 15, center = T, scale = F, cv = NULL, 
               x.test = NULL, y.test = NULL, method = 'simpls', alpha = 0.05, 
               coeffs.ci = NULL, coeffs.alpha = 0.1, info = '')
{
   # Calibrate and validate a PLS model.
   #
   # Arguments:
   #   x: a matrix with predictor values    
   #   y: a vector with response values
   #   ncomp: maximum number of components to calculate
   #   center: logical, center or not x and y data
   #   scale: logical, standardize or not x data
   #   cv: cross-validation settings (see \code{?crossval} for details)
   #   x.test: a matrix with predictor values for test set validation
   #   y.test: a vector with response values for test set validation
   #   method: a method to calculate PLS model
   #   alpha: a sigificance limit for Q2 values
   #   info: a short string with information about the model
   #
   # Returns:
   #   model: a PLS model (object of class pls) 

   x = as.matrix(x)
   y = as.matrix(y)
   
   if (is.null(colnames(y)))
      colnames(y) = paste('y', 1:ncol(y), sep = '')
   
   # set LOO cross-validation if jack.knife is selected
   jack.knife = F
   if (!is.null(coeffs.ci) && coeffs.ci == 'jk')
   {
      jack.knife = T
      if (is.null(cv))
      {   
         cv = 1      
      }   
      else if (cv < 10)
      {
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
   model = pls.cal(x, y, ncomp, center = center, scale = scale, method = method)
   model$alpha = alpha   
   model$calres = predict.pls(model, x, y)
   model = selectCompNum.pls(model, model$ncomp)
   # do cross-validation if needed
   if (!is.null(cv))
   {   
      res = pls.crossval(model, x, y, cv, center = center, scale = scale, 
                                 jack.knife = T)    
      if (jack.knife == T)
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
   if (!is.null(x.test) && !is.null(y.test))
   {
      x.test = as.matrix(x.test)
      y.test = as.matrix(y.test)
      
      if (is.null(colnames(y.test)))
         colnames(y.test) = paste('y', 1:ncol(y), sep = '')
      
      if (ncol(x.test) != ncol(x) || ncol(y.test) != ncol(y))
         stop('Calibration and data set should have the same number of variables!')
      
      if (nrow(x.test) != nrow(y.test))
         stop('Number of rows in x.test and y.test should be the same!')
      
      model$testres = predict.pls(model, x.test, y.test)
   }
   
   model$call = match.call()
   model$info = info
   
   class(model) = "pls"
   
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
#' algorithm for compiting PC space (only 'simpls' is supported so far)
#' @param cv
#' logical, does calibration for cross-validation or not
#' 
#' @return model
#' an object with calibrated PLS model
#' 
pls.cal = function(x, y, ncomp, center, scale, method = 'simpls', cv = FALSE)
{   
   # center and scale data according to arguments
   x = prep.autoscale(as.matrix(x), center = center, scale = scale)
   y = prep.autoscale(as.matrix(y), center = center, scale = scale)
   
   # do PLS
   if (method == 'simpls')
      res = pls.simpls(x, y, ncomp, cv = cv)
   else
      stop('Method with this name is not supported!')

   # return a list with model parameters
   model = list(
      xloadings = res$xloadings,
      yloadings = res$yloadings,
      weights = res$weights,
      coeffs = regcoeffs(res$coeffs),
      method = method,
      xtnorm = sqrt(colSums(res$xscores ^ 2)/(nrow(res$xscores) - 1)),   
      ytnorm = sqrt(colSums(res$yscores ^ 2)/(nrow(res$yscores) - 1)),   
      xcenter = attr(x, 'prep:center'),
      xscale = attr(x, 'prep:scale'),
      ycenter = attr(y, 'prep:center'),
      yscale = attr(y, 'prep:scale'),
      ncomp = res$ncomp
   )
   
   if (cv == FALSE)
   {   
      model$selratio = pls.calculateSelectivityRatio(model, x)
   }
   
   model
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
pls.simpls = function(x, y, ncomp, cv = FALSE)
{
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
   for (n in 1:ncomp)
   {
      # get the dominate eigenvector of A'A
      e = eigen(t(A) %*% A)
      q = e$vectors[1:nresp]
      
      
      # calculate and store weights
      w = A %*% q
      c = t(w) %*% M %*% w
      w = w/sqrt(as.numeric(c))
      W[, n] = w
      
      # calculate and store x loadings
      p = M %*% w
      P[, n] = p
      
      # calculate and store y loadings
      q = t(A) %*% w
      Q[, n] = q
      
      v = C %*% p
      v = v/sqrt(as.numeric(t(v) %*% v))
      
      # calculate and store regression coefficients
      B[, n, ] = W[, 1:n, drop = FALSE] %*% t(Q[, 1:n, drop = FALSE])
      
      # recalculate matrices for the next compnonent
      C = C - v %*% t(v)
      M = M - p %*% t(p)
      A = C %*% A      
      
      if (cv == F && e$value < 10^-12) {
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
   
   # set names for all results
   compnames = paste('Comp', 1:n)
   colnames(Q) = colnames(B) = colnames(P) = colnames(W) = compnames      
   colnames(TT) = colnames(U) = compnames      
   rownames(P) = rownames(B) = rownames(W) = prednames
   rownames(TT) = rownames(U) = objnames
   rownames(Q) = dimnames(B)[[3]] = respnames

   res = list(
      coeffs = B,
      weights = W,
      xloadings = P,
      xscores = TT,
      yloadings = Q,
      yscores = U,
      ncomp = n
   )   
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
#' @param jack.knife
#' logical, do jack-knifing or not
#' 
#' @return
#' object of class \code{plsres} with results of cross-validation
#'  
pls.crossval = function(model, x, y, cv, center = T, scale = F, jack.knife = T)
{
   x = as.matrix(x)
   y = as.matrix(y)
   
   ncomp = model$ncomp
   nobj = nrow(x)
   nvar = ncol(x)
   nresp = ncol(y)
   
   # get matrix with indices for cv segments
   idx = crossval(nobj, cv)
   seglen = ncol(idx);
   nseg = nrow(idx);
   nrep = dim(idx)[3]

   yp = array(0, dim = c(nobj, ncomp, nresp))
   Q2x = matrix(0, ncol = ncomp, nrow = nobj)   
   T2x = matrix(0, ncol = ncomp, nrow = nobj)   
   Q2y = matrix(0, ncol = ncomp, nrow = nobj)   
   T2y = matrix(0, ncol = ncomp, nrow = nobj)   
   jkcoeffs = array(0, dim = c(nvar, ncomp, ncol(y), nrow(idx)))
   
   # loop over segments and repetitions
   for (iRep in 1:nrep)
   {   
      for (iSeg in 1:nseg)
      {         
         ind = na.exclude(idx[iSeg, , iRep])
         
         if (length(ind) > 0)
         {   
            xc = x[-ind, , drop = F]
            yc = y[-ind, , drop = F]
            xt = x[ind, , drop = F]
            yt = y[ind, , drop = F]
            
            m = pls.cal(xc, yc, ncomp, center = center, scale = scale, cv = TRUE)            
            res = predict.pls(m, xt, yt, cv = T)
            
            xdist = ldecomp.getDistances(res$xscores, m$xloadings, res$xresiduals, 
                                         model$calres$xdecomp$tnorm)
            ydist = ldecomp.getDistances(res$xscores, m$yloadings, res$yresiduals, 
                                         model$calres$xdecomp$tnorm)
            
            dim(m$coeffs$values) = c(dim(m$coeffs$values), 1)
            
            yp[ind, , ] = yp[ind, , , drop = F]  + res$yp
            Q2x[ind, ]  = Q2x[ind, , drop = F] + xdist$Q2        
            T2x[ind, ]  = T2x[ind, , drop = F] + xdist$T2        
            Q2y[ind, ]  = ydist$Q2 + Q2y[ind, , drop = F]       
            T2y[ind, ]  = ydist$T2 + T2y[ind, , drop = F]        
            jkcoeffs[, , , iSeg] = jkcoeffs[, , , iSeg, drop = F] + m$coeffs$values
         }   
      }      
   }  
   
   # average results over repetitions
   yp = yp / nrep
   Q2x = Q2x / nrep
   T2x = T2x / nrep
   Q2y = Q2y / nrep
   T2y = T2y / nrep
   jkcoeffs = jkcoeffs / nrep
   
   dimnames(yp) = list(rownames(x), colnames(model$coeffs$values), colnames(model$calres$y.ref))         
   
   res = plsres(yp, y.ref = y, ncomp.selected = model$ncomp.selected,
                xdecomp = ldecomp(totvar = model$calres$xdecomp$totvar,
                                  tnorm = model$calres$xdecomp$tnorm,
                                  ncomp.selected = model$ncomp.selected,
                                  Q2 = Q2x, T2 = T2x),
                ydecomp = ldecomp(totvar = model$calres$ydecomp$totvar,
                                  tnorm = model$calres$ydecomp$tnorm,
                                  ncomp.selected = model$ncomp.selected,
                                  Q2 = Q2y, T2 = T2y)                   
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
#' See examples in help for \code{\link{pls}} function.
#' 
selectCompNum.pls = function(model, ncomp)
{
   if (ncomp > model$ncomp || ncomp < 0)
      stop('Wrong number of selected components!')
   
   model$ncomp.selected = ncomp      
   model$calres$ncomp.selected = ncomp
   
   if (!is.null(model$cvres)) 
      model$cvres$ncomp.selected = ncomp
   
   if (!is.null(model$testres)) 
      model$testres$ncomp.selected = ncomp

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
#' @param y.ref
#' a matrix with reference y values (responses)
#' @param cv
#' logical, are predictions for cross-validation or not
#' @param ...
#' other arguments
#' 
#' @return
#' PLS results (an object of class \code{plsres})
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'  
predict.pls = function(object, x, y.ref = NULL, cv = F, ...)
{   
   nresp = dim(object$coeffs$values)[3]
   nobj = nrow(x)
   
   # preprocess x and calculate scores, total and full variance
   x = as.matrix(x)
   x  = prep.autoscale(x, center = object$xcenter, scale = object$xscale)
   xscores = x %*% (object$weights %*% solve(t(object$xloadings) %*% object$weights))  
   xresiduals = x - xscores %*% t(object$xloadings)
   
   # make predictions
   yp = array(0, dim = c(nrow(x), object$ncomp, nresp))
   for (i in 1:nresp)
   {   
      yp[, , i] = x %*% object$coeffs$values[, , i]      
   }
   
   # if y is provided, calculate y residuals
   if (!is.null(y.ref))
   {   
      yy = prep.autoscale(y.ref, center = object$ycenter, scale = object$yscale)
      yscores = as.matrix(yy) %*% object$yloadings   

      ypp = yp[, ncol(yp), , drop = F]
      dim(ypp) = dim(yp)[c(1, 3)]
      yresiduals = ypp - yy
      
      dimnames(yp) = list(rownames(x), colnames(object$coeffs$values), colnames(y.ref))   
   }
   else
   {
      dimnames(yp) = list(rownames(x), colnames(object$coeffs$values), colnames(object$calres$y.ref))         
   }   
      
   # unscale predicted y values
   if (is.numeric(object$yscale))
      for (i in 1:nresp)
         yp[, , i] = sweep(yp[, , i, drop = F], 2L, object$yscale[i], '*', check.margin = F)
   
   # uncenter predicted y values
   if (length(object$ycenter) > 1 || object$ycenter != F)
      for (i in 1:nresp)
         yp[, , i] = sweep(yp[, , i, drop = F], 2L, object$ycenter[i], '+', check.margin = F)
   
   if (cv == F)
   {   
      # normal predictions      
      # calculate PLS decomposition for x and y data
      # and return all results
      rownames(xscores) = rownames(x)
      colnames(xscores) = paste("Comp", 1:object$ncomp)
      xdecomp = ldecomp(xscores, object$xloadings, xresiduals, 
                        totvar = sum(x^2),
                        ncomp.selected = object$ncomp.selected)
      if (!is.null(y.ref))
      {   
         dimnames(yscores) = dimnames(xscores)
         ydist = ldecomp.getDistances(xscores, object$yloadings, yresiduals)
         ydecomp = ldecomp(yscores, object$yloadings, yresiduals, sum(yy^2),
                           object$ytnorm, object$ncomp.selected,
                           ydist$T2, ydist$Q2)
      }
      else
      {
         ydecomp = NULL
      }      
      
      res = plsres(yp, y.ref = y.ref, ncomp.selected = object$ncomp.selected, 
                   xdecomp = xdecomp, ydecomp = ydecomp)    
   }      
   else
   {
      # predictions for cross-validation      
      # just return predictions, scores and residuals
      # decomposition will be calculated in crossval() function
      res = list(
         yp = yp,
         xscores = xscores,
         yscores = yscores,
         xresiduals = xresiduals,
         yresiduals = yresiduals
      )
   }         
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
pls.calculateSelectivityRatio = function(model, x)
{
   ny = dim(model$coeffs$values)[3]
   ncomp = dim(model$coeffs$values)[2]
   nvar = dim(model$coeffs$values)[1]
   
   selratio = array(0, dim = c(nvar, ncomp, ny))
   
   for (y in 1:ny)
   {   
      for (comp in 1:model$ncomp)
      {   
         b = model$coeffs$values[, comp, y]
         bnorm = sqrt(sum(b^2))
         w = b/bnorm

         ttp = x %*% w
         ptp = t(ttp) %*% x / as.numeric(t(ttp) %*% ttp)
   
         expvar = ttp %*% ptp
         resvar = apply(x - expvar, 2, var)
         expvar = apply(expvar, 2, var)
         
         selratio[, comp, y] = expvar / resvar        
      }
   }   
   
   dimnames(selratio) = dimnames(model$coeffs$values)

   selratio
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
pls.calculateVIPScores = function(object)
{
   ny = dim(object$coeffs$values)[3]
   nvar = dim(object$coeffs$values)[1]
   vipscores = matrix(0, nrow = nvar, ncol = ny)
   
   comp = object$ncomp.selected   

   w = object$weights[, 1:comp, drop = F]
   xloads = object$xloadings[, 1:comp, drop = F];
   xscores = object$calres$xdecomp$scores[, 1:comp, drop = F];

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
   
   for (y in 1:ny)
   {   
      b = object$coeffs$values[, comp, y, drop = F]
      dim(b) = c(dim(b)[1], 1) 
      
      bscores = ( t(wpw) %*% pinv(wpw %*% t(wpw)) ) %*% b
         
      TT = colSums(xscores^2)
      dim(TT) = c(1, comp)
      SS = bscores^2 * t(TT)
      vipscores[, y] = nvar * wnorm^2 %*% as.matrix(SS) / sum(SS)
   }   
   
   rownames(vipscores) = dimnames(object$coeffs$values)[[1]]
   colnames(vipscores) = dimnames(object$coeffs$values)[[3]]
   
   vipscores   
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
getSelectivityRatio.pls = function(obj, ncomp = NULL, ny = 1, ...)
{
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
      
   selratio = obj$selratio[, ncomp, ny]
   dim(selratio) = c(dim(obj$selratio)[1], 1)
   
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
getVIPScores.pls = function(obj, ny = 1, ...)
{
   vipscores = obj$vipscores[, ny, drop = F]
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
#' @param xlab
#' label for x axies
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @references
#' [1] Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), pp. 103-112.
#' 
plotVIPScores.pls = function(obj, ny = 1, type = 'l', main = NULL, 
                             xlab = 'Variables', ylab = '', ...)
{   
   main = getMainTitle(main, NULL, 'VIP scores')
   
   vipscores = getVIPScores.pls(obj, ny)
   mdaplot(cbind(1:nrow(vipscores), vipscores), type = type, main = main, xlab = xlab,
           ylab = ylab, ...)
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
#' @param xlab
#' label for x axies
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#' 
plotSelectivityRatio.pls = function(obj, ncomp = NULL, ny = 1, type = 'l', main = NULL, 
                                    xlab = 'Variables', ylab = '', ...)
{
   main = getMainTitle(main, ncomp, 'Seelctivity ratio')
   ncomp = getSelectedComponents(obj, ncomp)
   
   selratio = getSelectivityRatio(obj, ncomp, ny)
   mdaplot(cbind(1:nrow(selratio), selratio), type = type, main = main, xlab = xlab,
           ylab = ylab, ...)
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
#' @param show.legend
#' logical, show or not legend for the plot
#' @param show.labels
#' logical, show or not labels for the plot elements
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotRMSE.pls = function(obj, ny = 1, type = 'b',
                        main = 'RMSE', xlab = 'Components', ylab = NULL, 
                        show.legend = T, show.labels = F, ...)
{
   ncomp = obj$ncomp
   
   if (is.null(ylab))
   {   
      if (nrow(obj$calres$rmse) == 1)
         ylab = 'RMSE'
      else if (is.null(rownames(obj$calres$rmse)))
         ylab = sprintf('RMSE (#%d)', ny)
      else         
         ylab = sprintf('RMSE (%s)', rownames(obj$calres$rmse)[ny])
   }
   
   legend = c('cal')
   data = cbind(1:ncomp, obj$calres$rmse[ny, ])   
   labels = mdaplot.formatValues(obj$calres$rmse[ny, ])
   
   if (!is.null(obj$cvres)) 
   { 
      data = cbind(data, obj$cvres$rmse[ny, ])
      labels = cbind(labels, mdaplot.formatValues(obj$cvres$rmse[ny, ]))
      legend = c(legend, 'cv')
   }   
   
   if (!is.null(obj$testres)) 
   { 
      data = cbind(data, obj$testres$rmse[ny, ])
      labels = cbind(labels, mdaplot.formatValues(obj$testres$rmse[ny, ]))
      legend = c(legend, 'test')
   }     
   
   if (show.legend == F)
      legend = NULL
   
   if (show.labels == F)
      labels = NULL
   
   mdaplotg(data, type = type, legend = legend, labels = labels, 
            main = main, xlab = xlab, ylab = ylab, ...)
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
#' @param show.legend
#' logical, show or not legend for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotXVariance.pls = function(obj, type = 'b',
                             main = 'X variance', xlab = 'Components', 
                             ylab = 'Explained variance, %', 
                             show.legend = T, ...)
{
   plotVariance(obj, decomp = 'xdecomp', type = type, main = main,
                xlab = xlab, ylab = ylab, show.legend = show.legend, ...)
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
#' @param show.legend
#' logical, show or not legend for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotYVariance.pls = function(obj, type = 'b',
                             main = 'Y variance', xlab = 'Components', 
                             ylab = 'Explained variance, %', 
                             show.legend = T, ...)
{
   plotVariance(obj, decomp = 'ydecomp', type = type, main = main,
                xlab = xlab, ylab = ylab, show.legend = show.legend, ...)
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
#' @param show.legend
#' logical, show or not legend for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotXCumVariance.pls = function(obj, type = 'b',
                                main = 'X cumulative variance', 
                                xlab = 'Components', ylab = 'Explained variance, %', 
                                show.legend = T, ...)
{
   plotVariance(obj, decomp = 'xdecomp', variance = 'cumexpvar', 
                type = type, main = main, xlab = xlab, ylab = ylab, 
                show.legend = show.legend, ...)   
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
#' @param show.legend
#' logical, show or not legend for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotYCumVariance.pls = function(obj, type = 'b',
                                main = 'Y cumulative variance', 
                                xlab = 'Components', ylab = 'Explained variance, %', 
                                show.legend = T, ...)
{
   plotVariance(obj, decomp = 'ydecomp', variance = 'cumexpvar', 
                type = type, main = main, xlab = xlab, ylab = ylab, 
                show.legend = show.legend, ...)   
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
#' @param show.labels
#' logical, show or not labels for the plot elements
#' @param show.legend
#' logical, show or not legend for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotVariance.pls = function(obj, decomp = 'xdecomp', variance = 'expvar',
                            type = 'b',
                            main = 'X variance', xlab = 'Components', 
                            ylab = 'Explained variance, %', 
                            show.labels = F,
                            show.legend = T, ...)
{
   ncomp = obj$ncomp
   
   legend = c('cal')
   data = cbind(1:ncomp, obj$calres[[decomp]][[variance]])   
   labels = mdaplot.formatValues(obj$calres[[decomp]][[variance]])
   
   if (!is.null(obj$cvres)) 
   { 
      data = cbind(data, obj$cvres[[decomp]][[variance]])
      labels = cbind(labels, mdaplot.formatValues(obj$cvres[[decomp]][[variance]]))
      legend = c(legend, 'cv')
   }   
   
   if (!is.null(obj$testres)) 
   { 
      data = cbind(data, obj$testres[[decomp]][[variance]])
      labels = cbind(labels, mdaplot.formatValues(obj$testres[[decomp]][[variance]]))
      legend = c(legend, 'test')
   }     
   
   if (show.legend == F)
      legend = NULL
   
   if (show.labels == F)
      labels = NULL
   
   mdaplotg(data, legend = legend, type = type, main = main, 
            xlab = xlab, ylab = ylab, labels = labels, ...)
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
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param show.labels
#' logical, show or not labels for the plot elements
#' @param show.legend
#' logical, show or not legend for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotXScores.pls = function(obj, comp = c(1, 2), main = 'X scores',
                           xlab = NULL, ylab = NULL,
                           show.axes = T, 
                           show.labels = F, show.legend = T, ...)
{
   ncomp = length(comp)
   
   if (ncomp < 1 || ncomp > 2)
      stop('The plot can be made for one or two components only!')
      
   if (is.null(main))
      main = sprintf('XY scores (ncomp = %d)', ncomp)
   
   if (is.null(ylab))
   {   
      if (ncomp == 2)
         ylab = colnames(obj$calres$xdecomp$scores)[comp[2]]
      else
         ylab = colnames(obj$calres$xdecomp$scores)[comp]
   }      
   
   if (is.null(xlab))
   {   
      if (ncomp == 2)
         xlab = colnames(obj$calres$xdecomp$scores)[comp[1]]
      else
         xlab = 'Objects'
   }
   
   if (ncomp == 2)
   {   
      cdata = cbind(obj$calres$xdecomp$scores[, comp[1], drop = F], 
                    obj$calres$xdecomp$scores[, comp[2], drop = F])
   }
   else
   {
      ncrows = nrow(obj$calres$xdecomp$scores);
      cdata = cbind(1:ncrows, 
                    obj$calres$xdecomp$scores[, comp, drop = F])
   }   
   
   data = list(cdata = cdata)
   legend = c('cal')
   
   if (!is.null(obj$testres)) 
   { 
      if (ncomp == 2)
      {   
         tdata = cbind(obj$testres$xdecomp$scores[, comp[1], drop = F], 
                    obj$testres$xdecomp$scores[, comp[2], drop = F])
      }
      else
      {
         tdata = cbind(ncrows + (1:nrow(obj$testres$xdecomp$scores)), 
                       obj$testres$xdecomp$scores[, comp, drop = F])         
      }   
      data$tdata = tdata
      legend = c(legend, 'test')      
   }   
   
   if (show.legend == F)
      legend = NULL
   
   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, legend = legend, show.labels = show.labels, show.lines = show.lines,
            main = main, xlab = xlab, ylab = ylab, ...)   
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
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param show.labels
#' logical, show or not labels for the plot elements
#' @param show.legend
#' logical, show or not legend for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
plotXYScores.pls = function(obj, comp = 1, main = NULL,
                            xlab = 'X scores', ylab = 'Y scores',
                            show.axes = T,
                            show.labels = F, show.legend = T, ...)
{
   if (is.null(comp)) 
      comp = obj$ncomp.selected   
   else if (comp <= 0 || comp > obj$ncomp) 
      stop('Wrong component number!')
   
   if (is.null(main))
      main = sprintf('XY scores (comp %d)', comp)
   
   cdata = cbind(obj$calres$xdecomp$scores[, comp, drop = F], 
                 obj$calres$ydecomp$scores[, comp, drop = F])
   colnames(cdata) = c('X scores', 'Y scores')
   data = list(cdata = cdata)
   legend = c('cal')
   
   if (!is.null(obj$testres)) 
   { 
      tdata = cbind(obj$testres$xdecomp$scores[, comp, drop = F], 
                    obj$testres$ydecomp$scores[, comp, drop = F])
      data$tdata = tdata
      legend = c(legend, 'test')      
   }   
   
   if (show.legend == F)
      legend = NULL
   
   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, legend = legend, show.labels = show.labels, main = main,
            xlab = xlab, ylab = ylab, show.lines = show.lines, ...)   
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
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param legend.position
#' position of legend on the plot (if shown)
#' @param show.line
#' logical, show or not line fit for the plot points
#' @param show.labels
#' logical, show or not labels for the plot elements
#' @param show.legend
#' logical, show or not legend for the plot
#' @param colmap
#' a colormap to use for coloring the plot items
#' @param col
#' a vector with color values for plot items
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotPredictions.pls = function(obj, ncomp = NULL, ny = 1, main = NULL, xlab = NULL,
                               ylab = NULL, legend.position = 'topleft', show.line = T, show.labels = F,
                               show.legend = T, colmap = 'default', col = NULL, ...)
{
   # set default values for main title, x and y axis labels
   if (is.null(main))
   {   
      if (is.null(ncomp))
         main = 'Predictions'
      else
         main = sprintf('Predictions (ncomp = %d)', ncomp)
   }
   
   if (is.null(xlab))
   {
      if (!is.null(colnames(obj$calres$y.ref)))
         xlab = sprintf('%s, measured', colnames(obj$calres$y.ref)[ny])
      else
         xlab = 'y, measured'
   }  
   
   if (is.null(ylab))
   {
      if (!is.null(colnames(obj$calres$y.ref)))
         ylab = sprintf('%s, predicted', colnames(obj$calres$y.ref)[ny])
      else
         ylab = 'y, predicted'
   }  
   
   # check number of components
   if (is.null(ncomp)) 
      ncomp = obj$ncomp.selected
   else if (ncomp <= 0 || ncomp > obj$ncomp) 
      stop('Wrong value for number of components!')
   
   # make list with groups of measured and predicted values for each type of result
   cdata = cbind(obj$calres$y.ref[, ny], obj$calres$y.pred[, ncomp, ny])
   colnames(cdata) = c('y, measured', 'y, predicted')
   rownames(cdata) = rownames(obj$calres$y.ref)
   legend = c('cal')
   data = list(cdata = cdata)
   
   if (!is.null(obj$cvres)) 
   { 
      cvdata = cbind(obj$cvres$y.ref[, ny], obj$cvres$y.pred[, ncomp, ny])
      colnames(cvdata) = c('y, measured', 'y, predicted')
      rownames(cvdata) = rownames(obj$cvres$y.ref)
      legend = c(legend, 'cv')
      data$cvdata = cvdata
   }   
   
   if (!is.null(obj$testres)) 
   { 
      testdata = cbind(obj$testres$y.ref[, ny], obj$testres$y.pred[, ncomp, ny])
      colnames(testdata) = c('y, measured', 'y, predicted')
      rownames(testdata) = rownames(obj$testres$y.ref)
      legend = c(legend, 'test')
      data$testdata = testdata
   }   
   
   if (show.legend == F)
      legend = NULL
   
   mdaplotg(data, legend = legend, show.labels = show.labels, main = main, 
            colmap = colmap, xlab = xlab, ylab = ylab, col = col, 
            legend.position = legend.position, ...)  
   
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
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.line
#' logical, show or not line for 0 value
#' @param show.labels
#' logical, show or not labels for the plot elements
#' @param show.legend
#' logical, show or not legend for the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotYResiduals.pls = function(obj, ncomp = NULL, ny = 1, type = 'p', 
                              main = NULL, 
                              ylab = NULL,
                              xlab = NULL,
                              show.line = T,
                              show.labels = F,
                              show.legend = T,
                              ...)
{
   if (is.null(main))
   {   
      if (is.null(ncomp))
         main = 'Y residuals'
      else
         main = sprintf('Y residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(xlab))
   {
      if (type == 'p')
      {   
         if (ncol(obj$calres$y.ref) == 1)
            xlab = 'y, measured'
         else
            xlab = sprintf('%s, measured', colnames(obj$calres$y.ref)[ny])
      }
      else
      {
         xlab = 'Objects'
      }   
   }  

   if (is.null(ylab))
   {
      if (ncol(obj$calres$y.ref) == 1)
         ylab = 'y residuals'
      else
         ylab = sprintf('y residuals (%s)', colnames(obj$calres$y.ref)[ny])
   }  
   
   if (is.null(ncomp)) 
      ncomp = obj$ncomp.selected
   else if (ncomp <= 0 || ncomp > obj$ncomp) 
      stop('Wrong value for number of components!')
   
   if (type == 'p')
   {   
      cdata = cbind(obj$calres$y.ref[, ny], obj$calres$y.pred[, ncomp, ny] - obj$calres$y.ref[, ny])
      colnames(cdata) = c('y values', 'y residuals')
   }
   else
   {   
      nobj = nrow(obj$calres$y.ref)
      cdata = cbind(1:nrow(obj$calres$y.ref), obj$calres$y.pred[, ncomp, ny] - obj$calres$y.ref[, ny])
      colnames(cdata) = c('Objects', 'y residuals')
   }
   rownames(cdata) = rownames(obj$calres$y.pred)
   legend = c('cal')
   data = list(cdata = cdata)
   
   if (!is.null(obj$cvres)) 
   { 
      if (type == 'p')
      {   
         cvdata = cbind(obj$cvres$y.ref[, ny], 
                        obj$cvres$y.pred[, ncomp, ny] - obj$cvres$y.ref[, ny])
         colnames(cvdata) = c('y values', 'y residuals')
      }
      else
      {   
         cvdata = cbind(nobj + (1:nrow(obj$cvres$y.ref)), 
                        obj$cvres$y.pred[, ncomp, ny] - obj$cvres$y.ref[, ny])
         colnames(cvdata) = c('Objects', 'y residuals')
         nobj = nobj + nrow(obj$cvres$y.ref)
      }      
      rownames(cvdata) = rownames(obj$cvres$y.pred)
      legend = c(legend, 'cv')
      data$cvdata = cvdata
   }   
   
   if (!is.null(obj$testres)) 
   { 
      if (type == 'p')
      {   
         tdata = cbind(obj$testres$y.ref[, ny], 
                       obj$testres$y.pred[, ncomp, ny] - obj$testres$y.ref[, ny])
         colnames(tdata) = c('y values', 'y residuals')
      }
      else
      {   
         tdata = cbind(nobj + (1:nrow(obj$testres$y.ref)), 
                       obj$testres$y.pred[, ncomp, ny] - obj$testres$y.ref[, ny])
         colnames(tdata) = c('Objects', 'y residuals')
      }
      rownames(tdata) = rownames(obj$testres$y.pred)
      legend = c(legend, 'test')
      data$tdata = tdata
   }   
   
   if (show.legend == F)
      legend = NULL
   
   if (show.line == T)
      show.line = c(NA, 0)
   
   mdaplotg(data, legend = legend, type = type, show.labels = show.labels, xlab = xlab, ylab = ylab, 
            main = main, show.lines = show.line, ...)      
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
#' @param ny
#' number of response variable to make the plot for (if y is multivariate)
#' @param main
#' main plot title
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotRegcoeffs.pls = function(obj, ncomp = NULL, ny = 1, main = NULL, ylab = NULL, ...)
{
   if (is.null(main))
   {   
      if (is.null(ncomp))
         main = 'Regression coefficients'
      else
         main = sprintf('Regression coefficients (ncomp = %d)', ncomp)
   }
   
   if (is.null(ylab))
   {   
      if (ncol(obj$calres$y.ref) == 1)
         ylab = 'Coefficients'
      else
         ylab = sprintf('Coefficients (%s)', colnames(obj$calres$y.ref)[ny])
   }
   
   if (is.null(ncomp)) 
      ncomp = obj$ncomp.selected
   else if (ncomp <= 0 || ncomp > obj$ncomp) 
      stop('Wrong value for number of components!')
   
   plot(obj$coeffs, ncomp = ncomp, ny = ny, main = main, ylab = ylab, ...)
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
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotXLoadings.pls = function(obj, comp = c(1, 2), type = 'p', main = 'X loadings', 
                             ylab = NULL, xlab = NULL, show.axes = T, ...)
{
   ncomp = length(comp)
   
   if (min(comp) < 1 || max(comp) > obj$ncomp)
      stop('Wrong components number!')
   
   if (type == 'p' && ncomp != 2)
      stop('Scatter plot can be made only for two components!')

   if (type == 'p')
   {
      data = cbind(obj$xloadings[, comp[1], drop = F],
                   obj$xloadings[, comp[2], drop = F])      
      legend = NULL
   }   
   else
   {
      data = cbind(1:nrow(obj$xloadings),
                   obj$xloadings[, comp, drop = F])  
      legend = colnames(obj$xloadings[, comp, drop = F])
   }   
   
   if (is.null(ylab))
   {   
      if (type == 'p')
         ylab = colnames(obj$xloadings)[comp[2]]
      else   
         ylab = 'Loadings'
   }

   if (is.null(xlab))
   {   
      if (type == 'p')
         xlab = colnames(obj$xloadings)[comp[1]]
      else   
         xlab = 'Variables'
   }

   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, main = main, type = type, ylab = ylab, xlab = xlab, legend = legend, 
            show.lines = show.lines, ...)
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
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotXYLoadings.pls = function(obj, comp = c(1, 2), main = 'XY loadings', 
                             ylab = NULL, xlab = NULL, show.axes = F, ...)
{
   if (length(comp) != 2)
      stop('This plot can be made for only two components!')
   
   xdata = cbind(obj$xloadings[, comp[1], drop = F],
                 obj$xloadings[, comp[2], drop = F])      
   ydata = cbind(obj$yloadings[, comp[1], drop = F],
                 obj$yloadings[, comp[2], drop = F])      
   data = list(xdata = xdata, ydata = ydata)
   
   legend = c('X', 'Y')
   
   if (is.null(ylab))
      ylab = colnames(obj$xloadings)[comp[2]]
   
   if (is.null(xlab))
      xlab = colnames(obj$xloadings)[comp[1]]
   
   if (show.axes == T)
      show.lines = c(0, 0)
   else
      show.lines = F
   
   mdaplotg(data, main = main, type = 'p', ylab = ylab, xlab = xlab, legend = legend, 
            show.lines = show.lines, ...)
}


#' X residuals plot for PLS
#' 
#' @description
#' Shows a plot with Q2 residuals vs. Hotelling T2 values for PLS decomposition of x data.
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
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{pls}} function.
#' 
plotXResiduals.pls = function(obj, ncomp = NULL, 
                              main = NULL, xlab = 'T2', ylab = 'Q2',
                              show.labels = F, show.legend = T, ...)
{
   if (is.null(main))
   {   
      if (is.null(ncomp))
         main = 'X residuals'
      else
         main = sprintf('X residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp)) 
      ncomp = obj$ncomp.selected
   else if (ncomp <= 0 || ncomp > obj$ncomp) 
      stop('Wrong value for number of components!')
   
   cdata = cbind(obj$calres$xdecomp$T2[, ncomp], obj$calres$xdecomp$Q2[, ncomp])
   
   colnames(cdata) = c('T2', 'Q2')
   rownames(cdata) = rownames(obj$calres$xdecomp$scores)
   
   data = list(cdata = cdata)
   legend = NULL
   
   if (!is.null(obj$testres))
   {
      tdata = cbind(obj$testres$xdecomp$T2[, ncomp], obj$testres$xdecomp$Q2[, ncomp])      
      colnames(tdata) = c('T2', 'Q2')
      rownames(tdata) = rownames(obj$testres$xdecomp$scores)
      
      data$tdata = tdata
      legend = c('cal', 'test')
   }      
   
   mdaplotg(data, main = main, xlab = xlab, ylab = ylab,
            show.labels = show.labels, legend = legend, ...)
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
plot.pls = function(x, ncomp = NULL, ny = 1, show.legend = T, show.labels = F, ...)
{
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
#' @method summary pls
#' @S3method summary pls
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
summary.pls = function(object, ncomp = NULL, ny = NULL, ...)
{
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
   
   for (y in ny)
   {   
      if (ncol(obj$calres$y.ref) > 1)
         cat(sprintf('\nResponse variable #%d (%s)\n', y, colnames(obj$calres$y.ref)[y]))
      
      data = as.matrix(obj$calres, ncomp = ncomp, ny = y)
      rownames(data) = 'Cal'
      
      if (!is.null(obj$cvres))
      {
         data = rbind(data, as.matrix(obj$cvres, ncomp = ncomp, ny = y))      
         rownames(data)[2] = 'CV'
      }
      
      if (!is.null(obj$testres))
      {
         data = rbind(data, as.matrix(obj$testres, ncomp = ncomp, ny = y))
         rownames(data)[nrow(data)] = 'Test'
      }   
      
      data = data[, -c(1, 3), drop = F]
      
      data[, 1:2] = round(data[, 1:2], 2)      
      data[, 4:5] = round(data[, 4:5], 2)  
      data[, 3] = mdaplot.formatValues(data[, 3], round.only = T)
      data[, 6] = round(data[, 6], 4)      
      data[, 7] = round(data[, 7], 2)      
      
      print(data)
   }   
   cat('\n')
}

#' Print method for PLS model object
#' 
#' @method print pls
#' @S3method print pls
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a PLS model (object of class \code{pls})
#' @param ...
#' other arguments
#' 
print.pls = function(x, ...)
{
   obj = x
   
   cat('\nPLS model (class pls)\n')
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
