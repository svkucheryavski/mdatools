## class and methods for linear decomposition X = TP' + E ##

#' Linear decomposition of data
#' 
#' @description
#' Creates an object of ldecomp class.
#'
#' @param scores
#' matrix with score values (nobj x ncomp).
#' @param residuals
#' matrix with data residuals 
#' @param loadings
#' matrix with loading values (nvar x ncomp).
#' @param ncomp.selected
#' number of selected components
#' @param attrs
#' list with attributes of original dataset
#' @param totvar
#' full variance of original data, preprocessed and centered
#' @param dist
#' list with calculated T2 and Q values (e.g. for CV)
#' @param var
#' list with explained and cumulative explained variance (e.g. for CV)
#' @param cal
#' logical, true if data is for calibration of a LDECOMP based model
#' @param tnorm
#' singular values for score normalization
#'
#' @return
#' Returns an object (list) of \code{ldecomp} class with following fields:
#' \item{scores }{matrix with score values (nobj x ncomp).}
#' \item{residuals }{matrix with data residuals (nobj x nvar).}
#' \item{T2 }{matrix with T2 distances (nobj x ncomp).}
#' \item{Q }{matrix with Q statistic (nobj x ncomp).}
#' \item{tnorm }{vector with singular values used for scores normalization.}
#' \item{ncomp.selected }{selected number of components.}
#' \item{expvar }{explained variance for each component.}
#' \item{cumexpvar }{cumulative explained variance.}
#' \item{modpower}{modelling power of variables.}
#'
#' @details
#' \code{ldecomp} is a general class for decomposition X = TP' + E. Here, X is a data matrix, 
#' T - matrix with scores, P - matrix with loadings and E - matrix with residuals. It is used, 
#' for example, for PCA results (\code{\link{pcares}}), in PLS and other methods. The class also 
#' includes methods for calculation and plotting residuals, variances, and so on.
#'
#' There is no need to use the \code{ldecomp} manually. For example, when build PCA model 
#' with \code{\link{pca}} or apply it to a new data, the results will automatically inherit 
#' all methods of \code{ldecomp}.
#'
#' @importFrom methods show
#' @importFrom stats convolve cor lm na.exclude predict pt qf qnorm qt sd var
#'
#' @export
ldecomp = function(scores = NULL, residuals = NULL, loadings = NULL, ncomp.selected = NULL, 
                   attrs = NULL, tnorm = NULL, dist = NULL, var = NULL, cal = FALSE, 
                   totvar = NULL) {
  
   if (!is.null(scores) && !is.null(loadings)) {
      # there are scores and loadings, calculate the rest and add attributes
      ncomp = ncol(loadings)
      scores = mda.setattr(scores, attrs, type = 'row')
      residuals = mda.setattr(residuals, attrs)

      rownames(scores) = rownames(residuals) = attrs$dimnames[[1]]
      colnames(scores) = colnames(loadings)
      colnames(residuals) = attrs$dimnames[[2]]
      attr(scores, 'name') = 'Scores'
      attr(scores, 'xaxis.name') = 'Components'
      attr(residuals, 'name') = 'Residuals'
      
   } else {
      ncomp = ncol(dist$Q)
   }

   # calculate residual distances
   if (is.null(dist))
      dist = ldecomp.getDistances(scores, loadings, residuals, tnorm, cal = cal) 
      
   # set names and attributes for common results 
   colnames(dist$Q) = colnames(dist$T2) = colnames(loadings)
   rownames(dist$Q) = rownames(dist$T2) = attrs$dimnames[[1]]
   if (!is.null(dist$modpower)) {
      # set attributes for modpower and exclude rows if necessary
      rownames(dist$modpower) = rownames(loadings)
      colnames(dist$modpower) = colnames(loadings)
      attr(dist$modpower, 'name') = 'Modelling power'
      attr(dist$modpower, 'xaxis.name') = 'Components'
      attr(dist$modpower, 'yaxis.name') = attr(loadings, 'yaxis.name')
      attr(dist$modpower, 'yaxis.values') = attr(loadings, 'yaxis.values')
      if (length(attrs$exclcols > 0))
         dist$modpower = mda.exclrows(dist$modpower, attrs$exclcols)
   }     
      
   # set attributes for Q
   dist$Q = mda.setattr(dist$Q, attrs, type = 'row')
   attr(dist$Q, 'name') = 'Squared residual distance (Q)'
   attr(dist$Q, 'xaxis.name') = 'Components'
      
   # set attributes for T2 
   dist$T2 = mda.setattr(dist$T2, mda.getattr(dist$Q))
   attr(dist$T2, 'name') = 'T2 residuals'
  
   # calculate explained variance
   if (is.null(var))
      var = ldecomp.getVariances(dist$Q, totvar) 
   
   obj = list(
      scores = scores,
      residuals = residuals,
      ncomp = ncomp,
      ncomp.selected = ncomp.selected,
      T2 = dist$T2,
      Q = dist$Q,
      tnorm = dist$tnorm,
      modpower = dist$modpower,
      cumexpvar = var$cumexpvar,
      expvar = var$expvar,
      totvar = totvar,
      Qlim = NULL,
      T2lim = NULL,
      lim.type = NULL
   )
   
   if (is.null(ncomp.selected))
      obj$ncomp.selected = ncol(scores)
   else
      obj$ncomp.selected = ncomp.selected
   
   obj$call = match.call()
   class(obj) = "ldecomp"
   obj
}

#' Residuals distances for linear decomposition
#'
#' @description
#' Computes residual distances (Q and T2) and modelling power for a data decomposition X = TP' + E.
#' 
#' @param scores
#' matrix with scores (T).
#' @param loadings
#' matrix with loadings (P).
#' @param residuals
#' matrix with residuals (E).
#' @param tnorm
#' vector with singular values for scores normalisation
#' @param cal
#' if TRUE method will realize that these distances are calculated for calibration set
#' 
#' @details
#' The distances are calculated for every 1:n components, where n goes from 1 to ncomp 
#' (number of columns in scores and loadings). 
#' 
#' @return
#' Returns a list with Q, Qvar, T2 and modelling power values for each component.
#'  
ldecomp.getDistances = function(scores, loadings, residuals, tnorm = NULL, cal = FALSE) {
   # get attributes
   attrs.scores = mda.getattr(scores)
   attrs.loadings = mda.getattr(loadings)
   
   # get sizes
   ncomp = ncol(scores)
   nobj = nrow(scores)
   nvar = nrow(loadings)
   
   # prepare zero matrices for the distances
   T2 = matrix(0, nrow = nobj, ncol = ncomp)
   Q = matrix(0, nrow = nobj, ncol = ncomp)
   
   # calculate normalized scores
   if (is.null(tnorm) && nrow(scores) > 0)
      tnorm = sqrt(colSums(scores^2)/(nrow(scores) - 1));   
   
   scoresn = sweep(scores, 2L, tnorm, '/', check.margin = F);  

   # calculate variance for data columns
   data = tcrossprod(scores, loadings) + residuals;
  
   # correct data and loadings for excluded columns
   varnames = rownames(loadings)
   if (length(attrs.loadings$exclrows) > 0){
      loadings = loadings[-attrs.loadings$exclrows, , drop = F]
      data = data[, -attrs.loadings$exclrows, drop = F]
   }
   
   modpower = matrix(0, nrow = nrow(loadings), ncol = ncomp)
   
   # standard deviation for data values 
   if (nobj > 1 && cal == TRUE)
      datasd = sqrt(colSums(data^2)/(nobj - 1))
   
   # calculate distances for each set of components
   for (i in 1:ncomp) {
      exp = tcrossprod(scores[, 1:i, drop = F], loadings[, 1:i, drop = F]);
      res = data - exp;

      Q[, i] = rowSums(res^2)
      T2[, i] = rowSums(scoresn[, 1:i, drop = F]^2)
      
      if (nobj > i && cal == TRUE){
         modpower[, i] = 1 - sqrt(colSums(res^2)/(nobj - i - 1))/datasd
      }
   }  
   
   if (length(attrs.loadings$exclrows) > 0){
      out.modpower = matrix(0, nrow = nvar, ncol = ncomp)
      out.modpower[-attrs.loadings$exclrows, ] = modpower
   } else {
      out.modpower = modpower
   }
   
   # return the results
   res = list(
      Q = Q,
      T2 = T2,
      modpower = out.modpower,
      tnorm = tnorm
   )
   
   res
}


#' Explained variance for linear decomposition
#' 
#' @description
#' Computes explained variance and cumulative explained variance for a data decomposition 
#' X = TP' + E.
#'
#' @param Q
#' Q values (squared residuals distance from object to component space).
#' @param totvar
#' Total variance of the original data (after preprocessing).
#' 
#' @return
#' Returns a list with two vectors.
#' 
ldecomp.getVariances = function(Q, totvar) { 
   if (length(attr(Q, 'exclrows')) > 0)
      Q = Q[-attr(Q, 'exclrows'), , drop = F]
   cumresvar = colSums(Q) / totvar * 100
   cumexpvar = 100 - cumresvar
   expvar = c(cumexpvar[1], diff(cumexpvar))
  
   res = list(
      expvar = expvar,
      cumexpvar = cumexpvar
   )
   
   res
}

#' Calculates critical limits for T2-residuals using Hotelling T2 distribution
#' 
#' @description 
#' The method is based on n
#' 
#' @param ncomp 
#' number of components
#' @param T2
#' vector with T2-residuals for selected component
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#' @param T2lim
#' T2 limits from a PCA model (needed for probabilities)
#' @param return
#' what to return: \code{'limits'} or \code{'probability'}
#' 
#' @export
reslim.hotelling = function(ncomp, T2 = NULL, alpha = 0.05, gamma = 0.01, 
                            T2lim = NULL, return = 'limits') {
   if (return == 'limits') {
      nobj = length(T2)
      out = rep(0, 4)
      if (nobj > ncomp) {
         out[1:2] = 
            (ncomp * (nobj - 1) / (nobj - ncomp)) * qf(c(1 - alpha, 1 - gamma), ncomp, nobj - ncomp)
      } else {
         out[1:2] = 0
      }
      out[3] = mean(T2)
      out[4] = nobj - ncomp
   } else {
      nobj = T2lim[4] + ncomp
      out = pf(T2 * (nobj - ncomp) / (ncomp * (nobj - 1)), ncomp, nobj - ncomp)
   }
   
   out
}

#' Calculates critical limits or statistic values for Q-residuals using Chi-squared distribution
#' 
#' @description 
#' The method is based on Chi-squared distribution with DF = 2 * (m(Q)/s(Q)^2
#' 
#' @param Q
#' vector with Q-residuals for selected component
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#' @param Qlim
#' vector with Q limits for selected number of components (from model)
#' @param return
#' what to return: \code{'limits'} or \code{'probability'}
#' 
#' @export
reslim.chisq = function(Q, alpha = 0.05, gamma = 0.01, Qlim = NULL, return = 'limits') {
   nobj = length(Q)
   if (return == 'limits') {
      
      # calculate mean and DF for Q
      Q.mean = mean(Q)
      Q.DF = 2 * (Q.mean/sd(Q))^2      
      
      out = rep(0, 4)
      if (Q.mean == 0) {
         out = c(0, 0, 0, 1)
      } else {
         out[1] = qchisq(1 - alpha, floor(Q.DF)) * Q.mean / Q.DF
         out[2] = qchisq((1 - gamma)^(1/nobj), floor(Q.DF)) * Q.mean / Q.DF
         out[3] = Q.mean
         out[4] = Q.DF
      }
   } else {
      # get mean and DF from model
      Q.mean = Qlim[3]
      Q.DF = Qlim[4]
      out = pchisq(Q.DF * Q / Q.mean, floor(Q.DF))
   }

   out      
}

#' Calculates critical limits for Q-residuals using classic JM approach
#' 
#' @description 
#' The method is based on 
#' 
#' @param eigenvals
#' vector with eigenvalues for all variables
#' @param Q
#' vector with Q-residuals for selected component
#' @param ncomp 
#' number of components
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#' @param return
#' what to return: \code{'limits'} or \code{'probability'}
#' 
#' @export
reslim.jm = function(eigenvals, Q, ncomp, alpha = 0.05, gamma = 0.01, return = 'limits') {
   evals = eigenvals[(ncomp + 1):length(eigenvals)]   
   t1 = sum(evals)
   t2 = sum(evals^2)
   t3 = sum(evals^3)
   h0 = 1 - 2 * t1 * t3/3/(t2^2);
   if (h0 < 0.001) h0 = 0.001
   
   if (return == 'limits') {
      # inverse error function 
      erfinv = function (x) qnorm((1 + x)/2)/sqrt(2)
      
      out = rep(0, 4)
      ca = sqrt(2) * erfinv(c(1 - 2 * alpha, 1 - 2 * gamma))
      h1 = ca * sqrt(2 * t2 * h0^2)/t1
      h2 = t2 * h0 * (h0 - 1)/(t1^2)
      out[1:2] = t1 * (1 + h1 + h2)^(1/h0)
      out[3] = mean(Q)
      out[4] = 1
   } else {
      # error function 
      erf = function (x) 1 - pnorm(-x * sqrt(2)) * 2
      
      h1 = (Q/t1)^h0
      h2 = t2 * h0 * (h0 - 1)/t1^2
      d = t1 * (h1 - 1 - h2)/(sqrt(2 * t2) * h0)
      out = 0.5 * (1+erf(d/sqrt(2)))
   }
   
   out
}

#' Statistical limits for Q and T2 residuals using Data Driven approach
#' 
#' @description 
#' Method is based on paper by Pomerantsev, Rodionova (JChem, 2014)
#' 
#' @param Q
#' vector with Q-residuals for selected component
#' @param T2
#' vector with T2-residuals for selected component
#' @param type
#' which estimator to use: 'moments' or 'robust'
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#' @param Qlim
#' vector with Q limits for selected number of components (from model)
#' @param T2lim
#' vector with T2 limits for selected number of components (from model)
#' @param return
#' what to return: \code{'limits'} or \code{'probability'}
#' 
#' @import stats
#' 
#' @export
reslim.dd = function(Q, T2, type = 'ddmoments', alpha = 0.05, gamma = 0.01, Qlim = NULL, 
                     T2lim = NULL, return = 'limits') {
   
   # classic estimator of u0 and DF for DD method
   ddclassic = function(u) {
      u0 = mean(u)
      DF = round(2 * (u0/sd(u))^2)
      if (DF == 1) DF = 1
      res = c(u0, ifelse(DF > 250, 250, DF))
      res
   }

   # robust estimator of u0 and DF for DD method
   ddrobust = function(u) {
      M = median(u)
      R = quantile(u, 0.75) - quantile(u, 0.25)
      DF = R/M
      if (DF > 2.685592117) {
         DF = 1
      } else if (DF < 0.194565995) {
         DF = 100
      } else {
         DF = exp(1.380948 * log(2.68631 / DF))^1.185785
      }
      u0 = 0.5 * DF * (M/qchisq(0.5, DF) + R/(qchisq(0.75, DF) - qchisq(0.25, DF)))
      res = c(u0, ifelse(DF > 250, 250, DF))
      res
   }

   
   if (return == 'limits') {
      
      # estimate mean and DF
      if (type == 'ddmoments') {
         T2.res = ddclassic(T2)
         Q.res = ddclassic(Q)
      } else if (type == "ddrobust"){
         T2.res = ddrobust(T2)
         Q.res = ddrobust(Q)
      } else {
         stop('Wrong value for "type" parameter!')
      }
      
      T2.mean = T2.res[1]; T2.DF = T2.res[2]
      Q.mean = Q.res[1]; Q.DF = Q.res[2]
      
      nobj = length(Q)

      # calculate critical values
      Ocrit = qchisq((1 - gamma)^(1/nobj), T2.DF + Q.DF)
      Dcrit = qchisq(1 - alpha, T2.DF + Q.DF)
   
      # calculate critical limits
      Qlim = c(Dcrit * Q.mean / Q.DF, Ocrit * Q.mean / Q.DF, Q.mean, round(Q.DF))
      T2lim = c(rep(-(Q.mean / T2.mean) * (T2.DF / Q.DF), 2), T2.mean, round(T2.DF))
      out = list(Qlim = Qlim, T2lim = T2lim)
   } else {
      # get mean and DF from matrix with limits from model      
      Q.mean = Qlim[3]; Q.DF = Qlim[4]
      T2.mean = T2lim[3]; T2.DF = T2lim[4]
      out = pchisq(Q / Q.mean * Q.DF + T2 / T2.mean * T2.DF, round(Q.DF) + round(T2.DF))
   }
   
   out
}

#' Shows lines with critical limits on residuals plot
#' 
#' @description 
#' Shows lines with critical limits on residuals plot
#' 
#' @param lim
#' matrix with residual limits (2x2)
#' @param lim.type
#' type of limits
#' @param lim.col
#' vector with two values - line color for extreme and outlier borders 
#' @param lim.lwd
#' vector with two values - line width for extreme and outlier borders 
#' @param lim.lty
#' vector with two values - line type for extreme and outlier borders 
#' 
#' @export
ldecomp.plotLimits = function(lim, lim.type, lim.col, lim.lwd, lim.lty) {
   if (substr(lim.type, 1, 2) != 'dd') {
      # rectangular limits
      lines(c(0, lim[1, 1]), c(lim[1, 2], lim[1, 2]), lty = lim.lty[1], lwd = lim.lwd[1], 
            col = lim.col[1])
      lines(c(lim[1, 1], lim[1, 1]), c(0, lim[1, 2]), lty = lim.lty[1], lwd = lim.lwd[1], 
            col = lim.col[1])
      lines(c(0, lim[2, 1]), c(lim[2, 2], lim[2, 2]), lty = lim.lty[2], lwd = lim.lwd[2], 
            col = lim.col[2])
      lines(c(lim[2, 1], lim[2, 1]), c(0, lim[2, 2]), lty = lim.lty[2], lwd = lim.lwd[2], 
            col = lim.col[2])
   } else {
      # triangular limits
      abline(a = lim[1, 2], b = lim[1, 1], lty = lim.lty[1], lwd = lim.lwd[1], col = lim.col[1])
      abline(a = lim[2, 2], b = lim[2, 1], lty = lim.lty[2], lwd = lim.lwd[2], col = lim.col[2])
   }
}

#' Cumulative explained variance plot for linear decomposition
#' 
#' @description
#' Shows a plot with cumulative explained variance values vs. number of components.
#' 
#' @param obj
#' object of \code{ldecomp} class.
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param labels
#' what to show as labels for plot objects
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
#' @export
plotCumVariance.ldecomp = function(obj, type = 'b', main = 'Cumulative variance',
                                   xlab = 'Components', ylab = 'Explained variance, %',
                                   show.labels = F, labels = 'values', ...) {
   mdaplot(obj$cumexpvar, main = main, xticks = 1:obj$ncomp, xlab = xlab, ylab = ylab, type = type, 
           show.labels = show.labels, labels = labels, ...)
}

#' Explained variance plot for linear decomposition
#' 
#' @description
#' Shows a plot with explained variance values vs. number of components.
#' 
#' @param obj
#' object of \code{ldecomp} class.
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for plot objects.
#' @param labels
#' what to show as labels for plot objects.
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
#' @export
plotVariance.ldecomp = function(obj, type = 'b', main = 'Variance',
                                xlab = 'Components', ylab = 'Explained variance, %',
                                show.labels = F, labels = 'values', ...) {
   mdaplot(obj$expvar, main = main, xticks = 1:obj$ncomp, xlab = xlab, ylab = ylab, 
           show.labels = show.labels, labels = labels, type = type, ...)
}


#' Scores plot for linear decomposition
#' 
#' @description
#' Shows a plot with scores values for data objects.
#' 
#' @param obj
#' object of \code{ldecomp} class.
#' @param comp
#' which components to show the plot for (can be one value or vector with two values).
#' @param main
#' main title for the plot
#' @param type
#' type of the plot
#' @param xlab
#' label for x-axis.
#' @param ylab
#' label for y-axis.
#' @param show.legend
#' logical, show or not a legend on the plot.
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
#' @export
plotScores.ldecomp = function(obj, comp = c(1, 2), main = 'Scores', 
                              type = 'p', xlab = NULL, ylab = NULL,
                              show.labels = FALSE, show.legend = TRUE, 
                              show.axes = TRUE, ...) {
   if (is.null(obj$scores)) {
      warning('Scores values are not specified!')
   } else {   
      data = mda.subset(obj$scores, select = comp)
      colnames(data) = paste('Comp ', comp, ' (', round(obj$expvar[comp], 2) , '%)', sep = '')
      
      if (show.axes == T) {
         if (type == 'p') {
            if (length(comp) > 1)
               show.lines = c(0, 0)
            else
               show.lines = c(NA, 0)
         } else if (type != 'h') {
            show.lines = c(NA, 0)      
         } else {
            show.lines = FALSE
         }
      } else {
         show.lines = FALSE
      }
     
      if (is.null(xlab) && length(comp) == 1)
         xlab = ifelse(is.null(attr(data, 'yaxis.name')), 'Objects', attr(data, 'yaxis.name'))
      
      if (type != 'p') {
         data = mda.t(data)
        
         if (nrow(data) == 1) {
            if (is.null(ylab))
               ylab = rownames(data)[1]
            mdaplot(data, main = main, type = type, show.labels = show.labels, 
                    show.lines = show.lines, xlab = xlab, ylab = ylab, ...)
         } else {
            if (is.null(show.legend))
               show.legend = TRUE
            
            mdaplotg(data, main = main, type = type, show.labels = show.labels, 
                     show.lines = show.lines, xlab = xlab, ylab = ylab, 
                     show.legend = show.legend, ...)
         }
      } else {
         mdaplot(data, main = main, type = type, show.labels = show.labels, show.lines = show.lines, 
                 xlab = xlab, ylab = ylab, ...)
      }
   }  
}
#' Residuals plot for linear decomposition
#' 
#' @description
#' Shows a plot with T2 vs Q values for data objects.
#' 
#' @param obj
#' object of \code{ldecomp} class.
#' @param ncomp
#' what number of components to show the plot for (if NULL, model selected value will be used).
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
#' @export
plotResiduals.ldecomp = function(obj, ncomp = NULL, main = NULL, xlab = NULL, ylab = NULL, 
                                 show.labels = F, ...) {
   if (is.null(main)) {
      if (is.null(ncomp))
         main = 'Residuals'
      else
         main = sprintf('Residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   
   if (is.null(xlab))
      xlab = expression(paste('Hotelling ', T^2, ' distance'))
   if (is.null(ylab))
      ylab = 'Squared residual distance, Q'      
   
   data = mda.cbind(
      mda.subset(obj$T2, select = ncomp), 
      mda.subset(obj$Q, select = ncomp)
   )
   
   # show plot
   mdaplot(data, main = main, xlab = xlab, ylab = ylab, show.labels = show.labels, ...)
}  

#' Print method for linear decomposition
#'
#' @description
#' Generic \code{print} function for linear decomposition. Prints information about 
#' the \code{ldecomp} object.
#' 
#' @param x
#' object of class \code{ldecomp}
#' @param str
#' user specified text to show as a description of the object
#' @param ...
#' other arguments
#' 
#' @export
print.ldecomp = function(x, str = NULL, ...) {   
   if (is.null(str))
      str ='Results of data decomposition (class ldecomp)'
   
   if (nchar(str) > 0)   
      cat(sprintf('\n%s\n', str))
   
   cat('\nMajor fields:\n')   
   cat('$scores - matrix with score values\n')
   cat('$T2 - matrix with T2 distances\n')
   cat('$Q - matrix with Q residuals\n')
   cat('$ncomp.selected - selected number of components\n')
   cat('$expvar - explained variance for each component\n')
   cat('$cumexpvar - cumulative explained variance\n')
}

#' as.matrix method for ldecomp object
#' 
#' @description
#' Generic \code{as.matrix} function for linear decomposition. Returns a matrix with information 
#' about the decomposition.
#' 
#' @param x
#' object of class \code{ldecomp}
#' @param ...
#' other arguments
#' 
#' @export
as.matrix.ldecomp = function(x, ...) {
   data = cbind(x$expvar, x$cumexpvar)   
   rownames(data) = colnames(x$Q)
   colnames(data) = c('Expvar', 'Cumexpvar')   
   data
}  

#' Summary statistics for linear decomposition
#'
#' @description
#' Generic \code{summary} function for linear decomposition. Prints statistic about 
#' the decomposition.
#' 
#' @param object
#' object of class \code{ldecomp}
#' @param str
#' user specified text to show as a description of the object
#' @param ...
#' other arguments
#' 
#' @export
summary.ldecomp = function(object, str = NULL, ...) {
   if (is.null(str))
      str ='Summary for data decomposition (class ldecomp)'
   
   cat('\n')
   cat(str, '\n')
   cat(sprintf('\nSelected components: %d\n\n', object$ncomp.selected))      
   
   data = as.matrix(object)
   print(round(data, 2))   
}

#' Inverse error function 
#' 
#' @param x
#' a matrix or vector with data values
#' 
#' @export
erfinv = function (x) qnorm((1 + x)/2)/sqrt(2)


