#' Class for storing and visualising of data linear decomposition (X = TP' + E)
#'
#' @description
#' Creates an object of ldecomp class.
#'
#' @param scores
#' matrix with score values (I x A).
#' @param loadings
#' matrix with loading values (J x A).
#' @param residuals
#' matrix with data residuals (I x J)
#' @param tnorm
#' vector with singular values for score normalization
#' @param ncomp.selected
#' number of selected components
#' @param cal
#' logical, true if data is for calibration of a LDECOMP based model
#'
#' @return
#' Returns an object (list) of \code{ldecomp} class with following fields:
#' \item{scores }{matrix with score values (nobj x ncomp).}
#' \item{residuals }{matrix with data residuals (nobj x nvar).}
#' \item{tnorm }{vector with singular values used for scores normalization.}
#' \item{ncomp.selected }{selected number of components.}
#' \item{expvar }{explained variance for each component.}
#' \item{cumexpvar }{cumulative explained variance.}
#' \item{T2 }{matrix with T2 distances (nobj x ncomp).}
#' \item{Q }{matrix with Q statistic (nobj x ncomp).}
#'
#' @details
#' \code{ldecomp} is a general class for decomposition of data in form X = TP' + E. Here, X is a
#' data matrix, T - matrix with scores, P - matrix with loadings and E - matrix with residuals.
#' It is used, for example, for PCA results (\code{\link{pcares}}), in PLS and other methods.
#' The class also includes methods for calculation and plotting residuals, variances, and so on.
#'
#' There is no need to use the \code{ldecomp} manually. For example, when build PCA model
#' with \code{\link{pca}} or apply it to a new data, the results will automatically inherit
#' all methods of \code{ldecomp}.
#'
#' @importFrom methods show
#' @importFrom stats convolve cor lm na.exclude predict pt qf qnorm qt sd var
#'
#' @export
ldecomp <- function(scores, loadings, residuals, ncomp.selected = ncol(scores),
   Qlim = NULL, T2lim = NULL, tnorm = NULL) {

   ncomp <- ncol(scores)

   obj <- list(
      scores = scores,
      residuals = residuals,
      ncomp = ncomp,
      ncomp.selected = ncomp.selected
   )

   # get distances and add them to the object
   dist <- ldecomp.getDistances(scores, loadings, residuals, tnorm)
   obj <- c(obj, dist)

   # get variance and add it to the object
   var <- ldecomp.getVariances(scores, loadings, residuals, dist$Q)
   obj <- c(obj, var)

   # if limits are not provided - compute them

   obj$call <- match.call()
   class(obj) <- "ldecomp"

   return(obj)
}

#' Compute score and residual distances
#'
#' @description
#' Compute orthogonal Euclidean distance from object to PC space (Q, q) and Mahalanobis
#' squared distance between projection of the object to the space and its origin (T2, h).
#'
#' @param scores
#' matrix with scores (T).
#' @param loadings
#' matrix with loadings (P).
#' @param residuals
#' matrix with residuals (E).
#' @param tnorm
#' vector with singular values for scores normalisation (can be provided)
#'
#' @details
#' The distances are calculated for every 1:n components, where n goes from 1 to ncomp
#' (number of columns in scores and loadings).
#'
#' @return
#' Returns a list with Q, T2 and tnorm values for each component.
#'
ldecomp.getDistances = function(scores, loadings, residuals, tnorm = NULL) {

   # get names and attributes
   var_names <- rownames(loadings)
   obj_names <- rownames(scores)
   rows_excluded <- attr(scores, "exclrows")
   cols_excluded <- attr(scores, "exclcols")

   # get sizes
   ncomp <- ncol(scores)
   nobj <- nrow(scores)
   nvar <- nrow(loadings)

   # remove excluded variables from loadings and residuals
   if (length(cols_excluded) > 0) {
      loadings <- loadings[-cols_excluded, , drop = FALSE]
      residuals <- residuals[-cols_excluded, , drop = FALSE]
   }

   # get rid of hidden scores and residuals (needed for some calculations)
   scores_visible <- scores
   residuals_visible <- residuals
   nobj_visible <- nobj
   if (length(rows_excluded) > 0) {
      scores_visible <- scores_visible[-rows_excluded, , drop = FALSE]
      residuals_visible <- residuals_visible[-rows_excluded, , drop = FALSE]
      nobj_visible <- nobj - length(rows_excluded)
   }


   # calculate singular values for score normalization (if not provided)
   if (is.null(tnorm)) {
      tnorm <- sqrt(colSums(scores_visible^2)/(nobj_visible - 1));
   }

   # normalize the scores
   scoresn <- scale(scores, center = FALSE, scale = tnorm)

   # prepare zero matrices for the and model power
   T2 <- matrix(0, nrow = nobj, ncol = ncomp)
   Q <- matrix(0, nrow = nobj, ncol = ncomp)

   # calculate distances and model power for each possible number of components in model
   for (i in 1:ncomp) {
      res <- residuals
      if (i < ncomp) {
         res <- res +
            tcrossprod(
               scores[, (i + 1):ncomp, drop = F],
               loadings[, (i + 1):ncomp, drop = F]
            )
      }

      Q[, i] <- rowSums(res^2)
      T2[, i] <- rowSums(scoresn[, 1:i, drop = F]^2)
   }

   # set attributes for Q
   Q <- mda.setattr(Q, mda.getattr(scores), type = 'row')
   attr(Q, "name") <- "Squared residual distance (q)"
   attr(Q, "xaxis.name") <- "Components"

   # set attributes for T2
   T2 = mda.setattr(T2, mda.getattr(Q))
   attr(T2, 'name') = 'Score distance (h)'

   colnames(Q) <- colnames(T2) <- colnames(loadings)
   rownames(Q) <- rownames(T2) <- rownames(scores)

   # return the results
   return(list(Q = Q, T2 = T2, tnorm = tnorm))
}


#' Compute explained variance
#'
#' @description
#' Computes explained variance and cumulative explained variance for data decomposition.
#'
#' @param scores
#' matrix with scores (T).
#' @param loadings
#' matrix with loadings (P).
#' @param residuals
#' matrix with residuals (E).
#' @param Q
#' matrix with squared orthogonal distances.
#'
#' @return
#' Returns a list with two vectors.
#'
ldecomp.getVariances <- function(scores, loadings, residuals, Q) {

   # get names and attributes
   rows_excluded <- attr(scores, "exclrows")
   cols_excluded <- attr(scores, "exclcols")


   # remove excluded columns from loadings and residuals
   if (length(cols_excluded) > 0) {
      loadings <- loadings[-cols_excluded, , drop = FALSE]
      residuals <- residuals[-cols_excluded, , drop = FALSE]
   }

   # remove excluded rows from scores, residuals and Q
   if (length(rows_excluded) > 0) {
      scores <- scores[-rows_excluded, , drop = FALSE]
      residuals <- residuals[-rows_excluded, , drop = FALSE]
      Q <- Q[-rows_excluded, , drop = FALSE]
   }

   # compute total variance
   totvar <- sum(tcrossprod(scores, loadings)^2) + sum(residuals^2)

   # compute explained variance
   cumexpvar <- 100 * (1 - colSums(Q) / totvar)
   expvar <- c(cumexpvar[1], diff(cumexpvar))

   names(cumexpvar) <- names(expvar) <- colnames(Q)
   return(list(expvar = expvar, cumexpvar = cumexpvar))
}

ldecomp.getQLimits <- function(Q, Qlim, alpha, gamma) {

}

ldecomp.getT2Limits <- function() {

}

#' Calculate critical limits for score distance using Hotelling T2 distribution
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
reslim.hotelling <- function(ncomp, T2 = NULL, alpha = 0.05, gamma = 0.01,
                            T2lim = NULL, return = 'limits') {

   # return critical values
   if (return == "limits") {
      nobj <- T2lim[4] + ncomp
      out <- pf(T2 * (nobj - ncomp) / (ncomp * (nobj - 1)), ncomp, nobj - ncomp)
      return(out)
   }

   # return probabilities
   nobj <- length(T2)
   DoF <- nobj - ncomp
   out <- c(0, 0, mean(T2), DoF)
   if (nobj > ncomp) {
      out[1:2] <- (ncomp * (nobj - 1) / DoF) * qf(c(1 - alpha, 1 - gamma), ncomp, DoF)
   }

   return(out)
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


