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
   cols_excluded <- attr(loadings, "exclrows")

   # get sizes
   ncomp <- ncol(scores)
   nobj <- nrow(scores)
   nvar <- nrow(loadings)

   # remove excluded variables from loadings and residuals
   if (length(cols_excluded) > 0) {
      loadings <- loadings[-cols_excluded, , drop = FALSE]
      residuals <- residuals[, -cols_excluded, drop = FALSE]
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


#' Compute critical limits for orthogonal distance
#'
#' @param Q
#' matrix with distances
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#' @param lim.type
#' which method to use ("chisq", "ddrobust", "ddmoments")
#'
#' @description
#' Compute critical values for extremes and outliers based on the matrix with distances.
#'
#' If data driven approach is used ("ddrobust" or "ddmoments") then instead of the limit,
#' number of degrees of freedom (Nq) and scaling factor (q0) are computed instead.
#'
ldecomp.getLimits <- function(Q, T2, alpha, gamma, lim.type) {

   ncomp <- ncol(Q)

   # list of functions to compute critical values for Q
   fqcrit <- list(
      "chisq" = chisq.crit,
      "ddmoments" = ddmoments.crit,
      "ddrobust" = ddrobust.crit
   )

   # list of functions to compute critical values for T2
   fhcrit <- list(
      "chisq" = hotelling.crit,
      "ddmoments" = ddmoments.crit,
      "ddrobust" = ddrobust.crit
   )

   # remove excluded rows from Q
   rows_excluded <- attr(Q, "exclrows")
   if (length(rows_excluded) > 0) {
      Q <- Q[-rows_excluded, , drop = FALSE]
      T2 <- T2[-rows_excluded, , drop = FALSE]
   }

   # compute the values
   Qlim <- fqcrit[[lim.type]](Q, 1:ncomp, alpha, gamma)
   T2lim <- fhcrit[[lim.type]](T2, 1:ncomp, alpha, gamma)

   # in case of data driven approach we need to compute the critical values manually
   if (regexpr("dd", lim.type) > 0) {
      # get the parameters
      h0 <- T2lim[3, ]
      Nh <- T2lim[4, ]
      q0 <- Qlim[3, ]
      Nq <- Qlim[4, ]

      # calculate critical limits
      crit1 <- qchisq(1 - alpha, Nh + Nq)
      crit2 <- qchisq((1 - gamma)^(1/nrow(Q)), Nh + Nq)

      # save limits as slope and intercept of critical limit line
      Qlim[1, ] <- crit1 * q0 / Nq
      Qlim[2, ] <- crit2 * q0 / Nq
      T2lim[1, ] <- crit1 * h0 / Nh
      T2lim[2, ] <- crit2 * h0 / Nh
   }

   # set column and row names
   colnames(Qlim) <- colnames(T2lim) <- colnames(Q)
   rownames(Qlim) <- rownames(T2lim) <- c("Extremes limit", "Outliers limit", "Mean (u0)", "DoF")

   # set attributes
   attr(Qlim, "name") <- "Critical limits for Q"
   attr(T2lim, "name") <- "Critical limits for T2"
   attr(Qlim, "alpha") <- attr(T2lim, "alpha") <- alpha
   attr(Qlim, "gamma") <- attr(T2lim, "gamma") <- gamma

   return(list(Qlim = Qlim, T2lim = T2lim))
}

#' Calculate critical limits for distance values using Hotelling T2 distribution
#'
#' @param U
#' matrix or vector with distances
#' @param ncomp
#' number of components
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @return
#' vector with four values: critical limits for given alpha and gamma, mean distance and DoF.
#'
#' @export
hotelling.crit <- function(U, ncomp, alpha = 0.05, gamma = 0.01) {

   if (is.null(dim(U))) dim(U) <- c(length(U), 1)

   nobj <- nrow(U)
   u0 <- apply(U, 2, mean)
   DoF <- (nobj - ncomp)
   out <- rbind(
      (ncomp * (nobj - 1) / DoF) * qf(1 - alpha, ncomp, DoF),
      (ncomp * (nobj - 1) / DoF) * qf((1 - gamma)^(1/nobj), ncomp, DoF),
      u0, DoF
   )

   return(out)
}

#' Calculate probabilities for distance values and given parameters using Hotelling T2 distribution
#'
#' @param U
#' matrix or vector with distances
#' @param ncomp
#' number of components
#' @param lim
#' vector with parameters
#'
#' @export
hotelling.prob <- function(U, ncomp, lim){
      nobj <- lim[4] + ncomp
      return(pf(u * (DoF) / (ncomp * (nobj - 1)), ncomp, DoF))
}

#' Calculates critical limits for distance values using Chi-square distribution
#'
#' @description
#' The method is based on Chi-squared distribution with DF = 2 * (m(u)/s(u)^2
#'
#' @param U
#' matrix or vector with distance values
#' @param ncomp
#' number of components
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @export
chisq.crit <- function(U, ncomp, alpha = 0.05, gamma = 0.01) {

   if (is.null(dim(U))) dim(U) <- c(length(U), 1)

   nobj <- nrow(U)
   u0 <- apply(U, 2, mean)
   su <- apply(U, 2, sd)
   Nu <- 2 * (u0/su)^2

   DoF <- floor(Nu)
   DoF[DoF == 0] <- 1

   out <- rbind(
      qchisq(1 - alpha, DoF) * u0 / Nu,
      qchisq((1 - gamma)^(1/nobj), DoF) * u0 / Nu,
      u0, Nu
   )

   return(out)
}


#' Calculate probabilities for distance values using Chi-square distribution
#'
#' @param u
#' vector with distances
#' @param ncomp
#' number of components
#' @param lim
#' vector with parameters
#'
#' @export
chisq.prob <- function(u, ncomp, lim){
   u0 <- lim[3, ]
   Nu <- lim[4, ]
   DoF <- floor(Nu)
   DoF[DoF == 0] <- 1

   return(pchisq(Nu * u / u0, DoF))
}


#' Calculates critical limits for distance values using Data Driven moments approach
#'
#' @param U
#' matrix or vector with distance values
#' @param ncomp
#' number of components
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @export
ddmoments.crit <- function(U, ncomp, alpha = 0.05, gamma = 0.01) {

   if (is.null(dim(U))) dim(U) <- c(length(U), 1)

   u0 <- apply(U, 2, mean)
   su <- apply(U, 2, sd)

   Nu <- round(2 * (u0/su)^2)
   Nu <- ifelse(Nu > 250, 250, Nu)
   Nu <- ifelse(Nu < 1, 1, Nu)

   return(rbind(0, 0, u0, Nu))
}


#' Calculates critical limits for distance values using Data Driven robust approach
#'
#' @param U
#' matrix or vector with distance values
#' @param ncomp
#' number of components
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @export
ddrobust.crit <- function(U, ncomp, alpha, gamma) {

   if (is.null(dim(U))) dim(U) <- c(length(U), 1)

   Mu <- apply(U, 2, median)
   Su <- apply(U, 2, IQR)

   RM <- Su / Mu
   Nu <- round(exp((1.380948*log(2.68631 / RM)) ^ 1.185785))
   Nu[RM > 2.685592117] <- 1
   Nu[RM < 0.194565995] <- 100

   u0 <- 0.5 * Nu * (Mu/qchisq(0.50, Nu) + Su/(qchisq(0.75, Nu) - qchisq(0.25, Nu)))
   return(rbind(0, 0, u0, Nu))
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

