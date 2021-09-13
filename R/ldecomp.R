#' Class for storing and visualising linear decomposition of dataset (X = TP' + E)
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
#' @param eigenvals
#' vector with eigenvalues for the loadings
#' @param ncomp.selected
#' number of selected components
#'
#' @return
#' Returns an object (list) of \code{ldecomp} class with following fields:
#' \item{scores }{matrix with score values (I x A).}
#' \item{residuals }{matrix with data residuals (I x J).}
#' \item{T2 }{matrix with score distances (I x A).}
#' \item{Q }{matrix with orthogonal distances (I x A).}
#' \item{ncomp.selected }{selected number of components.}
#' \item{expvar }{explained variance for each component.}
#' \item{cumexpvar }{cumulative explained variance.}
#'
#' @details
#' \code{ldecomp} is a general class for storing results of decomposition of dataset in
#' form X = TP' + E. Here, X is a data matrix, T - matrix with scores, P - matrix with
#' loadings and E - matrix with residuals. It is used, for example, for PCA results
#' (\code{\link{pcares}}), in PLS and other methods. The class also includes methods for
#' calculation of residual distances and explained variance.
#'
#' There is no need to use the \code{ldecomp} manually. For example, when build PCA model
#' with \code{\link{pca}} or apply it to a new data, the results will automatically inherit
#' all methods of \code{ldecomp}.
#'
#' @importFrom methods show
#' @importFrom stats convolve cor lm na.exclude predict pt qf qnorm qt sd var
#'
#' @export
ldecomp <- function(scores, loadings, residuals, eigenvals, ncomp.selected = ncol(scores)) {

   ncomp <- ncol(scores)

   obj <- list(
      scores = scores,
      residuals = residuals,
      ncomp = ncomp,
      ncomp.selected = ncomp.selected,
      categories = NULL
   )

   # get distances and add them to the object
   dist <- ldecomp.getDistances(scores, loadings, residuals, eigenvals)
   obj <- c(obj, dist)

   # get variance and add it to the object
   var <- ldecomp.getVariances(scores, loadings, residuals, dist$Q)
   obj <- c(obj, var)

   obj$call <- match.call()
   class(obj) <- "ldecomp"

   return(obj)
}

#' Cumulative explained variance plot
#'
#' @description
#' Shows a plot with cumulative explained variance vs. number of components.
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param type
#' type of the plot
#' @param labels
#' what to show as labels for plot objects
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotCumVariance.ldecomp <- function(obj, type = "b", labels = "values", show.plot = TRUE, ...) {

   return(
      plotVariance(obj, variance = "cumexpvar", type = type, labels = labels,
         show.plot = show.plot, ...)
   )
}

#' Explained variance plot
#'
#' @description
#' Shows a plot with explained variance vs. number of components.
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param type
#' type of the plot
#' @param variance
#' string, which variance to make the plot for ("expvar", "cumexpvar")
#' @param labels
#' what to show as labels for plot objects.
#' @param xticks
#' vector with ticks for x-axis
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ylab
#' label for y-axis
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotVariance.ldecomp <- function(obj, type = "b", variance = "expvar", labels = "values",
   xticks = seq_len(obj$ncomp), show.plot = TRUE, ylab = "Explained variance, %", ...) {

   if (!show.plot) return(obj[[variance]])

   return(
      mdaplot(obj[[variance]], xticks = xticks, labels = labels, type = type, ylab = ylab, ...)
   )
}

#' Scores plot
#'
#' @description
#' Shows a plot with scores values for data objects.
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param comp
#' which components to show the plot for (can be one value or vector with two values).
#' @param type
#' type of the plot
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotScores.ldecomp <- function(obj, comp = c(1, 2), type = "p", show.axes = TRUE,
   show.plot = TRUE, ...) {

   # get scores for given components and generate column names with explained variance
   plot_data <- mda.subset(obj$scores, select = comp)
   colnames(plot_data) <- paste0("Comp ", comp, " (", round(obj$expvar[comp], 2), "%)")
   attr(plot_data, "name") <- "Scores"

   # if no plot required - return plot series object
   if (!show.plot) {
      return(plot_data)
   }

   # set up values for showing axes lines
   show.lines <- FALSE
   if (show.axes) {
      show.lines <- if (length(comp) == 2 && type == "p") c(0, 0) else c(NA, 0)
   }

   # scatter plot
   if (type == "p") {
      p <- mdaplot(plot_data, type = type, show.lines = show.lines, ...)
      return(invisible(p))
   }

   # line or bar plot
   plot_data <- mda.t(plot_data)
   attr(plot_data, "yaxis.name") <- "Score"
   return(mdaplotg(plot_data, type = type, show.lines = show.lines, ...))
}

#' Residual distance plot
#'
#' @description
#' Shows a plot with orthogonal (Q, q) vs. score (T2, h) distances for data objects.
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param ncomp
#' number of components to show the plot for (if NULL, selected by model value will be used).
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param labels
#' what to show as labels if necessary
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotResiduals.ldecomp <- function(obj, ncomp = obj$ncomp.selected, norm = FALSE, log = FALSE,
   show.labels = FALSE, labels = "names", show.plot = TRUE, ...) {

   attrs <- mda.getattr(obj$Q)

  # function for transforming distances
   transform <- function(u, u0, norm, log) {
      if (norm) u <- u / u0
      if (log) u <- log(1 + u)
      return(u)
   }

   # function for creating labels depending on transformation
   get_label <- function(lab, norm, log) {
      if (norm) lab <- paste0(lab, "/", lab, "0")
      if (log) lab <- paste0("log(1 + ", lab, ")")
      return(lab)
   }

   # get scale factors
   h0 <- if (!is.null(attr(obj$T2, "u0"))) attr(obj$T2, "u0")[[ncomp]]
   q0 <- if (!is.null(attr(obj$Q, "u0"))) attr(obj$Q, "u0")[[ncomp]]

   # check that scaling values exist
   if (norm && (is.null(h0) || is.null(q0))) {
      warning("Can not normalize distances as scaling values are absent.")
      norm <- FALSE
   }

   # prepare plot data
   h <- transform(obj$T2[, ncomp], h0, norm, log)
   q <- transform(obj$Q[, ncomp], q0, norm, log)

   # default values for local labels
   lxlab <- get_label("h", norm, log)
   lylab <- get_label("q", norm, log)

   # combine everything to dataset and assign attributes
   plot_data <- mda.cbind(h, q)
   plot_data <- mda.setattr(plot_data, attrs, "row")
   rownames(plot_data) <- rownames(obj$Q)
   colnames(plot_data) <- c(
      paste0("Score distance, ", lxlab),
      paste0("Orthogonal distance, ", lylab)
   )

   attr(plot_data, "name") <- sprintf("Distances (ncomp = %d)", ncomp)

   # if no plot required - return plot series object
   if (!show.plot) return(plot_data)

   # show plot
   return(mdaplot(plot_data, ...))
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
print.ldecomp <- function(x, str = NULL, ...) {
   if (is.null(str)) {
      str <- "Results of data decomposition (class ldecomp)."
   }

   if (nchar(str) > 0) {
      fprintf("\n%s\n", str)
   }

   cat("\nMajor fields:\n")
   cat("$scores - matrix with score values\n")
   cat("$T2 - matrix with T2 distances\n")
   cat("$Q - matrix with Q residuals\n")
   cat("$ncomp.selected - selected number of components\n")
   cat("$expvar - explained variance for each component\n")
   cat("$cumexpvar - cumulative explained variance\n")
}

#' as.matrix method for ldecomp object
#'
#' @description
#' Generic \code{as.matrix} function for linear decomposition. Returns a matrix with information
#' about the decomposition.
#'
#' @param x
#' object of class \code{ldecomp}
#' @param ncomp
#' number of components to get the result for (if NULL will return for each available)
#' @param ...
#' other arguments
#'
#' @export
as.matrix.ldecomp <- function(x, ncomp = NULL, ...) {

   out <- cbind(x$expvar, x$cumexpvar)
   rownames(out) <- colnames(x$Q)
   colnames(out) <- c("Expvar", "Cumexpvar")

   if (!is.null(ncomp)) {
      out <- out[ncomp, , drop = FALSE]
   }

   return(out)
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
summary.ldecomp <- function(object, str = NULL, ...) {
   if (is.null(str)) {
      str <- "Summary for data decomposition (class ldecomp)."
   }

   fprintf("\n%s\n", str)
   fprintf("\nSelected components: %d\n\n", object$ncomp.selected)

   print(round(as.matrix(object), 2))
}


##########################
# Static methods         #
##########################


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
   attr(expvar, "name") <- "Variance"
   attr(cumexpvar, "name") <- "Cumulative variance"
   attr(expvar, "xaxis.name") <- attr(cumexpvar, "xaxis.name") <- "Components"

   return(list(expvar = expvar, cumexpvar = cumexpvar))
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
#' @param eigenvals
#' vector with eigenvalues for the components
#'
#' @details
#' The distances are calculated for every 1:n components, where n goes from 1 to ncomp
#' (number of columns in scores and loadings).
#'
#' @return
#' Returns a list with Q, T2 and tnorm values for each component.
#'
ldecomp.getDistances <- function(scores, loadings, residuals, eigenvals) {

   # get names and attributes
   rows_excluded <- attr(scores, "exclrows")
   cols_excluded <- attr(loadings, "exclrows")

   # get sizes
   ncomp <- ncol(scores)
   nobj <- nrow(scores)

   # remove excluded variables from loadings and residuals
   if (length(cols_excluded) > 0) {
      loadings <- loadings[-cols_excluded, , drop = FALSE]
      residuals <- residuals[, -cols_excluded, drop = FALSE]
   }

   # get rid of hidden scores and residuals (needed for some calculations)
   scores_visible <- scores
   residuals_visible <- residuals
   if (length(rows_excluded) > 0) {
      scores_visible <- scores_visible[-rows_excluded, , drop = FALSE]
      residuals_visible <- residuals_visible[-rows_excluded, , drop = FALSE]
   }

   # normalize the scores
   scoresn <- scale(scores, center = FALSE, scale = sqrt(eigenvals))

   # prepare zero matrices for the and model power
   T2 <- matrix(0, nrow = nobj, ncol = ncomp)
   Q <- matrix(0, nrow = nobj, ncol = ncomp)

   # calculate distances and model power for each possible number of components in model
   for (i in seq_len(ncomp)) {
      res <- residuals
      if (i < ncomp) {
         res <- res +
            tcrossprod(
               scores[, (i + 1):ncomp, drop = F],
               loadings[, (i + 1):ncomp, drop = F]
            )
      }

      Q[, i] <- rowSums(res^2)
      T2[, i] <- rowSums(scoresn[, seq_len(i), drop = F]^2)
   }

   # set attributes for Q
   Q <- mda.setattr(Q, mda.getattr(scores), type = "row")
   attr(Q, "name") <- "Squared residual distance (q)"
   attr(Q, "xaxis.name") <- "Components"

   # set attributes for T2
   T2 <- mda.setattr(T2, mda.getattr(Q))
   attr(T2, "name") <- "Score distance (h)"

   colnames(Q) <- colnames(T2) <- colnames(loadings)
   rownames(Q) <- rownames(T2) <- rownames(scores)

   # return the results
   return(list(Q = Q, T2 = T2))
}

###############################
# Methods for critical limits #
###############################

#' Calculate critical limits for distance values using Jackson-Mudholkar approach
#'
#' @param residuals
#' matrix with PCA residuals
#' @param eigenvals
#' vector with eigenvalues for PCA components
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @return
#' vector with four values: critical limits for given alpha and gamma, mean distance and DoF.
#'
#' @export
jm.crit <- function(residuals, eigenvals, alpha = 0.05, gamma = 0.01) {

   # if not all eigenvalues available - ise residuals to compute the rest
   ncomp <- length(eigenvals)
   nobj <- nrow(residuals)
   max_ncomp <- min(nrow(residuals) - 1, ncol(residuals))
   if (length(eigenvals) < max_ncomp) {
      eigenvals <- c(eigenvals, svd(residuals)$d[seq_len(max_ncomp - ncomp)]^2 / (nobj - 1))
   }

   # since it is residuals we do not need eigenvalue for PC1
   eigenvals <- eigenvals[-1]
   t1 <- rev(cumsum(rev(eigenvals)))[seq_len(ncomp)]
   t2 <- rev(cumsum(rev(eigenvals)^2))[seq_len(ncomp)]
   t3 <- rev(cumsum(rev(eigenvals)^3))[seq_len(ncomp)]

   h0 <- 1 - 2 * t1 * t3 / 3 / (t2^2);
   ifelse(h0 < 0.001, h0 <- 0.001, h0)

   # inverse error function
   erfinv <- function(x) qnorm((1 + x) / 2) / sqrt(2)
   gcl <- 1 - (1 - gamma) ^ (1 / nobj)
   ca <- sqrt(2) * erfinv(c(1 - 2 * alpha, (1 - 2 * gcl)))

   # compute h1 for alpha and gamma
   h1a <- ca[1] * sqrt(2 * t2 * h0^2) / t1
   h1g <- ca[2] * sqrt(2 * t2 * h0^2) / t1
   h2 <- t2 * h0 * (h0 - 1) / (t1 ^ 2)

   out <- rbind(
      t1 * (1 + h1a + h2) ^ (1 / h0),
      t1 * (1 + h1g + h2) ^ (1 / h0)
   )

   if (ncomp == max_ncomp) {
      out[, max_ncomp] <- 0
   }

   attr(out, "eigenvals") <- eigenvals
   return(out)
}

#' Calculate probabilities for distance values and given parameters using Hotelling T2 distribution
#'
#' @param u
#' vector with distances
#' @param eigenvals
#' vector with eigenvalues for PCA components
#' @param ncomp
#' number of components
#'
#' @export
jm.prob <- function(u, eigenvals, ncomp) {

   erf <- function(x) 1 - pnorm(-x * sqrt(2)) * 2

   t1 <- rev(cumsum(rev(eigenvals)))[ncomp]
   t2 <- rev(cumsum(rev(eigenvals)^2))[ncomp]
   t3 <- rev(cumsum(rev(eigenvals)^3))[ncomp]

   h0 <- 1 - 2 * t1 * t3 / 3 / (t2^2);
   ifelse(h0 < 0.001, h0 <- 0.001, h0)

   h1 <- (u / t1)^h0
   h2 <- t2 * h0 * (h0 - 1) / t1^2
   d <- t1 * (h1 - 1 - h2) / (sqrt(2 * t2) * h0)

   return(0.5 * (1 + erf(d / sqrt(2))))
}

#' Calculate critical limits for distance values using Hotelling T2 distribution
#'
#' @param nobj
#' number of objects in calibration set
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
hotelling.crit <- function(nobj, ncomp, alpha = 0.05, gamma = 0.01) {
   return(
      rbind(
         (ncomp * (nobj - 1) / (nobj - ncomp)) * qf(1 - alpha, ncomp, (nobj - ncomp)),
         (ncomp * (nobj - 1) / (nobj - ncomp)) * qf((1 - gamma) ^ (1 / nobj), ncomp, (nobj - ncomp))
      )
   )
}

#' Calculate probabilities for distance values and given parameters using Hotelling T2 distribution
#'
#' @param u
#' vector with distances
#' @param ncomp
#' number of components
#' @param nobj
#' number of objects in calibration set
#'
#' @export
hotelling.prob <- function(u, ncomp, nobj) {
   return(pf(u * (nobj - ncomp) / (ncomp * (nobj - 1)), ncomp, (nobj - ncomp)))
}

#' Calculates critical limits for distance values using Chi-square distribution
#'
#' @description
#' The method is based on Chi-squared distribution with DF = 2 * (m(u)/s(u)^2
#'
#' @param param
#' matrix with distribution parameters
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @export
chisq.crit <- function(param, alpha = 0.05, gamma = 0.01) {

   u0 <- param$u0
   nobj <- param$nobj
   Nu <- param$Nu

   DoF <- floor(Nu)
   DoF <- ifelse(DoF < 1, 1, DoF)

   return(
      rbind(
         qchisq(1 - alpha, DoF) * u0 / Nu,
         qchisq((1 - gamma) ^ (1 / nobj), DoF) * u0 / Nu
      )
   )
}

#' Calculate probabilities for distance values using Chi-square distribution
#'
#' @param u
#' vector with distances
#' @param param
#' vector with distribution parameters
#'
#' @export
chisq.prob <- function(u, param) {
   u0 <- param[1]
   Nu <- param[2]

   DoF <- floor(Nu)
   DoF[DoF == 0] <- 1
   return(pchisq(Nu * u / u0, DoF))
}

#' Calculates critical limits for distance values using Data Driven moments approach
#'
#' @param paramQ
#' matrix with parameters for distribution of Q distances
#' @param paramT2
#' matrix with parameters for distribution of T2 distances
#' @param alpha
#' significance level for extreme objects
#' @param gamma
#' significance level for outliers
#'
#' @export
dd.crit <- function(paramQ, paramT2, alpha = 0.05, gamma = 0.01) {

   nobj <- paramQ$nobj
   Nq <- round(paramQ$Nu)
   Nq[Nq < 1] <- 1
   Nq[Nq > 250] <- 250

   Nh <- round(paramT2$Nu)
   Nh[Nh < 1] <- 1
   Nh[Nh > 250] <- 250

   return(
      rbind(
         qchisq(1 - alpha, Nq + Nh),
         qchisq((1 - gamma) ^ (1 / nobj), Nq + Nh)
      )
   )
}

#' Calculates critical limits for distance values using Data Driven moments approach
#'
#' @param U
#' matrix or vector with distance values
#'
#' @export
ddmoments.param <- function(U) {

   if (is.null(dim(U))) dim(U) <- c(length(U), 1)

   u0 <- apply(U, 2, mean)
   su <- apply(U, 2, sd)
   Nu <- 2 * (u0 / su)^2

   return(list(u0 = u0, Nu = Nu, nobj = nrow(U)))
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
ddrobust.param <- function(U, ncomp, alpha, gamma) {

   if (is.null(dim(U))) dim(U) <- c(length(U), 1)

   Mu <- apply(U, 2, median)
   Su <- apply(U, 2, IQR)

   RM <- Su / Mu
   Nu <- round(exp((1.380948 * log(2.68631 / RM)) ^ 1.185785))
   Nu[RM > 2.685592117] <- 1
   Nu[RM < 0.194565995] <- 100

   u0 <- 0.5 * Nu * (Mu / qchisq(0.50, Nu) + Su / (qchisq(0.75, Nu) - qchisq(0.25, Nu)))
   return(list(u0 = u0, Nu = Nu, nobj = nrow(U)))
}

#' Compute parameters for critical limits based on calibration results
#'
#' @param U
#' matrix with residual distances
#'
#' @export
ldecomp.getLimParams <- function(U) {

   U <- mda.purgeRows(U)

   return(
      list(
         "moments" = ddmoments.param(U),
         "robust" = ddrobust.param(U),
         "nobj" = nrow(U)
      )
   )
}

#' Compute critical limits for orthogonal distances (Q)
#'
#' @param lim.type
#' which method to use for calculation of critical limits for residuals
#' @param alpha
#' significance level for extreme limits.
#' @param gamma
#' significance level for outlier limits.
#' @param params
#' distribution parameters returned by ldecomp.getLimParams
#' @param residuals
#' matrix with residuals (E)
#' @param eigenvals
#' egenvalues for the components used to decompose the data
#'
#' @export
ldecomp.getQLimits <- function(lim.type, alpha, gamma, params, residuals, eigenvals) {

   pQ <- if (regexpr("robust", lim.type) > 0) params$Q$robust else params$Q$moments
   ncomp <- length(pQ$u0)

   if (lim.type == "jm") {
      # methods based on Jackson-Mudholkar approach
      residuals <- mda.purge(residuals)
      lim <- jm.crit(residuals, eigenvals, alpha, gamma)
      eigenvals <- attr(lim, "eigenvals")
      lim <- rbind(lim, pQ$u0, nrow(residuals))
      attr(lim, "eigenvals") <- eigenvals

   } else {
      # methods based on chi-square distribution
      pT2 <- if (regexpr("robust", lim.type) > 0) params$T2$robust else params$T2$moments
      DoF <- round(pQ$Nu)
      DoF[DoF < 1] <- 1
      DoF[DoF > 250] <- 250

      lim <- switch(lim.type,
         "chisq" = chisq.crit(pQ, alpha, gamma),
         "ddmoments" = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pQ$u0),
         "ddrobust"  = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pQ$u0),
         stop("Wrong value for 'lim.type' parameter.")
      )

      lim <- rbind(lim, pQ$u0, DoF)
   }

   colnames(lim) <- paste("Comp", seq_len(ncomp))
   rownames(lim) <- c("Extremes limits", "Outliers limits", "Mean", "DoF")
   attr(lim, "name") <- "Critical limits for orthogonal distances (Q)"
   attr(lim, "alpha") <- alpha
   attr(lim, "gamma") <- gamma
   attr(lim, "lim.type") <- lim.type

   return(lim)
}

#' Compute critical limits for score distances (T2)
#'
#' @param lim.type
#' which method to use for calculation ("chisq", "ddmoments", "ddrobust")
#' @param alpha
#' significance level for extreme limits.
#' @param gamma
#' significance level for outlier limits.
#' @param params
#' distribution parameters returned by ldecomp.getLimParams
#'
#' @export
ldecomp.getT2Limits <- function(lim.type, alpha, gamma, params) {

   pQ <- if (regexpr("robust", lim.type) > 0) params$Q$robust else params$Q$moments
   pT2 <- if (regexpr("robust", lim.type) > 0) params$T2$robust else params$T2$moments
   ncomp <- length(pT2$u0)

   DoF <- round(pT2$Nu)
   DoF[DoF < 1] <- 1
   DoF[DoF > 250] <- 250

   if (lim.type %in% c("jm", "chisq")) DoF <- pT2$nobj

   lim <- switch(lim.type,
      "jm" = hotelling.crit(pT2$nobj, seq_len(ncomp), alpha, gamma),
      "chisq" = hotelling.crit(pT2$nobj, seq_len(ncomp), alpha, gamma),
      "ddmoments" = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pT2$u0),
      "ddrobust"  = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pT2$u0),
      stop("Wrong value for 'lim.type' parameter.")
   )

   lim <- rbind(lim, pT2$u0, DoF)
   colnames(lim) <- paste("Comp", seq_len(ncomp))
   rownames(lim) <- c("Extremes limits", "Outliers limits", "Mean", "DoF")
   attr(lim, "name") <- "Critical limits for score distances (T2)"
   attr(lim, "alpha") <- alpha
   attr(lim, "gamma") <- gamma
   attr(lim, "lim.type") <- lim.type

   return(lim)
}

#' Compute coordinates of lines or curves with critical limits
#'
#' @param Qlim
#' matrix with critical limits for orthogonal distances
#' @param T2lim
#' matrix with critical limits for score distances
#' @param ncomp
#' number of components for computing the coordinates
#' @param norm
#' logical, shall distance values be normalized or not
#' @param log
#' logical, shall log transformation be applied or not
#' @param show.limits
#' vector with two logical values defining if limits for extreme and/or outliers must be shown
#'
#' @return
#' list with two matrices (x and y coordinates of corresponding limits)
#'
#' @export
ldecomp.getLimitsCoordinates <- function(Qlim, T2lim, ncomp, norm, log,
   show.limits = c(TRUE, TRUE)) {

   # get parameters
   h0 <- T2lim[3, ncomp]
   q0 <- Qlim[3, ncomp]
   Nh <- T2lim[4, ncomp]
   Nq <- Qlim[4, ncomp]
   lim.type <- attr(Qlim, "lim.type")

   # check that show.limits is logical
   if (!all(is.logical(show.limits))) {
      stop("Parameter 'show.limits' must have logical value(s).")
   }

   # if show.limits has only one value - duplicate it
   if (length(show.limits) == 1) {
      show.limits <- rep(show.limits, 2)
   }

   # compute the limits
   if (lim.type %in% c("jm", "chisq")) {

      # quadratic limits
      hE <- c(0, T2lim[1, ncomp], T2lim[1, ncomp])
      hO <- c(0, T2lim[2, ncomp], T2lim[2, ncomp])

      qE <- c(Qlim[1, ncomp], Qlim[1, ncomp], 0)
      qO <- c(Qlim[2, ncomp], Qlim[2, ncomp], 0)

   } else {

      ## slope and intercepts
      eB <- Qlim[1, ncomp]
      oB <- Qlim[2, ncomp]
      eA <- oA <- -1 * (q0 / h0) * (Nh / Nq)

      hE <- seq(-0.95, -eB / eA, length.out = 100)
      hO <- seq(-0.95, -oB / oA, length.out = 100)
      qE <- eA * hE + eB
      qO <- oA * hO + oB
   }

   if (norm) {
      hE <- hE / h0
      qE <- qE / q0
      hO <- hO / h0
      qO <- qO / q0
   }

   if (log) {
      hE <- log(1 + hE)
      qE <- log(1 + qE)
      hO <- log(1 + hO)
      qO <- log(1 + qO)
   }

   return(list(
      extremes = if (show.limits[[1]]) cbind(hE, qE),
      outliers = if (show.limits[[2]]) cbind(hO, qO)
   ))
}

#' Residuals distance plot for a set of ldecomp objects
#'
#' @description
#' Shows a plot with score (T2, h) vs orthogonal (Q, q) distances and corresponding critical
#' limits for given number of components.
#'
#' @param res
#' list with result objects to show the plot for
#' @param Qlim
#' matrix with critical limits for orthogonal distance
#' @param T2lim
#' matrix with critical limits for score distance
#' @param ncomp
#' how many components to use (by default optimal value selected for the model will be used)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param cgroup
#' color grouping of plot points (works only if one result object is available)
#' @param xlim
#' limits for x-axis (if NULL will be computed automatically)
#' @param ylim
#' limits for y-axis (if NULL will be computed automatically)
#' @param show.legend
#' logical, show or not a legend on the plot (needed if several result objects are available)
#' @param show.limits
#' vector with two logical values defining if limits for extreme and/or outliers must be shown
#' @param lim.col
#' vector with two values - line color for extreme and outlier limits
#' @param lim.lwd
#' vector with two values - line width for extreme and outlier limits
#' @param lim.lty
#' vector with two values - line type for extreme and outlier limits
#' @param show.legend
#' logical, show or not legend on the plot (if more than one result object)
#' @param legend.position
#' if legend must be shown, where it should be
#' @param show.excluded
#' logical, show or hide rows marked as excluded (attribute `exclrows`).
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' The function is a bit more advanced version of \code{\link{plotResiduals.ldecomp}}. It allows to
#' show distance values for several result objects (e.g. calibration and test set or calibration
#' and new prediction set) as well as display the correspondng critical limits in form of lines
#' or curves.
#'
#' Depending on how many result objects your model has or how many you specified manually,
#' using the \code{res} parameter, the plot behaves in a bit different way.
#'
#' If only one result object is provided, then it allows to colorise the points using \code{cgroup}
#' parameter. If two or more result objects are provided, then the function show
#' distances in groups, and adds corresponding legend.
#'
#' The function can show distance values normalised (h/h0 and q/q0) as well as with log
#' transformation (log(1 + h/h0), log(1 + q/q0)). The latter is useful if distribution of the
#' points is skewed and most of them are densely located around bottom left corner.
#'
#' @export
ldecomp.plotResiduals <- function(res, Qlim, T2lim, ncomp, log = FALSE, norm = FALSE,
   cgroup = NULL, xlim = NULL, ylim = NULL, show.limits = c(TRUE, TRUE),
   lim.col = c("darkgray", "darkgray"), lim.lwd = c(1, 1), lim.lty = c(2, 3),
   show.legend = TRUE, legend.position = "topright", show.excluded = FALSE, ...) {

   # return column with values either with or without excluded outliers
   getValues <- function(x, dim) {
      return(if (show.excluded) x[, dim] else mda.purgeRows(x)[, dim])
   }

   # compute limits fo axis depending on values and position of critical limits
   getPlotLim <- function(lim, pd, ld, dim) {
      if (!is.null(lim) || all(!show.limits)) return(lim)
      limits <- if (show.limits[[2]]) max(ld$outliers[, dim]) else max(ld$extremes[, dim])
      return(c(0, max(sapply(pd, function(x) max(c(getValues(x, dim), limits)) * 1.05))))
   }

   # check that show.limits is logical
   if (!all(is.logical(show.limits))) {
      stop("Parameter 'show.limits' must have logical value(s).")
   }

   # if show.limits has only one value - duplicate it
   if (length(show.limits) == 1) {
      show.limits <- rep(show.limits, 2)
   }

   # keep only ojects of class "ldecomp" in result list
   res <- getRes(res, "ldecomp")

   # compute plot data for each result object
   plot_data <- lapply(res, plotResiduals.ldecomp, ncomp = ncomp, norm = norm, log = log,
      show.plot = FALSE)

   # get coordinates for critical limits
   lim_data <- ldecomp.getLimitsCoordinates(Qlim, T2lim, ncomp = ncomp, norm = norm, log = log)
   xlim <- getPlotLim(xlim, plot_data, lim_data, 1)
   ylim <- getPlotLim(ylim, plot_data, lim_data, 2)

   # make plot
   if (length(plot_data) == 1) {
      mdaplot(plot_data[[1]], type = "p", xlim = xlim, ylim = ylim, cgroup = cgroup,
         show.excluded = show.excluded, ...)
   } else {
      mdaplotg(plot_data, type = "p", xlim = xlim, ylim = ylim, show.legend = show.legend,
         show.excluded = show.excluded, legend.position = legend.position, ...)
   }

   # show critical limits
   if (show.limits[[1]]) {
      lines(lim_data$extremes[, 1], lim_data$extremes[, 2],
         col = lim.col[1], lty = lim.lty[1], lwd = lim.lwd[1])
   }

   if (show.limits[[2]]) {
      lines(lim_data$outliers[, 1], lim_data$outliers[, 2],
         col = lim.col[2], lty = lim.lty[2], lwd = lim.lwd[2])
   }
}
