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
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotCumVariance.ldecomp <- function(obj, type = "b", main = "Cumulative variance",
   xlab = "Components", ylab = "Explained variance, %", show.labels = FALSE,
   labels = "values", show.plot = TRUE, ...) {

   return(
      plotVariance(obj, variance = "cumexpvar", main = main, xlab = xlab, ylab = ylab,
            type = type, show.labels = show.labels, labels = labels, show.plot = show.plot, ...)
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
#' @param variance
#' string, which variance to make the plot for ("expvar", "cumexpvar")
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotVariance.ldecomp <- function(obj, type = "b", variance = "expvar", main = "Variance",
   xlab = "Components", ylab = "Explained variance, %", show.labels = FALSE, labels = "values",
   show.plot = TRUE, ...) {

   if (!show.plot) {
      return(obj[[variance]])
   }

   return(
      mdaplot(obj[[variance]], main = main, xticks = seq_len(obj$ncomp), xlab = xlab, ylab = ylab,
           show.labels = show.labels, labels = labels, type = type, ...)
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
#' @param main
#' main title for the plot
#' @param type
#' type of the plot
#' @param xlab
#' label for x-axis.
#' @param ylab
#' label for y-axis.
#' @param show.legend
#' logical, show or not a legend on the plot (needed in case of line or bar plot).
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param labels
#' what to show as labels
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotScores.ldecomp <- function(obj, comp = c(1, 2), main = "Scores", type = "p", xlab = NULL,
   ylab = NULL, show.labels = FALSE, labels = "names", show.legend = TRUE, show.axes = TRUE,
   show.plot = TRUE, ...) {

   # get scores for given components and generate column names with explained variance
   plot_data <- mda.subset(obj$scores, select = comp)
   colnames(plot_data) <- paste0("Comp ", comp, " (", round(obj$expvar[comp], 2), "%)")

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
      p <- mdaplot(plot_data, type = type, show.labels = show.labels, labels = labels,
         show.lines = show.lines, main = main, xlab = xlab, ylab = ylab, ...)

      return(invisible(p))
   }

   # line or bar plot
   if (is.null(ylab)) ylab <- "Score value"
   if (type == "h") show.lines <- FALSE

   return(mdaplotg(mda.t(plot_data), type = type, show.labels = show.labels, labels = labels,
      show.lines = show.lines, show.legend = show.legend, main = main, xlab = xlab,
      ylab = ylab, ...))
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
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param labels
#' what to show as labels if necessary
#' @param main
#' main title for the plot
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotResiduals.ldecomp <- function(obj, ncomp = obj$ncomp.selected,
   show.labels = FALSE, labels = "names", main = NULL, show.plot = TRUE, ...) {

   attrs <- mda.getattr(obj$Q)

   # prepare plot data
   h <- obj$T2[, ncomp]
   q <- obj$Q[, ncomp]

   # combine everything to dataset and assign attributes
   plot_data <- mda.cbind(h, q)
   plot_data <- mda.setattr(plot_data, attrs, "row")
   rownames(plot_data) <- rownames(obj$Q)
   colnames(plot_data) <- c("Score distance, h", "Orthogonal distance, q")

   # if no plot required - return plot series object
   if (!show.plot) {
      return(plot_data)
   }

   # set up main title for the plot
   if (is.null(main)) {
      main <- if (is.null(ncomp)) "Residuals" else sprintf("Residuals (ncomp = %d)", ncomp)
   }

   # show plot
   return(mdaplot(plot_data, main = main, show.labels = show.labels, labels = labels, ...))
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
#' @param ...
#' other arguments
#'
#' @export
as.matrix.ldecomp <- function(x, ...) {
   data <- cbind(x$expvar, x$cumexpvar)
   rownames(data) <- colnames(x$Q)
   colnames(data) <- c("Expvar", "Cumexpvar")
   return(data)
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
      T2[, i] <- rowSums(scoresn[, 1:i, drop = F]^2)
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
