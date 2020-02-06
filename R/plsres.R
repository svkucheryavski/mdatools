#' PLS results
#'
#' @description
#' \code{plsres} is used to store and visualize results of applying a PLS model to a new data.
#'
#' @param y.pred
#' predicted y values.
#' @param y.ref
#' reference (measured) y values.
#' @param ncomp.selected
#' selected (optimal) number of components.
#' @param xdecomp
#' PLS decomposition of X data (object of class \code{ldecomp}).
#' @param ydecomp
#' PLS decomposition of Y data (object of class \code{ldecomp}).
#' @param info
#' information about the object.
#'
#' @details
#' Do not use \code{plsres} manually, the object is created automatically when one applies a PLS
#' model to a new data set, e.g. when calibrate and validate a PLS model (all calibration and
#' validation results in PLS model are stored as objects of \code{plsres} class) or use function
#' \code{\link{predict.pls}}.
#'
#' The object gives access to all PLS results as well as to the plotting methods for visualisation
#' of the results. The \code{plsres} class also inherits all properties and methods of \code{regres}
#'  - general class for regression results.
#'
#' If no reference values provided, regression statistics will not be calculated and most of the
#' plots not available. The class is also used for cross-validation results, in this case some of
#' the values and methods are not available (e.g. scores and scores plot, etc.).
#'
#' All plots are based on \code{\link{mdaplot}} function, so most of its options can be used (e.g.
#' color grouping, etc.).
#'
#' RPD is ratio of standard deviation of response values to standard error of prediction (SDy/SEP).
#'
#' @return
#' Returns an object of \code{plsres} class with following fields:
#' \item{ncomp }{number of components included to the model.}
#' \item{ncomp.selected }{selected (optimal) number of components.}
#' \item{y.ref }{a matrix with reference values for responses.}
#' \item{y.pred }{a matrix with predicted values for responses.}
#' \item{rmse }{a matrix with root mean squared error values for each response and component.}
#' \item{slope }{a matrix with slope values for each response and component.}
#' \item{r2 }{a matrix with determination coefficients for each response and component.}
#' \item{bias }{a matrix with bias values for each response and component.}
#' \item{sep }{a matrix with standard error values for each response and component.}
#' \item{rpd }{a matrix with RPD values for each response and component.}
#' \item{xdecomp }{decomposition of predictors (object of class \code{ldecomp}).}
#' \item{ydecomp }{decomposition of responses (object of class \code{ldecomp}).}
#' \item{info }{information about the object.}
#'
#' @seealso
#'
#' Methods for \code{plsres} objects:
#' \tabular{ll}{
#'    \code{print} \tab prints information about a \code{plsres} object.\cr
#'    \code{\link{summary.plsres}} \tab shows performance statistics for the results.\cr
#'    \code{\link{plot.plsres}} \tab shows plot overview of the results.\cr
#'    \code{\link{plotXScores.plsres}} \tab shows scores plot for x decomposition.\cr
#'    \code{\link{plotXYScores.plsres}} \tab shows scores plot for x and y decomposition.\cr
#'    \code{\link{plotXVariance.plsres}} \tab shows explained variance plot for x decomposition.\cr
#'    \code{\link{plotYVariance.plsres}} \tab shows explained variance plot for y decomposition.\cr
#'    \code{\link{plotXCumVariance.plsres}} \tab shows cumulative explained variance plot for y
#'    decomposition.\cr
#'    \code{\link{plotYCumVariance.plsres}} \tab shows cumulative explained variance plot for y
#'    decomposition.\cr
#'    \code{\link{plotXResiduals.plsres}} \tab shows T2 vs. Q plot for x decomposition.\cr
#'    \code{\link{plotYResiduals.plsres}} \tab shows residuals plot for y values.\cr
#' }
#'
#' Methods inherited from \code{regres} class (parent class for \code{plsres}):
#' \tabular{ll}{
#'    \code{\link{plotPredictions.regres}} \tab shows predicted vs. measured plot.\cr
#'    \code{\link{plotRMSE.regres}} \tab shows RMSE plot.\cr
#' }
#'
#' See also \code{\link{pls}} - a class for PLS models.
#'
#' @examples
#' ### Examples of using PLS result class
#' library(mdatools)

#' ## 1. Make a PLS model for concentration of first component
#' ## using full-cross validation and get calibration results
#'
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#'
#' model = pls(x, y, ncomp = 8, cv = 1)
#' model = selectCompNum(model, 2)
#' res = model$calres
#'
#' summary(res)
#' plot(res)
#'
#' ## 2. Make a PLS model for concentration of first component
#' ## and apply model to a new dataset
#'
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#'
#' model = pls(x, y, ncomp = 6, cv = 1)
#' model = selectCompNum(model, 2)
#'
#' x.new = simdata$spectra.t
#' y.new = simdata$conc.t[, 1]
#' res = predict(model, x.new, y.new)
#'
#' summary(res)
#' plot(res)
#'
#' ## 3. Show variance and error plots for PLS results
#' par(mfrow = c(2, 2))
#' plotXCumVariance(res, type = 'h')
#' plotYCumVariance(res, type = 'b', show.labels = TRUE, legend.position = 'bottomright')
#' plotRMSE(res)
#' plotRMSE(res, type = 'h', show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' ## 4. Show scores plots for PLS results
#' ## (for results plot we can use color grouping)
#' par(mfrow = c(2, 2))
#' plotXScores(res)
#' plotXScores(res, show.labels = TRUE, cgroup = y.new)
#' plotXYScores(res)
#' plotXYScores(res, comp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' ## 5. Show predictions and residuals plots for PLS results
#' par(mfrow = c(2, 2))
#' plotXResiduals(res, show.label = TRUE, cgroup = y.new)
#' plotYResiduals(res, show.label = TRUE)
#' plotPredictions(res)
#' plotPredictions(res, ncomp = 4, xlab = 'C, reference', ylab = 'C, predictions')
#' par(mfrow = c(1, 1))
#'
#' @export
plsres <- function(y.pred, y.ref = NULL, ncomp.selected = dim(y.pred)[2], xdecomp = NULL,
   ydecomp = NULL, info = "") {

   obj <- regres(y.pred, y.ref = y.ref, ncomp.selected = ncomp.selected)
   obj$ncomp <- dim(y.pred)[2]
   obj$xdecomp <- xdecomp
   obj$ydecomp <- ydecomp
   obj$info <- info
   obj$ncomp.selected <- ncomp.selected

   obj$call <- match.call()
   class(obj) <- c("plsres", "regres")

   return(obj)
}

#' as.matrix method for PLS results
#'
#' @description
#' Returns a matrix with model performance statistics for PLS results
#'
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' number of components to calculate the statistics for
#' @param ny
#' for which response variable calculate the statistics for
#' @param ...
#' other arguments
#'
#' @export
as.matrix.plsres <- function(x, ncomp = NULL, ny = 1, ...) {

   xdecomp_res <- if (!is.null(x$xdecomp)) as.matrix(x$xdecomp) else matrix(NA, x$ncomp, 2)
   ydecomp_res <- if (!is.null(x$ydecomp)) as.matrix(x$ydecomp) else matrix(NA, x$ncomp, 2)

   out <- cbind(
      xdecomp_res,
      ydecomp_res,
      as.matrix.regres(x, ny = ny)
   )

   rownames(out) <- paste("Comp", seq_len(x$ncomp))
   colnames(out)[1:4] <- c("X expvar", "X cumexpvar", "Y expvar", "Y cumexpvar")

   if (!is.null(ncomp)) {
      out <- out[ncomp, , drop = FALSE]
   }

   return(out)
}

#' summary method for PLS results object
#'
#' @description
#' Shows performance statistics for the results.
#'
#' @param object
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' for which response variable show the summary for
#' @param ...
#' other arguments
#'
#' @export
summary.plsres <- function(object, ny = seq_len(object$nresp), ncomp = NULL, ...) {

   cat("\nPLS regression results (class plsres) summary\n")
   fprintf("Info: %s\n", object$info)
   fprintf("Number of selected components: %d\n", object$ncomp.selected)

   if (is.null(object$y.ref)) {
      cat("No reference data provided to calculate prediction performance.")
      return()
   }

   if (length(ncomp) == 1) {
      fprintf("\nNumber of selected components: %d\n", ncomp)
   }

   for (y in ny) {
      fprintf("\nResponse variable %s:\n", dimnames(object$y.pred)[[3]][y])
      out <- as.matrix.plsres(object, ny = y, ncomp = ncomp)
      if (!any(is.na(out[, 1:4]))) out[, 1:4] <- round(out[, 1:4], 3)
      out[, 5] <- round(out[, 5], 3)
      out[, 6] <- mdaplot.formatValues(out[, 6], round.only = T)
      out[, 7] <- round(out[, 7], 3)
      out[, 8] <- round(out[, 8], 4)
      out[, 9] <- round(out[, 9], 2)
      print(out)
   }
}

#' print method for PLS results object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ...
#' other arguments
#'
#' @export
print.plsres <- function(x, ...) {
   cat("\nPLS results (class plsres)\n")
   cat("\nCall:\n")
   print(x$call)

   cat("\nMajor fields:\n")
   cat("$ncomp.selected - number of selected components\n")
   cat("$y.pred - array with predicted y values\n")

   if (!is.null(x$y.ref)) {
      cat("$y.ref - matrix with reference y values\n")
      cat("$rmse - root mean squared error\n")
      cat("$r2 - coefficient of determination\n")
      cat("$slope - slope for predicted vs. measured values\n")
      cat("$bias - bias for prediction vs. measured values\n")
      cat("$ydecomp - decomposition of y values (ldecomp object)\n")
   }

   cat("$xdecomp - decomposition of x values (ldecomp object)\n")
}

################################
#  Plotting methods            #
################################

#' Explained X variance plot for PLS results
#'
#' @description
#' Shows plot with explained X variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param decomp
#' which dcomposition to use ("xdecomp" or "ydecomp")
#' @param variance
#' which variance to use ("expvar", "cumexpvar")
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotVariance.plsres <- function(obj, decomp = "xdecomp", variance = "expvar", ...) {
   if (is.null(obj[[decomp]])) return(NULL)
   return(plotVariance.ldecomp(obj[[decomp]], variance = variance, ...))
}

#' Explained X variance plot for PLS results
#'
#' @description
#' Shows plot with explained X variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXVariance.plsres <- function(obj, main = "Variance (X)", ...) {
   return(plotVariance.plsres(obj, decomp = "xdecomp", main = main, ...))
}

#' Explained Y variance plot for PLS results
#'
#' @description
#' Shows plot with explained Y variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotYVariance.plsres <- function(obj, main = "Variance (Y)", ...) {
   return(plotVariance.plsres(obj, decomp = "ydecomp", main = main, ...))
}

#' Explained cumulative X variance plot for PLS results
#'
#' @description
#' Shows plot with cumulative explained X variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXCumVariance.plsres <- function(obj, main = "Cumulative variance (X)", ...) {
   return(plotVariance.plsres(obj, decomp = "xdecomp", variance = "cumexpvar", main = main, ...))
}

#' Explained cumulative Y variance plot for PLS results
#'
#' @description
#' Shows plot with cumulative explained Y variance vs. number of components.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotYCumVariance.plsres <- function(obj, main = "Cumulative variance (Y)", ...) {
   return(plotVariance.plsres(obj, decomp = "ydecomp", variance = "cumexpvar", main = main, ...))
}

#' X scores plot for PLS results
#'
#' @description
#' Shows plot with scores values for PLS decomposition of x data.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXScores.plsres <- function(obj, comp = c(1, 2), main = "Scores (X)", ...) {
   if (is.null(obj$xdecomp)) return(invisible(NULL))
   return(plotScores.ldecomp(obj$xdecomp, comp = comp, main = main, ...))
}

#' XY scores plot for PLS results
#'
#' @description
#' Shows plot with X vs. Y scores values for PLS results.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' which component to show the plot for
#' @param show.plot
#' logical, show plot or just return plot data
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXYScores.plsres <- function(obj, ncomp = 1, show.plot = TRUE, ...) {

   if (is.null(obj$xdecomp) || is.null(obj$ydecomp)) return(invisible(NULL))

   if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
      stop("Wrong value for ncomp argument.")
   }

   plot_data <- cbind(
      obj$xdecomp$scores[, ncomp, drop = FALSE],
      obj$ydecomp$scores[, ncomp, drop = FALSE]
   )

   plot_data <- mda.setattr(plot_data, mda.getattr(obj$xdecomp$scores))
   rownames(plot_data) <- rownames(obj$xdecomp$scores)
   colnames(plot_data) <- c(
      sprintf("X-scores (Comp %d, %.2f%%)", ncomp, obj$xdecomp$expvar[ncomp]),
      sprintf("Y-scores (Comp %d, %.2f%%)", ncomp, obj$ydecomp$expvar[ncomp])
   )

   attr(plot_data, "name") <- "Scores (XY)"

   if (!show.plot) {
      return(plot_data)
   }

   return(mdaplot(plot_data, type = "p", ...))
}

#' X residuals plot for PLS results
#'
#' @description
#' Shows a plot with Q residuals vs. Hotelling T2 values for PLS decomposition of x data.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param main
#' main title for the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plotXResiduals.plsres <- function(obj, ncomp = obj$ncomp.selected, norm = TRUE, log = FALSE,
   main = sprintf("X-residuals (ncomp = %d)", ncomp), ...) {

   if (is.null(obj$xdecomp)) return(invisible(NULL))

   if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
      stop("Wrong value for ncomp argument.")
   }

   return(plotResiduals.ldecomp(obj$xdecomp, ncomp = ncomp, main = main,
      norm = norm, log = log, ...))
}

#' Y residuals plot for PLS results
#'
#' @description
#' Shows a plot with Y residuals vs reference Y values for selected component.
#'
#' @param obj
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param main
#' main title for the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' Proxy for \code{\link{plotResiduals.regres}} function.
#'
#' @export
plotYResiduals.plsres <- function(obj, ncomp = obj$ncomp.selected,
   main = sprintf("Y-residuals (ncomp = %d)", ncomp), ...) {

   if (is.null(obj$y.ref)) return(invisible(NULL))
   return(plotResiduals.regres(obj, ncomp = ncomp, main = main, ...))
}


#' Overview plot for PLS results
#'
#' @description
#' Shows a set of plots for PLS results.
#'
#' @param x
#' PLS results (object of class \code{plsres})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' which y variable to show the summary for (if NULL, will be shown for all)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{plsres}} function.
#'
#' @export
plot.plsres <- function(x, ncomp = x$ncomp.selected, ny = 1, show.labels = FALSE, ...) {

   if (is.null(x$y.ref)) {
      par(mfrow = c(1, 2))
      plotXResiduals(x, ...)
      plotPredictions.regres(x, ncomp = ncomp, ny = ny, ...)
      par(mfrow = c(1, 1))
      return()
   }

   par(mfrow = c(2, 2))
   plotXResiduals(x, ncomp = ncomp, ...)
   plotYVariance(x, ...)
   plotRMSE(x, ny = ny, ...)
   plotPredictions.regres(x, ncomp = ncomp, ny = ny, ...)
   par(mfrow = c(1, 1))
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
plotXYResiduals.plsres <- function(obj, ncomp = obj$ncomp.selected, norm = TRUE, log = FALSE,
   show.labels = FALSE, labels = "names", show.plot = TRUE, ...) {

   attrs <- mda.getattr(obj$xdecomp$Q)

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
   h0 <- if (!is.null(attr(obj$xdecomp$T2, "u0"))) attr(obj$xdecomp$T2, "u0")[[ncomp]]
   q0 <- if (!is.null(attr(obj$xdecomp$Q, "u0"))) attr(obj$xdecomp$Q, "u0")[[ncomp]]
   z0 <- if (!is.null(attr(obj$ydecomp$Q, "u0"))) attr(obj$ydecomp$Q, "u0")[[ncomp]]

   # get DoF factors
   Nh <- if (!is.null(attr(obj$xdecomp$T2, "Nu"))) attr(obj$xdecomp$T2, "Nu")[[ncomp]]
   Nq <- if (!is.null(attr(obj$xdecomp$Q, "Nu"))) attr(obj$xdecomp$Q, "Nu")[[ncomp]]

   # get distances
   h <- obj$xdecomp$T2[, ncomp]
   q <- obj$xdecomp$Q[, ncomp]
   z <- obj$ydecomp$Q[, ncomp]

   # compute full distance for X
   f <- Nh * h / h0 + Nq * q / q0
   f0 <- Nh + Nq

   # prepare plot data
   f <- transform(f, f0, norm, log)
   z <- transform(z, z0, norm, log)

   # default values for local labels
   lxlab <- get_label("f", norm, log)
   lylab <- get_label("z", norm, log)

   # combine everything to dataset and assign attributes
   plot_data <- mda.cbind(f, z)
   plot_data <- mda.setattr(plot_data, attrs, "row")

   rownames(plot_data) <- rownames(obj$xdecomp$Q)
   colnames(plot_data) <- c(
      paste0("Full X-distance, ", lxlab),
      paste0("Y-distance, ", lylab)
   )

   attr(plot_data, "name") <- sprintf("XY-distances (ncomp = %d)", ncomp)

   # if no plot required - return plot series object
   if (!show.plot) return(plot_data)

   # show plot
   return(mdaplot(plot_data, ...))
}
