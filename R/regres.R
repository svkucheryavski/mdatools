#' Regression results
#'
#' @description
#' Class for storing and visualisation of regression predictions
#'
#' @param y.pred
#' vector or matrix with y predicted values
#' @param y.ref
#' vector with reference (measured) y values
#' @param ncomp.selected
#' if y.pred calculated for different components, which to use as default
#'
#' @return
#' a list (object of \code{regres} class) with fields, including:
#' \tabular{ll}{
#'    \code{y.pred} \tab a matrix with predicted values \cr
#'    \code{y.ref} \tab a vector with reference (measured) values \cr
#'    \code{ncomp.selected} \tab selected column/number of components for predictions \cr
#'    \code{rmse} \tab root mean squared error for predicted vs measured values \cr
#'    \code{slope} \tab slope for predicted vs measured values \cr
#'    \code{r2} \tab coefficient of determination for predicted vs measured values \cr
#'    \code{bias} \tab bias for predicted vs measured values \cr
#'    \code{rpd} \tab RPD values \cr
#' }
#'
#' @export
regres <- function(y.pred, y.ref = NULL, ncomp.selected = 1) {

   if (is.null(y.pred) || length(dim(y.pred)) != 3) {
      stop("Parameter 'y.pred' should be a 3-way array.")
   }

   if (ncomp.selected > dim(y.pred)[2]) {
      stop("Wrong value for 'ncomp.selected' parameter.")
   }

   if (!is.null(y.ref)) y.ref <- as.matrix(y.ref)

   obj <- list()
   obj$y.pred <- y.pred
   obj$y.ref <- y.ref
   obj$ncomp <- dim(y.pred)[2]
   obj$ncomp.selected <- ncomp.selected
   obj$nresp <- dim(y.pred)[3]
   obj$respnames <- dimnames(y.pred)[[3]]
   if (is.null(obj$respnames)) obj$respnames <- paste0("y", seq_len(obj$nresp))

   obj <- c(obj, regres.getPerformanceStats(y.pred, y.ref))
   obj$call <- match.call()
   class(obj) <- "regres"

   return(obj)
}

#' as.matrix method for regression results
#'
#' @description
#' Returns a matrix with model performance statistics for regression results
#'
#' @param x
#' regression results (object of class \code{regres})
#' @param ncomp
#' model complexity (number of components) to calculate the statistics for (can be a vector)
#' @param ny
#' for which response variable calculate the statistics for
#' @param ...
#' other arguments
#'
#' @export
as.matrix.regres <- function(x, ncomp = NULL, ny = 1, ...) {
   if (is.null(x$y.ref)) return()

   out <- cbind(x$r2[ny, ], x$rmse[ny, ], x$slope[ny, ],
      x$bias[ny, ], x$rpd[ny, ])

   colnames(out) <- c("R2", "RMSE", "Slope", "Bias", "RPD")
   rownames(out) <- dimnames(x$y.pred)[[2]]

   if (!is.null(ncomp)) {
      out <- out[ncomp, , drop = FALSE]
   }

   return(out)
}

#' summary method for regression results object
#'
#' @description
#' Shows performance statistics for the regression results.
#'
#' @param object
#' regression results (object of class \code{regres})
#' @param ncomp
#' model complexity to show the summary for (if NULL - shows for all available values)
#' @param ny
#' for which response variable show the summary for (one value or a vector)
#' @param ...
#' other arguments
#'
#' @export
summary.regres <- function(object, ncomp = object$ncomp.selected, ny = seq_len(object$nresp), ...) {

   cat("\nRegression results (class regres) summary\n")
   if (is.null(object$y.ref)) {
      cat("No reference data provided to calculate prediction performance.")
      return()
   }

   if (!is.null(ncomp)) {
      fprintf("\nNumber of selected components: %d\n\n", ncomp)
   }

   for (i in ny) {
      fprintf("\nResponse variable %s:\n", object$respnames[i])
      print(as.matrix.regres(object, ny = i, ncomp = ncomp))
   }
}

#' print method for regression results object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' regression results (object of class \code{regres})
#' @param ...
#' other arguments
#'
#' @export
print.regres <- function(x, ...) {
   cat("\nRegression results (class regres)\n")
   cat("\nCall:\n")
   print(x$call)

   cat("\nMajor fields:\n")
   cat("$y.pred - matrix or vector with predicted y values\n")
   if (!is.null(x$y.ref)) {
      cat("$y.ref - vector with reference y values\n")
      cat("$rmse - root mean squared error\n")
      cat("$r2 - coefficient of determination\n")
      cat("$slope - slope for predicted vs. measured values\n")
      cat("$bias - bias for prediction vs. measured values\n")
   }

   if (ncol(x$y.pred) > 1) {
      cat("$ncomp.selected - number of selected components\n")
   }
}


################################
#  Static methods              #
################################


regres.getPerformanceStats <- function(y.pred, y.ref) {
   if (is.null(y.ref)) return(NULL)

   attrs <- mda.getattr(y.pred)

   # remove excluded rows so they are not counted
   # when calculating statistics
   if (length(attrs$exclrows) > 0) {
      y.pred <- y.pred[-attrs$exclrows, , , drop = F]
      y.ref <- y.ref[-attrs$exclrows, , drop = F]
   }

   # residuals (errors) based statistics
   err <- regres.err(y.pred, y.ref)
   ytot <- colSums(scale(y.ref, center = TRUE, scale = FALSE)^2)

   stats <- list(
      "r2" = regres.r2(err, ytot),
      "bias" = regres.bias(err),
      "rmse" = regres.rmse(err)
   )

   stats$slope <- regress.addattrs(regres.slope(y.pred, y.ref), attributes(err), "Slope")
   stats$sep <- regress.addattrs(sqrt(stats$rmse^2 - stats$bias^2), attributes(err), "SEP")
   stats$rpd <- regress.addattrs(apply(y.ref, 2, sd) / stats$sep, attributes(err), "RPD")

   return(stats)
}

#' Add names and attributes to matrix with statistics
#'
#' @param stat
#' matrix with statistics
#' @param attrs
#' attributes from error matrix
#' @param name
#' name of statistic
#'
regress.addattrs <- function(stat, attrs, name) {

   attr(stat, "name") <- name
   dimnames(stat) <- attrs$dimnames[c(3, 2)]
   attr(stat, "xaxis.name") <- attrs$xaxis.name
   attr(stat, "yaxis.name") <- attrs$yaxis.name

   return(stat)
}

#' Error of prediction
#'
#' @description
#' Calculates array of differences between predicted and reference values.
#'
#' @param y.ref
#' vector with reference values
#' @param y.pred
#' matrix with predicted values
#'
regres.err <- function(y.pred, y.ref) {
   ncomp <- dim(y.pred)[2]
   err <- array(repmat(y.ref, ncomp, 1), dim = dim(y.pred)) - y.pred

   attr(err, "name") <- "Error of prediction"
   attr(err, "xaxis.name") <- "Components"
   attr(err, "yaxis.name") <- "Predictors"

   return(err)
}

#' Determination coefficient
#'
#' @description
#' Calculates matrix with coeffient of determination for every response and components
#'
#' @param err
#' vector with difference between reference and predicted y-values
#' @param ytot
#' total variance for y-values
#'
#' @export
regres.r2 <- function(err, ytot) {
   r2 <- t(1 - scale(colSums(err^2), center = F, scale = ytot))
   return(regress.addattrs(r2, attributes(err), "Coefficient of determination"))
}

#' Prediction bias
#'
#' @description
#' Calculates matrix with bias (average prediction error) for every response and components
#'
#' @param err
#' vector with difference between reference and predicted y-values
#'
regres.bias <- function(err) {
   bias <- t(colSums(err) / nrow(err))
   return(regress.addattrs(bias, attributes(err), "Bias"))
}

#' RMSE
#'
#' @description
#' Calculates matrix with root mean squared error of prediction for every response and components.
#'
#' @param err
#' vector with difference between reference and predicted y-values
#'
regres.rmse <- function(err) {
   rmse <- t(sqrt(colSums(err^2) / nrow(err)))
   return(regress.addattrs(rmse, attributes(err), "RMSE"))
}

#' Slope
#'
#' @description
#' Calculates matrix with slope of predicted and measured values for every response and components.
#'
#' @param y.ref
#' vector with reference values
#' @param y.pred
#' matrix with predicted values
#'
regres.slope <- function(y.pred, y.ref) {
   nresp <- ncol(y.ref)
   ncomp <- ncol(y.pred)
   slope <- matrix(0, nrow = nresp, ncol = ncomp)
   for (i in seq_len(nresp)) {
      slope[i, ] <- matrix(coefficients(lm(y.pred[, , i] ~ y.ref[, i])), nrow = 2)[2, ]
   }

   return(slope)
}


################################
#  Plotting methods            #
################################


#' Predictions plot for regression results
#'
#' @description
#' Shows plot with predicted y values.
#'
#' @param obj
#' regression results (object of class \code{regres})
#' @param ny
#' number of predictor to show the plot for (if y is multivariate)
#' @param ncomp
#' complexity of model (e.g. number of components) to show the plot for
#' @param show.line
#' logical, show or not line fit for the plot points
#' @param show.stat
#' logical, show or not legend with statistics on the plot
#' @param stat.col
#' color of text in legend with statistics
#' @param stat.cex
#' size of text in legend with statistics
#' @param xlim
#' limits for x-axis (if NULL will be computed automatically)
#' @param ylim
#' limits for y-axis (if NULL will be computed automatically)
#' @param axes.equal
#' logical, make limits for x and y axes equal or not
#' @param show.plot
#' logical, show plot or just return plot data
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' If reference values are available, the function shows a scatter plot with predicted vs.
#' reference values, otherwise predicted values are shown vs. object numbers.
#'
#' @export
plotPredictions.regres <- function(obj, ny = 1, ncomp = obj$ncomp.selected, show.line = TRUE,
   show.stat = FALSE, stat.col = "#606060", stat.cex = 0.85, xlim = NULL, ylim = NULL,
   axes.equal = TRUE, show.plot = TRUE, ...) {

   if (length(ny) != 1) {
      stop("You can show prediction plot only for one selected response variable.")
   }

   if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
      stop("Wrong value for ncomp argument.")
   }

   if (is.null(obj$y.ref)) {
      plot_data <- matrix(obj$y.pred[, ncomp, ny], ncol = 1)
      attr(plot_data, "xaxis.name") <- colnames(plot_data) <-
         paste0(obj$respnames[ny], ", predicted")
      attr(plot_data, "yaxis.name") <- attr(obj$y.pred, "yaxis.name")
      attr(plot_data, "yaxis.values") <- attr(obj$y.pred, "yaxis.values")
   } else {
      plot_data <- cbind(obj$y.ref[, ny], obj$y.pred[, ncomp, ny])
      colnames(plot_data) <- c(
         paste0(obj$respnames[ny], ", reference"),
         paste0(obj$respnames[ny], ", predicted")
      )
   }

   plot_data <- mda.setattr(plot_data, mda.getattr(obj$y.pred))
   rownames(plot_data) <- rownames(obj$y.pred)
   attr(plot_data, "name") <- paste0("Predictions (ncomp = ", ncomp, ")")

   if (!show.plot) {
      return(plot_data)
   }

   if (axes.equal && !is.null(obj$yref)) {
      xlim <- ylim <- range(plot_data)
   }

   p <- mdaplot(plot_data, type = "p", xlim = xlim, ylim = ylim, ...)

   if (is.null(obj$y.ref)) {
      return(invisible(p))
   }

   if (show.stat) {
      stat.text <- sprintf("nLV = %d\nRMSE = %.3f\nR2 = %.3f\n",
         ncomp, obj$rmse[ny, ncomp], obj$r2[ny, ncomp])
      text(p$xlim[1], p$ylim[2], stat.text, adj = c(0, 1), col = stat.col, cex = stat.cex)
   }

   if (show.line) {
      plotRegressionLine(p)
   }

   return(invisible(p))
}

#' Residuals plot for regression results
#'
#' @description
#' Shows plot with Y residuals (difference between predicted and reference values) for selected
#' response variable and complexity (number of components).
#'
#' @param obj
#' regression results (object of class \code{regres})
#' @param ny
#' number of predictor to show the plot for (if y is multivariate)
#' @param ncomp
#' complexity of model (e.g. number of components) to show the plot for
#' @param show.lines
#' allows to show the horisontal line at y = 0
#' @param show.plot
#' logical, show plot or just return plot data
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @export
plotResiduals.regres <- function(obj, ny = 1, ncomp = obj$ncomp.selected,
   show.lines = c(NA, 0), show.plot = TRUE, ...) {

   if (is.null(obj$y.ref)) {
      stop("Y-residuals can not be plotted without reference values.")
   }

   if (length(ny) != 1) {
      stop("You can make residuals plot only for one selected response variable.")
   }

   if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
      stop("Wrong value for ncomp argument.")
   }

   plot_data <- cbind(obj$y.ref[, ny], obj$y.ref[, ny] - obj$y.pred[, ncomp, ny])
   plot_data <- mda.setattr(plot_data, mda.getattr(obj$y.pred))
   colnames(plot_data) <- c(
      sprintf("%s, reference", obj$respnames[ny]),
      sprintf("%s, residuals", obj$respnames[ny])
   )
   attr(plot_data, "name") <- sprintf("Y-residuals (ncomp = %d)", ncomp)

   if (!show.plot) {
      return(plot_data)
   }

   return(mdaplot(plot_data, type = "p", show.lines = show.lines, ...))
}

#' RMSE plot for regression results
#'
#' @description
#' Shows plot with RMSE values vs. model complexity (e.g. number of components).
#'
#' @param obj
#' regression results (object of class \code{regres})
#' @param ny
#' number of predictor to show the plot for (if y is multivariate)
#' @param type
#' type of the plot
#' @param xticks
#' vector with ticks for x-axis
#' @param labels
#' what to use as labels ("names", "values" or "indices")
#' @param show.plot
#' logical, show plot or just return plot data
#' @param ylab
#' label for y-axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @export
plotRMSE.regres <- function(obj, ny = 1, type = "b", xticks = seq_len(obj$ncomp),
   labels = "values", show.plot = TRUE, ylab = paste0("RMSE (", obj$respnames[ny], ")"), ...) {

   if (is.null(obj$rmse)) {
      stop("RMSE values are not available.")
   }

   if (length(ny) != 1) {
      stop("You can make residuals plot only for one selected response variable.")
   }

   plot_data <- mda.subset(obj$rmse, ny)
   attr(plot_data, "name") <- "RMSE"

   if (!show.plot) {
      return(plot_data)
   }

   return(mdaplot(plot_data, type = type, xticks = xticks, labels = labels, ylab = ylab, ...))
}

#' Plot method for regression results
#'
#' @details
#' This is a shortcut for \code{\link{plotPredictions.regres}}
#'
#' @param x
#' regression results (object of class \code{regres})
#' @param ...
#' other arguments
#'
#' @export
plot.regres <- function(x, ...) {
   plotPredictions.regres(x, ...)
}
