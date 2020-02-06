#' Regression coefficients
#'
#' @description
#' class for storing and visualisation of regression coefficients
#' for regression models
#'
#' @param coeffs
#' array (npred x ncomp x nresp) with regression coefficients
#' @param ci.coeffs
#' array (npred x ncomp x nresp x cv) with regression coefficients for
#' computing confidence intervals (e.g. from cross-validation) using Jack-Knifing method
#' @param use.mean
#' logical, tells how to compute standard error for regression coefficients. If \code{TRUE}
#' mean values for ci.coeffs is computed first. If \code{FALSE}, \code{values} (coefficients
#' computed for global model) are used as mean.
#'
#' @return
#' a list (object of \code{regcoeffs} class) with fields, including:
#' \tabular{ll}{
#'    \code{values} \tab an array (nvar x ncomp x ny) with regression coefficients \cr
#'    \code{se} \tab an array (nvar x ncomp x ny) with standard errors for the coefficients \cr
#'    \code{t.values} \tab an array (nvar x ncomp x ny) with t-values for the coefficients \cr
#'    \code{p.values} \tab an array (nvar x ncomp x ny) with p-values for coefficients \cr
#' }
#'
#' last three fields are available if parameter \code{ci.coeffs} was provided.
#'
#' Check also \code{\link{confint.regcoeffs}}, \code{\link{summary.regcoeffs}} and
#' \code{\link{plot.regcoeffs}}.
#'
#' @export
regcoeffs <- function(coeffs, ci.coeffs = NULL, use.mean = TRUE) {

   if (is.null(dim(coeffs)) || length(dim(coeffs)) != 3) {
      stop("Coefficients must be provided as 3-way array.")
   }

   obj <- list()
   obj$values <- coeffs
   obj$nvar <- dim(coeffs)[1]

   # assign response number and names
   obj$nresp <- dim(coeffs)[3]
   obj$respnames <- dimnames(coeffs)[[3]]
   if (is.null(obj$respnames)) obj$respnames <- paste0("y", seq_len(obj$nresp))

   # add statistics and class name
   obj <- c(obj, regcoeffs.getStats(coeffs, ci.coeffs, use.mean))
   obj$call <- match.call()
   class(obj) <- "regcoeffs"

   return(obj)
}

#' Confidence intervals for regression coefficients
#'
#' @description
#' returns matrix with confidence intervals for regression coeffocoents
#' for given response number and number of components.
#'
#' @param object
#' regression coefficients object (class \code{regcoeffs})
#' @param parm
#' not used, needed for compatiility with general method
#' @param level
#' confidence level
#' @param ncomp
#' number of components (one value)
#' @param ny
#' index of response variable (one value)
#' @param ...
#' other arguments
#'
#' @export
confint.regcoeffs <- function(object, parm = NULL, level = 0.95, ncomp = 1, ny = 1, ...) {

   if (length(ncomp) != 1) {
      stop("Parameter 'ncomp' should be just one value.")
   }

   if (ncomp < 1 || ncomp > dim(object$values)[2]) {
      stop("Wrong value for parameter 'ncomp'.")
   }

   if (length(ny) != 1) {
      stop("Parameter 'ny' should be just one value.")
   }

   if (ny < 1 || ny > dim(object$values)[3]) {
      stop("Wrong value for parameter 'ny'.")
   }

   alpha <- 1 - level
   t <- -qt(alpha / 2, object$DoF)
   ci <- repmat(object$se[, ncomp, ny], 1, 2) %*% diag(c(-t, t)) +
      repmat(object$values[, ncomp, ny], 1, 2)

   if (length(attr(object$values, "exclrows"))) {
      ci[attr(object$values, "exclrows"), ] <- 0
   }

   ci <- mda.setattr(ci, mda.getattr(object$values))
   rownames(ci) <- dimnames(object$values)[[1]]
   colnames(ci) <- paste0(round(c(alpha / 2, alpha / 2 + level) * 100, 1), "%")
   attr(ci, "name") <- "Confidence interval"

   return(ci)
}

#' as.matrix method for regression coefficients class
#'
#' @description
#' returns matrix with regression coeffocoents for given response number and amount of components
#'
#' @param x
#' regression coefficients object (class \code{regcoeffs})
#' @param ncomp
#' number of components to return the coefficients for
#' @param ny
#' number of response variable to return the coefficients for
#' @param ...
#' other arguments
#'
#' @export
as.matrix.regcoeffs <- function(x, ncomp = 1, ny = 1, ...) {
   return(x$values[, ncomp, ny])
}

#' Summary method for regcoeffs object
#'
#' @description
#' Shows estimated coefficients and statistics (if available).
#'
#' @param object
#' object of class \code{regcoeffs}
#' @param ncomp
#' how many components to use
#' @param ny
#' which y variable to show the summary for
#' @param alpha
#' significance level for confidence interval (if statistics available)
#' @param ...
#' other arguments
#'
#' @details
#' Statistcs are shown if Jack-Knifing was used when model is calibrated.
#'
#' @export
summary.regcoeffs <- function(object, ncomp = 1, ny = 1, alpha = 0.05, ...) {

   if (length(ncomp) != 1) {
      stop("Parameter 'ncomp' should be just one value.")
   }

   if (ncomp < 1 || ncomp > dim(object$values)[2]) {
      stop("Wrong value for parameter 'ncomp'.")
   }

   if (length(ny) != 1) {
      stop("Parameter 'ny' should be just one value.")
   }

   if (ny < 1 || ny > dim(object$values)[3]) {
      stop("Wrong value for parameter 'ny'.")
   }

   attrs <- mda.getattr(object$values)
   coeffs <- object$values[, ncomp, ny, drop = F]
   dim(coeffs) <- c(dim(object$values)[1], 1)
   colnames(coeffs)[1] <- "Coeffs"
   if (!is.null(object$se)) {
      coeffs <- cbind(
         coeffs,
         object$se[, ncomp, ny],
         round(object$t.values[, ncomp, ny], 2),
         round(object$p.values[, ncomp, ny], 3),
         confint(object, ncomp = ncomp, ny = ny, level = 1 - alpha)
      )
      colnames(coeffs)[1:4] <- c("Coeffs", "Std. err.", "t-value", "p-value")
   }

   rownames(coeffs) <- dimnames(object$values)[[1]]
   attr(coeffs, "exclrows") <- attrs$exclrows
   attr(coeffs, "name") <- paste0(
      "Regression coefficients for ", object$respnames[ny], " (ncomp = ", ncomp, ")"
   )

   cat("\n")
   mda.show(coeffs)
   if (ncol(coeffs) > 1) {
      cat(sprintf("\nDegrees of freedom (Jack-Knifing): %d\n", object$DoF))
   }
   cat("\n")
}

#' print method for regression coefficients class
#'
#' @description
#' prints regression coeffocoent values for given response number and amount of components
#'
#' @param x
#' regression coefficients object (class \code{regcoeffs})
#' @param ...
#' other arguments
#'
#' @export
print.regcoeffs <- function(x, ...) {
   cat("\nRegression coefficients (class regcoeffs)\n")
   cat("\nCall:\n")
   print(x$call)
   cat("\nMajor fields:\n")
   cat("$values - array with regression coefficients\n")
   cat("$se - array with standard errors\n")
   cat("$t.values - array with t-values\n")
   cat("$p.values - array with p-values\n")
   cat("$DoF - degrees of freedom for Jack-Knifing\n")
   cat("\nThe last four fields available only if Jack-Knifing was used.\n\n")
}


################################
#  Static methods              #
################################


#' Distribution statistics for regression coeffificents
#'
#' @description
#' calculates standard error, t-values and p-values for
#' regression coefficients based on Jack-Knifing method.
#'
#' @param coeffs
#' array (npred x ncomp x nresp) with regression coefficients
#' @param ci.coeffs
#' array (npred x ncomp x nresp x cv) with regression coefficients for
#' computing confidence intervals (e.g. from cross-validation) using Jack-Knifing method
#' @param use.mean
#' logical, tells how to compute standard error for regression coefficients. If \code{TRUE}
#' mean values for ci.coeffs is computed first. If \code{FALSE}, \code{values} (coefficients
#' computed for global model) are used as mean.
#'
#' @return
#' a list with statistics three arrays: srandard error, t-values and p-values computed for
#' each regression coefficient.
#'
#' @export
regcoeffs.getStats <- function(coeffs, ci.coeffs = NULL, use.mean = TRUE) {

   if (is.null(ci.coeffs)) return()

   if (is.null(dim(ci.coeffs)) || length(dim(ci.coeffs)) != 4) {
      stop("Coefficients for distribution statistics must be provided as 4-way array.")
   }

   # get attributes and prepare arrays
   nseg <- dim(ci.coeffs)[4]
   dim_names <- dimnames(coeffs)
   DoF <- dim(ci.coeffs)[4] - 1
   attrs <- mda.getattr(coeffs)
   se <- p.values <- t.values <- array(0, dim = dim(coeffs))

   # to make sure p-values for excluded predictors are 1 (not important)
   p.values <- p.values + 1

   # prepare correct indices for predictors taking into accound excluded variables
   row_ind <- seq_len(dim(coeffs)[1])
   if (length(attrs$exclrows) > 0) {
      row_ind <- row_ind[-attrs$exclrows]
      coeffs <- coeffs[-attrs$exclrows, , , drop = FALSE]
   }

   # check the dimension
   if (any(dim(coeffs) != dim(ci.coeffs)[1:3])) {
      stop("Dimension of coefficients for distribution statistics does not much the 'coeffs'.")
   }

   # compute mean if needed
   if (use.mean) {
      coeffs <- apply(ci.coeffs, 1:3, mean)
   }

   # compute main statistics
   err <- sweep(ci.coeffs, 1:3, coeffs, "-")
   ssq <- apply(err^2, 1:3, sum)
   se[row_ind, , ] <- sqrt((nseg - 1) / nseg * ssq)
   t.values[row_ind, , ] <- coeffs / se[row_ind, , , drop = FALSE]
   p.values[row_ind, , ] <- 2 * (1 - pt(abs(t.values[row_ind, , , drop = FALSE]), nseg - 1))

   # add names and attributes
   dimnames(se) <- dimnames(t.values) <- dimnames(p.values) <- dim_names
   t.values <- mda.setattr(t.values, attrs)
   p.values <- mda.setattr(p.values, attrs)
   se <- mda.setattr(se, attrs)
   attr(se, "name") <- "Standard error (Jack-knifing)"
   attr(t.values, "name") <- "t-values (Jack-knifing)"
   attr(p.values, "name") <- "p-values (Jack-knifing)"

   return(
      list(
         se = se,
         t.values = t.values,
         p.values = p.values,
         DoF = DoF
      )
   )
}


################################
#  Plotting methods            #
################################


#' Regression coefficients plot
#'
#' @description
#' Shows plot with regression coefficient values for every predictor variable (x)
#'
#' @param x
#' regression coefficients object (class \code{regcoeffs})
#' @param ncomp
#' number of components to use for creating the plot
#' @param ny
#' index of response variable to make the plot for
#' @param type
#' type of the plot
#' @param col
#' vector with two colors for the plot (one is used to show real coefficient and another one to
#' show confidence intervals)
#' @param show.lines
#' allows to show horizontal line at c(NA, 0)
#' @param show.ci
#' logical, show or not confidence intervals if they are available
#' @param alpha
#' significance level for confidence intervals (a number between 0 and 1, e.g.
#' for 95\% alpha = 0.05)
#' @param ...
#' other arguments for plotting methods (e.g. main, xlab, etc)
#'
#' @export
plot.regcoeffs <- function(x, ncomp = 1, ny = 1, type = (if (x$nvar > 30) "l" else "h"),
   col = c(mdaplot.getColors(1), "lightgray"), show.lines = c(NA, 0), show.ci = FALSE,
   alpha = 0.05, ...) {

   if (length(ncomp) != 1) {
      stop("Parameter 'ncomp' should be just one value.")
   }

   if (ncomp < 1 || ncomp > dim(x$values)[2]) {
      stop("Wrong value for parameter 'ncomp'.")
   }

   if (length(ny) != 1) {
      stop("Parameter 'ny' should be just one value.")
   }

   if (ny < 1 || ny > dim(x$values)[3]) {
      stop("Wrong value for parameter 'ny'.")
   }

   plot_data <- matrix(x$values[, ncomp, ny], nrow = 1)
   attr(plot_data, "exclcols") <- attr(x$values, "exclrows")
   attr(plot_data, "xaxis.name") <- attr(x$values, "yaxis.name")
   attr(plot_data, "xaxis.values") <- attr(x$values, "yaxis.values")

   colnames(plot_data) <- rownames(x$values)
   attr(plot_data, "name") <- paste0("Regression coefficients (ncomp = ", ncomp, ")")
   attr(plot_data, "yaxis.name") <- paste0("Coefficients (", x$respnames[ny], ")")

   if (!(show.ci && !is.null(x$se))) {
      return(mdaplot(plot_data, col = col[1], type = type, show.lines = show.lines, ...))
   }

   ci <- confint(x, ncomp = ncomp, ny = ny, level = 1 - alpha)
   if (type == "l") {
      return(
         mdaplot(mda.rbind(plot_data, mda.t(ci)), type = "l", col = col[c(2, 1, 1)],
            show.lines = show.lines, ...)
      )
   }

   err <- (ci[, 2] - ci[, 1]) / 2
   mdaplotg(list(plot_data, mda.rbind(plot_data, err)), type = c(type, "e"), show.legend = FALSE,
      col = col[c(2, 1)], show.grid = T, show.lines = show.lines, ...)
}
