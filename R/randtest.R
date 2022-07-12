#' Randomization test for PLS regression
#'
#' @description
#' \code{randtest} is used to carry out randomization/permutation test for a PLS regression model
#'
#' @param x
#' matrix with predictors.
#' @param y
#' vector or one-column matrix with response.
#' @param ncomp
#' maximum number of components to test.
#' @param center
#' logical, center or not predictors and response values.
#' @param scale
#' logical, scale (standardize) or not predictors and response values.
#' @param nperm
#' number of permutations.
#' @param sig.level
#' significance level.
#' @param silent
#' logical, show or not test progress.
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#'
#' @return
#' Returns an object of \code{randtest} class with following fields:
#' \item{nperm }{number of permutations used for the test.}
#' \item{stat }{statistic values calculated for each component.}
#' \item{alpha }{alpha values calculated for each component.}
#' \item{statperm }{matrix with statistic values for each permutation.}
#' \item{corrperm }{matrix with correlation between predicted and reference y-vales for each
#' permutation.}
#' \item{ncomp.selected }{suggested number of components.}
#'
#' @details
#' The class implements a method for selection of optimal number of components in PLS1 regression
#' based on the randomization test [1]. The basic idea is that for each component from 1 to
#' \code{ncomp} a statistic T, which is a covariance between t-score (X score, derived from a PLS
#' model) and the reference Y values, is calculated. By repeating this for randomly permuted
#' Y-values a distribution of the statistic is obtained. A parameter \code{alpha} is computed to
#' show how often the statistic T, calculated for permuted Y-values, is the same or higher than
#' the same statistic, calculated for original data without permutations.
#'
#' If a component is important, then the covariance for unpermuted data should be larger than the
#' covariance for permuted data and therefore the value for \code{alpha} will be quie small (there
#' is still a small chance to get similar covariance). This makes \code{alpha} very similar to
#' p-value in a statistical test.
#'
#' The \code{randtest} procedure calculates alpha for each component, the values can be observed
#' using \code{summary} or \code{plot} functions. There are also several function, allowing e.g.
#' to show distribution of statistics and the critical value for each component.
#'
#' @references
#' S. Wiklund et al. Journal of Chemometrics 21 (2007) 427-439.
#'
#' @seealso
#' Methods for \code{randtest} objects:
#' \tabular{ll}{
#'  \code{print.randtest} \tab prints information about a \code{randtest} object.\cr
#'  \code{\link{summary.randtest}} \tab shows summary statistics for the test.\cr
#'  \code{\link{plot.randtest}} \tab shows bar plot for alpha values.\cr
#'  \code{\link{plotHist.randtest}} \tab shows distribution of statistic plot.\cr
#'  \code{\link{plotCorr.randtest}} \tab shows determination coefficient plot.\cr
#' }
#'
#' @examples
#' ### Examples of using the test
#'
#' ## Get the spectral data from Simdata set and apply SNV transformation
#'
#' data(simdata)
#'
#' y = simdata$conc.c[, 3]
#' x = simdata$spectra.c
#' x = prep.snv(x)
#'
#' ## Run the test and show summary
#' ## (normally use higher nperm values > 1000)
#' r = randtest(x, y, ncomp = 4, nperm = 200, silent = FALSE)
#' summary(r)
#'
#' ## Show plots
#'
#' par( mfrow = c(3, 2))
#' plot(r)
#' plotHist(r, ncomp = 3)
#' plotHist(r, ncomp = 4)
#' plotCorr(r, 3)
#' plotCorr(r, 4)
#' par( mfrow = c(1, 1))
#'
#' @export
randtest <- function(x, y, ncomp = 15, center = TRUE, scale = FALSE, nperm = 1000, sig.level = 0.05,
                    silent = TRUE, exclcols = NULL, exclrows = NULL) {
   x <- as.matrix(x)
   y <- as.matrix(y)

   # remove excluded columns and rows
   if (length(exclcols) > 0) {
      x <- mda.exclcols(x, exclcols)
      x <- x[, -attr(x, "exclcols"), drop = FALSE]
   }

   if (length(exclrows) > 0) {
      x <- mda.exclrows(x, exclrows)
      exclrows <- attr(x, "exclrows")
      x <- x[-exclrows, , drop = FALSE]
      y <- y[-exclrows, , drop = FALSE]
   }

   nobj <- nrow(x)

   x <- prep.autoscale(as.matrix(x), center = center, scale = scale)
   y <- prep.autoscale(as.matrix(y), center = center, scale = scale)

   stat <- matrix(0, ncol = ncomp, nrow = 1)
   alpha <- matrix(0, ncol = ncomp, nrow = 1)
   statperm <- matrix(0, ncol = ncomp, nrow = nperm)
   corrperm <- matrix(0, ncol = ncomp, nrow = nperm)

   m <- NULL
   for (icomp in seq_len(ncomp)) {

      if (!silent) fprintf("Permutations for component #%d...\n", icomp)

      if (icomp > 1) {
         x <- x - xscores %*% t(m$xloadings)
         y <- y - xscores %*% t(m$yloadings)
      }

      m <- pls.simpls(x, y, 1)
      xscores <- x %*% (m$weights %*% solve(crossprod(m$xloadings, m$weights)))

      stat[icomp] <- (t(xscores) %*% y) / nobj
      for (iperm in seq_len(nperm)) {
         yp <- y[sample(1:nobj)]
         mp <- pls.simpls(x, yp, 1)
         pxscores <- x %*% (mp$weights %*% solve(crossprod(mp$xloadings, mp$weights)))
         statperm[iperm, icomp] <- crossprod(pxscores, yp) / nobj
         corrperm[iperm, icomp] <- cor(y, yp)
      }

      alpha[icomp] <- sum(statperm[, icomp] > stat[icomp]) / nperm
   }

   ncomp.selected <- max(which(alpha <= sig.level))
   colnames(alpha) <- colnames(stat) <- paste("Comp", seq_len(ncomp))
   colnames(statperm) <- colnames(corrperm) <- paste("Comp", seq_len(ncomp))
   rownames(statperm) <- rownames(corrperm) <- seq_len(nperm)
   rownames(alpha) <- "Alpha"
   rownames(stat) <- "Statistic"

   res <- list(
      nperm = nperm,
      stat = stat,
      alpha = alpha,
      statperm = statperm,
      corrperm = corrperm,
      ncomp.selected = ncomp.selected
   )

   res$call <- match.call()
   class(res) <- "randtest"

   return(res)
}


#' Histogram plot for randomization test results
#'
#' @description
#' Makes a histogram for statistic values distribution for particular component, also
#' show critical value as a vertical line.
#'
#' @param obj
#' results of randomization test (object of class `randtest`)
#' @param ncomp
#' number of component to make the plot for
#' @param bwd
#' width of bars (between 0 and 1)
#' @param ...
#' other optional arguments
#'
#' @details
#' See examples in help for \code{\link{randtest}} function.
#'
#' @export
plotHist.randtest <- function(obj, ncomp = obj$ncomp.selected, bwd = 0.9, ...) {

   h <- hist(obj$statperm[, ncomp], plot = FALSE)
   plot_data <- h$counts
   attr(plot_data, "xaxis.values") <- h$mids
   attr(plot_data, "xaxis.name") <- "Test statistic"
   attr(plot_data, "yaxis.name") <- "Frequency"
   attr(plot_data, "name") <- sprintf("Distribution for permutations (ncomp = %d)", ncomp)

   mdaplot(plot_data, type = "h", show.lines = c(obj$stat[ncomp], NA), bwd = bwd, ...)
}

#' Correlation plot for randomization test results
#'
#' @description
#' Makes a plot with statistic values vs. coefficient of determination between permuted
#' and reference y-values.
#'
#' @param obj
#' results of randomization test (object of class `randtest`)
#' @param ncomp
#' number of component to make the plot for
#' @param ylim
#' limits for y axis
#' @param xlab
#' label for x-axis
#' @param ylab
#' label for y-axis
#' @param ...
#' other optional arguments
#'
#' @details
#' See examples in help for \code{\link{randtest}} function.
#'
#' @export
plotCorr.randtest <- function(obj, ncomp = obj$ncomp.selected, ylim = NULL,
   xlab = expression(r^2), ylab = "Test statistic", ...) {

   plot_data <- list(
      "perm" = cbind(obj$corrperm[, ncomp]^2, obj$statperm[, ncomp]),
      "est" = cbind(1, obj$stat[, ncomp])
   )


   if (is.null(ylim)) {
      ylim <- c(
         min(plot_data[[1]][, 2], plot_data[[2]][, 2]),
         max(plot_data[[1]][, 2], plot_data[[2]][, 2])
      )
   }

   attr(plot_data[[1]], "name") <- sprintf("Permutations (ncomp = %d)", ncomp)
   colnames(plot_data[[1]]) <- c(expression(r^2), "Test statistic")
   mdaplotg(plot_data, type = "p", ylim = ylim, xlab = xlab, ylab = ylab,
      legend.position = "bottomright", ...)

   fit_data <- rbind(apply(plot_data[[1]], 2, mean), plot_data[[2]])
   lines(fit_data[, 1], fit_data[, 2], col = rgb(0.6, 0.6, 0.6), lty = 2, lwd = 0.75)
}

#' Plot for randomization test results
#'
#' @description
#' Makes a bar plot with alpha values for each component.
#'
#' @param x
#' results of randomization test (object of class `randtest`)
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other optional arguments
#'
#' @details
#' See examples in help for \code{\link{randtest}} function.
#'
#' @export
plot.randtest <- function(x, main = "Alpha", xlab = "Components", ylab = "", ...) {
   mdaplot(x$alpha, show.lines = c(NA, 0.05), type = "h", main = main, xlab = xlab,
      ylab = ylab, ...)
}

#' Summary method for randtest object
#'
#' @description
#' Shows summary for randomization test results.
#'
#' @param object
#' randomization test results (object of class \code{randtest})
#' @param ...
#' other arguments
#'
#' @export
summary.randtest <- function(object, ...) {
   data <- rbind(object$alpha, object$stat)
   cat("\nSummary for permutation test results\n")
   fprintf("Number of permutations: %d\n", object$nperm)
   fprintf("Suggested number of components: %d\n", object$ncomp.selected)
   cat("\nStatistics and alpha values:\n")
   show(data)
   cat("\n")
}

#' Print method for randtest object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' a randomization test results (object of class \code{randtest})
#' @param ...
#' other arguments
#'
#' @export
print.randtest <- function(x, ...) {
   cat("\nRandomization test results (class randtest)\n")
   cat("\nCall:\n")
   print(x$call)
   cat("\nMajor fields:\n")
   cat("$nperm - number of permutations\n")
   cat("$ncomp.selected - number of selected components (suggested)\n")
   cat("$alpha - vector with alpha values calculated for each component.\n")
   cat("$stat - vector with statistic values calculated for each component.\n")
   cat("$statperm - matrix with statistic values for each permutation.\n")
   cat("$corrperm - correlations between predicted and reference y-vales for permutations.\n")
   cat("\nTry summary(obj) and plot(obj) to see the test results.\n")
}
