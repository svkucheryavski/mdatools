## class and methods for iPLS ##

#' Variable selection with interval PLS
#'
#' @description
#' Applies iPLS alrogithm to find variable intervals most important for
#' prediction
#'
#' @param x
#' a matrix with predictor values
#' @param y
#' a vector with response values
#' @param glob.ncomp
#' maximum number of components for a global PLS model
#' @param center
#' logical, center or not the data values
#' @param scale
#' logical, standardize or not the data values
#' @param cv
#' cross-validation settings (see details)
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param int.ncomp
#' maximum number of components for interval PLS models
#' @param int.num
#' number of intervals
#' @param int.width
#' width of intervals
#' @param int.limits
#' a two column matrix with manual intervals specification
#' @param int.niter
#' maximum number of iterations (if NULL it will be the same as number of intervals)
#' @param ncomp.selcrit
#' criterion for selecting optimal number of components ('min' for minimum of RMSECV)
#' @param method
#' iPLS method (\code{'forward'} or \code{'backward'})
#' @param x.test
#' matrix with predictors for test set (by default is NULL, if specified, is used instead of cv).
#' @param y.test
#' matrix with responses for test set.
#' @param silent
#' logical, show or not information about selection process
#'
#' @return
#' object of 'ipls' class with several fields, including:
#'    \item{var.selected}{a vector with indices of selected variables}
#'    \item{int.selected}{a vector with indices of selected intervals }
#'    \item{int.num}{total number of intervals}
#'    \item{int.width}{width of the intervals}
#'    \item{int.limits}{a matrix with limits for each interval}
#'    \item{int.stat}{a data frame with statistics for the selection algorithm}
#'    \item{glob.stat}{a data frame with statistics for the first step (individual intervals)}
#'    \item{gm}{global PLS model with all variables included}
#'    \item{om}{optimized PLS model with selected variables}
#'
#' @details
#' The algorithm splits the predictors into several intervals and tries to find a combination
#' of the intervals, which gives best prediction performance. There are two selection methods:
#' "forward" when the intervals are successively included, and "backward" when the intervals
#' are successively excluded from a model. On the first step the algorithm finds the best
#' (forward) or the worst (backward) individual interval. Then it tests the others to find the
#' one which gives the best model in a combination with the already selected/excluded one. The
#' procedure continues until the maximum number of iteration is reached.
#'
#' There are several ways to specify the intervals. First of all either number of intervals
#' (\code{int.num}) or width of the intervals (\code{int.width}) can be provided. Alternatively
#' one can specify the limits (first and last variable number) of the intervals manually
#' with \code{int.limits}.
#'
#' Cross-validation settings, \code{cv}, can be a number or a list. If \code{cv} is a number, it
#' will be used as a number of segments for random cross-validation (if \code{cv = 1}, full
#' cross-validation will be preformed). If it is a list, the following syntax can be used:
#' \code{cv = list('rand', nseg, nrep)} for random repeated cross-validation with \code{nseg}
#' segments and \code{nrep} repetitions or \code{cv = list('ven', nseg)} for systematic splits
#' to \code{nseg} segments ('venetian blinds').
#'
#' @references
#' [1] Lars Noergaard at al.  Interval partial least-squares regression (iPLS): a
#' comparative chemometric study with an example from near-infrared spectroscopy.
#' Appl.Spec. 2000; 54: 413-419
#'
#' @examples
#' library(mdatools)
#'
#' ## forward selection for simdata
#'
#' data(simdata)
#' Xc = simdata$spectra.c
#' yc = simdata$conc.c[, 3, drop = FALSE]
#'
#' # run iPLS and show results
#' im = ipls(Xc, yc, int.ncomp = 5, int.num = 10, cv = 4, method = "forward")
#' summary(im)
#' plot(im)
#'
#' # show "developing" of RMSECV during the algorithm execution
#' plotRMSE(im)
#'
#' # plot predictions before and after selection
#' par(mfrow = c(1, 2))
#' plotPredictions(im$gm)
#' plotPredictions(im$om)
#'
#' # show selected intervals on spectral plot
#' ind = im$var.selected
#' mspectrum = apply(Xc, 2, mean)
#' plot(simdata$wavelength, mspectrum, type = 'l', col = 'lightblue')
#' points(simdata$wavelength[ind], mspectrum[ind], pch = 16, col = 'blue')
#'
#' @export
ipls <- function(x, y, glob.ncomp = 10, center = TRUE, scale = FALSE, cv = list("ven", 10),
   exclcols = NULL, exclrows = NULL,  int.ncomp = glob.ncomp, int.num = NULL, int.width = NULL,
   int.limits = NULL, int.niter = NULL, ncomp.selcrit = "min", method = "forward",
   x.test = NULL, y.test = NULL, silent = FALSE) {

   # process names and values for xaxis
   x <- mda.df2mat(x)
   xaxis.name <- attr(x, "xaxis.name")
   xaxis.values <- attr(x, "xaxis.values")

   if (is.null(xaxis.name)) xaxis.name <- "Variables"
   if (is.null(xaxis.values)) xaxis.values <- seq_len(ncol(x))
   if (is.null(dim(y))) dim(y) <- c(length(y), 1)

   # if test set is provided check it and set cv to NULL
   if (!is.null(x.test) && !is.null(y.test)) {
      if (is.null(dim(y.test))) dim(y.test) <- c(length(y.test), 1)
      if (nrow(x.test) != nrow(y.test)) stop("Number of rows in 'x.test' and 'y.test' should be the same")
      if (ncol(x.test) != ncol(x)) stop("Number of columns in 'x.test' and 'x' should be the same")
      if (ncol(y.test) != ncol(y)) stop("Number of columns in 'y.test' and 'y' should be the same")
      cv <- NULL
   }

   # remove excluded columns and rows
   if (length(exclcols) > 0) {
      x <- mda.exclcols(x, exclcols)
      exclcols <- attr(x, "exclcols")
      x <- x[, -exclcols, drop = FALSE]
      xaxis.values <- xaxis.values[-exclcols]
      attr(x, "exclcols") <- NULL

      if (!is.null(x.test)) {
         x.test <- x.test[, -exclcols, drop = FALSE]
         attr(x.test, "exclcols") <- NULL
      }
   }

   if (length(exclrows) > 0) {
      x <- mda.exclrows(x, exclrows)
      exclrows <- attr(x, "exclrows")
      x <- x[-exclrows, , drop = FALSE]
      y <- y[-exclrows, , drop = FALSE]
      attr(x, "exclrows") <- NULL
      attr(y, "exclrows") <- NULL
   }

   attr(x, "xaxis.values") <- xaxis.values
   attr(x, "xaxis.name") <- xaxis.name

   if (length(dim(y)) != 2) {
      stop("Response variable (y) should be a matrix or a sequence.")
   }

   if (ncol(y) > 1) {
      warning("iPLS can work with one y-variable at time, selecting first column.")
      y <- y[, 1, drop = FALSE]
   }

   # add name to the single column
   if (is.null(colnames(y))) colnames(y) <- "y1"

   # get number of predictors
   npred <- ncol(x)

   # get matrix with limits
   int.limits <- ipls.getintlimits(int.limits, int.width, int.num, npred)

   # check limits
   if (ncol(int.limits) != 2 || nrow(int.limits) < 2) {
      stop("Interval limits shall be provided as a matrix with two columns.")
   }

   # get difference (size of intervals)
   df <- apply(int.limits, 1, diff)

   # check that difference and interval values are also fine
   if (min(int.limits) < 1 || max(int.limits) > npred || any(df < 0)) {
      stop("Wrong values for interval limits.")
   }

   # recompute interval numbers and width (average)
   int.num <- nrow(int.limits)
   int.width <- max(df) + 1

   # add names
   rownames(int.limits) <- seq_len(int.num)
   colnames(int.limits) <- c("Left", "Right")

   # define number of iterations
   int.niter <- if (!is.null(int.niter) || int.num < 30) int.num else 30

   # build an object
   obj <- list(
      cv = cv,
      glob.ncomp = glob.ncomp,
      int.ncomp = int.ncomp,
      xmean = apply(x, 2, "mean"),
      int.width = int.width,
      int.num = int.num,
      int.niter = int.niter,
      int.limits = int.limits,
      center = center,
      scale = scale,
      method = method,
      ncomp.selcrit = ncomp.selcrit,
      silent = silent,
      xaxis.name = xaxis.name,
      xaxis.values = xaxis.values,
      x.test = x.test,
      y.test = y.test
   )

   # make a global model
   obj$gm <- pls(x, y, ncomp = obj$glob.ncomp, center = obj$center, scale = obj$scale, cv = obj$cv,
      ncomp.selcrit = obj$ncomp.selcrit, x.test = obj$x.test, y.test = obj$y.test)


   # get statistics for global model
   gmres <- if (is.null(obj$cv)) obj$gm$res$test else obj$gm$res$cv

   glob.stat <- data.frame(
      "n" = 0,
      "start" = 1,
      "end" = ncol(x),
      "nComp" = obj$gm$ncomp.selected,
      "RMSE" = gmres$rmse[1, obj$gm$ncomp.selected],
      "R2" = gmres$r2[1, obj$gm$ncomp.selected]
   )

   # initialize statistics for local models
   int.stat <- data.frame(
      "n" = 0,
      "start" = 1,
      "end" = ncol(x),
      "selected" = FALSE,
      "nComp" = obj$gm$ncomp.selected,
      "RMSE" = gmres$rmse[1, obj$gm$ncomp.selected],
      "R2" = gmres$r2[1, obj$gm$ncomp.selected]
   )

   # do iPLS
   obj <- switch(method,
      "forward" = ipls.forward(x, y, obj, int.stat, glob.stat),
      "backward" = ipls.backward(x, y, obj, int.stat, glob.stat),
      stop("Wrong value for parameter 'method'.")
   )

   if (!is.null(obj$int.stat)) {
      rownames(obj$int.stat) <- seq_len(nrow(obj$int.stat))
      obj$int.stat$R2 <- round(obj$int.stat$R2, 3)
   }


   # prepare vector with selected variables
   int.limits <- obj$int.limits[obj$int.selected, , drop = FALSE]
   obj$var.selected <- apply(int.limits, 1, function(x) seq(x[1], x[2]))
   if (is.list(obj$var.selected)) obj$var.selected <- do.call(c, obj$var.selected)
   names(obj$var.selected) <- dim(obj$var.selected) <- NULL
   obj$var.selected <- sort(obj$var.selected)

   obj$om <- pls(x[, obj$var.selected, drop = FALSE], y, ncomp = obj$int.ncomp, center = obj$center,
      scale = obj$scale, cv = obj$cv, ncomp.selcrit = obj$ncomp.selcrit,
      x.test = obj$x.test[, obj$var.selected, drop = FALSE], y.test = obj$y.test)

   obj$call <- match.call()
   class(obj) <- "ipls"

   return(obj)
}

ipls.getintlimits <- function(int.limits, int.width, int.num, npred) {

   if (is.null(c(int.limits, int.width, int.num))) {
      stop("Specify either interval width or number of intervals.")
   }

   if (!is.null(int.limits)) {
      return(int.limits)
   }

   int.width <- if (is.null(int.num)) int.width else round(npred / int.num)
   int.num <- round(npred / int.width)

   # generate interval limits similar to the way from [1]
   if (int.num == npred) {
      return(cbind(seq_len(npred), seq_len(npred)))
   }

   varrem <- npred %% int.num
   nvarrem <- trunc(npred / int.num)

   int1 <- if (varrem == 0) seq(1, npred, by = nvarrem) else
         c(
            seq(1, (varrem - 1) * (nvarrem + 1) + 1, by = nvarrem + 1),
            seq((varrem - 1) * (nvarrem + 1) + 2 + nvarrem, npred, by = nvarrem)
         )
   int2 <- c(int1[2:int.num] - 1, npred)

   return(cbind(int1, int2))
}


#' Runs the forward iPLS algorithm
#'
#' @param x
#' a matrix with predictor values
#' @param y
#' a vector with response values
#' @param obj
#' object with initial settings for iPLS algorithm
#' @param int.stat
#' data frame with initial interval statistics
#' @param glob.stat
#' data frame with initial global statistics
#'
ipls.forward <- function(x, y, obj, int.stat, glob.stat) {

   # define vectors with status, selected and non-selected intervals
   int.nonselected <- seq_len(obj$int.num)
   int.selected <- NULL


   if (!obj$silent) {
      fprintf("\nModel with all intervals: RMSE = %f, nLV = %d\n",
         int.stat$RMSE[1], int.stat$nComp[1])
   }

   # do loop for max number of intervals
   selind <- NULL
   rmse <- Inf
   for (i in seq_len(obj$int.niter)) {
      if (!obj$silent) fprintf("Iteration %3d/%3d... ", i, obj$int.niter)

      sel <- NULL
      for (l in int.nonselected) {

         # combine already selected intervals with the current
         ind <- obj$int.limits[l, 1]:obj$int.limits[l, 2]
         xc <- x[, c(selind, ind), drop = FALSE]
         xt <- if (!is.null(obj$x.test)) obj$x.test[, c(selind, ind), drop = FALSE] else NULL

         # build a model
         m <- pls(xc, y, ncomp = obj$int.ncomp, center = obj$center, scale = obj$scale, cv = obj$cv,
            ncomp.selcrit = obj$ncomp.selcrit, x.test = xt, y.test = obj$y.test)

         lres <- if (is.null(obj$cv)) m$res$test else m$res$cv

         # if first round, build a data frame with statistics for each interval
         if (i == 1) {
            glob.stat <- rbind(glob.stat, data.frame(
               "n" = l,
               "start" = obj$int.limits[l, 1],
               "end" = obj$int.limits[l, 2],
               "nComp" = m$ncomp.selected,
               "RMSE" = lres$rmse[1, m$ncomp.selected],
               "R2" = lres$r2[1, m$ncomp.selected]
            ))
         }

         # else check if rmse has been improved
         if (rmse > lres$rmse[1, m$ncomp.selected]) {
            ncomp <- m$ncomp.selected
            rmse <- lres$rmse[1, m$ncomp.selected]
            r2 <- lres$r2[1, m$ncomp.selected]
            sel <- l
         }
      }


      if (is.null(sel)) {
         if (!obj$silent) cat("no improvements, stop.\n\n")
         break
      }

      selind <- c(selind, obj$int.limits[sel, 1]:obj$int.limits[sel, 2])
      int.nonselected <- int.nonselected[int.nonselected != sel]
      int.selected <- c(int.selected, sel)
      int.stat <- rbind(int.stat, data.frame(
         "n" = sel,
         "start" = obj$int.limits[sel, 1],
         "end" = obj$int.limits[sel, 2],
         "selected" = TRUE,
         "nComp" = ncomp,
         "RMSE" = rmse,
         "R2" = r2
      ))

      if (!obj$silent) {
         fprintf("selected interval %3d (RMSE = %f, nLV = %d)\n", sel, rmse, m$ncomp.selected)
      }

   }

   # return the selection results
   obj$glob.stat <- glob.stat
   obj$int.stat <- int.stat
   obj$int.selected <- int.selected

   return(obj)
}

#' Runs the backward iPLS algorithm
#'
#' @param x
#' a matrix with predictor values
#' @param y
#' a vector with response values
#' @param obj
#' object with initial settings for iPLS algorithm
#' @param int.stat
#' data frame with initial interval statistics
#' @param glob.stat
#' data frame with initial global statistics
#'
ipls.backward <- function(x, y, obj, int.stat, glob.stat) {

   # define vectors with status, selected and non-selected intervals
   int.selected <- seq_len(obj$int.num)
   int.nonselected <- NULL

   if (!obj$silent) {
      fprintf("\nModel with all intervals: RMSE = %f, nLV = %d\n",
         int.stat$RMSE[1], int.stat$nComp[1])
   }

   # do loop for max number of intervals
   unselind <- NULL
   rmse <- Inf
   for (i in seq_len(obj$int.niter)) {

      if (length(int.selected) == 1) break

      if (!obj$silent) {
         fprintf("Iteration %3d/%3d... ", i, obj$int.niter)
      }

      # do loop to select an interval
      unsel <- NULL
      for (l in int.selected) {
         ind <- obj$int.limits[l, 1]:obj$int.limits[l, 2]

         # combine already selected intervals with the current
         xc <- x[, -c(unselind, ind), drop = FALSE]
         xt <- if (!is.null(obj$x.test)) obj$x.test[, -c(unselind, ind), drop = FALSE] else NULL

         # build a model
         m <- pls(xc, y, ncomp = obj$int.ncomp, center = obj$center, scale = obj$scale, cv = obj$cv,
            ncomp.selcrit = obj$ncomp.selcrit, x.test = xt, y.test = obj$y.test)

         lres <- if (is.null(obj$cv)) m$res$test else m$res$cv

         # if first round, build a data frame with statistics for each interval
         if (i == 1) {
            glob.stat <- rbind(glob.stat, data.frame(
               "n" = l,
               "start" = obj$int.limits[l, 1],
               "end" = obj$int.limits[l, 2],
               "nComp" = m$ncomp.selected,
               "RMSE" = lres$rmse[1, m$ncomp.selected],
               "R2" = lres$r2[1, m$ncomp.selected]
            ))
         }

         # if last two intervals are left keep them both
         if (length(int.selected) == 2) {
            int.stat <- rbind(int.stat, data.frame(
               "n" = l,
               "start" = obj$int.limits[l, 1],
               "end" = obj$int.limits[l, 2],
               "selected" = FALSE,
               "nComp" = m$ncomp.selected,
               "RMSE" = lres$rmse[1, m$ncomp.selected],
               "R2" = lres$r2[1, m$ncomp.selected]
            ))
            unsel <- NULL
            break
         }

         # else check if rmse has been improved
         if (rmse > lres$rmse[1, m$ncomp.selected]) {
            ncomp <- m$ncomp.selected
            rmse <- lres$rmse[1, m$ncomp.selected]
            r2 <- lres$r2[1, m$ncomp.selected]
            unsel <- l
         }
      }

      if (is.null(unsel)) {
         if (!obj$silent) cat("no improvements, stop.\n\n")
         break
      }

      unselind <- c(unselind, obj$int.limits[unsel, 1]:obj$int.limits[unsel, 2])
      int.selected <- int.selected[int.selected != unsel]
      int.nonselected <- c(int.nonselected, unsel)

      int.stat <- rbind(int.stat, data.frame(
         "n" = unsel,
         "start" = obj$int.limits[unsel, 1],
         "end" = obj$int.limits[unsel, 2],
         "selected" = FALSE,
         "nComp" = ncomp,
         "RMSE" = rmse,
         "R2" = r2
      ))

      if (!obj$silent) {
         fprintf("excluded interval %3d (RMSE = %f, nLV = %d)\n", unsel, rmse, m$ncomp.selected)
      }
   }

   # return the selection results
   obj$glob.stat <- glob.stat
   obj$int.stat <- int.stat
   obj$int.selected <- int.selected

   return(obj)
}

#' iPLS performance plot
#'
#' @description
#' Shows PLS performance for each selected or excluded intervals at the
#' first iteration
#'
#' @param obj
#' iPLS results (object of class ipls)
#' @param glob.ncomp
#' number of components for global PLS model with all intervals
#' @param main
#' main title for the plot
#' @param xlab
#' label for x-axis
#' @param ylab
#' label for y-axis
#' @param xlim
#' limits for x-axis
#' @param ylim
#' limits for y-axis
#' @param ...
#' other arguments
#'
#' @details
#' The plot shows intervals as bars, which height corresponds to RMSECV obtained when particular
#' interval was selected (forward) or excluded (backward) from a model at the first iteration.
#' The intervals found optimal after backward/forward iPLS selection are shown with green color
#' while the other intervals are gray.
#'
#' See examples in help for \code{\link{ipls}} function.
#'
#'  @seealso
#'  \code{\link{summary.ipls}}, \code{\link{plotRMSE.ipls}}
#'
plotSelection.ipls <- function(obj, glob.ncomp = obj$gm$ncomp.selected, main = "iPLS results",
   xlab = obj$xaxis.name, ylab = if (is.null(obj$cv)) "RMSEP" else "RMSECV", xlim = NULL, ylim = NULL, ...) {

   if (glob.ncomp < 1 || glob.ncomp > obj$gm$ncomp) {
      stop("Wrong value for number of components.")
   }

   int <- obj$int.limits
   xlabels <- if (is.null(obj$xaxis.values)) seq_len(max(int)) else obj$xaxis.values
   if (!is.numeric(xlabels)) {
      stop("Parameter 'xlabels' should be a vector with numbers.")
   }

   if (length(xlabels) != (max(int) - min(int) + 1)) {
      stop("Wrong values for 'xlabels' parameter.")
   }

   # redefine intervals using xlabels
   int[, 1] <- xlabels[int[, 1]]
   int[, 2] <- xlabels[int[, 2]]

   # compute values for barplot
   mids <- apply(int, 1, mean)
   rmse <- obj$glob.stat$RMSE[2:nrow(obj$glob.stat)]
   ncomp <- obj$glob.stat$nComp[2:nrow(obj$glob.stat)]

   # adjust axis limits
   if (is.null(xlim)) xlim <- c(min(int), max(int))
   if (is.null(ylim)) ylim <- c(0, max(rmse) * 1.1)

   # rescale mean X values to fit the plot
   xmean <- (obj$xmean - min(obj$xmean)) / (max(obj$xmean) - min(obj$xmean)) * (ylim[2] - ylim[1])

   # make plot
   plot(0, 0, type = "n", main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)

   # show intervals as gray bars
   dim(rmse) <- c(1, length(rmse))
   plotBars(
      list(x_values = mids, y_values = rmse),
      col = rgb(0.9, 0.9, 0.9),
      bwd = 1,
      border = rgb(0.8, 0.8, 0.8))

   # show selected intervals as green bars
   rmse[, -obj$int.selected] <- 0
   plotBars(
      list(x_values = mids, y_values = rmse),
      col = rgb(0.5, 1.0, 0.6),
      bwd = 1,
      border = rgb(0.75, 0.8, 0.75)
   )

   # show mean signal
   lines(xlabels, xmean, col = rgb(1.0, 0.7, 0.7), lwd = 2)

   # number of components for each interval
   text(mids,  matrix(0.05 * ylim[2], ncol = length(mids)), ncomp,
      col = rgb(0.4, 0.4, 0.4), cex = 0.85)

   dx <- diff(xlim) / 50
   abline(h = obj$gm$cvres$rmse[1, glob.ncomp], lty = 2, col = rgb(0.5, 0.5, 0.5))
   text(xlim[2] + dx, obj$gm$cvres$rmse[1, glob.ncomp], glob.ncomp, cex = 0.85,
        col = rgb(0.3, 0.3, 0.3), font = 2, pos = 3)
}

#' RMSE development plot
#'
#' @description
#' Shows how RMSE develops for each iteration of iPLS selection
#' algorithm
#'
#' @param obj
#' iPLS results (object of class ipls)
#' @param glob.ncomp
#' number of components for global PLS model with all intervals
#' @param main
#' main title for the plot
#' @param xlab
#' label for x-axis
#' @param ylab
#' label for y-axis
#' @param xlim
#' limits for x-axis
#' @param ylim
#' limits for y-axis
#' @param ...
#' other arguments
#'
#' @details
#' The plot shows RMSE values obtained at each iteration of the iPLS algorithm as bars. The first
#' bar correspond to the global model with all variables included, second - to the model obtained
#' at the first iteration and so on. Number at the bottom of each bar corresponds to the interval
#' included or excluded at the particular iteration.
#'
#' @seealso
#' \code{\link{summary.ipls}}, \code{\link{plotSelection.ipls}}
#'
#' @export
plotRMSE.ipls <- function(obj, glob.ncomp = obj$gm$ncomp.selected, main = "RMSE development",
   xlab = "Iterations", ylab = if (is.null(obj$cv)) "RMSEP" else "RMSECV", xlim = NULL, ylim = NULL, ...) {

   if (glob.ncomp < 1 || glob.ncomp > obj$gm$ncomp) {
      stop("Wrong value for number of components.")
   }

   rmse <- obj$int.stat$RMSE
   n <- obj$int.stat$n
   mids <- 0:(length(n) - 1)

   if (is.null(xlim)) xlim <- c(min(mids) - 0.5, max(mids) + 0.5)
   if (is.null(ylim)) ylim <- c(0, max(rmse) * 1.1)

   # make plot
   plot(0, 0, type = "n", main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
   abline(h = rmse[1], lty = 2, col = "gray")
   lines(mids, rmse, type = "b", col = mdaplot.getColors(1), pch = 16)
}

#' Overview plot for iPLS results
#'
#' @description
#' Shows a plot for iPLS results.
#'
#' @param x
#' a  (object of class \code{pca})
#' @param ...
#' other arguments
#'
#' @details
#' See details for \code{\link{plotSelection.ipls}}.
#'
#' @export
plot.ipls <- function(x, ...) {
   plotSelection(x, ...)
}

#' Print method for iPLS
#'
#' @description
#' Prints information about the iPLS object structure
#'
#' @param x
#' a iPLS (object of class \code{ipls})
#' @param ...
#' other arguments
#'
#' @export
print.ipls <- function(x, ...) {
   cat("\niPLS results (class ipls)\n")
   cat("\nCall:\n")
   print(x$call)
   cat("\nMajor fields:\n")
   cat("$var.selected - vector with selected variables\n")
   cat("$int.selected - vector with selected intervals\n")
   cat("$int.num - number of intervals\n")
   cat("$int.width - width of the intervals\n")
   cat("$int.limits - limits for the intervals\n")
   cat("$int.stat - table with statistics for the interval selection results\n")
   cat("$glob.stat - table with statistics for the first iteration of the algorithm\n")
   cat("\nTry summary(obj) and plot(obj) to see details.\n")
}

#' Summary for iPLS results
#'
#' @description
#' Shows statistics and algorithm parameters for iPLS results.
#'
#' @param object
#' a iPLS (object of class \code{ipls})
#' @param glob.ncomp
#' number of components for global PLS model with all intervals
#' @param ...
#' other arguments
#'
#' @details
#' The method shows information on the algorithm parameters as well as a table with selected or
#' excluded interval. The table has the following columns: 'step' showing on which iteration
#' an interval was selected or excluded, 'start and 'end' show variable indices for the interval,
#' 'nComp' is a number of components used in a model, 'RMSE' is RMSECV for the model and 'R2' is
#' coefficient of determination for the same model.
#'
#' @export
summary.ipls <- function(object, glob.ncomp = object$gm$ncomp.selected, ...) {

   if (glob.ncomp < 1 || glob.ncomp > object$gm$ncomp) {
      stop("Wrong value for number of components.")
   }

   glob.rmse <- object$gm$cvres$rmse[1, glob.ncomp]
   opt.ncomp <- object$om$ncomp.selected
   opt.rmse <- object$om$cvres$rmse[1, opt.ncomp]
   rmse.suffix <- if (is.null(object$cv)) "P" else "CV"

   cat("\niPLS variable selection results\n")
   fprintf("  Method: %s\n", object$method)
   fprintf("  Validation: %s\n", if (is.null(object$cv)) "test set" else crossval.str(object$cv))
   fprintf("  Number of intervals: %d\n", object$int.num)
   fprintf("  Number of selected intervals: %d\n", length(object$int.selected))
   fprintf("  RMSE%s for global model: %f (%d LVs)\n", rmse.suffix, glob.rmse, glob.ncomp)
   fprintf("  RMSE%s for optimized model: %f (%d LVs)\n", rmse.suffix, opt.rmse, opt.ncomp)
   cat("\nSummary for selection procedure:\n")
   show(object$int.stat)
   cat("\n")
}
