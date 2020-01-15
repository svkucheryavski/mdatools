#' Partial Least Squares regression
#'
#' @description
#' \code{pls} is used to calibrate, validate and use of partial least squares (PLS)
#' regression model.
#'
#' @param x
#' matrix with predictors.
#' @param y
#' matrix with responses.
#' @param ncomp
#' maximum number of components to calculate.
#' @param center
#' logical, center or not predictors and response values.
#' @param scale
#' logical, scale (standardize) or not predictors and response values.
#' @param cv
#' cross-validation settings (see details).
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param x.test
#' matrix with predictors for test set.
#' @param y.test
#' matrix with responses for test set.
#' @param method
#' algorithm for computing PLS model (only 'simpls' is supported so far)
#' @param alpha
#' significance level for calculating statistical limits for residuals.
#' @param coeffs.ci
#' method to calculate p-values and confidence intervals for regression coefficients (so far only
#' jack-knifing is availavle: \code{='jk'}).
#' @param coeffs.alpha
#' significance level for calculating confidence intervals for regression coefficients.
#' @param info
#' short text with information about the model.
#' @param light
#' run normal or light (faster) version of PLS without calculationg some performance statistics.
#' @param ncomp.selcrit
#' criterion for selecting optimal number of components (\code{'min'} for
#' first local minimum of RMSECV and \code{'wold'} for Wold's rule.)
#'
#' @return
#' Returns an object of \code{pls} class with following fields:
#' \item{ncomp }{number of components included to the model.}
#' \item{ncomp.selected }{selected (optimal) number of components.}
#' \item{xloadings }{matrix with loading values for x decomposition.}
#' \item{yloadings }{matrix with loading values for y decomposition.}
#' \item{weights }{matrix with PLS weights.}
#' \item{selratio }{array with selectivity ratio values.}
#' \item{vipscores }{matrix with VIP scores values.}
#' \item{coeffs }{object of class \code{\link{regcoeffs}} with regression coefficients calculated
#' for each component.}
#' \item{info }{information about the model, provided by user when build the model.}
#' \item{calres }{an object of class \code{\link{plsres}} with PLS results for a calibration data.}
#' \item{testres }{an object of class \code{\link{plsres}} with PLS results for a test data, if it
#' was provided.}
#' \item{cvres }{an object of class \code{\link{plsres}} with PLS results for cross-validation, if
#' this option was chosen.}
#'
#' @details
#' So far only SIMPLS method [1] is available, more coming soon. Implementation works both with one
#' and multiple response variables.
#'
#' Like in \code{\link{pca}}, \code{pls} uses number of components (\code{ncomp}) as a minimum of
#' number of objects - 1, number of x variables and the default or provided value. Regression
#' coefficients, predictions and other results are calculated for each set of components from 1
#' to \code{ncomp}: 1, 1:2, 1:3, etc. The optimal number of components, (\code{ncomp.selected}),
#' is found using Wold's R criterion, but can be adjusted by user using function
#' (\code{\link{selectCompNum.pls}}). The selected optimal number of components is used for all
#' default operations - predictions, plots, etc.
#'
#' Cross-validation settings, \code{cv}, can be a number or a list. If \code{cv} is a number, it
#' will be used as a number of segments for random cross-validation (if \code{cv = 1}, full
#' cross-validation will be preformed). If it is a list, the following syntax can be used:
#' \code{cv = list('rand', nseg, nrep)} for random repeated cross-validation with \code{nseg}
#' segments and \code{nrep} repetitions or \code{cv = list('ven', nseg)} for systematic splits
#' to \code{nseg} segments ('venetian blinds').
#'
#' Selectivity ratio [2] and VIP scores [3] are calculated for any PLS model authomatically, however
#' while selectivity ratio values are calculated for all computed components, the VIP scores are
#' computed only for selected components (to save calculation time) and recalculated every time when
#' \code{selectCompNum()} is called for the model.
#'
#' Calculation of confidence intervals and p-values for regression coefficients are available
#' only by jack-knifing so far. See help for \code{\link{regcoeffs}} objects for details.
#'
#' @references
#' 1. S. de Jong, Chemometrics and Intelligent Laboratory Systems 18 (1993) 251-263.
#' 2. Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), 35-48.
#' 3. Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), 103-112.
#'
#' @seealso
#' Methods for \code{pls} objects:
#' \tabular{ll}{
#'  \code{print} \tab prints information about a \code{pls} object.\cr
#'  \code{\link{summary.pls}} \tab shows performance statistics for the model.\cr
#'  \code{\link{plot.pls}} \tab shows plot overview of the model.\cr
#'  \code{\link{pls.simpls}} \tab implementation of SIMPLS algorithm.\cr
#'  \code{\link{predict.pls}} \tab applies PLS model to a new data.\cr
#'  \code{\link{selectCompNum.pls}} \tab set number of optimal components in the model.\cr
#'  \code{\link{plotPredictions.pls}} \tab shows predicted vs. measured plot.\cr
#'  \code{\link{plotRegcoeffs.pls}} \tab shows regression coefficients plot.\cr
#'  \code{\link{plotXScores.pls}} \tab shows scores plot for x decomposition.\cr
#'  \code{\link{plotXYScores.pls}} \tab shows scores plot for x and y decomposition.\cr
#'  \code{\link{plotXLoadings.pls}} \tab shows loadings plot for x decomposition.\cr
#'  \code{\link{plotXYLoadings.pls}} \tab shows loadings plot for x and y decomposition.\cr
#'  \code{\link{plotRMSE.pls}} \tab shows RMSE plot.\cr
#'  \code{\link{plotXVariance.pls}} \tab shows explained variance plot for x decomposition.\cr
#'  \code{\link{plotYVariance.pls}} \tab shows explained variance plot for y decomposition.\cr
#'  \code{\link{plotXCumVariance.pls}} \tab shows cumulative explained variance plot for y
#'  decomposition.\cr
#'  \code{\link{plotYCumVariance.pls}} \tab shows cumulative explained variance plot for y
#'  decomposition.\cr
#'  \code{\link{plotXResiduals.pls}} \tab shows T2 vs. Q plot for x decomposition.\cr
#'  \code{\link{plotYResiduals.pls}} \tab shows residuals plot for y values.\cr
#'  \code{\link{plotSelectivityRatio.pls}} \tab shows plot with selectivity ratio values.\cr
#'  \code{\link{plotVIPScores.pls}} \tab shows plot with VIP scores values.\cr
#'  \code{\link{getSelectivityRatio.pls}} \tab returns vector with selectivity ratio values.\cr
#'  \code{\link{getVIPScores.pls}} \tab returns vector with VIP scores values.\cr
#'  \code{\link{getRegcoeffs.pls}} \tab returns matrix with regression coefficients.\cr
#' }
#'
#' Most of the methods for plotting data (except loadings and regression coefficients) are also
#' available for PLS results
#' (\code{\link{plsres}}) objects. There is also a randomization test for PLS-regression
#' (\code{\link{randtest}}).
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @examples
#' ### Examples of using PLS model class
#' library(mdatools)
#'
#' ## 1. Make a PLS model for concentration of first component
#' ## using full-cross validation and automatic detection of
#' ## optimal number of components and show an overview
#'
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#'
#' model = pls(x, y, ncomp = 8, cv = 1)
#' summary(model)
#' plot(model)
#'
#' ## 2. Make a PLS model for concentration of first component
#' ## using test set and 10 segment cross-validation and show overview
#'
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#' x.t = simdata$spectra.t
#' y.t = simdata$conc.t[, 1]
#'
#' model = pls(x, y, ncomp = 8, cv = 10, x.test = x.t, y.test = y.t)
#' model = selectCompNum(model, 2)
#' summary(model)
#' plot(model)
#'
#' ## 3. Make a PLS model for concentration of first component
#' ## using only test set validation and show overview
#'
#' data(simdata)
#' x = simdata$spectra.c
#' y = simdata$conc.c[, 1]
#' x.t = simdata$spectra.t
#' y.t = simdata$conc.t[, 1]
#'
#' model = pls(x, y, ncomp = 6, x.test = x.t, y.test = y.t)
#' model = selectCompNum(model, 2)
#' summary(model)
#' plot(model)
#'
#' ## 4. Show variance and error plots for a PLS model
#' par(mfrow = c(2, 2))
#' plotXCumVariance(model, type = 'h')
#' plotYCumVariance(model, type = 'b', show.labels = TRUE, legend.position = 'bottomright')
#' plotRMSE(model)
#' plotRMSE(model, type = 'h', show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' ## 5. Show scores plots for a PLS model
#' par(mfrow = c(2, 2))
#' plotXScores(model)
#' plotXScores(model, comp = c(1, 3), show.labels = TRUE)
#' plotXYScores(model)
#' plotXYScores(model, comp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' ## 6. Show loadings and coefficients plots for a PLS model
#' par(mfrow = c(2, 2))
#' plotXLoadings(model)
#' plotXLoadings(model, comp = c(1, 2), type = 'l')
#' plotXYLoadings(model, comp = c(1, 2), legend.position = 'topleft')
#' plotRegcoeffs(model)
#' par(mfrow = c(1, 1))
#'
#' ## 7. Show predictions and residuals plots for a PLS model
#' par(mfrow = c(2, 2))
#' plotXResiduals(model, show.label = TRUE)
#' plotYResiduals(model, show.label = TRUE)
#' plotPredictions(model)
#' plotPredictions(model, ncomp = 4, xlab = 'C, reference', ylab = 'C, predictions')
#' par(mfrow = c(1, 1))
#'
#' ## 8. Selectivity ratio and VIP scores plots
#' par(mfrow = c(2, 2))
#' plotSelectivityRatio(model)
#' plotSelectivityRatio(model, ncomp = 1)
#' par(mfrow = c(1, 1))
#'
#' ## 9. Variable selection with selectivity ratio
#' selratio = getSelectivityRatio(model)
#' selvar = !(selratio < 8)
#'
#' xsel = x[, selvar]
#' modelsel = pls(xsel, y, ncomp = 6, cv = 1)
#' modelsel = selectCompNum(modelsel, 3)
#'
#' summary(model)
#' summary(modelsel)
#'
#' ## 10. Calculate average spectrum and show the selected variables
#' i = 1:ncol(x)
#' ms = apply(x, 2, mean)
#'
#' par(mfrow = c(2, 2))
#'
#' plot(i, ms, type = 'p', pch = 16, col = 'red', main = 'Original variables')
#' plotPredictions(model)
#'
#' plot(i, ms, type = 'p', pch = 16, col = 'lightgray', main = 'Selected variables')
#' points(i[selvar], ms[selvar], col = 'red', pch = 16)
#' plotPredictions(modelsel)
#'
#' par(mfrow = c(1, 1))
#'
#' @export
pls <- function(x, y, ncomp = min(nrow(x) - 1, ncol(x), 20), center = TRUE, scale = FALSE,
   cv = NULL, exclcols = NULL, exclrows = NULL, x.test = NULL, y.test = NULL, method = "simpls",
   alpha = 0.05, coeffs.ci = "jk", coeffs.alpha = 0.05, info = "", light = (ncol(x) > 300),
   ncomp.selcrit = "min") {

   # exclude columns if "exclcols" is provided
   if (length(exclcols) > 0) {
      x <- mda.exclcols(x, exclcols)
   }

   # exclude rows if "exclrows" is provided
   if (length(exclrows) > 0) {
      x <- mda.exclrows(x, exclrows)
      y <- mda.exclrows(y, exclrows)
   }

   # build a model and apply to calibration set
   model <- pls.cal(x, y, ncomp, center = center, scale = scale, method = method)
   model$info <- info
   model$call <- match.call()

   # get calibration results
   model$res <- list()
   model$res[["cal"]] <- predict.pls(model, x, y)
   model$res[["cal"]]$info <- "Calibration results"
   model$calres <- model$res[["cal"]]

   # do cross-validation if needed
   if (!is.null(cv)) {
      cvres <- crossval.regmodel(model, x, y, cv, cal.fun = pls.cal)
      model$res[["cv"]] <- plsres(cvres$y.pred, cvres$y.ref, ncomp.selected = ncomp)
      model$res[["cv"]]$info <- "Cross-validation results"
      model$cvres <- model$res[["cv"]]

      if (coeffs.ci == "jk") {
         model$coeffs <- regcoeffs(model$coeffs$values, cvres$jk.coeffs)
      }
   }

   # do test set validation if provided
   if (!is.null(x.test) && !is.null(y.test)){
      model$res[["test"]] <- predict.pls(model, x.test, y.test)
      model$res[["test"]]$info <- "Test set validation results"
      model$testres <- model$res[["test"]]
   }

   if (!light) {
      # we calculate both for data without excluded rows and columns
      # model$selratio <- pls.calculateSelectivityRatio(model, x)
      # model$vipscores <- pls.calculateVIPScores(model)
   }

   # select optimal number of components
   model$ncomp.selcrit <- ncomp.selcrit
   model <- selectCompNum(model, selcrit = ncomp.selcrit)
   return(model)
}

#' Select optimal number of components for PLS model
#'
#' @description
#' Allows user to select optimal number of components for PLS model
#'
#' @param model
#' PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to select
#' @param selcrit
#' criterion for selecting optimal number of components (\code{'min'} for
#' first local minimum of RMSECV and \code{'wold'} for Wold's rule.)
#' @param ...
#' other parameters if any
#'
#' @return
#' the same model with selected number of components
#'
#' @details
#'
#' The method sets \code{ncomp.selected} parameter for the model and return it back. The parameter
#' points out to the optimal number of components in the model. You can either specify it manually
#' as argument \code{ncomp} or use one of the algorithms for automatic selection.
#'
#' Automatic selection by default based on cross-validation statistics. If no cross-validation
#' results are found in the model, the method will use test set validation results. If they are
#' not available as well, the model will use calibration results and give a warning as in this case
#' the selected number of components will lead to overfitted model.
#'
#' There are two algorithms for automatic selection you can chose between: either first local
#' minimum of RMSE (`selcrit='min'`) or Wold's rule (`selcrit='wold'`).
#'
#' The first local minimum criterion finds at which component, A, error of prediction starts
#' raising and select (A - 1) as the optimal number. The Wold's criterion finds which component A
#' does not make error smaller at least by 5% comparing to previous value and selects (A - 1) as
#' the optimal number.
#'
#' If model is PLS2 model (has several response variables) the method computes optimal number of
#' components for each response and returns the smallest value. For example, if for first
#' response 2 components give the smallest error and for the second response this number is 3,
#' A = 2 will be selected as a final result.
#'
#' It is not recommended to use automatic selection for real applications, always investigate
#' your model (RMSE and Y-variance plot, regression coefficients) to make correct decision.
#'
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
selectCompNum.pls <- function(obj, ncomp = NULL, selcrit = obj$ncomp.selcrit, ...) {

   # returns index based on Wold's R criterion
   # finds first number of components which does not make error smaller by 5%
   # comparing to the previous value
   fWold <- function(press, ncomp) {
      r <- press[, 2:ncomp, drop = FALSE] / press[, 1:(ncomp - 1), drop = FALSE]
      return(which(r > 0.95, arr.ind = TRUE))
   }

   # returns index based on first local minimum
   fMin <- function(press, ncomp) {
      r <- press[, 2:ncomp, drop = FALSE] - press[, 1:(ncomp - 1), drop = FALSE]
      return(which(r > 0, arr.ind = TRUE))
   }

   # returns number of components based on PRESS and index values
   f <- function(res, indFun) {
      press <- res$rmse^2 * nrow(res$y.pred)
      ncomp <- ncol(press)

      # if number of components is 2 - for every row find which component
      # gives smallest error and take the smallest number
      if (ncomp < 3) return(min(apply(press, 1, which.min)))

      # otherwise use dedicated function
      ind <- indFun(press, ncomp)

      return((if (length(ind) == 0) ncomp else min(ind[, 2])))
   }

   # if only one component in the model do nothing
   if (obj$ncomp == 1) return(obj)

   # if user provided ncomp use it
   if (!is.null(ncomp)) selcrit <- ""

   # get name of first result available in the sequence
   name <- intersect(c("cv", "test", "cal"), names(obj$res))[1]
   res <- obj$res[[name]]

   if (name == "cal") {
      warning("No validation results were found.")
   }

   # estimate number of optimal components
   if (is.null(selcrit)) selcrit = ""
   ncomp <- switch(selcrit,
      "wold" = f(res, fWold),
      "min" = f(res, fMin),
      ncomp
   )

   # if NULL - somthing went wrong
   if (is.null(ncomp)) {
      stop("Can not estimate correct number of PLS components.")
   }

   # if not, check that value is meaningful
   if (ncomp > obj$ncomp || ncomp < 0) {
      stop("Wrong number of selected components.")
   }

   # correct number of model and calibration results
   obj$ncomp.selected <- ncomp
   obj$res[["cal"]]$ncomp.selected <- ncomp
   obj$calres <- obj$res[["cal"]]

   # correct number of components for cross-validation results
   if (!is.null(obj$res[["cv"]])) {
      obj$res[["cv"]]$ncomp.selected <- ncomp
      obj$cvres <- obj$res[["cv"]]
   }

   # correct number of components for test set results
   if (!is.null(obj$res[["test"]])) {
      obj$res[["test"]]$ncomp.selected <- ncomp
      obj$testres <- obj$res[["test"]]
   }

   obj$call <- match.call()
   return(obj)
}

#' PLS predictions
#'
#' @description
#' Applies PLS model to a new data set
#'
#' @param object
#' a PLS model (object of class \code{pls})
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with reference y values (responses)
#' @param cv
#' logical, shall predictions be made for cross-validation procedure or not
#' @param ...
#' other arguments
#'
#' @return
#' PLS results (an object of class \code{plsres})
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
predict.pls <- function(object, x, y = NULL, cv = FALSE, ...) {

   # get names
   prednames <- rownames(object$xloadings)
   respnames <- rownames(object$yloadings)
   objnames <- rownames(x)

   # preprocess x and calculate scores, total and full variance
   x.attrs <- mda.getattr(x)
   y.attrs <- mda.getattr(y)

   # set names for y-axis (rows if it is empty)
   if (is.null(x.attrs$yaxis.name)) {
      x.attrs$yaxis.name <- "Objects"
   }

   # correct dimension for y if it is a vector
   if (!is.null(y) && is.null(dim(y))) {
      dim(y) <- c(length(y), 1)
   }

   # convert to matrices
   x <- mda.df2mat(x)
   y <- mda.df2mat(y)

   # get dimensions
   nresp <- dim(object$coeffs$values)[3]
   ncomp <- dim(object$coeffs$values)[2]

   # check dimensions of predictors
   if (ncol(x) != dim(object$coeffs$values)[1]) {
      stop("Wrong number of columns in matrix with predictors (x).")
   }

   # check dimensions of responses
   if (!is.null(y)) {
      if (nrow(x) != nrow(y))
         stop("Matrices with predictors (x) and response (y) should have the same number of rows.")
      if (ncol(y) != nresp)
         stop("Wrong number of columns in matrix with response values (y).")
   }

   # autoscale x
   x <- prep.autoscale(x, center = object$xcenter, scale = object$xscale)

   # compute x scores and residuals
   xscores <- x %*% (object$weights %*% solve(crossprod(object$xloadings, object$weights)))
   xresiduals <- x - tcrossprod(xscores, object$xloadings)

   # set attributes
   xscores <- mda.setattr(xscores, x.attrs, "row")
   xresiduals <- mda.setattr(xresiduals, x.attrs)
   attr(xscores, "name") <- "X-scores"
   attr(xscores, "xaxis.name") <- "Components"
   attr(xresiduals, "name") <- "Residuals"

   # set names
   rownames(xscores) <- rownames(xresiduals) <- objnames
   colnames(xscores) <- colnames(object$xloadings)
   colnames(xresiduals) <- prednames

   # make predictions
   yp <- apply(object$coeffs$values, 3, function(x, y)(y %*% x), x)
   dim(yp) <- c(nrow(x), ncomp, dim(object$coeffs$values)[3])

   # if reference values are provided calculate and set up names for ydecomp
   y.ref <- y
   ydecomp <- NULL
   if (!is.null(y.ref)) {

      # autoscale y-values
      y <- prep.autoscale(y, center = object$ycenter, scale = object$yscale)

      # compute and orthogonalize y-scores
      yscores <- as.matrix(y) %*% object$yloadings
      for (a in seq_len(ncomp)) {
         for (n in 1:2) {
            for (j in seq_len(a - 1)) {
               yscores[, a] <- yscores[, a] -
                  tcrossprod(xscores[, j], yscores[, a]) %*% xscores[, j]
            }
         }
      }

      # compute y-residuals
      yresiduals <- y - yp[, ncomp, ]

      # set names
      rownames(yscores) <- rownames(yresiduals) <- objnames
      colnames(yscores) <- colnames(object$yloadings)
      colnames(yresiduals) <- respnames

      # set attributes
      yscores <- mda.setattr(yscores, x.attrs, "row")
      yresiduals <- mda.setattr(yresiduals, y.attrs)
      attr(yscores, "exclrows") <- attr(yresiduals, "exclrows") <- x.attrs$exclrows
      attr(yscores, "name") <- "Y-scores"
      attr(yscores, "xaxis.name") <- "Components"
      attr(yresiduals, "name") <- "Residuals"

      # create ydecomp object (we use xscores as residuals for different components are computed
      # as xscores %*% t(yloadings)), but then we assign correct residuals
      ydecomp <- ldecomp(scores = xscores, loadings = object$yloadings, residuals = yresiduals,
            eigenvals = object$yeigenvals, ncomp.selected = object$ncomp.selected)
      ydecomp$scores <- yscores
   }

   # unscale predicted y values
   if (is.numeric(object$yscale)) {
      yp <- sweep(yp, 3, object$yscale, '*')
   }

   # uncenter predicted y values
   if (is.numeric(object$ycenter)) {
      yp <- sweep(yp, 3, object$ycenter, '+')
   }

   # if predictions for cross-validation - return
   if (cv) {
      return(list(y.pred = yp))
   }

   # set up all attributes and names
   yp <- mda.setattr(yp, x.attrs, "row")
   attr(yp, "exclrows") <- x.attrs$exclrows
   attr(yp, "name") <- "Response values, predicted"
   dimnames(yp) <- c(list(rownames(x)), dimnames(object$coeffs$values)[2:3])

   # create xdecomp object
   xdecomp <- ldecomp(scores = xscores, residuals = xresiduals, loadings = object$xloadings,
      eigenvals = object$xeigenvals, ncomp.selected = object$ncomp.selected)

   return(
      plsres(yp, y.ref = y.ref, ncomp.selected = object$ncomp.selected,
         xdecomp = xdecomp, ydecomp = ydecomp)
   )
}


################################
#  Plotting methods            #
################################


#' Explained X variance plot for PLS
#'
#' @description
#' Shows plot with explained X variance vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot("b", "l" or "h")
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXVariance.pls <- function(obj, type = "b", main = "Variance (X)", ...) {
   plotVariance(obj, decomp = "xdecomp", type = type, main = main, ...)
}

#' Explained Y variance plot for PLS
#'
#' @description
#' Shows plot with explained Y variance vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot("b", "l" or "h")
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotYVariance.pls <- function(obj, type = "b", main = "Variance (Y)", ...) {
   plotVariance(obj, decomp = "ydecomp", type = type, main = main, ...)
}

#' Cumulative explained X variance plot for PLS
#'
#' @description
#' Shows plot with cumulative explained X variance vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot("b", "l" or "h")
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXCumVariance.pls <- function(obj, type = "b", main = "Cumulative variance (X)", ...) {
   plotVariance(obj, decomp = "xdecomp", variance = "cumexpvar", type = type, main = main, ...)
}

#' Cumulative explained Y variance plot for PLS
#'
#' @description
#' Shows plot with cumulative explained Y variance vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param type
#' type of the plot("b", "l" or "h")
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotYCumVariance.pls <- function(obj, type = "b", main = "Cumulative variance (Y)", ...) {
   plotVariance(obj, decomp = "ydecomp", variance = "cumexpvar", type = type, main = main, ...)
}

#' Variance plot for PLS
#'
#' @description
#' Shows plot with variance values vs. number of components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param decomp
#' which decomposition to use ('xdecomp' for x or 'ydecomp' for y)
#' @param variance
#' which variance to use ('expvar', 'cumexpvar)
#' @param type
#' type of the plot('b', 'l' or 'h')
#' @param labels
#' what to show as labels for plot objects.
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotVariance.pls <- function(obj, decomp = "xdecomp", variance = "expvar", type = "b",
   labels = "values", res = obj$res, ...) {

   plot_data <- lapply(res, plotVariance, decomp = decomp, variance = variance, show.plot = FALSE)
   mdaplotg(plot_data, labels = labels, type = type, ...)
}

#' X scores plot for PLS
#'
#' @description
#' Shows plot with X scores values for selected components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param main
#' main plot title
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXScores.pls <- function(obj, comp = c(1, 2), show.axes = T,  main = "Scores (X)",
   res = obj$res, ...) {


   # set up values for showing axes lines
   show.lines <- FALSE
   if (show.axes) {
      show.lines <- if (length(comp) == 2) c(0, 0) else c(NA, 0)
   }

   plot_data <- lapply(res, plotXScores, comp = comp, type = "p", show.plot = FALSE)
   mdaplotg(plot_data, show.lines = show.lines, type = "p", main = main, ...)
}

#' XY scores plot for PLS
#'
#' @description
#' Shows plot with X vs. Y scores values for selected component.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which component to show the plot for
#' @param main
#' main plot title
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXYScores.pls <- function(obj, ncomp = 1, show.axes = T,  res = obj$res, ...) {

   show.lines <- if (show.axes) c(0, 0) else FALSE
   plot_data <- lapply(res, plotXYScores, ncomp = ncomp, type = "p", show.plot = FALSE)
   mdaplotg(plot_data, show.lines = show.lines, type = "p", ...)
}

#' X residuals plot for PLS
#'
#' @description
#' Shows a plot with Q residuals vs. Hotelling T2 values for PLS decomposition of x data.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXResiduals.pls <- function(obj, ncomp = 1, res = obj$res,
   main = sprintf("X-residuals (ncomp = %d)", ncomp), ...) {

   # TODO: implement norm and log attributes support (u0 attribute for Q and T2)
   # TODO: implement showing residual limits
   plot_data <- lapply(res, plotXResiduals, ncomp = ncomp, show.plot = FALSE)
   mdaplotg(plot_data, type = "p", main = main, ...)
}

#' X loadings plot for PLS
#'
#' @description
#' Shows plot with X loading values for selected components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXLoadings.pls <- function(obj, comp = c(1, 2), type = "p", show.axes = TRUE,
   show.legend = TRUE, ...) {

   plot_data <- mda.subset(obj$xloadings, select = comp)
   colnames(plot_data) <- sprintf("Comp %d (%.2f%%)", comp, obj$res[["cal"]]$xdecomp$expvar[comp])
   attr(plot_data, "name") <- "Loadings (X)"

   # set up values for showing axes lines
   show.lines <- FALSE
   if (show.axes) {
      show.lines <- if (length(comp) == 2 && type == "p") c(0, 0) else c(NA, 0)
   }

   if (type == "p") {
      return(mdaplot(plot_data, type = type, show.lines = show.lines, ...))
   }

   plot_data <- mda.t(plot_data)
   attr(plot_data, "yaxis.name") <- "Loading"
   mdaplotg(plot_data, show.legend = show.legend, type = type, show.lines = show.lines, ...)
}

#' X loadings plot for PLS
#'
#' @description
#' Shows plot with X loading values for selected components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param type
#' type of the plot
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotWeights.pls <- function(obj, comp = 1, type = (if (nrow(obj$weights) < 20) "h" else "l"),
   show.axes = TRUE, show.legend = TRUE, ...) {

   plot_data <- mda.subset(obj$weights, select = comp)
   colnames(plot_data) <- sprintf("Comp %d (%.2f%%)", comp, obj$res[["cal"]]$xdecomp$expvar[comp])
   attr(plot_data, "name") <- "Weights"

   # set up values for showing axes lines
   show.lines <- FALSE
   if (show.axes) {
      show.lines <- if (length(comp) == 2 && type == "p") c(0, 0) else c(NA, 0)
   }

   if (type == "p") {
      return(mdaplot(plot_data, type = type, show.lines = show.lines, ...))
   }

   plot_data <- mda.t(plot_data)
   attr(plot_data, "yaxis.name") <- "Weight"
   mdaplotg(plot_data, show.legend = show.legend, type = type, show.lines = show.lines, ...)
}

#' XY loadings plot for PLS
#'
#' @description
#' Shows plot with X and Y loading values for selected components.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param comp
#' which components to show the plot for (one or vector with several values)
#' @param main
#' main plot title
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plotXYLoadings.pls <- function(obj, comp = c(1, 2), show.axes = TRUE, ...) {

   if (length(comp) != 2) {
      stop("This plot can be made for only two components.")
   }

   plot_data <- list(
      "X" = mda.subset(obj$xloadings, select = comp),
      "Y" = mda.subset(obj$yloadings, select = comp)
   )
   colnames(plot_data[[1]]) <- sprintf("Comp %d (%.2f%%)", comp,
      obj$res[["cal"]]$xdecomp$expvar[comp])

   attr(plot_data, "name") <- "Loadings (XY)"
   show.lines <- if (show.axes) c(0, 0) else FALSE
   mdaplotg(plot_data, type = "p", show.lines = show.lines, ...)
}

#' Model overview plot for PLS
#'
#' @description
#' Shows a set of plots (x residuals, regression coefficients, RMSE and predictions) for PLS model.
#'
#' @param x
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' which y variable to show the summary for (if NULL, will be shown for all)
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{pls}} function.
#'
#' @export
plot.pls <- function(x, ncomp = x$ncomp.selected, ny = 1, show.legend = TRUE, ...) {

   if (!is.null(ncomp) && (ncomp <= 0 || ncomp > x$ncomp)) {
      stop("Wrong value for number of components.")
   }

   par(mfrow = c(2, 2))
   plotXResiduals(x, ncomp = ncomp, show.legend = show.legend)
   plotRegcoeffs(x, ncomp = ncomp, ny = ny)
   plotRMSE(x, ny = ny, show.legend = show.legend)
   plotPredictions(x, ncomp = ncomp, ny = ny, show.legend = show.legend)
   par(mfrow = c(1, 1))
}


################################
#  Static methods              #
################################


#' Runs selected PLS algorithm
#'
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param y
#' a matrix with y values (responses from calibration set)
#' @param ncomp
#' how many components to compute
#' @param method
#' algorithm for computing PLS model
#' @param cv
#' logical, is this for CV or not
#'
#' @export
pls.run <- function(x, y, ncomp = min(nrow(x) - 1, ncol(x)), method = "simpls", cv = FALSE) {

   if (ncomp < 1 || ncomp > min(nrow(x) - 1, ncol(x))) {
      stop("Wrong value for 'ncomp' parameter.")
   }

   methods <- list("simpls" = pls.simpls)

   if (!(method %in% names(methods))) {
      stop("Method with this name is not supported.")
   }

   return(methods[[method]](x, y, ncomp, cv = cv))
}

#' SIMPLS algorithm
#'
#' @description
#' SIMPLS algorithm for calibration of PLS model
#'
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with y values (responses)
#' @param ncomp
#' number of components to calculate
#' @param cv
#' logical, is model calibrated during cross-validation or not
#'
#' @return
#' a list with computed regression coefficients, loadings and scores for x and y matrices,
#' and weights.
#'
#' @references
#' [1]. S. de Jong. SIMPLS: An Alternative approach to partial least squares regression.
#' Chemometrics and Intelligent Laboratory Systems, 18, 1993 (251-263).
#'
pls.simpls <- function(x, y, ncomp, cv = FALSE) {

   x <- as.matrix(x)
   y <- as.matrix(y)

   npred <- ncol(x)
   nresp <- ncol(y)

   # initial estimation
   A <- crossprod(x, y)
   M <- crossprod(x, x)
   C <- diag(npred)

   # prepare space for results
   B <- array(0, dim = c(npred, ncomp, nresp))
   W <- matrix(0, nrow = npred, ncol = ncomp)
   P <- matrix(0, nrow = npred, ncol = ncomp)
   Q <- matrix(0, nrow = nresp, ncol = ncomp)

   # loop for each components
   for (n in seq_len(ncomp)) {
      # get the dominate eigenvector of A'A
      e <- eigen(crossprod(A))
      q <- e$vectors[seq_len(nresp)]

      # calculate and store weights
      w <- A %*% q
      c <- crossprod(w, (M %*% w))
      w <- w / sqrt(as.numeric(c))
      W[, n] <- w

      # calculate and store x loadings
      p <- M %*% w
      P[, n] <- p

      # calculate and store y loadings
      q <- crossprod(A, w)
      Q[, n] <- q

      v <- C %*% p
      v <- v / sqrt(as.numeric(crossprod(v)))

      # compute coefficients for current component
      B[, n, ] <- tcrossprod(W[, seq_len(n), drop = FALSE], Q[, seq_len(n), drop = FALSE])

      # recalculate matrices for the next compnonent
      C <- C - tcrossprod(v)
      M <- M - tcrossprod(p)
      A <- C %*% A

      if (!cv && max(e$values) < 10 * .Machine$double.eps) {
         # stop cycle is egienvalue is almost zero
         break
      }
   }

   # truncate results if n is smaller than ncomp
   W <- W[, seq_len(n), drop = FALSE]
   P <- P[, seq_len(n), drop = FALSE]
   Q <- Q[, seq_len(n), drop = FALSE]
   B <- B[, seq_len(n), , drop = FALSE]

   return(list(coeffs = B, weights = W, xloadings = P, yloadings = Q, ncomp = n))
}

#' PLS model calibration
#'
#' @description
#' Calibrates (builds) a PLS model for given data and parameters
#'
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with y values (responses)
#' @param ncomp
#' number of components to calculate
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#' @param method
#' algorithm for computing PLS model (only 'simpls' is supported so far)
#'
#' @return model
#' an object with calibrated PLS model
#'
pls.cal <- function(x, y, ncomp, center, scale, method = "simpls", cv = FALSE) {

   # prepare empty list for model object and assign
   # several properties, which do not depend on calculations below
   model <- list()
   model$center <- center
   model$scale <- scale
   model$method <- method
   class(model) <- c("pls", "regmodel")

   # get attributes
   x.attrs <- mda.getattr(x)
   y.attrs <- mda.getattr(y)

   # get names of variables
   prednames <- colnames(x)
   respnames <- colnames(y)

   # if y is a vector convert it to a matrix
   if (is.null(dim(y))) dim(y) <- c(length(y), 1)

   # check dimensions
   if (nrow(x) != nrow(y)) {
      stop("Number of rows for predictors and responses should be the same.")
   }

   # convert data to a matrix
   x <- mda.df2mat(x)
   y <- mda.df2mat(y)

   # get dimension of original data
   x.nrows <- nrow(x)
   x.ncols <- ncol(x)
   y.ncols <- ncol(y)

   # check if data has missing values
   if (any(is.na(x))) {
      stop("Predictors have missing values, try to fix this using pca.mvreplace.")
   }

   if (any(is.na(y))) {
      stop("Responses have missing values, try to fix this using pca.mvreplace.")
   }

   # set column names for predictors if missing
   if (is.null(colnames(y))) {
      colnames(y) <- paste0("y", seq_len(ncol(y)))
   }

   # correct x-axis name
   if (is.null(x.attrs$xaxis.name)) {
     x.attrs$xaxis.name <- "Predictors"
   }

   if (is.null(y.attrs$xaxis.name)) {
     y.attrs$xaxis.name <- "Responses"
   }

   # remove excluded rows
   if (length(x.attrs$exclrows) > 0) {
      x <- x[-x.attrs$exclrows, , drop = FALSE]
      y <- y[-x.attrs$exclrows, , drop = FALSE]
   }

   # autoscale and save the mean and std values for predictors
   x <- prep.autoscale(x, center = center, scale = scale)
   model$xcenter <- attr(x, "prep:center")
   model$xscale <- attr(x, "prep:scale")

   # autoscale and save the mean and std values for responses
   y <- prep.autoscale(y, center = center, scale = scale)
   model$ycenter <- attr(y, "prep:center")
   model$yscale <- attr(y, "prep:scale")

   # remove excluded columns for predictors
   if (length(x.attrs$exclcols) > 0) {
      x <- x[, -x.attrs$exclcols, drop = FALSE]
   }

   # remove excluded columns for responses
   if (length(y.attrs$exclcols) > 0) {
      y <- y[, -y.attrs$exclcols, drop = FALSE]
   }

   # get dimensions of reduced datasets
   xc.nrows <- nrow(x)
   xc.ncols <- ncol(x)
   yc.ncols <- ncol(y)

   # find maximum number of objects in a segment
   nobj.cv <- 0
   if (!is.logical(cv)) {
      nseg <- if (is.numeric(cv)) cv else cv[[2]]
      nobj.cv <- if (nseg == 1) 1 else ceiling(xc.nrows/nseg)
   }

   # correct maximum number of components
   ncomp <- min(xc.ncols, xc.nrows - 1 - nobj.cv, ncomp)

   # fit model

   fit <- pls.run(x, y, method = method, ncomp = ncomp)
   ncomp <- fit$ncomp
   model$ncomp <- ncomp

   # if it is for cross-validation return the results as is
   if (cv) {
      model$coeffs <- regcoeffs(fit$coeffs)
      model$xloadings <- fit$xloadings
      model$weights <- fit$weights
      return(model)
   }

   # compute eigenvalues
   xscores <- x %*% (fit$weights %*% solve(crossprod(fit$xloadings, fit$weights)))
   yscores <- as.matrix(y) %*% fit$yloadings
   xeigenvals <- colSums(xscores^2) / (xc.nrows - 1)
   yeigenvals <- colSums(yscores^2) / (xc.nrows - 1)

   # correct results related to predictors for missing columns in x
   # corresponding rows will be set to 0 and excluded
   xloadings <- matrix(0, nrow = x.ncols, ncol = ncomp)
   yloadings = matrix(0, nrow = y.ncols, ncol = ncomp)
   weights <- matrix(0, nrow = x.ncols, ncol = ncomp)
   coeffs <- array(0, dim = c(x.ncols, ncomp, yc.ncols))

   pred_ind <- seq_len(x.ncols)
   if (length(x.attrs$exclcols) > 0) pred_ind <- pred_ind[-x.attrs$exclcols]

   resp_ind <- seq_len(y.ncols)
   if (length(y.attrs$exclcols) > 0) resp_ind <- resp_ind[-y.attrs$exclcols]

   # x-loadings
   xloadings[pred_ind, ] <- fit$xloadings
   xloadings <- mda.exclrows(xloadings, x.attrs$exclcols)

   # y-loadings
   yloadings[resp_ind, ] <- fit$yloadings
   yloadings <- mda.exclrows(yloadings, y.attrs$exclcols)

   # weights
   weights[pred_ind, ] <- fit$weights
   weights <- mda.exclrows(weights, x.attrs$exclcols)

   # coeffs
   coeffs[pred_ind, , ] <- fit$coeffs
   coeffs <- mda.exclrows(coeffs, x.attrs$exclcols)

   # set names and attributes
   compnames <- paste("Comp", seq_len(ncomp))

   # x-loadings and weights
   rownames(xloadings) <- rownames(weights) <- prednames
   colnames(xloadings) <- colnames(weights) <- compnames
   attr(xloadings, "name") <- "X loadings"
   attr(weights, "name") <- "Weights"
   attr(xloadings, "xaxis.name") <- attr(weights, "xaxis.name") <- "Components"
   attr(xloadings, "yaxis.name") <- attr(weights, "yaxis.name") <- x.attrs$xaxis.name
   attr(xloadings, "yaxis.values") <- attr(weights, "yaxis.values") <- x.attrs$xaxis.values

   # coefficients
   dimnames(coeffs) <- list(prednames, compnames, colnames(y))
   attr(coeffs, "yaxis.name") <- x.attrs$xaxis.name
   attr(coeffs, "yaxis.values") <- x.attrs$xaxis.values

   # y-loadings
   rownames(yloadings) <- respnames
   colnames(yloadings) <- colnames(xloadings)
   attr(yloadings, "name") <- "Y loadings"
   attr(yloadings, "xaxis.name") <- "Components"
   attr(yloadings, "yaxis.name") <- y.attrs$xaxis.name
   attr(yloadings, "yaxis.values") <- y.attrs$xaxis.values

   # set up and return model parameters
   model$xloadings <- xloadings
   model$yloadings <- yloadings
   model$weights <- weights
   model$coeffs <- regcoeffs(coeffs)

   model$ncomp.selected <- ncomp
   model$exclrows <- x.attrs$exclrows
   model$exclcols <- x.attrs$exclcols
   model$xeigenvals <- xeigenvals
   model$yeigenvals <- yeigenvals

   return(model)
}

# ! Stopped refactoring here


#' Selectivity ratio calculation
#'
#' @description
#' Calculates selectivity ration for each component and response variable in
#' the PLS model
#'
#' @param model
#' a PLS model (object of class \code{pls})
#' @param x
#' predictor values from calibration set, preprocessed, centered and scaled
#'
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#'
#' @return
#' array \code{nvar x ncomp x ny} with selectivity ratio values
#'
pls.calculateSelectivityRatio = function(model, x) {
   ny = dim(model$coeffs$values)[3]
   ncomp = dim(model$coeffs$values)[2]
   nvar = ncol(x)

   # get coefficients
   coeffs = model$coeffs$values
   attrs = mda.getattr(coeffs)
   exclvars = attrs$exclrows
   nexclvar = length(exclvars)
   if (nexclvar > 0) {
      coeffs = coeffs[-exclvars, , , drop = FALSE]
      attrs$exclrow = NULL
   }

   selratio = array(0, dim = c(nvar, ncomp, ny))
   for (y in 1:ny) {
      for (comp in 1:model$ncomp) {
         b = coeffs[, comp, y]
         bnorm = sqrt(as.numeric(crossprod(b)))
         w = b/bnorm

         ttp = x %*% w
         ptp = crossprod(ttp, x) / as.numeric(crossprod(ttp))

         expvar = ttp %*% ptp
         resvar = apply(x - expvar, 2, var)
         expvar = apply(expvar, 2, var)

         selratio[, comp, y] = expvar / resvar
      }
   }

   if (nexclvar > 0) {
      selratio.out = array(0, dim = c(nvar + nexclvar, ncomp, ny))
      selratio.out[-exclvars, , ] = selratio
   } else {
      selratio.out = selratio
   }

   dimnames(selratio.out) = dimnames(model$coeffs$values)
   selratio.out = mda.setattr(selratio.out, attrs)
   attr(selratio.out, 'name') = 'Selectivity ratio'
   selratio.out
}

#' VIP scores calculation for PLS model
#'
#' @description
#' Calculates VIP (Variable Importance in Projection) scores for each component and
#' response variable in the PLS model
#'
#' @param object
#' a PLS model (object of class \code{pls})
#'
#' @return
#' matrix \code{nvar x ny} with VIP score values for selected number of components
#'
pls.calculateVIPScores = function(object) {

   # get values
   comp = object$ncomp.selected
   coeffs = object$coeffs$values
   w = object$weights[, 1:comp, drop = F]
   xloads = object$xloadings[, 1:comp, drop = F];
   xscores = object$calres$xdecomp$scores[, 1:comp, drop = F];

   # remove hidden variables
   attrs = mda.getattr(coeffs)
   exclvars = attrs$exclrows
   nexclvar = length(exclvars)
   if (nexclvar > 0) {
      coeffs = coeffs[-exclvars, , , drop = FALSE]
      w = w[-exclvars, , drop = FALSE]
      xloads = xloads[-exclvars, , drop = FALSE]
   }

   # remove hidden objects
   if (length(attr(xscores, 'exclrows')) > 0)
      xscores = xscores[-attr(xscores, 'exclrows'), , drop = FALSE]

   ny = dim(coeffs)[3]
   nvar = dim(coeffs)[1]
   vipscores = matrix(0, nrow = nvar, ncol = ny)

   # regression coefficients for working with scores instead of x
   # T = X * WPW
   # T * WPW' = X * WPW * WPW'
   # T * WPW' * (WPW * WPW')^-1 = X
   # YP = X * b
   # YP = T * WPW' * (WPW * WPW')^-1 * b
   # YP = T * bT, where bT = WPW' * (WPW * WPW)^-1 * b
   wpw = (w %*% solve(t(xloads) %*% w))

   # normalise weights
   n = 1/sqrt(colSums(w^2))
   if (comp == 1)
      dim(n) = c(1, 1)
   else
      n = diag(n)
   wnorm = w %*% n

   for (y in 1:ny) {
      b = coeffs[, comp, y, drop = F]
      dim(b) = c(dim(b)[1], 1)

      bscores = ( t(wpw) %*% pinv(wpw %*% t(wpw)) ) %*% b

      TT = colSums(xscores^2)
      dim(TT) = c(1, comp)
      SS = bscores^2 * t(TT)
      vipscores[, y] = nvar * wnorm^2 %*% as.matrix(SS) / sum(SS)
   }

   if (nexclvar > 0) {
      vipscores.out = matrix(0, nrow = nvar + nexclvar, ncol = ny)
      vipscores.out[-exclvars, ] = vipscores
   } else {
      vipscores.out = vipscores
   }

   rownames(vipscores.out) = dimnames(object$coeffs$values)[[1]]
   colnames(vipscores.out) = dimnames(object$coeffs$values)[[3]]
   vipscores.out = mda.setattr(vipscores.out, attrs)
   attr(vipscores.out, 'name') = 'VIP scores'
   vipscores.out
}

#' Selectivity ratio for PLS model
#'
#' @description
#' Returns vector with selectivity ratio values for given number of components
#' and response variable
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to get the values for (if NULL user selected as optimal will be
#' used)
#' @param ny
#' which response to get the values for (if y is multivariate)
#' @param ...
#' other parameters
#'
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#'
#' @return
#' vector with selectivity ratio values
#'
#' @export
getSelectivityRatio.pls = function(obj, ncomp = NULL, ny = 1, ...) {
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected

   attrs = mda.getattr(obj$selratio)
   selratio = obj$selratio[, ncomp, ny]
   dim(selratio) = c(dim(obj$selratio)[1], 1)
   selratio = mda.setattr(selratio, attrs)
   rownames(selratio) = dimnames(obj$selratio)[[1]]
   colnames(selratio) = dimnames(obj$selratio)[[3]][[ny]]

   selratio
}

#' VIP scores for PLS model
#'
#' @description
#' Returns vector with VIP scores values for given number of components
#' and response variable
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ny
#' which response to get the values for (if y is multivariate)
#' @param ...
#' other parameters
#'
#' @references
#' [1] Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), pp. 103-112.
#'
#' @return
#' vector with VIP scores values
#'
#' @export
getVIPScores.pls = function(obj, ny = 1, ...) {
   return(mda.subset(obj$vipscores, select = ny))
}

#' VIP scores plot for PLS model
#'
#' @description
#' Shows a plot with VIP scores values for given number of components
#' and response variable
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ny
#' which response to get the values for (if y is multivariate)
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @references
#' [1] Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), pp. 103-112.
#'
#' @export
plotVIPScores.pls = function(obj, ny = 1, type = 'l', main = NULL, ylab = '', ...) {

   main = getMainTitle(main, NULL, 'VIP scores')
   vipscores = getVIPScores.pls(obj, ny)
   mdaplot(mda.t(vipscores), type = type, main = main, ylab = ylab, ...)
}

#' Selectivity ratio plot for PLS model
#'
#' @description
#' Shows a plot with selectivity ratio values for given number of components
#' and response variable
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to get the values for (if NULL user selected as optimal will be used)
#' @param ny
#' which response to get the values for (if y is multivariate)
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#'
#' @export
plotSelectivityRatio.pls = function(obj, ncomp = NULL, ny = 1, type = 'l',
                                    main = NULL, ylab = '', ...) {
   main = getMainTitle(main, ncomp, 'Selectivity ratio')
   ncomp = getSelectedComponents(obj, ncomp)

   selratio = getSelectivityRatio(obj, ncomp, ny)
   mdaplot(mda.t(selratio), type = type, main = main, ylab = ylab, ...)
}


#' Summary method for PLS model object
#'
#' @description
#' Shows performance statistics for the model.
#'
#' @param object
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ny
#' which y variable to show the summary for (if NULL, will be shown for all)
#' @param ...
#' other arguments
#'
#' @export
summary.pls = function(object, ncomp = NULL, ny = NULL, ...) {
   obj = object

   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   else if (ncomp <= 0 || ncomp > obj$ncomp)
      stop('Wrong value for number of components!')

   if (is.null(ny))
      ny = 1:ncol(obj$calres$y.ref)

   cat('\nPLS model (class pls) summary\n')
   cat('\nPerformance and validation:\n')
   cat(sprintf('Number of selected components: %d\n', ncomp))

   for (y in ny) {
      if (ncol(obj$calres$y.ref) > 1)
         cat(sprintf('\nResponse variable #%d (%s)\n', y, colnames(obj$calres$y.ref)[y]))

      data = as.matrix(obj$calres, ncomp = ncomp, ny = y)
      rownames(data) = 'Cal'

      if (!is.null(obj$cvres)) {
         data = rbind(data, as.matrix(obj$cvres, ncomp = ncomp, ny = y))
         rownames(data)[nrow(data)] = 'CV'
      }

      if (!is.null(obj$testres)) {
         data = rbind(data, as.matrix(obj$testres, ncomp = ncomp, ny = y))
         rownames(data)[nrow(data)] = 'Test'
      }

      data = data[, -c(1, 3), drop = F]
      data[, 1:2] = round(data[, 1:2], 2)
      data[, 3] = round(data[, 3], 3)
      data[, 4] = mdaplot.formatValues(data[, 4], round.only = T)
      data[, 5] = round(data[, 5], 3)
      data[, 6] = round(data[, 6], 4)
      data[, 7] = round(data[, 7], 1)

      print(data)
   }
   cat('\n')
}

#' Print method for PLS model object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' a PLS model (object of class \code{pls})
#' @param ...
#' other arguments
#'
#' @export
print.pls = function(x, ...) {
   obj = x

   cat('\nPLS model (class pls)\n')
   cat('\nCall:\n')
   print(obj$call)
   cat('\nMajor fields:\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$coeffs - object (regcoeffs) with regression coefficients\n')
   cat('$xloadings - vector with x loadings\n')
   cat('$yloadings - vector with y loadings\n')
   cat('$weights - vector with weights\n')
   cat('$calres - results for calibration set\n')

   if (!is.null(obj$cvres)) {
      cat('$cvres - results for cross-validation\n')
   }

   if (!is.null(obj$testres)) {
      cat('$testres - results for test set\n')
   }

   cat('\nTry summary(model) and plot(model) to see the model performance.\n')
}
