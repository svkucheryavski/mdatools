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
#' @param lim.type
#' which method to use for calculation of critical limits for residual distances (see details)
#' @param alpha
#' significance level for extreme limits for T2 and Q disances.
#' @param gamma
#' significance level for outlier limits for T2 and Q distances.
#' @param info
#' short text with information about the model.
#' @param ncomp.selcrit
#' criterion for selecting optimal number of components (\code{'min'} for
#' first local minimum of RMSECV and \code{'wold'} for Wold's rule.)
#'
#' @return
#' Returns an object of \code{pls} class with following fields:
#' \item{ncomp }{number of components included to the model.}
#' \item{ncomp.selected }{selected (optimal) number of components.}
#' \item{xcenter }{vector with values used to center the predictors (x).}
#' \item{ycenter }{vector with values used to center the responses (y).}
#' \item{xscale }{vector with values used to scale the predictors (x).}
#' \item{yscale }{vector with values used to scale the responses (y).}
#' \item{xloadings }{matrix with loading values for x decomposition.}
#' \item{yloadings }{matrix with loading values for y decomposition.}
#' \item{xeigenvals }{vector with eigenvalues of components (variance of x-scores).}
#' \item{yeigenvals }{vector with eigenvalues of components (variance of y-scores).}
#' \item{weights }{matrix with PLS weights.}
#' \item{coeffs }{object of class \code{\link{regcoeffs}} with regression coefficients calculated
#' for each component.}
#' \item{info }{information about the model, provided by user when build the model.}
#' \item{cv }{information cross-validation method used (if any).}
#' \item{res }{a list with result objects (e.g. calibration, cv, etc.)}
#'
#' @details
#' So far only SIMPLS method [1] is available. Implementation works both with one
#' and multiple response variables.
#'
#' Like in \code{\link{pca}}, \code{pls} uses number of components (\code{ncomp}) as a minimum of
#' number of objects - 1, number of x variables and the default or provided value. Regression
#' coefficients, predictions and other results are calculated for each set of components from 1
#' to \code{ncomp}: 1, 1:2, 1:3, etc. The optimal number of components, (\code{ncomp.selected}),
#' is found using first local minumum, but can be also forced to user defined value using function
#' (\code{\link{selectCompNum.pls}}). The selected optimal number of components is used for all
#' default operations - predictions, plots, etc.
#'
#' Cross-validation settings, \code{cv}, can be a number or a list. If \code{cv} is a number, it
#' will be used as a number of segments for random cross-validation (if \code{cv = 1}, full
#' cross-validation will be preformed). If it is a list, the following syntax can be used:
#' \code{cv = list("rand", nseg, nrep)} for random repeated cross-validation with \code{nseg}
#' segments and \code{nrep} repetitions or \code{cv = list("ven", nseg)} for systematic splits
#' to \code{nseg} segments ('venetian blinds').
#'
#' Calculation of confidence intervals and p-values for regression coefficients can by done
#' based on Jack-Knifing resampling. This is done automatically if cross-validation is used.
#' However it is recommended to use at least 10 segments for stable JK result. See help for
#' \code{\link{regcoeffs}} objects for more details.
#'
#' @references
#' 1. S. de Jong, Chemometrics and Intelligent Laboratory Systems 18 (1993) 251-263.
#' 2. Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), 35-48.
#' 3. Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), 103-112.
#'
#' @seealso
#' Main methods for \code{pls} objects:
#' \tabular{ll}{
#'  \code{print} \tab prints information about a \code{pls} object.\cr
#'  \code{\link{summary.pls}} \tab shows performance statistics for the model.\cr
#'  \code{\link{plot.pls}} \tab shows plot overview of the model.\cr
#'  \code{\link{pls.simpls}} \tab implementation of SIMPLS algorithm.\cr
#'  \code{\link{predict.pls}} \tab applies PLS model to a new data.\cr
#'  \code{\link{selectCompNum.pls}} \tab set number of optimal components in the model.\cr
#'  \code{\link{setDistanceLimits.pls}} \tab allows to change parameters for critical limits.\cr
#'  \code{\link{categorize.pls}} \tab categorize data rows similar to
#'    \code{\link{categorize.pca}}.\cr
#'  \code{\link{selratio}} \tab computes matrix with selectivity ratio values.\cr
#'  \code{\link{vipscores}} \tab computes matrix with VIP scores values.\cr
#' }
#'
#' Plotting methods for \code{pls} objects:
#' \tabular{ll}{
#'  \code{\link{plotXScores.pls}} \tab shows scores plot for x decomposition.\cr
#'  \code{\link{plotXYScores.pls}} \tab shows scores plot for x and y decomposition.\cr
#'  \code{\link{plotXLoadings.pls}} \tab shows loadings plot for x decomposition.\cr
#'  \code{\link{plotXYLoadings.pls}} \tab shows loadings plot for x and y decomposition.\cr
#'  \code{\link{plotXVariance.pls}} \tab shows explained variance plot for x decomposition.\cr
#'  \code{\link{plotYVariance.pls}} \tab shows explained variance plot for y decomposition.\cr
#'  \code{\link{plotXCumVariance.pls}} \tab shows cumulative explained variance plot for y
#'  decomposition.\cr
#'  \code{\link{plotYCumVariance.pls}} \tab shows cumulative explained variance plot for y
#'  decomposition.\cr
#'  \code{\link{plotXResiduals.pls}} \tab shows distance/residuals plot for x decomposition.\cr
#'  \code{\link{plotXYResiduals.pls}} \tab shows joint distance plot for x and y decomposition.\cr
#'  \code{\link{plotWeights.pls}} \tab shows plot with weights.\cr
#'  \code{\link{plotSelectivityRatio.pls}} \tab shows plot with selectivity ratio values.\cr
#'  \code{\link{plotVIPScores.pls}} \tab shows plot with VIP scores values.\cr
#' }
#'
#' Methods inherited from \code{regmodel} object (parent class for \code{pls}):
#' \tabular{ll}{
#'  \code{\link{plotPredictions.regmodel}} \tab shows predicted vs. measured plot.\cr
#'  \code{\link{plotRMSE.regmodel}} \tab shows RMSE plot.\cr
#'  \code{\link{plotYResiduals.regmodel}} \tab shows residuals plot for y values.\cr
#'  \code{\link{getRegcoeffs.regmodel}} \tab returns matrix with regression coefficients.\cr
#' }
#'
#' Most of the methods for plotting data (except loadings and regression coefficients) are also
#' available for PLS results (\code{\link{plsres}}) objects. There is also a randomization test
#' for PLS-regression (\code{\link{randtest}}) and implementation of interval PLS algorithm
#' for variable selection (\code{\link{ipls}})
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
   info = "", ncomp.selcrit = "min", lim.type = "ddmoments", alpha = 0.05, gamma = 0.01) {

   # if y is a vector, convert it to matrix
   if (is.null(dim(y))) {
      dim(y) <- c(length(y), 1)
   }

   # check calibration data and process excluded rows and columns
   x <- prepCalData(x, exclrows = exclrows, exclcols = exclcols, min.nrows = 2, min.ncols = 1)
   y <- prepCalData(y, exclrows = exclrows, exclcols = NULL, min.nrows = 2, min.ncols = 1)

   # build a model and apply to calibration set
   model <- pls.cal(x, y, ncomp, center = center, scale = scale, method = method, cv = cv)
   model$info <- info
   model$call <- match.call()

   # get calibration results
   model$res <- list()
   model$res[["cal"]] <- predict.pls(model, x, y)
   model$res[["cal"]]$info <- "calibration results"
   model$calres <- model$res[["cal"]]

   # compute critical limit parameters
   model$limParams <- list(
      "Q" = ldecomp.getLimParams(model$res[["cal"]]$xdecomp$Q),
      "T2" = ldecomp.getLimParams(model$res[["cal"]]$xdecomp$T2),
      "Z" = ldecomp.getLimParams(model$res[["cal"]]$ydecomp$Q)
   )

   # do cross-validation if needed
   if (!is.null(cv)) {
      cvres <- crossval.regmodel(model, x, y, cv, cal.fun = pls.cal)
      model$res[["cv"]] <- plsres(cvres$y.pred, cvres$y.ref, ncomp.selected = model$ncomp)
      model$res[["cv"]]$info <- "cross-validation results"
      model$cvres <- model$res[["cv"]]
      model$coeffs <- regcoeffs(model$coeffs$values, cvres$jk.coeffs)
   }

   # do test set validation if provided
   if (!is.null(x.test) && !is.null(y.test)) {
      model$res[["test"]] <- predict.pls(model, x.test, y.test)
      model$res[["test"]]$info <- "test set validation results"
      model$testres <- model$res[["test"]]
   }

   # select optimal number of components
   model$cv <- cv
   model$ncomp.selcrit <- ncomp.selcrit
   model <- selectCompNum(model, selcrit = ncomp.selcrit)

   # set distance limits
   model <- setDistanceLimits(model, lim.type = lim.type, alpha = alpha, gamma = gamma)

   return(model)
}

#' Select optimal number of components for PLS model
#'
#' @description
#' Allows user to select optimal number of components for PLS model
#'
#' @param obj
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
#' points out to the optimal number of components in the model. You can either specify it manually,
#' as argument \code{ncomp}, or use one of the algorithms for automatic selection.
#'
#' Automatic selection by default based on cross-validation statistics. If no cross-validation
#' results are found in the model, the method will use test set validation results. If they are
#' not available as well, the model will use calibration results and give a warning as in this case
#' the selected number of components will lead to overfitted model.
#'
#' There are two algorithms for automatic selection you can chose between: either first local
#' minimum of RMSE (`selcrit="min"`) or Wold's rule (`selcrit="wold"`).
#'
#' The first local minimum criterion finds at which component, A, error of prediction starts
#' raising and selects (A - 1) as the optimal number. The Wold's criterion finds which component A
#' does not make error smaller at least by 5% comparing to the previous value and selects (A - 1)
#' as the optimal number.
#'
#' If model is PLS2 model (has several response variables) the method computes optimal number of
#' components for each response and returns the smallest value. For example, if for the first
#' response 2 components give the smallest error and for the second response this number is 3,
#' A = 2 will be selected as a final result.
#'
#' It is not recommended to use automatic selection for real applications, always investigate
#' your model (via RMSE, Y-variance plot, regression coefficients) to make correct decision.
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
   if (is.null(selcrit)) selcrit <- ""
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

#' Compute and set statistical limits for residual distances.
#'
#' @description
#' Computes statisticsl limits for orthogonal and score distances (x-decomposition) and
#' orthogonal distance (y-decomposition) based on calibration set and assign the calculated
#' values as model properties.
#'
#' @param obj
#' object with PLS model
#' @param lim.type
#' type of limits ("jm", "chisq", "ddmoments", "ddrobust")
#' @param alpha
#' significance level for detection of extreme objects
#' @param gamma
#' significance level for detection of outliers (for data driven approach)
#' @param ...
#' other arguments
#'
#' @details
#'
#' The limits can be accessed as fields of model objects: \code{$Qlim}, \code{$T2lim}, and
#' \code{$Zlim}. Each is a matrix with four rows and \code{ncomp} columns. In case of limits
#' for x-decomposition, first row contains critical limits for extremes, second row - for outliers,
#' third row contains mean value for corresponding distances (or its robust estimate in case of
#' \code{lim.type = "ddrobust"}) and last row contains the degrees of freedom.
#'
#' @return
#' Object models with the three fields updated.
#'
#' @export
setDistanceLimits.pls <- function(obj, lim.type = obj$lim.type, alpha = obj$alpha,
   gamma = obj$gamma, ...) {

   obj$T2lim <- ldecomp.getT2Limits(lim.type, alpha, gamma, obj$limParams)
   obj$Qlim <- ldecomp.getQLimits(lim.type, alpha, gamma, obj$limParams,
      obj$res[["cal"]]$xdecomp$residuals, obj$xeigenvals)
   obj$Zlim <- pls.getZLimits(lim.type, alpha, gamma, obj$limParams)

   obj$alpha <- alpha
   obj$gamma <- gamma
   obj$lim.type <- lim.type

   attr(obj$res[["cal"]]$xdecomp$Q, "u0") <- obj$Qlim[3, ]
   attr(obj$res[["cal"]]$xdecomp$Q, "Nu") <- obj$Qlim[4, ]

   attr(obj$res[["cal"]]$xdecomp$T2, "u0") <- obj$T2lim[3, ]
   attr(obj$res[["cal"]]$xdecomp$T2, "Nu") <- obj$T2lim[4, ]

   attr(obj$res[["cal"]]$ydecomp$Q, "u0") <- obj$Zlim[3, ]
   attr(obj$res[["cal"]]$ydecomp$Q, "Nu") <- obj$Zlim[4, ]

   obj$calres <- obj$res[["cal"]]

   if (!is.null(obj$res$test)) {
      attr(obj$res[["test"]]$xdecomp$Q, "u0") <- obj$Qlim[3, ]
      attr(obj$res[["test"]]$xdecomp$Q, "Nu") <- obj$Qlim[4, ]

      attr(obj$res[["test"]]$xdecomp$T2, "u0") <- obj$T2lim[3, ]
      attr(obj$res[["test"]]$xdecomp$T2, "Nu") <- obj$T2lim[4, ]

      attr(obj$res[["test"]]$ydecomp$Q, "u0") <- obj$Zlim[3, ]
      attr(obj$res[["test"]]$ydecomp$Q, "Nu") <- obj$Zlim[4, ]

      obj$testres <- obj$res[["test"]]
   }

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
#' PLS results (an object of class \code{\link{plsres}})
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

   # check datasets and convert to matrix if needed
   x <- prepCalData(x, min.nrows = 1, min.ncols = nrow(object$xloadings) - length(x.attrs$exclcols))

   # get dimensions
   nresp <- dim(object$coeffs$values)[3]
   ncomp <- dim(object$coeffs$values)[2]

   # check dimensions of predictors
   if (ncol(x) != dim(object$coeffs$values)[1]) {
      stop("Wrong number of columns in matrix with predictors (x).")
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
   y.ref <- NULL
   ydecomp <- NULL
   if (!is.null(y)) {

      if (is.null(dim(y))) dim(y) <- c(length(y), 1)

      if (nrow(x) != nrow(y)) {
         stop("Matrices with predictors (x) and response (y) should have the same number of rows.")
      }

      if (ncol(y) != nresp) {
         stop("Wrong number of columns in matrix with response values (y).")
      }

      # keep the original y values as reference
      y.ref <- y

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
      attr(ydecomp$Q, "u0") <- object$Zlim[3, ]
      attr(ydecomp$Q, "Nu") <- object$Zlim[4, ]
   }

   # unscale predicted y values
   yp <- if (is.numeric(object$yscale)) sweep(yp, 3, object$yscale, "*") else yp

   # uncenter predicted y values
   yp <- if (is.numeric(object$ycenter)) sweep(yp, 3, object$ycenter, "+") else yp

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

   # add u0 and Nu parameters as arguments, so the residuals can be normalized
   attr(xdecomp$Q, "u0") <- object$Qlim[3, ]
   attr(xdecomp$Q, "Nu") <- object$Qlim[4, ]

   attr(xdecomp$T2, "u0") <- object$T2lim[3, ]
   attr(xdecomp$T2, "Nu") <- object$T2lim[4, ]

   return(
      plsres(yp, y.ref = y.ref, ncomp.selected = object$ncomp.selected,
         xdecomp = xdecomp, ydecomp = ydecomp)
   )
}

#' Categorize data rows based on PLS results and critical limits for total distance.
#'
#' @description
#' The method uses full distance for decomposition of X-data and squared Y-residuals of PLS results
#' from \code{res} with critical limits computed for the PLS model and categorizes the
#' corresponding objects as "regular", "extreme" or "outlier".
#'
#' @param obj
#' object with PCA model
#' @param res
#' object with PCA results
#' @param ncomp
#' number of components to use for the categorization
#' @param ...
#' other parameters
#'
#' @details
#' The method does not categorize hidden values if any. It is based on the approach described in
#' [1] and works only if data driven approach is used for computing critical limits.
#'
#' @return
#' vector (factor) with results of categorization.
#'
#' @references
#' 1. Rodionova O. Ye., Pomerantsev A. L. Detection of Outliers in Projection-Based Modeling.
#' Analytical Chemistry (2020, in publish). doi: 10.1021/acs.analchem.9b04611
#'
#' @export
categorize.pls <- function(obj, res = obj$res$cal, ncomp = obj$ncomp.selected, ...) {

   create_categories <- function(extremes_ind, outliers_ind) {
      categories <- rep(1, length(extremes_ind))
      categories[extremes_ind] <- 2
      categories[outliers_ind] <- 3
      return(factor(categories, levels = 1:3, labels = c("regular", "extreme", "outlier")))
   }

   # if not data driven - quit
   if (!(obj$lim.type %in% c("ddmoments", "ddrobust"))) {
      stop("categorize.pls() works only with data driven limit types ('ddoments' or 'ddrobust').")
   }

   # get distance values for selected number of components
   h <- res$xdecomp$T2[, ncomp]
   q <- res$xdecomp$Q[, ncomp]
   z <- res$ydecomp$Q[, ncomp]

   # remove excluded values if any
   rows_excluded <- attr(res$xdecomp$Q, "exclrows")
   if (length(rows_excluded) > 0) {
      h <- h[-rows_excluded]
      q <- q[-rows_excluded]
      z <- z[-rows_excluded]
   }

   # get DoF
   Nh <- obj$T2lim[4, ncomp]
   Nq <- obj$Qlim[4, ncomp]
   Nz <- obj$Zlim[4, ncomp]
   Nf <- Nq + Nh

   # get scale factor
   h0 <- obj$T2lim[3, ncomp]
   q0 <- obj$Qlim[3, ncomp]
   z0 <- obj$Zlim[3, ncomp]

   # process degrees of freedom for (Z)
   Nz <- round(Nz)
   Nz[Nz < 1] <- 1
   Nz[Nz > 250] <- 250

   # process degrees of freedom for (F)
   Nf <- round(Nf)
   Nf[Nf < 1] <- 1
   Nf[Nf > 250] <- 250

   # compute total distance and DoF for it
   g <- Nh * h / h0 + Nq * q / q0 + Nz * z / z0
   Ng <- Nh + Nq + Nz
   nobj <- nrow(obj$res$cal$xdecomp$scores)

   # compute limits for total distance
   ext_lim <- qchisq(1 - obj$alpha, Ng)
   out_lim <- qchisq((1 - obj$gamma) ^ (1 / nobj), Ng)

   outliers_ind <- g > out_lim
   extremes_ind <- g > ext_lim & g < out_lim

   return(create_categories(extremes_ind, outliers_ind))
}

#' Summary method for PLS model object
#'
#' @description
#' Shows performance statistics for the model.
#'
#' @param object
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to count.
#' @param ny
#' which y variables to show the summary for (can be a vector)
#' @param ...
#' other arguments
#'
#' @export
summary.pls <- function(object, ncomp = object$ncomp.selected,
   ny = seq_len(nrow(object$yloadings)), ...) {

   if (length(ncomp) != 1 || ncomp < 0 || ncomp > object$ncomp) {
      stop("Wrong value for the 'ncomp' parameter.")
   }

   cat("\nPLS model (class pls) summary\n")
   cat("-------------------------------\n")
   fprintf("Info: %s\n", object$info)
   fprintf("Number of selected components: %d\n", ncomp)
   fprintf("Cross-validation: %s\n", crossval.str(object$cv))

   cat("\n")
   for (y in ny) {
      fprintf("Response variable: %s\n", rownames(object$yloadings)[y])
      out <- do.call(rbind, lapply(object$res, as.matrix, ncomp = ncomp, ny = y))
      rownames(out) <- capitalize(names(object$res))

      if (!any(is.na(out[, 1:4]))) out[, 1:4] <- round(out[, 1:4], 3)
      out[, 5] <- round(out[, 5], 3)
      out[, 6] <- mdaplot.formatValues(out[, 6], round.only = T)
      out[, 7] <- round(out[, 7], 3)
      out[, 8] <- round(out[, 8], 4)
      out[, 9] <- round(out[, 9], 2)

      print(out[, -c(1, 3), drop = FALSE])
      cat("\n")
   }
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
print.pls <- function(x, ...) {
   cat("\nPLS model (class pls)\n")
   cat("\nCall:\n")
   print(x$call)

   cat("\nMajor fields:\n")
   cat("$ncomp - number of calculated components\n")
   cat("$ncomp.selected - number of selected components\n")
   cat("$coeffs - object (regcoeffs) with regression coefficients\n")
   cat("$xloadings - vector with x loadings\n")
   cat("$yloadings - vector with y loadings\n")
   cat("$weights - vector with weights\n")
   cat("$res - list with results (calibration, cv, etc)\n")

   cat("\nTry summary(model) and plot(model) to see the model performance.\n")
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
#' @param main
#' title for the plot
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
#' @param main
#' title for the plot
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
#' @param main
#' title for the plot
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
#' @param main
#' title for the plot
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
#' which decomposition to use ("xdecomp" for x or "ydecomp" for y)
#' @param variance
#' which variance to use ("expvar", "cumexpvar")
#' @param type
#' type of the plot("b", "l" or "h")
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
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param main
#' main plot title
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
#' @param ncomp
#' which component to show the plot for
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

#' Residual distance plot for decomposition of X data
#'
#' @description
#' Shows a plot with orthogonal distance vs score distance for PLS decomposition of X data.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (by default optimal value selected for the model will be used)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param cgroup
#' color grouping of plot points (works only if one result object is available)
#' @param xlim
#' limits for x-axis
#' @param ylim
#' limits for y-axis
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
#' @param main
#' title for the plot
#' @param legend.position
#' position of legend (if shown)
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' The function is almost identical to \code{\link{plotResiduals.pca}}.
#'
#' @export
plotXResiduals.pls <- function(obj, ncomp = obj$ncomp.selected, norm = TRUE, log = FALSE,
   main = sprintf("X-distances (ncomp = %d)", ncomp), cgroup = NULL, xlim = NULL, ylim = NULL,
   show.limits = c(TRUE, TRUE), lim.col = c("darkgray", "darkgray"), lim.lwd = c(1, 1),
   lim.lty = c(2, 3), show.legend = TRUE, legend.position = "topright", res = obj$res, ...) {

   # get xdecomp from list with result objects
   res <- lapply(res, function(x) if ("ldecomp" %in% class(x$xdecomp)) x$xdecomp)
   res <- res[!sapply(res, is.null)]

   ldecomp.plotResiduals(res, obj$Qlim, obj$T2lim, ncomp = ncomp, log = log, norm = norm,
      cgroup = cgroup, xlim = xlim, ylim = ylim, show.limits = show.limits, lim.col = lim.col,
      lim.lwd = lim.lwd, show.legend = show.legend, main = main, ...)
}

#' Residual XY-distance plot
#'
#' @description
#' Shows a plot with full X-distance (f) vs. orthogonal Y-distance (z) for PLS model results.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' how many components to use (by default optimal value selected for the model will be used)
#' @param log
#' logical, apply log tranformation to the distances or not (see details)
#' @param norm
#' logical, normalize distance values or not (see details)
#' @param cgroup
#' color grouping of plot points (works only if one result object is available)
#' @param xlim
#' limits for x-axis
#' @param ylim
#' limits for y-axis
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
#' @param main
#' title for the plot
#' @param legend.position
#' position of legend (if shown)
#' @param res
#' list with result objects to show the plot for (by defaul, model results are used)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' The function presents a way to identify extreme objects and outliers based on both full distance
#' for X-decomposition (known as f) and squared residual distance for Y-decomposition (z). The
#' approach has been proposed in [1].
#'
#' The plot is available only if data driven methods (classic or robust) have been used for
#' computing of critical limits.
#'
#' @references
#' 1. Rodionova O. Ye., Pomerantsev A. L. Detection of Outliers in Projection-Based Modeling.
#' Analytical Chemistry (2020, in publish). doi: 10.1021/acs.analchem.9b04611
#'
#' @export
plotXYResiduals.pls <- function(obj, ncomp = obj$ncomp.selected, norm = TRUE, log = FALSE,
   main = sprintf("XY-distances (ncomp = %d)", ncomp), cgroup = NULL, xlim = NULL, ylim = NULL,
   show.limits = c(TRUE, TRUE), lim.col = c("darkgray", "darkgray"), lim.lwd = c(1, 1),
   lim.lty = c(2, 3), show.legend = TRUE, legend.position = "topright", res = obj$res, ...) {

   if (!(obj$lim.type %in% c("ddmoments", "ddrobust"))) {
      stop("plotXYResiduals() works only with data driven limit types ('ddoments' or 'ddrobust')")
   }

   # generate values for cgroup if categories should be used
   if (length(cgroup) == 1 && cgroup == "categories") {
      cgroup <- categorize(obj, res[[1]], ncomp = ncomp)
   }

   # get xdecomp from list with result objects
   res <- lapply(res, function(x) if ("ldecomp" %in% class(x$xdecomp)) x)
   res <- res[!sapply(res, is.null)]

   # function to compute plot limits
   getPlotLim <- function(lim, pd, ld, dim, show.limits) {
      if (!is.null(lim) || all(!show.limits)) return(lim)
      limits <- if (show.limits[[2]]) ld$outliers else ld$extremes
      return(c(0, max(sapply(pd, function(x) max(x[, dim])), limits[, dim])) * 1.05)
   }

   # check that show.limits is logical
   if (!all(is.logical(show.limits))) {
      stop("Parameter 'show.limits' must have logical value(s).")
   }

   # if show.limits has only one value - duplicate it
   if (length(show.limits) == 1) {
      show.limits <- rep(show.limits, 2)
   }

   # compute plot data for each result object
   plot_data <- lapply(res, plotXYResiduals.plsres, ncomp = ncomp, norm = norm, log = log,
      show.plot = FALSE)

   # get coordinates for critical limits
   lim_data <- pls.getLimitsCoordinates(obj$Qlim, obj$T2lim, obj$Zlim,
      ncomp = ncomp, nobj = obj$limParams$Q$moments$nobj, norm = norm, log = log)

   xlim <- getPlotLim(xlim, plot_data, lim_data, 1, show.limits)
   ylim <- getPlotLim(ylim, plot_data, lim_data, 2, show.limits)

   # make plot
   if (length(plot_data) == 1) {
      mdaplot(plot_data[[1]], type = "p", xlim = xlim, ylim = ylim, cgroup = cgroup, ...)
   } else {
      mdaplotg(plot_data, type = "p", xlim = xlim, ylim = ylim, show.legend = show.legend,
         legend.position = legend.position, ...)
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
#' @param show.legend
#' logical, show or not legend on the plot (when it is available)
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
#' @param show.legend
#' logical, show or not a legend
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

#' VIP scores plot for PLS model
#'
#' @description
#' Shows a plot with VIP scores values for given number of components
#' and response variable
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ny
#' which response to plot the values for (if y is multivariate), can be a vector.
#' @param ncomp
#' number of components to count
#' @param type
#' type of the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See \code{\link{vipscores}} for more details.
#'
#' @export
plotVIPScores.pls <- function(obj, ny = 1, ncomp = obj$ncomp.selected,
   type = "l", ...) {

   vipscores <- vipscores(obj, ncomp = ncomp)
   mdaplotg(mda.t(mda.subset(vipscores, select = ny)), type = type, ...)
}

#' Selectivity ratio plot for PLS model
#'
#' @description
#' Computes and shows a plot for Selectivity ratio values for given number of components
#' and response variable
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ny
#' which response to plot the values for (if y is multivariate), can be a vector.
#' @param ncomp
#' number of components to count
#' @param type
#' type of the plot
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#'
#' @details
#' See \code{\link{vipscores}} for more details.
#'
#' @export
plotSelectivityRatio.pls <- function(obj, ny = 1,
   ncomp = obj$ncomp.selected, type = "l", ...) {

   selratio <- selratio(obj, ncomp = ncomp)
   mdaplotg(mda.t(mda.subset(selratio, select = ny)), type = type, ...)
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
      c <- as.numeric(crossprod(w, (M %*% w)))

      # stop cycle since c-value is very small and can result in singular matrix
      if (c < .Machine$double.eps) {
         n <- n - 1
         warning(paste0(
            "PLS can not compute more than ", n, " components (eigenvalues are too small). "
         ), call. = FALSE)
         break
      }

      w <- w / sqrt(c)
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
#' @param cv
#' logical, is model calibrated during cross-validation or not (or cv settings for calibration)
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
   if (!is.logical(cv) && !is.null(cv)) {
      nseg <- if (is.numeric(cv)) cv else cv[[2]]
      nobj.cv <- if (nseg == 1) 1 else ceiling(xc.nrows / nseg)

      # we set cv to FALSE so fitting knows that it is not a part of cross-validation
      cv <- FALSE
   }

   # set cv to FALSE also if it was null (needed for correct call of pls.run() method)
   if (is.null(cv)) cv <- FALSE

   # correct maximum number of components
   ncomp <- min(xc.ncols, xc.nrows - 1 - nobj.cv, ncomp)

   # fit model
   fit <- pls.run(x, y, method = method, ncomp = ncomp, cv = cv)
   model$ncomp <- ncomp <- fit$ncomp

   # if it is for cross-validation return the results as is
   if (is.logical(cv) && cv) {
      model$coeffs <- regcoeffs(fit$coeffs)
      model$xloadings <- fit$xloadings
      model$weights <- fit$weights
      return(model)
   }

   # compute x-scores and residuals
   xscores <- x %*% (fit$weights %*% solve(crossprod(fit$xloadings, fit$weights)))
   yscores <- as.matrix(y) %*% fit$yloadings

   # compute eigenvalues
   xeigenvals <- colSums(xscores^2) / (xc.nrows - 1)
   attr(xeigenvals, "DoF") <- (xc.nrows - 1)
   yeigenvals <- colSums(yscores^2) / (xc.nrows - 1)
   attr(yeigenvals, "DoF") <- (xc.nrows - 1)

   # correct results related to predictors for missing columns in x
   # corresponding rows will be set to 0 and excluded
   xloadings <- matrix(0, nrow = x.ncols, ncol = ncomp)
   yloadings <- matrix(0, nrow = y.ncols, ncol = ncomp)
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

#' VIP scores for PLS model
#'
#' @description
#' Calculates VIP (Variable Importance in Projection) scores for predictors for given number
#' of components and response variable.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to count
#'
#' @return
#' matrix \code{nvar x ny} with VIP score values (columns correspond to responses).
#'
#' @details
#' May take some time in case of large number of predictors Returns results as a column-vector,
#' with all necessary attributes inherited (e.g. xaxis.values, excluded variables, etc.). If you
#' want to make a plot use for example: \code{mdaplot(mda.t(v), type = "l")}, where \code{v} is
#' a vector with computed VIP scores. Or just try \code{\link{plotVIPScores.pls}}.
#'
#' @references
#' [1] Il-Gyo Chong, Chi-Hyuck Jun. Chemometrics and Laboratory Systems, 78 (2005), pp. 103-112.
#'
#' @export
vipscores <- function(obj, ncomp = obj$ncomp.selected) {

   if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
      stop("Wrong value for the 'ncomp' parameter.")
   }

   # subset needed model parameters
   comp <- seq_len(ncomp)
   weights <- obj$weights[, comp, drop = FALSE]
   yloads <- obj$yloadings[, comp, drop = FALSE];

   # get eigenvalues and multiply them to degrees of freedom
   xeigenvals <- obj$xeigenvals[comp]

   # get number and indices of variables and adjust dimension for regcoeffs
   nvar <- nrow(weights)
   var_ind <- seq_len(nvar)

   # remove hidden variables
   if (length(obj$exclcols) > 0) {
      weights <- weights[-obj$exclcols, , drop = FALSE]
      var_ind <- var_ind[-obj$exclcols]
   }

   # prepare matrix for vipscores
   vipscores <- matrix(0, nrow = nvar, ncol = nrow(yloads))

   # normalize scores
   wnorm <- weights %*% diag(1 / sqrt(colSums(weights^2)), nrow = ncomp, ncol = ncomp)

   # compute sum of squares for explained y variance and normalize it
   ssq <- yloads^2 %*% diag(xeigenvals, nrow = ncomp, ncol = ncomp)
   ssq <- ssq %*% diag(1 / rowSums(ssq), nrow = ncomp, ncol = ncomp)

   # compute VIP scores
   vipscores[var_ind, ] <- sqrt(nvar * wnorm^2 %*% t(ssq))

   rownames(vipscores) <- rownames(obj$xloadings)
   colnames(vipscores) <- rownames(obj$yloadings)

   attr(vipscores, "exclrows") <- obj$exclcols
   attr(vipscores, "yaxis.values") <- attr(obj$xloadings, "yaxis.values")
   attr(vipscores, "yaxis.name") <- attr(obj$xloadings, "yaxis.name")
   attr(vipscores, "xaxis.name") <- ""
   attr(vipscores, "name") <- sprintf("VIP scores (ncomp = %d)", ncomp)

   return(vipscores)
}

#' VIP scores for PLS model
#'
#' @description
#' Returns vector with VIP scores values. This function is a proxy for \code{\link{vipscores}}
#' and will be removed in future releases.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to count
#' @param ...
#' other parameters
#'
#' @return
#' matrix \code{nvar x 1} with VIP score values
#'
#' @export
getVIPScores.pls <- function(obj, ncomp = obj$ncomp.selected, ...) {

   warning("This function is deprecated and will be removed in future. Use 'vipscores()' insted.")
   return(vipscores(obj, ncomp = ncomp))
}

#' Selectivity ratio calculation
#'
#' @description
#' Calculates selectivity ratio for each component and response variable in
#' the PLS model
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to count
#'
#' @references
#' [1] Tarja Rajalahti et al. Chemometrics and Laboratory Systems, 95 (2009), pp. 35-48.
#'
#' @return
#' array \code{nvar x ncomp x ny} with selectivity ratio values
#'
#' @export
selratio <- function(obj, ncomp = obj$ncomp.selected) {

   if (length(ncomp) != 1 || ncomp < 1 || ncomp > obj$ncomp) {
      stop("Wrong value for the 'ncomp' parameter.")
   }

   # get number and indices of variables and adjust dimension for regcoeffs
   nvar <- nrow(obj$weights)
   nresp <- nrow(obj$yloadings)
   var_ind <- seq_len(nvar)

   # reproduce x values
   xresiduals <- obj$res[["cal"]]$xdecomp$residuals
   xscores <- obj$res[["cal"]]$xdecomp$scores
   x <- xresiduals + tcrossprod(xscores, obj$xloadings)

   # remove excluded rows
   if (length(obj$exclrows) > 0) {
      x <- x[-obj$exclrows, , drop = FALSE]
   }

   # subset needed model parameters
   coeffs <- obj$coeffs$values[, ncomp, , drop = FALSE]

   # correct dimension for coefficients
   dim(coeffs) <- c(nvar, nresp)

   # remove hidden variables
   if (length(obj$exclcols) > 0) {
      x <- x[, -obj$exclcols, drop = FALSE]
      coeffs <- coeffs[-obj$exclcols, , drop = FALSE]
      var_ind <- var_ind[-obj$exclcols]
   }

   # prepare matrix for vipscores
   selratio <- matrix(0, nrow = nvar, ncol = nresp)

   # get norm value for regression coefficients
   bnorm <- sqrt(colSums(coeffs^2))

   # compute target projections
   ttp <-  x %*% (coeffs %*% diag(1 / bnorm, nrow = nresp, ncol = nresp))
   ptp <- t(crossprod(x, ttp) %*% diag(1 / colSums(ttp^2), nrow = nresp, ncol = nresp))

   # compute selectivity ratio
   for (y in seq_len(nresp)) {
      expvar <- ttp[, y, drop = FALSE] %*% ptp[y, , drop = FALSE]
      selratio[var_ind, y] <- colSums(expvar^2) / colSums((x - expvar)^2)
   }

   rownames(selratio) <- rownames(obj$xloadings)
   colnames(selratio) <- rownames(obj$yloadings)

   attr(selratio, "exclrows") <- obj$exclcols
   attr(selratio, "yaxis.values") <- attr(obj$xloadings, "yaxis.values")
   attr(selratio, "yaxis.name") <- attr(obj$xloadings, "yaxis.name")
   attr(selratio, "xaxis.name") <- ""
   attr(selratio, "name") <- sprintf("Selectivity ratio (ncomp = %d)", ncomp)

   return(selratio)
}

#' Selectivity ratio for PLS model
#'
#' @description
#' Returns vector with Selectivity ratio values. This function is a proxy for \code{\link{selratio}}
#' and will be removed in future releases.
#'
#' @param obj
#' a PLS model (object of class \code{pls})
#' @param ncomp
#' number of components to get the values for (if NULL user selected as optimal will be used)
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
getSelectivityRatio.pls <- function(obj, ncomp = obj$ncomp.selected, ...) {
   warning("This function is deprecated and will be removed in future. Use 'selratio()' insted.")
   return(selratio(obj, ncomp = ncomp))
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
#'
#' @export
pls.getZLimits <- function(lim.type, alpha, gamma, params) {

   if (!(lim.type %in% c("ddmoments", "ddrobust"))) {
      return(NULL)
   }

   pZ <- if (regexpr("robust", lim.type) > 0) params$Z$robust else params$Z$moments
   DoF <- round(pZ$Nu)
   DoF[DoF < 1] <- 1
   DoF[DoF > 250] <- 250

   ncomp <- length(pZ$u0)
   lim <- rbind(0, 0, pZ$u0, DoF)

   colnames(lim) <- paste("Comp", seq_len(ncomp))
   rownames(lim) <- c("Extremes limits", "Outliers limits", "Mean", "DoF")
   attr(lim, "name") <- "Critical limits for orthogonal distance (Z)"
   attr(lim, "alpha") <- alpha
   attr(lim, "gamma") <- gamma
   attr(lim, "lim.type") <- lim.type

   return(lim)
}

#' Compute coordinates of lines or curves with critical limits
#'
#' @param Qlim
#' matrix with critical limits for orthogonal distances (X)
#' @param T2lim
#' matrix with critical limits for score distances (X)
#' @param Zlim
#' matrix with critical limits for orthogonal distances (Y)
#' @param nobj
#' number of objects to compute the limits for
#' @param ncomp
#' number of components for computing the coordinates
#' @param norm
#' logical, shall distance values be normalized or not
#' @param log
#' logical, shall log transformation be applied or not
#'
#' @return
#' list with two matrices (x and y coordinates of corresponding limits)
#'
#' @export
pls.getLimitsCoordinates <- function(Qlim, T2lim, Zlim, nobj, ncomp, norm, log) {

   # get DoF
   Nh <- T2lim[4, ncomp]
   Nq <- Qlim[4, ncomp]
   Nz <- Zlim[4, ncomp]
   Nf <- Nq + Nh

   # get scaling factor
   z0 <- Zlim[3, ncomp]
   f0 <- Nf

   # process degrees of freedom for (Z)
   Nz <- round(Nz)
   Nz[Nz < 1] <- 1
   Nz[Nz > 250] <- 250

   # process degrees of freedom for (F)
   Nf <- round(Nf)
   Nf[Nf < 1] <- 1
   Nf[Nf > 250] <- 250

   # get limit parameters
   alpha <- attr(Qlim, "alpha")
   gamma <- attr(Qlim, "gamma")

   ## slope and intercepts
   eB <- qchisq(1 - alpha, Nf + Nz) / Nz * z0
   oB <- qchisq((1 - gamma) ^ (1 / nobj), Nf + Nz) / Nz * z0
   eA <- oA <- -1 * (z0 / f0) * (Nf / Nz)

   fE <- seq(-0.95, -eB / eA, length.out = 100)
   fO <- seq(-0.95, -oB / oA, length.out = 100)
   zE <- eA * fE + eB
   zO <- oA * fO + oB

   if (norm) {
      fE <- fE / f0
      zE <- zE / z0
      fO <- fO / f0
      zO <- zO / z0
   }

   if (log) {
      fE <- log(1 + fE)
      zE <- log(1 + zE)
      fO <- log(1 + fO)
      zO <- log(1 + zO)
   }

   return(list(
      extremes = cbind(fE, zE),
      outliers = cbind(fO, zO)
   ))
}
