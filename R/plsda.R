#' Partial Least Squares Discriminant Analysis
#'
#' @description
#' \code{plsda} is used to calibrate, validate and use of partial least squares discrimination
#' analysis (PLS-DA) model.
#'
#' @param x
#' matrix with predictors.
#' @param c
#' vector with class membership (should be either a factor with class names/numbers in case of
#' multiple classes or a vector with logical values in case of one class model).
#' @param ncomp
#' maximum number of components to calculate.
#' @param center
#' logical, center or not predictors and response values.
#' @param scale
#' logical, scale (standardize) or not predictors and response values.
#' @param cv
#' cross-validation settings (see details).
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
#' @param x.test
#' matrix with predictors for test set.
#' @param c.test
#' vector with reference class values for test set (same format as calibration values).
#' @param method
#' method for calculating PLS model.
#' @param lim.type
#' which method to use for calculation of critical limits for residual distances (see details)
#' @param alpha
#' significance level for extreme limits for T2 and Q disances.
#' @param gamma
#' significance level for outlier limits for T2 and Q distances.
#' @param info
#' short text with information about the model.
#' @param ncomp.selcrit
#' criterion for selecting optimal number of components (\code{'min'} for first local minimum of
#' RMSECV and \code{'wold'} for Wold's rule.)
#' @param classname
#' name (label) of class in case if PLS-DA is used for one-class discrimination model. In this case
#' it is expected that parameter `c` will be a vector with logical values.
#'
#' @return
#' Returns an object of \code{plsda} class with following fields (most inherited from class
#' \code{pls}):
#' \item{ncomp }{number of components included to the model.}
#' \item{ncomp.selected }{selected (optimal) number of components.}
#' \item{xloadings }{matrix with loading values for x decomposition.}
#' \item{yloadings }{matrix with loading values for y (c)  decomposition.}
#' \item{weights }{matrix with PLS weights.}
#' \item{coeffs }{matrix with regression coefficients calculated for each component.}
#' \item{info }{information about the model, provided by user when build the model.}
#' \item{calres }{an object of class \code{\link{plsdares}} with PLS-DA results for a calibration
#' data.}
#' \item{testres }{an object of class \code{\link{plsdares}} with PLS-DA results for a test data,
#' if it was provided.}
#' \item{cvres }{an object of class \code{\link{plsdares}} with PLS-DA results for cross-validation,
#' if this option was chosen.}
#'
#' @details
#' The \code{plsda} class is based on \code{pls} with extra functions and plots covering
#' classification functionality. All plots for \code{pls} can be used. E.g. of you want to see the
#' real predicted values (y in PLS) instead of classes use \code{plotPredictions.pls(model)} instead
#' of \code{plotPredictions(model)}.
#'
#' Cross-validation settings, \code{cv}, can be a number or a list. If \code{cv} is a number, it
#' will be used as a number of segments for random cross-validation (if \code{cv = 1}, full
#' cross-validation will be preformed). If it is a list, the following syntax can be used:
#' \code{cv = list('rand', nseg, nrep)} for random repeated cross-validation with \code{nseg}
#' segments and \code{nrep} repetitions or \code{cv = list('ven', nseg)} for systematic splits
#' to \code{nseg} segments ('venetian blinds').
#'
#' Calculation of confidence intervals and p-values for regression coefficients are available
#' only by jack-knifing so far. See help for \code{\link{regcoeffs}} objects for details.
#'
#' @seealso
#' Specific methods for \code{plsda} class:
#' \tabular{ll}{
#'  \code{print.plsda} \tab prints information about a \code{pls} object.\cr
#'  \code{summary.plsda} \tab shows performance statistics for the model.\cr
#'  \code{plot.plsda} \tab shows plot overview of the model.\cr
#'  \code{\link{predict.plsda}} \tab applies PLS-DA model to a new data.\cr
#' }
#'
#' Methods, inherited from \code{classmodel} class:
#' \tabular{ll}{
#'  \code{\link{plotPredictions.classmodel}} \tab shows plot with predicted values.\cr
#'  \code{\link{plotSensitivity.classmodel}} \tab shows sensitivity plot.\cr
#'  \code{\link{plotSpecificity.classmodel}} \tab shows specificity plot.\cr
#'  \code{\link{plotMisclassified.classmodel}} \tab shows misclassified ratio plot.\cr
#' }
#'
#' See also methods for class \code{\link{pls}}.
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @examples
#' ### Examples for PLS-DA model class
#'
#' library(mdatools)
#'
#' ## 1. Make a PLS-DA model with full cross-validation and show model overview
#'
#' # make a calibration set from iris data (3 classes)
#' # use names of classes as class vector
#' x.cal = iris[seq(1, nrow(iris), 2), 1:4]
#' c.cal = iris[seq(1, nrow(iris), 2), 5]
#'
#' model = plsda(x.cal, c.cal, ncomp = 3, cv = 1, info = 'IRIS data example')
#' model = selectCompNum(model, 1)
#'
#' # show summary and basic model plots
#' # misclassification will be shown only for first class
#' summary(model)
#' plot(model)
#'
#' # summary and model plots for second class
#' summary(model, nc = 2)
#' plot(model, nc = 2)
#'
#' # summary and model plot for specific class and number of components
#' summary(model, nc = 3, ncomp = 3)
#' plot(model, nc = 3, ncomp = 3)
#'
#' ## 2. Show performance plots for a model
#' par(mfrow = c(2, 2))
#' plotSpecificity(model)
#' plotSensitivity(model)
#' plotMisclassified(model)
#' plotMisclassified(model, nc = 2)
#' par(mfrow = c(1, 1))
#'
#' ## 3. Show both class and y values predictions
#' par(mfrow = c(2, 2))
#' plotPredictions(model)
#' plotPredictions(model, res = "cal", ncomp = 2, nc = 2)
#' plotPredictions(structure(model, class = "regmodel"))
#' plotPredictions(structure(model, class = "regmodel"), ncomp = 2, ny = 2)
#' par(mfrow = c(1, 1))
#'
#' ## 4. All plots from ordinary PLS can be used, e.g.:
#' par(mfrow = c(2, 2))
#' plotXYScores(model)
#' plotYVariance(model)
#' plotXResiduals(model)
#' plotRegcoeffs(model, ny = 2)
#' par(mfrow = c(1, 1))
#'
#' @export
plsda <- function(x, c, ncomp = min(nrow(x) - 1, ncol(x), 20), center = TRUE, scale = FALSE,
   cv = NULL, exclcols = NULL, exclrows = NULL, x.test = NULL, c.test = NULL, method = "simpls",
   lim.type = "ddmoments", alpha = 0.05, gamma = 0.01, info = "", ncomp.selcrit = "min",
   classname = NULL) {

   # check vector with class values and adjust if necessary
   c <- classmodel.processRefValues(c, classname)
   classnames <- if (!is.null(classname)) classname else levels(as.factor(c))
   classnames <- classnames[classnames != "None"]

   # get reference y-values for regression
   y <- mda.df2mat(as.factor(c), full = TRUE)
   y <- y[, colnames(y) != "None", drop = FALSE]
   y[y == 0] <- -1
   colnames(y) <- classnames
   rownames(y) <- rownames(x)

   # exclude columns if "exclcols" is provided
   if (length(exclcols) > 0) {
      x <- mda.exclcols(x, exclcols)
   }

   # exclude rows if "exclrows" is provided
   if (length(exclrows) > 0) {
      x <- mda.exclrows(x, exclrows)
   }

   # build a model and apply to calibration set
   model <- pls.cal(x, y, ncomp, center = center, scale = scale, method = method)
   model$classnames <- classnames
   model$nclasses <- length(model$classnames)
   class(model) <- c("pls", "regmodel")
   model$res <- list()

   # do cross-validation first (if specified)
   cvres <- NULL
   if (!is.null(cv)) {

      # cross-validation for regression model
      res <- crossval.regmodel(model, x, y, cv, pls.cal)

      # classification results
      c.pred <- classify.plsda(model, res$y.pred)
      c.ref <- if (length(attr(x, "exclrows") > 0)) c[-attr(x, "exclrows")] else c
      class_res <- classres(c.pred, c.ref = c.ref, p.pred = res$y.pred, ncomp.selected = ncomp)

      # regression results
      pls_res <- plsres(res$y.pred, res$y.ref, ncomp.selected = ncomp)

      # redefine regression coefficients
      model$coeffs <- regcoeffs(model$coeffs$values, res$jk.coeffs)

      # redefine
      cvres <- plsdares(pls_res, class_res)
      cvres$info <- "cross-validation results"
   }

   # add additional class names, info and call
   class(model) <- c("plsda", "classmodel", class(model))
   model$info <- info
   model$call <- match.call()

   # get calibration results
   model$res[["cal"]] <- predict.plsda(model, x, c)
   model$res[["cal"]]$info <- "calibration results"
   model$calres <- model$res[["cal"]]

   # compute critical limit parameters
   model$limParams <- list(
      "Q" = ldecomp.getLimParams(model$res[["cal"]]$xdecomp$Q),
      "T2" = ldecomp.getLimParams(model$res[["cal"]]$xdecomp$T2),
      "Z" = ldecomp.getLimParams(model$res[["cal"]]$ydecomp$Q)
   )

   # assign cross-validation results to the model (so they are under calibration)
   model$res[["cv"]] <- cvres
   model$cvres <- model$res[["cv"]]

   # do test set validation if provided
   if (!is.null(x.test) && !is.null(c.test)) {
      model$res[["test"]] <- predict.plsda(model, x.test, c.test)
      model$res[["test"]]$info <- "test set validation results"
      model$testres <- model$res[["test"]]
   }

   model$ncomp.selcrit <- ncomp.selcrit
   model <- selectCompNum(model, selcrit = ncomp.selcrit)

   # set distance limits
   model <- setDistanceLimits(model, lim.type = lim.type, alpha = alpha, gamma = gamma)

   model$cv <- cv
   return(model)
}

#' PLS-DA predictions
#'
#' @description
#' Applies PLS-DA model to a new data set
#'
#' @param object
#' a PLS-DA model (object of class \code{plsda})
#' @param x
#' a matrix with x values (predictors)
#' @param c.ref
#' a vector with reference class values (should be a factor)
#' @param ...
#' other arguments
#'
#' @return
#' PLS-DA results (an object of class \code{plsdares})
#'
#' @details
#' See examples in help for \code{\link{plsda}} function.
#'
#' @export
predict.plsda <- function(object, x, c.ref = NULL, ...) {

   y.ref <- NULL
   # prepare matrix with y-reference values according to classes used in the model
   if (!is.null(c.ref)) {
      attrs <- mda.getattr(c.ref)
      c.ref <- classmodel.processRefValues(c.ref, object$classnames)
      y.ref <- sapply(object$classnames, function(c) c == c.ref) * 2 - 1
      y.ref <- mda.setattr(y.ref, attrs)
      colnames(y.ref) <- object$classnames
      rownames(y.ref) <- rownames(x)
   }

   # do PLS predictions
   plsres <- predict.pls(object, x, y.ref)

   # classify objects and set attributes
   c.pred <- classify.plsda(object, plsres$y.pred)

   # combine everything to plsdares object
   cres <- classres(c.pred, c.ref = c.ref, p.pred = plsres$y.pred,
      ncomp.selected = object$ncomp.selected)
   return(plsdares(plsres, cres))
}

#' PLS-DA classification
#'
#' @description
#' Converts PLS predictions of y values to predictions of classes
#'
#' @param model
#' a PLS-DA model (object of class \code{plsda})
#' @param y
#' a matrix with predicted y values
#'
#' @return
#' Classification results (an object of class \code{classres})
#'
#' @details
#' This is a service function for PLS-DA class, do not use it manually.
#'
classify.plsda <- function(model, y) {
   c.pred <- ifelse(y < 0, -1, 1)
   c.pred <- mda.setattr(c.pred, mda.getattr(y))
   attr(c.pred, "name") <- "Class, predicted values"
   dimnames(c.pred) <- dimnames(y)

   return(c.pred)
}

#' Model overview plot for PLS-DA
#'
#' @description
#' Shows a set of plots (x residuals, regression coefficients, misclassification ratio and
#' predictions) for PLS-DA model.
#'
#' @param x
#' a PLS-DA model (object of class \code{plsda})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param nc
#' which class to show the plots
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{plsda}} function.
#'
#' @export
plot.plsda <- function(x, ncomp = x$ncomp.selected, nc = 1, show.legend = TRUE, ...) {

   if (!is.null(ncomp) && (ncomp <= 0 || ncomp > x$ncomp)) {
      stop("Wrong value for number of components.")
   }

   par(mfrow = c(2, 2))
   plotXResiduals(x, ncomp = ncomp, show.legend = show.legend, ...)
   plotRegcoeffs(x, ncomp = ncomp, ny = nc, ...)
   plotMisclassified(x, nc = nc, show.legend = show.legend, ...)
   plotPredictions(x, ncomp = ncomp, show.colorbar = show.legend, ...)
   par(mfrow = c(1, 1))
}

#' Summary method for PLS-DA model object
#'
#' @description
#' Shows some statistics for the model.
#'
#' @param object
#' a PLS-DA model (object of class \code{plsda})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param nc
#' which class to show the summary for (if NULL, will be shown for all)
#' @param ...
#' other arguments
#'
#' @export
summary.plsda <- function(object, ncomp = object$ncomp.selected,
   nc = seq_len(object$nclasses), ...) {

   if (length(ncomp) != 1 || ncomp < 0 || ncomp > object$ncomp) {
      stop("Wrong value for the 'ncomp' parameter.")
   }

   cat("\nPLS-DA model (class plsda) summary\n")
   cat("------------------------------------\n")
   fprintf("Info: %s\n", object$info)
   fprintf("Number of selected components: %d\n", ncomp)
   fprintf("Cross-validation: %s\n", crossval.str(object$cv))

   cat("\n")
   for (n in nc) {
      fprintf("Class #%d (%s)\n", n, object$classnames[n])
      out <- do.call(rbind, lapply(object$res, as.matrix, nc = n, ncomp = ncomp))
      rownames(out) <- capitalize(names(object$res))

      if (!any(is.na(out[, 1:4]))) out[, 1:4] <- round(out[, 1:4], 3)
      out[, 1:4] <- round(out[, 1:4], 2)
      print(out[, -c(1, 3), drop = FALSE])
      cat("\n")
   }
}

#' Print method for PLS-DA model object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' a PLS-DA model (object of class \code{plsda})
#' @param ...
#' other arguments
#'
#' @export
print.plsda <- function(x, ...) {
   cat("\nPLS-DA model (class plsda)\n")
   cat("\nCall:\n")
   print(x$call)

   cat("\nMajor fields:\n")
   cat("$ncomp - number of calculated components\n")
   cat("$ncomp.selected - number of selected components\n")
   cat("$coeffs - vector with regression coefficients\n")
   cat("$xloadings - vector with x loadings\n")
   cat("$yloadings - vector with Y loadings\n")
   cat("$weights - vector with weights\n")
   cat("$calres - results for calibration set\n")

   cat("\nTry summary(model) and plot(model) to see the model performance.\n")
}
