#' SIMCA one-class classification
#'
#' @description
#' \code{simca} is used to make SIMCA (Soft Independent Modelling of Class Analogies) model for
#' one-class classification.
#'
#' @param x
#' a numerical matrix with data values.
#' @param classname
#' short text (up to 20 symbols) with class name.
#' @param ncomp
#' maximum number of components to calculate.
#' @param center
#' logical, do mean centering of data or not.
#' @param scale
#' logical, do sdandardization of data or not.
#' @param cv
#' cross-validation settings (see details).
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values)
#' @param x.test
#' a numerical matrix with test data.
#' @param c.test
#' a vector with classes of test data objects (can be text with names of classes or logical).
#' @param method
#' method to compute principal components.
#' @param rand
#' vector with parameters for randomized PCA methods (if NULL, conventional PCA is used instead)
#' @param lim.type
#' which method to use for calculation of critical limits for residuals (see details)
#' @param alpha
#' significance level for calculating critical limits for T2 and Q residuals.
#' @param gamma
#' significance level for calculating outlier limits for T2 and Q residuals.
#' @param info
#' text with information about the model.
#'
#' @details
#' SIMCA is in fact PCA model with additional functionality, so \code{simca} class inherits most
#' of the functionality of \code{\link{pca}} class. It uses critical limits calculated for Q and T2
#' residuals calculated for PCA model for making classification decistion.
#'
#' Cross-validation settings, \code{cv}, can be a number or a list. If \code{cv} is a number, it
#' will be used as a number of segments for random cross-validation (if \code{cv = 1}, full
#' cross-validation will be preformed). If it is a list, the following syntax can be used:
#' \code{cv = list('rand', nseg, nrep)} for random repeated cross-validation with \code{nseg}
#' segments and \code{nrep} repetitions or \code{cv = list('ven', nseg)} for systematic splits
#' to \code{nseg} segments ('venetian blinds').
#'
#' @return
#' Returns an object of \code{simca} class with following fields:
#' \item{classname }{a short text with class name.}
#' \item{modpower }{a matrix with modelling power of variables.}
#' \item{calres }{an object of class \code{\link{simcares}} with classification results for a
#' calibration data.}
#' \item{testres }{an object of class \code{\link{simcares}} with classification results for a test
#' data, if it was provided.}
#' \item{cvres }{an object of class \code{\link{simcares}} with classification results for
#' cross-validation, if this option was chosen.}
#'
#' Fields, inherited from \code{\link{pca}} class:
#' \item{ncomp }{number of components included to the model.}
#' \item{ncomp.selected }{selected (optimal) number of components.}
#' \item{loadings }{matrix with loading values (nvar x ncomp).}
#' \item{eigenvals }{vector with eigenvalues for all existent components.}
#' \item{expvar }{vector with explained variance for each component (in percent).}
#' \item{cumexpvar }{vector with cumulative explained variance for each component (in percent).}
#' \item{T2lim }{statistical limit for T2 distance.}
#' \item{Qlim }{statistical limit for Q residuals.}
#' \item{info }{information about the model, provided by user when build the model.}
#'
#' @references
#' S. Wold, M. Sjostrom. "SIMCA: A method for analyzing chemical data in terms of similarity and
#' analogy" in B.R. Kowalski (ed.), Chemometrics Theory and Application, American Chemical Society
#' Symposium Series 52, Wash., D.C., American Chemical Society, p. 243-282.
#'
#' @seealso
#' Methods for \code{simca} objects:
#' \tabular{ll}{
#'  \code{print.simca} \tab shows information about the object.\cr
#'  \code{summary.simca} \tab shows summary statistics for the model.\cr
#'  \code{plot.simca} \tab makes an overview of SIMCA model with four plots.\cr
#'  \code{\link{predict.simca}} \tab applies SIMCA model to a new data.\cr
#'  \code{\link{plotModellingPower.simca}} \tab shows plot with modelling power of variables.\cr
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
#' Methods, inherited from \code{\link{pca}} class:
#' \tabular{ll}{
#'  \code{\link{selectCompNum.pca}} \tab set number of optimal components in the model\cr
#'  \code{\link{plotScores.pca}} \tab shows scores plot.\cr
#'  \code{\link{plotLoadings.pca}} \tab shows loadings plot.\cr
#'  \code{\link{plotVariance.pca}} \tab shows explained variance plot.\cr
#'  \code{\link{plotCumVariance.pca}} \tab shows cumulative explained variance plot.\cr
#'  \code{\link{plotResiduals.pca}} \tab shows Q vs. T2 residuals plot.\cr
#' }
#'
#' @examples
#' ## make a SIMCA model for Iris setosa class with full cross-validation
#' library(mdatools)
#'
#' data = iris[, 1:4]
#' class = iris[, 5]
#'
#' # take first 20 objects of setosa as calibration set
#' se = data[1:20, ]
#'
#' # make SIMCA model and apply to test set
#' model = simca(se, 'setosa', cv = 1)
#' model = selectCompNum(model, 1)
#'
#' # show infromation, summary and plot overview
#' print(model)
#' summary(model)
#' plot(model)
#'
#' # show predictions
#' par(mfrow = c(2, 1))
#' plotPredictions(model, show.labels = TRUE)
#' plotPredictions(model, res = 'calres', ncomp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' # show performance, modelling power and residuals for ncomp = 2
#' par(mfrow = c(2, 2))
#' plotSensitivity(model)
#' plotMisclassified(model)
#' plotModellingPower(model, ncomp = 2, show.labels = TRUE)
#' plotResiduals(model, ncomp = 2)
#' par(mfrow = c(1, 1))
#'
#' @export
simca <- function(x, classname, ncomp = min(nrow(x) - 1, ncol(x), 20),
   x.test = NULL, c.test = NULL, cv = NULL, ...) {

   if (!is.character(classname)) {
      stop("Argument 'classname' must be a text.")
   }

   if (length(classname) > 20) {
      stop("Argument 'classname' must have up to 20 symbols.")
   }

   # calibrate model
   model <- pca(x, ncomp = ncomp, ...)
   model$nclasses <- 1
   model$classnames <- c(classname)
   model$call <- match.call()
   class(model) <- c("simca", "classmodel", "pca")

   # apply model to calibration set
   model$res[["cal"]] <- predict(model, x, c.ref = rep(classname, nrow(x)))
   model$calres <- model$res[["cal"]]

   # do cross-validation if needed
   if (!is.null(cv)) {
      model$res[["cv"]] <- crossval.simca(model, x, cv)
      model$cvres <- model$res[["cv"]]
   }

   # apply model to test set if provided
   if (!is.null(x.test)){
      # if classes are not provided we assume the object are from the same class
      if (is.null(c.test)) c.test <- as.factor(rep(classname, nrow(x.test)))
      model$res[["test"]] <- predict.simca(model, x.test, c.ref = c.test)
      model$testres <- model$res[["test"]]
   }

   model
}


#' Probabilities of class belonging for PCA/SIMCA results
#'
#' @details
#' Computes p-value for every object being from the same populaion as calibration set
#' based on its orthogonal and score distances.
#'
#' @param obj
#' object with PCA model
#' @param ncomp
#' number of components to compute the probability for
#' @param q
#' vector with squared orthogonal distances for given number of components
#' @param h
#' vector with score distances for given number of components
#'
#' @export
getProbabilities.simca <- function(obj, ncomp, q, h) {

   p <- function(x, alpha){ ifelse((p <- (1 - x) / alpha * 0.5) > 1, 1, p) }
   return(p(getProbabilities.pca(obj, ncomp, q, h), obj$alpha))
}

#' SIMCA predictions
#'
#' @description
#' Applies SIMCA model to a new data set
#'
#' @param object
#' a SIMCA model (object of class \code{simca})
#' @param x
#' a matrix with x values (predictors)
#' @param c.ref
#' a vector with reference class names (same as class names for models)
#' @param cal
#' logical, are predictions for calibration set or not
#' @param ...
#' other arguments
#'
#' @return
#' SIMCA results (an object of class \code{simcares})
#'
#' @details
#' See examples in help for \code{\link{simca}} function.
#'
#' @export
predict.simca = function(object, x, c.ref = NULL, cal = FALSE, ...) {

   # get PCA results
   pca.res <- predict.pca(object, x)

   # do classification and set attributes
   class.res <- simca.classify(object, pca.res, c.ref)
   return(simcares(pca.res, class.res))
}

#' SIMCA classification
#'
#' @description
#' Make classification based on calculated T2 and Q values and corresponding limits
#'
#' @param obj
#' a SIMCA model (object of class \code{simca})
#' @param res
#' results of projection data to PCA space
#'
#' @return
#' vector with predicted class values (\code{c.pred})
#'
#' @details
#' This is a service function for SIMCA class, do not use it manually.
#'
simca.classify <- function(obj, pca.res, c.ref = NULL) {

   # check reference values
   c.ref <- classmodel.processRefValues(c.ref, object$classnames)

   ncomp <- obj$ncomp
   nobj <- nrow(pca.res$Q)
   c.pred <- array(0, dim = c(nobj, ncomp, 1))
   p.pred <- array(0, dim = c(nobj, ncomp, 1))

   for (i in seq_len(ncomp)) {
      p.pred[, i, 1] <- getProbabilities.simca(obj, i, pca.res$Q[, i], pca.res$T2[, i])
   }

   c.pred <- (p.pred >= 0.5) * 2 - 1
   dimnames(c.pred) <- dimnames(p.pred) <- list(
      rownames(pca.res$scores),
      colnames(obj$loadings),
      obj$classnames
   )

   c.pred <- mda.setattr(c.pred, mda.getattr(pca.res$scores), "row")
   attr(c.pred, "name") <- "Class, predicted"
   p.pred <- mda.setattr(p.pred, mda.getattr(pca.res$scores), "row")
   attr(p.pred, "name") <- "Class, probabilities"

   # combine everything to classres object
   return(
      classres(
         c.pred = c.pred,
         p.pred = p.pred,
         c.ref = c.ref,
         ncomp.selected = obj$ncomp.selected
      )
   )
}

#' Cross-validation of a SIMCA model
#'
#' @description
#' Does the cross-validation of a SIMCA model
#'
#' @param obj
#' a SIMCA model (object of class \code{simca})
#' @param x
#' a matrix with x values (predictors from calibration set)
#' @param cv
#' number of segments (if cv = 1, full cross-validation will be used)
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#'
#' @return
#' object of class \code{simcares} with results of cross-validation
#'
crossval.simca = function(obj, x, cv) {
   ncomp <- obj$ncomp

   # convert data to a matrix
   attrs <- mda.getattr(x)
   x <- mda.df2mat(x)

   # remove excluded rows
   if (length(attrs$exclrows) > 0) {
      x <- x[-attrs$exclrows, , drop = F]
   }

   # remove excluded columns
   if (length(attrs$exclcols) > 0) {
      x <- x[, -attrs$exclcols, drop = F]
   }

   # get matrix with indices for cv segments
   nobj <- nrow(x)
   idx <- crossval(nobj, cv)
   nseg <- nrow(idx);
   nrep <- dim(idx)[3]

   p.pred <- array(0, dim = c(nobj, ncomp, 1))

   # loop over segments
   for (iRep in seq_len(nrep)) {
      for (iSeg in seq_len(nseg)) {
         ind <- na.exclude(idx[iSeg, ,iRep])

         if (length(ind) > 0) {
            x.cal <- x[-ind, , drop = F]
            x.val <- x[ind, , drop = F]

            # calibrate PCA model and set distance limits
            m.loc <- pca.cal(x.cal, obj$ncomp, center = obj$center, scale = obj$scale,
               method = obj$method, rand = obj$rand)

            # set distance limits
            m.loc <- setDistanceLimits(m.loc, lim.type = obj$lim.type, alpha = obj$alpha,
               gamma = obj$gamma)

            # make prediction for validation subset
            res.loc <- predict.pca(m.loc, x.val)

            # compute and save probabilities
            for (i in seq_len(obj$ncomp)) {
               p.pred[ind, i, 1] <- p.pred[ind, i, 1] +
                  getProbabilities.simca(m.loc, i, res.loc$Q[, i], res.loc$T2[, i])
            }
         }
      }
   }

   # prepare predicted and reference values
   p.pred <- p.pred / nrep
   c.pred <- (p.pred >= 0.5) * 2 - 1
   c.ref <- rep(obj$classnames[1], nobj)

   dimnames(c.pred) <- dimnames(p.pred) <- list(
      rownames(x),
      colnames(obj$loadings),
      obj$classnames
   )

   attr(c.pred, "name") <- "Class, predicted"
   attr(p.pred, "name") <- "Class, probabilities"
   attr(p.pred, "yaxis.name") <- attr(c.pred, "yaxis.name") <- attr(x, "yaxis.name")
   attr(p.pred, "yaxis.values") <- attr(c.pred, "yaxis.values") <- attr(x, "yaxis.values")

   # combine everything to simcares object and return
   return(
      simcares(
         pca.res = NULL,
         class.res = classres(
            c.pred = c.pred,
            p.pred = p.pred,
            c.ref = c.ref,
            ncomp.selected = obj$ncomp.selected
         )
      )
   )
}

#' Model overview plot for SIMCA
#'
#' @description
#' Shows a set of plots for SIMCA model.
#'
#' @param x
#' a SIMCA model (object of class \code{simca})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{simcam}} function.
#'
#' @export
plot.simca <- function(x, comp = c(1, 2), ncomp = x$ncomp.selected, ...) {
   par(mfrow = c(2, 2))
   plotScores(x, comp, ...)
   plotLoadings(x, ncomp = ncomp, ...)
   plotResiduals(x, ncomp = ncomp, ...)
   plotCumVariance(x, ...)
   par(mfrow = c(1, 1))
}

#' Summary method for SIMCA model object
#'
#' @description
#' Shows performance statistics for the model.
#'
#' @param object
#' a SIMCA model (object of class \code{simca})
#' @param ...
#' other arguments
#'
#' @export
summary.simca = function(object, ...) {

   fprintf("\nSIMCA model for class '%s' summary\n\n", object$classname)

   if (!is.null(object$info) && nchar(object$info)) {
      fprintf("Info: %s\n", object$info)
   }

   if (!is.null(object$rand)) {
      fprintf("\nParameters for randomized algorithm: q = %d, p = %d\n",
         object$rand[1], object$rand[2])
   }

   if (length(object$exclrows) > 0) {
      fprintf("Excluded rows: %d\n", length(object$exclrows))
   }

   if (length(object$exclcols) > 0) {
      fprintf("Excluded coumns: %d\n", length(object$exclcols))
   }

   fprintf("\nNumber of selected components: %d\n", object$ncomp.selected)
   fprintf("Type of limits: %s\n", object$lim.type)
   fprintf("Alpha: %s\n", object$alpha)
   fprintf("Gamma: %s\n", object$gamma)
   cat("\n")

   summary_data <- NULL
   for (r in object$res) {
      summary_data <- rbind(summary_data, as.matrix.simcares(r, ncomp = object$ncomp.selected))
   }

   rownames(summary_data) <- capitalize(names(object$res))
   print(summary_data)
   cat("\n")
}

#' Print method for SIMCA model object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' a SIMCA model (object of class \code{simca})
#' @param ...
#' other arguments
#'
#' @export
print.simca = function(x, ...) {
   obj = x

   cat('\nSIMCA one class model (class simca)\n')

   cat('\nCall:\n')
   print(obj$call)

   cat('\nMajor fields:\n')
   cat('$classname - name of the class\n')
   cat('$alpha - significance level\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$loadings - matrix with loadings\n')
   cat('$eigenvals - eigenvalues for components\n')
   cat('$center - values for centering data\n')
   cat('$scale - values for scaling data\n')
   cat('$info - information about the model\n')
   cat('$cv - number of segments for cross-validation\n')
   cat('$calres - results (scores, etc) for calibration set\n')

   if (!is.null(obj$cvres)) {
      cat('$cvres - results for cross-validation\n')
   }

   if (!is.null(obj$testres)) {
      cat('$testres - results for test set\n')
   }
}

