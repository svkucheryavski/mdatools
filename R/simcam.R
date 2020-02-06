#' SIMCA multiclass classification
#'
#' @description
#' \code{simcam} is used to combine several one-class SIMCA models for multiclass classification.
#'
#' @param models
#' list with SIMCA models (\code{simca} objects).
#' @param info
#' optional text with information about the the object.
#'
#' @details
#' Besides the possibility for multiclass classification, SIMCAM also provides tools for
#' investigation of relationship among individual models (classes), such as discrimination power of
#' variables, Cooman's plot, model distance, etc.
#'
#' When create \code{simcam} object, the calibration data from all individual SIMCA models is
#' extracted and combined for making predictions and calculate performance of the multi-class model.
#' The results are stored in \code{$calres} field of the model object.
#'
#' @return
#' Returns an object of \code{simcam} class with following fields:
#' \item{models }{a list with provided SIMCA models.}
#' \item{dispower }{an array with discrimination power of variables for each pair of individual
#' models.}
#' \item{moddist }{a matrix with distance between each each pair of individual models.}
#' \item{classnames }{vector with names of individual classes.}
#' \item{nclasses }{number of classes in the object.}
#' \item{info }{information provided by user when create the object.}
#' \item{calres }{an object of class \code{\link{simcamres}} with classification results for a
#' calibration data.}
#'
#' @seealso
#' Methods for \code{simca} objects:
#' \tabular{ll}{
#'  \code{print.simcam} \tab shows information about the object.\cr
#'  \code{summary.simcam} \tab shows summary statistics for the models.\cr
#'  \code{plot.simcam} \tab makes an overview of SIMCAM model with two plots.\cr
#'  \code{\link{predict.simcam}} \tab applies SIMCAM model to a new data.\cr
#'  \code{\link{plotModelDistance.simcam}} \tab shows plot with distance between individual
#'  models.\cr
#'  \code{\link{plotDiscriminationPower.simcam}} \tab shows plot with discrimination power.\cr
#'  \code{\link{plotCooman.simcam}} \tab shows Cooman's plot for calibration data.\cr
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
#' Since SIMCAM objects and results are calculated only for optimal number of components, there is
#' no sense to show such plots like sensitivity or specificity vs. number of components. However
#' they are available as for any other classification model.
#'
#' @examples
#' ## make a multiclass SIMCA model for Iris data
#' library(mdatools)
#'
#' # split data
#' caldata = iris[seq(1, nrow(iris), 2), 1:4]
#' x.se = caldata[1:25, ]
#' x.ve = caldata[26:50, ]
#' x.vi = caldata[51:75, ]
#'
#' x.test = iris[seq(2, nrow(iris), 2), 1:4]
#' c.test = iris[seq(2, nrow(iris), 2), 5]
#'
#' # create individual models
#' m.se = simca(x.se, classname = "setosa")
#' m.se = selectCompNum(m.se, 1)
#'
#' m.vi = simca(x.vi, classname = "virginica")
#' m.vi = selectCompNum(m.vi, 2)
#'
#' m.ve = simca(x.ve, classname = "versicolor")
#' m.ve = selectCompNum(m.ve, 1)
#'
#' # combine models into SIMCAM objects, show statistics and plots
#' m = simcam(list(m.se, m.vi, m.ve), info = "simcam model for Iris data")
#' summary(m)
#'
#' # show predictions and residuals for calibration data
#' par(mfrow = c(2, 2))
#' plotPredictions(m)
#' plotCooman(m, nc = c(1, 2))
#' plotModelDistance(m, nc = 1)
#' plotDiscriminationPower(m, nc = c(1, 2))
#' par(mfrow = c(1, 1))
#'
#' # apply the SIMCAM model to test set and show statistics and plots
#' res = predict(m, x.test, c.test)
#' summary(res)
#' plotPredictions(res)
#'
#' @export
simcam <- function(models, info = "") {
   nclasses <- length(models)
   classnames <- unlist(lapply(models, function(x) x$classnames[[1]]))

   model <- list()
   model$models <- models
   model$classnames <- classnames
   model$nclasses <- nclasses
   model$info <- info

   # calculate statistics
   stat <- simcam.getPerformanceStats(models, classnames)
   model <- c(model, stat)
   model$call <- match.call()
   class(model) <- c("simcam", "classmodel")

   # make predictions for calibration set
   model$res <- list()
   cal_data <- getCalibrationData.simcam(model)
   model$res[["cal"]] <- predict(model, cal_data$x, cal_data$c.ref)
   model$calres <- model$res[["cal"]]

   model
}

#' SIMCA multiple classes predictions
#'
#' @description
#' Applies SIMCAM model (SIMCA for multiple classes) to a new data set
#'
#' @param object
#' a SIMCAM model (object of class \code{simcam})
#' @param x
#' a matrix with x values (predictors)
#' @param c.ref
#' a vector with reference class names (same as class names in models)
#' @param ...
#' other arguments
#'
#' @return
#' SIMCAM results (an object of class \code{simcamres})
#'
#' @details
#' See examples in help for \code{\link{simcam}} function.
#'
#' @export
predict.simcam <- function(object, x, c.ref = NULL, ...) {

   attrs <- mda.getattr(x)
   c.pred <- array(0, dim = c(nrow(x), 1, object$nclasses))
   simca.res <- list()
   for (i in seq_len(object$nclasses)) {
      simca.res[[i]] <- predict.simca(object$models[[i]], x, c.ref)
      c.pred[, , i] <- simca.res[[i]]$c.pred[, object$models[[i]]$ncomp.selected, ]
   }

   c.pred <- mda.setattr(c.pred, attrs, "row")

   dimnames(c.pred) <- list(rownames(x), paste("Comp"), object$classnames)
   attr(c.pred, "name") <- "SIMCAM predictions"
   names(simca.res) <- object$classnames
   class.res <- classres(c.pred, c.ref)

   return(simcamres(class.res, simca.res))
}

#' Get calibration data
#'
#' @description
#' Get data, used for calibration of the SIMCAM individual models and combine to one dataset.
#'
#' @param obj
#' SIMCAM model (object of class \code{simcam})
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{simcam}} function.
#'
#' @export
getCalibrationData.simcam <- function(obj, ...) {
   x <- NULL
   c.ref <- NULL

   for (i in seq_len(obj$nclasses)) {
      x.loc <- getCalibrationData.pca(obj$models[[i]])
      c.ref <- c(c.ref, rep(obj$models[[i]]$classnames[1], nrow(x.loc)))
      x <- mda.rbind(x, x.loc)
   }

   return(list(x = x, c.ref = as.factor(c.ref)))
}

#' Summary method for SIMCAM model object
#'
#' @description
#' Shows performance statistics for the model.
#'
#' @param object
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' number of class to show summary for (can be vector)
#' @param ...
#' other arguments
#'
#' @export
summary.simcam <- function(object, nc = seq_len(object$nclasses), ...) {
   cat("\nSIMCA multiple classes classification (class simcam)\n")
   fprintf("\nNumber of classes: %d\n", length(nc))

   if (!is.null(object$info)) {
      cat(sprintf("Info: %s\n", object$info))
   }

   cat("\nSummary for calibration results\n")
   print(as.matrix.simcamres(object$res[["cal"]], nc = nc))
   cat("\n")
}

#' Print method for SIMCAM model object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' a SIMCAM model (object of class \code{simcam})
#' @param ...
#' other arguments
#'
#' @export
print.simcam <- function(x, ...) {
   cat("\nSIMCA multiple classes classification (class simcam)\n")

   cat("\nCall:\n")
   print(x$call)

   cat("\nMajor fields:\n")
   cat("$models - list wth individual SIMCA models for each class\n")
   cat("$classnames - vector with names of classes\n")
   cat("$moddist - matrix with distance between the models\n")
   cat("$dispower - matrix with discrimination power values\n")
   cat("$info - information about the object\n")
   cat("$res - list with calibration results\n")
}


################################
#  Static methods              #
################################


#' Performance statistics for SIMCAM model
#'
#' @description
#' Calculates discrimination power and distance between individual SIMCA models.
#'
#' @param models
#' list with SIMCA models (as provided to simcam class)
#' @param classnames
#' names of the classes for each model
#'
simcam.getPerformanceStats <- function(models, classnames) {

   # set up
   nc <- length(models)
   nvar <- nrow(models[[1]]$loadings)
   varnames <- rownames(models[[1]]$loadings)

   ## function to get calibration data and remove excluded rows
   getModelData <- function(m) {
      d <- getCalibrationData.pca(m)
      if (length(attr(d, "exclrows")) > 0) d <- d[-attr(d, "exclrows"), , drop = FALSE]
      return(d)
   }

   ## get datasets
   datasets <- lapply(models, getModelData)

   # loop through all combinations of classes
   dispower <- array(0, dim = c(nc, nc, nvar))
   moddist <- array(0, dim = c(nc, nc))
   class_comb <- expand.grid(seq_len(nc), seq_len(nc))

   for (i in seq_len(nrow(class_comb))) {

      # class indics
      nc1 <- class_comb[i, 1]
      nc2 <- class_comb[i, 2]

      # select two models
      m1 <- models[[nc1]]
      m2 <- models[[nc2]]

      # select two datasets
      d1 <- datasets[[nc1]]
      d2 <- datasets[[nc2]]

      # get calibration results
      r11 <- predict.pca(m1, d1)
      r22 <- predict.pca(m2, d2)

      # apply model 1 to data 2 and vice versa
      r12 <- predict.pca(m1, d2)
      r21 <- predict.pca(m2, d1)

      # calculate residuals for projections (first model)
      e11 <- r11$residuals
      e12 <- r12$residuals
      if (m1$ncomp.selected < m1$ncomp) {
         comp_ind <- (m1$ncomp.selected + 1):m1$ncomp
         e11 <- e11 + r11$scores[, comp_ind, drop = FALSE] %*%
            t(m1$loadings[, comp_ind, drop = FALSE])
         e12 <- e12 + r12$scores[, comp_ind, drop = FALSE] %*%
            t(m1$loadings[, comp_ind, drop = FALSE])
      }

      # calculate residuals for projections (second model model)
      e22 <- r22$residuals
      e21 <- r21$residuals
      if (m2$ncomp.selected < m2$ncomp) {
         comp_ind <- (m2$ncomp.selected + 1):m2$ncomp
         e22 <- e22 + r22$scores[, comp_ind, drop = FALSE] %*%
            t(m2$loadings[, comp_ind, drop = FALSE])
         e21 <- e21 + r21$scores[, comp_ind, drop = FALSE] %*%
            t(m2$loadings[, comp_ind, drop = FALSE])
      }

      # calculate variance for the residuals
      s11 <- colSums(e11^2) / (nrow(d1) - m1$ncomp.selected - 1)
      s22 <- colSums(e12^2) / (nrow(d2) - m2$ncomp.selected - 1)
      s12 <- colSums(e12^2) / (nrow(d1))
      s21 <- colSums(e21^2) / (nrow(d2))

      # calculate model distance and discrimination power
      dispower[nc1, nc2, ] <- sqrt((s12 + s21) / (s11 + s22))
      moddist[nc1, nc2] <- sqrt(sum(s12 + s21) / sum(s11 + s22))
   }

   dimnames(dispower) <- list(classnames, classnames, varnames)
   dimnames(moddist) <- list(classnames, classnames)

   attrs <- mda.getattr(models[[1]]$loadings)
   attr(dispower, "name") <- "Discrimination power"
   attr(dispower, "xaxis.name") <- attrs$yaxis.name
   attr(dispower, "xaxis.values") <- attrs$yaxis.values

   return(
      list(
         dispower = dispower,
         moddist = moddist
      )
   )
}


################################
#  Plotting methods            #
################################


#' Model distance plot for SIMCAM model
#'
#' @description
#' Shows a plot with distance between one SIMCA model to others.
#'
#' @param obj
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' one value - number of class (SIMCA model) to show the plot for
#' @param type
#' type of the plot ("h", "l" or "b")
#' @param xticks
#' vector with tick values for x-axis
#' @param xticklabels
#' vector with tick labels for x-axis
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' The plot shows similarity between a selected model and the others as a ratio of residual
#' variance using the following algorithm. Let's take two SIMCA/PCA models, m1 and m2, which
#' have optimal number of components A1 and A2. The models have been calibrated
#' using calibration sets X1 and X2 with number of rows n1 and n2.
#' Then we do the following:
#'
#' 1. Project X2 to model m1 and compute residuals, E12
#' 2. Compute variance of the residuals as s12 = sum(E12^2) / n1
#' 3. Project X1 to model m2 and compute residuals, E21
#' 4. Compute variance of the residuals as s21 = sum(E21^2) / n2
#' 5. Compute variance of residuals for m1 as s1 = sum(E1^2) / (n1 - A1 - 1)
#' 6. Compute variance of residuals for m2 as s2 = sum(E2^2) / (n2 - A2 - 1)
#'
#' The model distance then can be computed as: d = sqrt((s12 + s21) / (s1 + s2))
#'
#' As one can see, if the two models and corresponding calibration sets are identical, then the
#' distance will be sqrt((n - A - 1) / n). For example, if n = 25 and A = 2, then the distance
#' between the model and itself is sqrt(25/22) = sqrt(0.88) = 0.938. This case is demonstrated
#' in the example section.
#'
#' In general, if distance between models is below one classes are overlapping. If it is above 3
#' the classes are well separated.
#'
#'
#' @examples
#' # create two calibration sets with n = 25 objects in each
#' data(iris)
#' x1 <- iris[1:25, 1:4]
#' x2 <- iris[51:75, 1:4]
#'
#' # create to SIMCA models with A = 2
#' m1 <- simca(x1, 'setosa', ncomp = 2)
#' m2 <- simca(x2, 'versicolor', ncomp = 2)
#'
#' # combine the models into SIMCAM class
#' m <- simcam(list(m1, m2))
#'
#' # show the model distance plot with distance values as labels
#' # note, that distance between setosa and setosa is 0.938
#' plotModelDistance(m, show.labels = TRUE, labels = "values")
#'
#' @export
plotModelDistance.simcam <- function(obj, nc = 1, type = "h", xticks = seq_len(obj$nclasses),
   xticklabels = obj$classnames, main = paste0("Model distance (", obj$classnames[nc], ")"),
   xlab = "Models", ylab = "", ...) {

   if (length(nc) != 1 || nc < 1 || nc > obj$nclasses) {
      stop("Wrong values for 'nc' parameter.")
   }

   mdaplot(mda.t(obj$moddist[, nc, drop = F]), type = type, xticks = xticks,
      xticklabels = xticklabels, main = main, xlab = xlab, ylab = ylab, ...)
}

#' Discrimination power plot for SIMCAM model
#'
#' @description
#' Shows a plot with discrimination power of predictors for a pair of SIMCA models
#'
#' @param obj
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' vector with two values - classes (SIMCA models) to show the plot for
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' Discrimination power shows an ability of variables to separate classes. The power is computed
#' similar to model distance, using variance of residuals. However in this case instead of sum the
#' variance across all variables, we take the ratio separately for individual variables.
#'
#' Discrimination power equal or above 3 is considered as high.
#'
#' @export
plotDiscriminationPower.simcam <- function(obj, nc = c(1, 2), type = "h",
   main = paste0("Discrimination power: ", obj$classnames[nc[1]], " vs. ", obj$classname[nc[2]]),
   xlab = attr(obj$dispower, "xaxis.name"), ylab = "", ...) {

   if (length(nc) != 2 || min(nc) < 1 || max(nc) > obj$nclasses) {
      stop("Wrong values for 'nc' parameter.")
   }

   attrs <- mda.getattr(obj$dispower)
   varnames <- dimnames(obj$dispower)[[3]]
   nvar <- dim(obj$dispower)[3]

   plot_data <- obj$dispower[nc[1], nc[2], ]
   dim(plot_data) <- c(1, nvar)
   colnames(plot_data) <- varnames
   plot_data <- mda.setattr(plot_data, attrs)

   mdaplot(plot_data, type = type, main = main, xlab = xlab, ylab = ylab, ...)
}

#' Cooman's plot for SIMCAM model
#'
#' @description
#' Shows a Cooman's plot for a pair of SIMCA models
#'
#' @param obj
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' vector with two values - classes (SIMCA models) to show the plot for
#' @param res
#' list with results to show the plot for
#' @param groupby
#' factor to use for grouping points on the plot
#' @param main
#' title of the plot
#' @param show.limits
#' logical, show or not critical limits
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' Cooman's plot shows squared orthogonal distance from data points to two selected SIMCA models
#' as well as critical limits for the distance (optional). In case if critical limits must be shown
#' they are computed using chi-square distribution regardless which type of limits is employed for
#' classification.
#'
#' If only one result object is provided (e.g. results for calibration set or new predictions),
#' then the points can be color grouped using 'groupby' parameter (by default reference class values
#' are used to make the groups). In case of multiple result objects, the points are color grouped
#' according to the objects (e.g. calibration set and test set).
#'
#' @export
plotCooman.simcam <- function(obj, nc = c(1, 2), res = list("cal" = obj$res[["cal"]]),
   groupby = res[[1]]$c.ref, main = "Cooman's plot", show.limits = TRUE, ...) {

   if (!is.list(res)) {
      stop("Parameter 'res' should be list with 'simcamres' objects.")
   }

   if (length(nc) != 2 || min(nc) < 1 || max(nc) > obj$nclasses) {
      stop("Wrong values for 'nc' parameter.")
   }

   plot_data <- list()
   for (i in seq_along(res)) {
      plot_data[[names(res)[i]]] <- plotCooman.simcamres(res[[i]], nc = nc, show.plot = FALSE)
   }

   # compute limits using chi-square distribution
   show.lines <- FALSE
   if (show.limits) {
      m1 <- obj$models[[nc[1]]]
      m2 <- obj$models[[nc[2]]]
      show.lines <- c(
         chisq.crit(
            list(
               "u0" = m1$Qlim[3, m1$ncomp.selected],
               "Nu" = m1$Qlim[4, m1$ncomp.selected],
               "nobj" = 1 # we do not need outliers limit
            ),
            m1$alpha,
            m1$gamma
         )[1],
         chisq.crit(
            list(
               "u0" = m2$Qlim[3, m2$ncomp.selected],
               "Nu" = m2$Qlim[4, m2$ncomp.selected],
               "nobj" = 1 # we do not need outliers limit
            ),
            m1$alpha,
            m1$gamma
         )[1]
      )
   }

   if (length(plot_data) == 1) {
      mdaplotg(plot_data[[1]], groupby = groupby, type = "p", main = main,
         show.lines = show.lines, ...)
   } else {
      mdaplotg(plot_data, type = "p", main = main, show.lines = show.lines, ...)
   }
}

#' Predictions plot for SIMCAM model
#'
#' @description
#' Makes a plot with class predictions for calibration dataset.
#'
#' @param obj
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' vector with class numbers to make the plot for.
#' @param main
#' plot title.
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#'
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#'
#' @export
plotPredictions.simcam <- function(obj, nc = seq_len(obj$nclasses),
   main = "SIMCAM Predictions (cal)", ...) {

   plotPredictions.classmodel(obj, res.name = "cal", nc = nc, main = main, ncomp = 1, ...)
}

#' Model overview plot for SIMCAM
#'
#' @description
#' Shows a set of plots for SIMCAM model.
#'
#' @param x
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' vector with two values - classes (SIMCA models) to show the plot for
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{simcam}} function.
#'
#' @export
plot.simcam <- function(x, nc = c(1, 2), ...) {
   par(mfrow = c(2, 1))
   plotDiscriminationPower(x, nc)
   plotModelDistance(x, nc[1])
   par(mfrow = c(1, 1))
}
