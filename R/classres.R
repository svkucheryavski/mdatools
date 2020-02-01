#' Results of classification
#'
#' @description
#' \code{classres} is used to store results classification for one or multiple classes.
#'
#' @param c.pred
#' matrix with predicted values (+1 or -1) for each class.
#' @param c.ref
#' matrix with reference values for each class.
#' @param p.pred
#' matrix with probability values for each class.
#' @param ncomp.selected
#' vector with selected number of components for each class.
#'
#' @details
#' There is no need to create a \code{classres} object manually, it is created automatically when
#' build a classification model (e.g. using \code{\link{simca}} or \code{\link{plsda}}) or apply
#' the model to new data. For any classification method from \code{mdatools}, a class using to
#' represent results of classification (e.g. \code{\link{simcares}}) inherits fields and methods of
#' \code{classres}.
#'
#' @return
#' \item{c.pred}{predicted class values (+1 or -1).}
#' \item{p.pred}{predicted class probabilities.}
#' \item{c.ref}{reference (true) class values if provided.}
#'
#' The following fields are available only if reference values were provided.
#' \item{tp}{number of true positives.}
#' \item{tn}{number of true negatives.}
#' \item{fp}{nmber of false positives.}
#' \item{fn}{number of false negatives.}
#' \item{specificity}{specificity of predictions.}
#' \item{sensitivity}{sensitivity of predictions.}
#' \item{misclassified}{ratio of misclassified objects.}
#'
#' @seealso
#' Methods \code{classres} class:
#' \tabular{ll}{
#'  \code{\link{showPredictions.classres}} \tab shows table with predicted values.\cr
#'  \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
#'  \code{\link{plotSensitivity.classres}} \tab makes sn plot.\cr
#'  \code{\link{plotSpecificity.classres}} \tab makes specificity plot.\cr
#'  \code{\link{plotMisclassified.classres}} \tab makes ms ratio plot.\cr
#'  \code{\link{plotPerformance.classres}} \tab makes plot with misclassified ratio, specificity
#'  and sensitivity values.\cr
#' }
#'
#' @export
classres <- function(c.pred, c.ref = NULL, p.pred = NULL, ncomp.selected = 1) {

   if (length(dim(c.pred)) != 3) {
      stop("Wrong number of dimensions for 'c.pred' array (should be 3-way array).")
   }

   obj <- list()
   obj$c.ref <- if (!is.null(c.ref)) as.factor(c.ref)
   obj$c.pred <- c.pred
   obj$p.pred <- p.pred
   obj$nclasses <- dim(c.pred)[3]
   obj$ncomp <- dim(c.pred)[2]
   obj$classnames <- dimnames(c.pred)[[3]]

   # check that ncomp.selected is correct
   if (is.null(ncomp.selected)) ncomp.selected <- obj$ncomp
   if (ncomp.selected < 1 || ncomp.selected > obj$ncomp) {
      stop("Wrong value for 'ncomp.selected' parameer.")
   }

   obj$ncomp.selected <- ncomp.selected

   if (!is.null(c.ref)) {
      obj <- c(obj, classres.getPerformance(c.ref, c.pred))
   }


   obj$call <- match.call()
   class(obj) <- "classres"

   return(obj)
}

#' Confusion matrix for classification results
#'
#' @details
#' Returns confusion matrix for classification results represented by the object.
#'
#' @param obj
#' classification results (object of class \code{simcares}, \code{simcamres}, etc)
#' @param ncomp
#' number of components to make the matrix for (NULL - use selected for a model).
#' @param ...
#' other arguments
#'
#' @description
#' The columns of the matrix correspond to classification results, rows - to the real classes. In
#' case of soft classification with multiple classes (e.g. SIMCAM) sum of values for every row
#' will not correspond to the total number of class members as the same object can be classified
#' as a member of several classes or non of them.
#'
#' @export
getConfusionMatrix.classres <- function(obj, ncomp = obj$ncomp.selected, ...) {

   if (is.null(obj$c.ref)) {
      stop("Reference classes are not available!")
   }

   attrs <- mda.getattr(obj$c.pred)
   c.pred <- obj$c.pred[, ncomp, , drop = FALSE]
   c.ref <- obj$c.ref

   # remove excluded rows
   if (length(attrs$exclrows) > 0) {
      c.pred <- c.pred[-attrs$exclrows, , , drop = FALSE]
      c.ref <- c.ref[-attrs$exclrows]
   }

   # get class names and numbers
   ref.classes <- levels(c.ref)
   ref.nclasses <- length(ref.classes)

   # compute the confusion matrix
   out <- matrix(0, nrow = ref.nclasses, ncol = obj$nclasses + 1)
   none <- rep(TRUE, length(c.ref))
   for (i in seq_len(obj$nclasses)) {
      ind <- c.pred[, , i] > 0
      out[, i] <- table(c.ref[ind], exclude = FALSE)
      none[ind] <- FALSE
   }

   # find ones that were not classified as member of any class and names
   out[, obj$nclasses + 1] <- table(c.ref[none], exclude = FALSE)
   rownames(out) <- ref.classes
   colnames(out) <- c(obj$classnames, "None")

   # reorder the table to match class name order in results
   ind1 <- match(colnames(out), rownames(out))
   ind1 <- ind1[!is.na(ind1)]
   ind2 <- seq_len(nrow(out))[-ind1]
   out <- out[c(ind1, ind2), , drop = FALSE]

   return(out)
}

#' Show predicted class values
#'
#' @description
#' Shows a table with predicted class values for classification result.
#'
#' @param obj
#' object with classification results (e.g. \code{plsdares} or \code{simcamres}).
#' @param ncomp
#' number of components to show the predictions for (NULL - use selected for a model).
#' @param ...
#' other parameters
#'
#' @details
#' The function prints a matrix where every column is a class and every row is an data object.
#' The matrix has either -1 (does not belong to the class) or +1 (belongs to the class) values.
#'
#' @export
showPredictions.classres <- function(obj, ncomp = obj$ncomp.selected, ...) {

   pred <- obj$c.pred[, ncomp, ]
   dim(pred) <- dim(obj$c.pred)[c(1, 3)]
   dimnames(pred) <- dimnames(obj$c.pred)[c(1, 3)]

   print(pred)
   cat("\n")
}

#' as.matrix method for classification results
#'
#' @description
#' Generic \code{as.matrix} function for classification results. Returns matrix with performance
#' values for specific class.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param ncomp
#' model complexity (number of components) to show the parameters for.
#' @param nc
#' if there are several classes, which class to show the parameters for.
#' @param ...
#' other arguments
#'
#' @export
as.matrix.classres <- function(x, ncomp = NULL, nc = 1, ...) {

   if (is.null(x$c.ref)) return()

   if (length(nc) != 1) {
      stop("Wrong value for 'nc' parameter.")
   }

   specificity <- if (is.null(x$specificity)) matrix(NA, x$nclasses, x$ncomp) else x$specificity
   out <- cbind(
      x$tp[nc, ], x$fp[nc, ], x$tn[nc, ], x$fn[nc, ],
      round(specificity[nc, ], 3),
      round(x$sensitivity[nc, ], 3),
      round(1 - x$misclassified[nc, ], 3)
   )

   out[is.nan(out)] <- NA
   colnames(out) <- c("TP", "FP", "TN", "FN", "Spec.", "Sens.", "Accuracy")
   rownames(out) <- dimnames(x$c.pred)[[2]]

   if (!is.null(ncomp)) {
      out <- out[ncomp, , drop = FALSE]
   }

   return(out)
}

#' Print information about classification result object
#'
#' @description
#' Generic \code{print} function for classification results. Prints information about major fields
#' of the object.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param str
#' User specified text (e.g. to be used for particular method, like PLS-DA, etc).
#' @param ...
#' other arguments
#'
#' @export
print.classres <- function(x, str = "Classification results (class classres)\nMajor fields:", ...) {

   if (nchar(str) > 0) fprintf("\n%s\n", str)

   cat("$c.pred - predicted class values\n")

   if (!is.null(x$c.ref)) {
      cat("$c.ref - reference (true) class values\n")
      cat("$tp - number of true positives\n")
      cat("$tn - number of true negatives\n")
      cat("$fp - number of false positives\n")
      cat("$fn - number of false negatives\n")
      cat("$specificity - specificity of predictions\n")
      cat("$sensitivity - sn of predictions\n")
      cat("$misclassified - misclassification ratio for predictions\n")
   }
}

#' Summary statistics about classification result object
#'
#' @description
#' Generic \code{summary} function for classification results. Prints performance values for the
#' results.
#'
#' @param object
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param ncomp
#' which number of components to make the plot for (use NULL to show results for all available).
#' @param nc
#' vector with class numbers to show the summary for.
#' @param ...
#' other arguments
#'
#' @export
summary.classres <- function(object, ncomp = object$ncomp.selected,
   nc = seq_len(object$nclasses), ...) {

   cat("\nClassification results (class classres) summary\n")

   if (is.null(object$c.ref)) {
      cat("No reference data provided to calculate prediction performance.")
      return()
   }

   fprintf("\nNumber of selected components: %d", ncomp)
   fprintf("\nNumber of classes: %d\n", object$nclasses)

   # detailed results for several components
   for (i in nc) {
      fprintf("\nClass '%s':\n", object$classnames[[i]])
      print(as.matrix.classres(object, nc = i, ncomp = ncomp))
   }

}

################################
#  Static methods              #
################################

#' Calculation of  classification performance parameters
#'
#' @description
#' Calculates and returns performance parameters for classification result (e.g. number of false
#' negatives, false positives, sn, specificity, etc.).
#'
#' @param c.ref
#' reference class values for objects (vector with numeric or text values)
#' @param c.pred
#' predicted class values for objects (array nobj x ncomponents x nclasses)
#'
#' @return
#' Returns a list with following fields:
#' \tabular{ll}{
#'    \code{$fn} \tab number of false negatives (nclasses x ncomponents) \cr
#'    \code{$fp} \tab number of false positives (nclasses x ncomponents) \cr
#'    \code{$tp} \tab number of true positives (nclasses x ncomponents) \cr
#'    \code{$sensitivity} \tab sn values (nclasses x ncomponents) \cr
#'    \code{$specificity} \tab specificity values (nclasses x ncomponents) \cr
#'    \code{$specificity} \tab ms ratio values (nclasses x ncomponents) \cr
#' }
#'
#' @details
#' The function is called automatically when a classification result with reference values is
#' created, for example when applying a \code{plsda} or \code{simca} models.
#'
classres.getPerformance <- function(c.ref, c.pred) {

   if (is.null(c.ref) || is.null(c.pred)) {
      stop("Both reference and predicted class values are required.")
   }

   if (length(c.ref) != dim(c.pred)[1]) {
      stop("Number of objects in reference and predicted results should be the same.")
   }

   # remove excluded rows for correct calculation of performance
   dim(c.ref) <- NULL
   attrs <- mda.getattr(c.pred)
   if (length(attrs$exclrows) > 0) {
      c.pred <- c.pred[-attrs$exclrows, , , drop = F]
      c.ref <- c.ref[-attrs$exclrows]
   }

   ncomp <- dim(c.pred)[2]
   nclasses <- dim(c.pred)[3]

   tp <- matrix(0, nrow = nclasses, ncol = ncomp)
   fp <- matrix(0, nrow = nclasses, ncol = ncomp)
   fn <- matrix(0, nrow = nclasses, ncol = ncomp)
   tn <- matrix(0, nrow = nclasses, ncol = ncomp)

   # compute main performance indicators
   classnames <- dimnames(c.pred)[[3]]
   for (i in seq_len(nclasses)) {
      fn[i, ] <- colSums((c.ref == classnames[i]) & (c.pred[, , i, drop = F] == -1))
      fp[i, ] <- colSums((c.ref != classnames[i]) & (c.pred[, , i, drop = F] == 1))
      tp[i, ] <- colSums((c.ref == classnames[i]) & (c.pred[, , i, drop = F] == 1))
      tn[i, ] <- colSums((c.ref != classnames[i]) & (c.pred[, , i, drop = F] == -1))
   }

   # compute main statistics
   sn <- tp / (tp + fn)
   sp <- tn / (tn + fp)
   ms <- (fp + fn) / (tp + tn + fp + fn)

   # add row with summary for all classes
   sn <- rbind(sn, colSums(tp) / colSums(tp + fn))
   sp <- rbind(sp, colSums(tn) / colSums(tn + fp))
   ms <- rbind(ms, colSums(fp + fn) / colSums(tp + tn + fp + fn))

   # add names
   row_names <- dimnames(c.pred)[[3]]
   col_names <- dimnames(c.pred)[[2]]
   rownames(fn) <- rownames(fp) <- rownames(tp) <- rownames(tn) <- row_names
   colnames(fn) <- colnames(fp) <- colnames(tp) <- colnames(sn) <- colnames(sp) <- col_names
   rownames(sn) <- rownames(sp) <- rownames(ms) <- c(row_names, "Total")

   # in case of one class classifier set sensitivity NULL
   if (all(is.na(sp))) sp <- NULL

   return(
      list(
         "fn" = fn,
         "fp" = fp,
         "tp" = tp,
         "tn" = tn,
         "sensitivity" = sn,
         "specificity" = sp,
         "misclassified" = ms
      )
   )
}


################################
#  Plotting methods            #
################################

#' Plot for class belonging probability
#'
#' @description
#' Makes a plot with class belonging probabilities for each object of the classification results.
#' Works only with classification methods, which compute this probability (e.g. SIMCA).
#'
#' @param obj
#' classification results (e.g. object of class \code{simcamres}).
#' @param ncomp
#' number of components to use the probabilities for.
#' @param nc
#' if there are several classes, which class to make the plot for.
#' @param type
#' type of the plot
#' @param ylim
#' vector with limits for y-axis
#' @param show.lines
#' shows a horizontal line at p = 0.5
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotProbabilities.classres <- function(obj, ncomp = obj$ncomp.selected, nc = 1, type = "h",
   ylim = c(0, 1.1), show.lines = c(NA, 0.5), ...) {

   if (is.null(obj$p.pred)) {
      stop("No probability values are available.")
   }

   if (nc > obj$nclasses || nc < 1) {
      stop("Wrong value for argument 'nc'.")
   }

   plot_data <- obj$p.pred[, ncomp, nc]
   cname <- obj$classnames[[nc]]
   attr(plot_data, "name") <- sprintf("Class probabilities, %s (ncomp = %d)", cname, ncomp)
   attr(plot_data, "yaxis.name") <- "Probability"
   attr(plot_data, "xaxis.name") <- attr(obj$c.pred, "yaxis.name")

   return(mdaplot(plot_data, show.lines = show.lines, type = type, ylim = ylim, ...))
}

#' Sensitivity plot for classification results
#'
#' @description
#' Makes a plot with sn values vs. model complexity (e.g. number of components) for
#' classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param legend.position
#' position of the legend (as in \code{mdaplotg}).
#' @param ...
#' other parameters for \code{\link{plotPerformance.classres}}
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotSensitivity.classres <- function(obj, legend.position = "bottomright", ...) {
   return(plotPerformance(obj, param = "sensitivity", legend.position = legend.position, ...))
}

#' Specificity plot for classification results
#'
#' @description
#' Makes a plot with specificity values vs. model complexity (e.g. number of components) for
#' classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param legend.position
#' position of the legend (as in \code{mdaplotg}).
#' @param ...
#' other parameters for \code{\link{plotPerformance.classres}}
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotSpecificity.classres <- function(obj, legend.position = "bottomright", ...) {
   return(plotPerformance(obj, param = "specificity", legend.position = legend.position, ...))
}

#' Misclassified ratio plot for classification results
#'
#' @description
#' Makes a plot with ms ratio values vs. model complexity (e.g. number of components) for
#' classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param ...
#' other parameters for \code{\link{plotPerformance.classres}}
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotMisclassified.classres <- function(obj, ...) {
   return(plotPerformance(obj, param = "misclassified", ...))
}

#' Performance plot for classification results
#'
#' @description
#' Makes a plot with classification performance parameters vs. model complexity (e.g. number of
#' components) for classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for.
#' @param type
#' type of the plot
#' @param param
#' which performance parameter to make the plot for (can be a vector with several values).
#' @param labels
#' what to show as labels for plot objects.
#' @param ylab
#' label for y axis
#' @param ylim
#' vector with two values - limits for y axis
#' @param xticks
#' vector with x-axis tick values
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotPerformance.classres <- function(obj, nc = 1, type = "b",
   param = c("sensitivity", "specificity", "misclassified"), labels = "values",
   ylab = "", ylim = c(0, 1.1), xticks = seq_len(obj$ncomp), show.plot = TRUE, ...) {

   if (is.null(obj$c.ref)) {
      stop("No reference data available")
   }

   # check if parameters requested are not NULL
   param <- param[param %in% sapply(names(obj), function(x) if (!is.null(obj[[x]])) x)]
   if (length(param) == 0) {
      stop("Performance parameteres you requested are not available in this result object.")
   }

   # prepare plot data
   plot_data <- do.call(rbind, lapply(obj[param], function(x) x[nc, , drop = FALSE]))

   attr(plot_data, "name") <- sprintf(
      if (length(param) == 1) capitalize(param) else "Classification performance (%s)",
      obj$classnames[[nc]]
   )

   attr(plot_data, "xaxis.name") <- "Components"
   rownames(plot_data) <- param

   # if no plot needed return the plat data
   if (!show.plot) {
      return(plot_data)
   }

   mdaplotg(plot_data, type = type, xticks = xticks, ylim = ylim, ylab = ylab, labels = labels, ...)
}

#' Prediction plot for classification results
#'
#' @description
#' Makes a plot with predicted class values for classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' vector with classes to show predictions for.
#' @param ncomp
#' model complexity (number of components) to make the plot for.
#' @param ylab
#' label for y axis
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} or \code{\link{mdaplot}} function
#' can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotPredictions.classres <- function(obj, nc = seq_len(obj$nclasses), ncomp = obj$ncomp.selected,
   ylab = "", show.plot = TRUE, ...) {

   # prepare data and attributes
   attrs <- mda.getattr(obj$c.pred)
   c.pred <- as.matrix(obj$c.pred[, ncomp, nc])
   row_ind <- seq_len(nrow(c.pred))
   class_names <- obj$classnames[nc]
   class_numbers <- seq_along(nc) + 1

   # multiply classes to integers starting from 2 (1 will be for none)
   plot_data <- (c.pred > 0) %*% diag(class_numbers, length(nc), length(nc))

   # fine those which were not classified as members of any class (none)
   plot_data <- cbind(rowSums(plot_data) == 0, plot_data)

   # unfold matrix with class numbers and merge with row indices
   plot_data <- cbind(row_ind, as.numeric(plot_data))

   # remove rows with zeros as class number
   plot_data <- plot_data[plot_data[, 2] > 0, , drop = FALSE]

   # add row names and exclude hidden rows
   if (is.null(attrs$yaxis.name)) attrs$yaxis.name <- "Objects"
   plot_data <- mda.exclrows(plot_data, plot_data[, 1] %in% attrs$exclrows)
   rownames(plot_data) <- rownames(c.pred)[plot_data[, 1]]
   colnames(plot_data) <- c(attrs$yaxis.name, "Classes")
   attr(plot_data, "name") <- sprintf("Predictions (ncomp = %d)", ncomp)

   if (!show.plot) {
      return(plot_data)
   }

   yticks <- c(1, class_numbers)
   yticklabels <- c("None", class_names)

   cgroup <- if (!is.null(obj$c.ref)) as.factor(obj$c.ref[plot_data[, 1]])
   mdaplot(plot_data, type = "p", ylab = ylab, yticks = yticks,
      cgroup = cgroup, yticklabels = yticklabels, ...)
}

#' Plot function for classification results
#'
#' @description
#' Generic plot function for classification results.
#' Alias for \code{\link{plotPredictions.classres}}.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param ...
#' other arguments for \code{plotPredictions()} method.
#'
#' @export
plot.classres <- function(x, ...) {
   plotPredictions.classres(x, ...)
}
