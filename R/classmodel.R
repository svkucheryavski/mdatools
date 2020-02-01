#' Check reference class values and convert it to a factor if necessary
#'
#' @param c.ref
#' class reference values provided by user
#' @param classnames
#' text with class name in case of logical reference values
#'
#' @export
classmodel.processRefValues <- function(c.ref, classnames = NULL) {

   if (is.null(c.ref) || is.factor(c.ref)) return(c.ref)
   attrs <- mda.getattr(c.ref)

   if (is.logical(c.ref)) {
      if (length(classnames) > 1) stop("Logical class values can't be used with multiclass model.")
      if (length(classnames) == 0) stop("Classname must be specified with logical class values.")
      return(mda.setattr(as.factor(ifelse(c.ref, classnames, "None")), attrs))
   }

   if (is.character(c.ref)) {
      return(mda.setattr(as.factor(c.ref), attrs))
   }

   stop("Parameter c.ref should be either a factor or vector with logical or text values.")
}

#' Predictions plot for classification model
#'
#' @description
#' Makes a plot with class predictions for a classification model.
#'
#' @param obj
#' a classification model (object of class \code{simca}, \code{plsda}, etc.). if \code{NULL} value
#' is specified, the result will be selected automatically by checking the nearest available from
#' test, cv and calibration results.
#' @param res.name
#' name of result object to make the plot for ("test", "cv" or "cal").
#' @param nc
#' vector with class numbers to make the plot for.
#' @param ncomp
#' what number of components to make the plot for.
#' @param main
#' title of the plot (if NULL will be set automatically)
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#'
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#'
#' @export
plotPredictions.classmodel <- function(obj, res.name = NULL, nc = seq_len(obj$nclasses),
   ncomp = NULL, main = NULL, ...) {

   if (is.null(res.name)) {
      res_names <- c("test", "cv", "cal")
      name_ind <- which(res_names %in% names(obj$res[!sapply(obj$res, is.null)]))[1]
      res.name <- res_names[name_ind]
   }

   res <- obj$res[[res.name]]
   if (is.null(res)) {
      stop("Wong value for 'res.name' parameter.")
   }

   if (is.null(ncomp)) {
      ncomp <- res$ncomp.selected
   }

   if (is.null(main)) main <- sprintf("Predictions (%s, ncomp = %d)", res.name, ncomp)
   plotPredictions.classres(res, nc = nc, ncomp = ncomp, main = main, ...)
}

#' Specificity plot for classification model
#'
#' @description
#' Makes a plot with specificity values vs. model complexity (e.g. number of components)
#'
#' @param obj
#' classification model (object of class \code{plsda}, \code{simca}, etc.).
#' @param legend.position
#' position of the legend (as in \code{mdaplotg}).
#' @param ...
#' parameters for \code{\link{plotPerformance.classmodel}} function.
#'
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#'
#' @export
plotSpecificity.classmodel <- function(obj, legend.position = "bottomright", ...) {
   plotPerformance(obj, param = "specificity", legend.position = legend.position, ...)
}

#' Sensitivity plot for classification model
#'
#' @description
#' Makes a plot with sensitivity values vs. model complexity (e.g. number of components)
#'
#' @param obj
#' classification model (object of class \code{plsda}, \code{simca}, etc.).
#' @param legend.position
#' position of the legend (as in \code{mdaplotg}).
#' @param ...
#' parameters for \code{\link{plotPerformance.classmodel}} function.
#'
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#'
#' @export
plotSensitivity.classmodel <- function(obj, legend.position = "bottomright", ...) {
   plotPerformance(obj, param = "sensitivity", legend.position = legend.position, ...)
}

#' Misclassified ratio plot for classification model
#'
#' @description
#' Makes a plot with misclassified ratio values vs. model complexity (e.g. number of components)
#'
#' @param obj
#' classification model (object of class \code{plsda}, \code{simca}, etc.).
#' @param ...
#' parameters for \code{\link{plotPerformance.classmodel}} function.
#'
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#'
#' @export
plotMisclassified.classmodel <- function(obj, ...) {
   plotPerformance(obj, param = "misclassified", ...)
}

#' Performance plot for classification model
#'
#' @description
#' Makes a plot with sensitivity values vs. model complexity (e.g. number of components)
#'
#' @param obj
#' classification model (object of class \code{plsda}, \code{simca}, etc.).
#' @param nc
#' class number to make the plot for.
#' @param param
#' which parameter to make the plot for (\code{"specificity"}, \code{"sensitivity"},
#' or \code{"misclassified"})
#' @param type
#' type of the plot
#' @param labels
#' what to show as labels for plot objects.
#' @param ylab
#' label for y axis
#' @param ylim
#' vector with two values - limits for y axis
#' @param xticks
#' vector with tick values for x-axis
#' @param res
#' list with result objects to show the plot for
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#'
#' @export
plotPerformance.classmodel <- function(obj, nc = 1, param = "misclassified", type = "b",
   labels = "values", ylab = "", ylim = c(0, 1.15),
   xticks = seq_len(dim(obj$res$cal$c.pred)[2]), res = obj$res, ...) {

   if (length(param) != 1) {
      stop("Specify which paramete you want to make the plot for.")
   }

   plot_data <- lapply(res, plotPerformance, nc = nc, param = param, show.plot = FALSE)
   mdaplotg(plot_data, type = type, xticks = xticks, ylab = ylab, ylim = ylim, labels = labels, ...)
}
