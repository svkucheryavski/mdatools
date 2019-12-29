#' Check reference class values and convert it to a factor if necessary
#'
#' @param c.ref
#' class reference values provided by user
#' @param classname
#' text with class name in case of logical reference values
#'
#' @export
classmodel.processReferenceValues <- function(c.ref, classname = NULL) {

   if (is.null(c.ref) || is.factor(c.ref)) return(c.ref)
   attrs <- mda.getattr(c.ref)

   if (is.logical(c.ref)) {
      return(mda.setattr(as.factor(ifelse(c.ref, model$classname, "None")), attrs))
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
#' @param res
#' which result to make the plot for (\code{"cal"}, \code{"cv"} or \code{"test"}).
#' @param nc
#' vector with class numbers to make the plot for.
#' @param ncomp
#' what number of components to make the plot for.
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#'
#' @details
#' See examples in description of \code{\link{plsda}}, \code{\link{simca}} or \code{\link{simcam}}.
#'
#' @export
plotPredictions.classmodel <- function(obj, res = NULL, nc = 1:obj$nclasses,
   ncomp = obj$ncomp.selected, ...) {

   if (!is.null(res)) {
      plotPredictions.cassres(res, nc = nc, ncomp = ncomp, ...)
      return()
   }

   for (resname in c("test", "cv", "cal")) {
      res <- obj$res[[resname]]
      if (is.null(res)) continue
      plotPredictions.classres(res, nc = nc, ncomp = ncomp, ...)
      return()
   }

   stop("Wong value for 'res' parameter.")
}

#' Specificity plot for classification model
#'
#' @description
#' Makes a plot with specificity values vs. model complexity (e.g. number of components)
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
plotSpecificity.classmodel <- function(obj, ...) {
   plotPerformance(obj, param = "specificity", ...)
}

#' Sensitivity plot for classification model
#'
#' @description
#' Makes a plot with sensitivity values vs. model complexity (e.g. number of components)
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
plotSensitivity.classmodel <- function(obj, ...) {
   plotPerformance(obj, param = "sensitivity", ...)
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
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ylim
#' vector with two values - limits for y axis
#' @param main
#' main title for the plot
#' @param xticks
#' vector with tick values for x-axis
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} function can be used.
#'
#' @export
plotPerformance.classmodel <- function(obj, nc = 1, param = "specificity", type = "h",
   xlab = "Components", ylab = "", ylim = c(0, 1.15),
   main = sprintf("%s for %s (ncomp = %d)", param, obj$classnames[[nc]], ncomp),
   xticks = seq_len(dim(obj$c.pred)[2]), ...) {

   res_names <- names(obj$res)
   plot_data <- list()
   for (i in seq_along(obj$res)) {
      plot_data[[res_names[i]]] <- plotPerformance.classres(obj$res[[i]], nc = nc, type = type,
         param = param, show.plot = FALSE)
   }

   mdaplotg(plot_data, type = type, main = main, xticks = xticks, xlab = xlab, ylab = ylab,
      ylim = ylim, ...)
}

