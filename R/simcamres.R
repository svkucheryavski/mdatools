#' Results of SIMCA multiclass classification
#'
#' @description
#' \code{simcamres} is used to store results for SIMCA multiclass classification.
#'
#' @param cres
#' results of classification (class \code{classres}).
#' @param pred.res
#' list with prediction results from each model (pcares objects)
#'
#' @details
#' Class \code{simcamres} inherits all properties and methods of class \code{\link{classres}}, plus
#' store values necessary to visualise prediction decisions (e.g. Cooman's plot or Residuals plot).
#'
#' In cotrast to \code{simcares} here only values for optimal (selected) number of components in
#' each individual SIMCA models are presented.
#'
#' There is no need to create a \code{simcamres} object manually, it is created automatically when
#' make a SIMCAM model (see \code{\link{simcam}}) or apply the model to a new data (see
#' \code{\link{predict.simcam}}). The object can be used to show summary and plots for the results.
#'
#' @return
#' Returns an object (list) of class \code{simcamres} with the same fields as \code{\link{classres}}
#' plus extra fields for Q and T2 values and limits:
#'
#' \item{c.pred}{predicted class values.}
#' \item{c.ref}{reference (true) class values if provided.}
#' \item{T2}{matrix with T2 values for each object and class.}
#' \item{Q}{matrix with Q values for each object and class.}
#' \item{T2lim}{vector with T2 statistical limits for each class.}
#' \item{Qlim}{vector with Q statistical limits for each class.}
#'
#' The following fields are available only if reference values were provided.
#' \item{tp}{number of true positives.}
#' \item{fp}{nmber of false positives.}
#' \item{fn}{number of false negatives.}
#' \item{specificity}{specificity of predictions.}
#' \item{sensitivity}{sensitivity of predictions.}
#'
#' @seealso
#' Methods for \code{simcamres} objects:
#' \tabular{ll}{
#'  \code{print.simcamres} \tab shows information about the object.\cr
#'  \code{summary.simcamres} \tab shows statistics for results of classification.\cr
#'  \code{\link{plotCooman.simcamres}} \tab makes Cooman's plot.\cr
#' }
#'
#' Methods, inherited from \code{\link{classres}} class:
#' \tabular{ll}{
#'  \code{\link{showPredictions.classres}} \tab show table with predicted values.\cr
#'  \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
#' }
#'
#' Check also \code{\link{simcam}}.
#'
#' @examples
#' ## see examples for simcam method.
#'
#' @export
simcamres <- function(cres, pred.res) {
   res <- cres
   res$pred.res <- pred.res
   res$classnames <- dimnames(cres$c.pred)[[3]]
   class(res) <- c("simcamres", "classres")

   return(res)
}

#' as.matrix method for SIMCAM results
#'
#' @description
#' Generic \code{as.matrix} function for SIMCAM results. Returns matrix with performance
#' values for specific class.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' vector with classes to use.
#' @param ...
#' other arguments
#'
#' @export
as.matrix.simcamres <- function(x, nc = seq_len(x$nclasses), ...) {
   comp <- sapply(x$pred.res, function(r) r$ncomp.selected)

   out <- do.call(rbind, lapply(nc, function(n) as.matrix.classres(x, nc = n)))
   out <- cbind(comp, out)

   rownames(out) <- x$classnames
   colnames(out)[1] <- "Ncomp"
   return(out)
}

#' Summary method for SIMCAM results object
#'
#' @description
#' Shows performance statistics for the results.
#'
#' @param object
#' SIMCAM results (object of class \code{simcamres})
#' @param nc
#' number of class to show summary for (can be vector)
#' @param ...
#' other arguments
#'
#' @export
summary.simcamres <- function(object, nc = seq_len(object$nclasses), ...) {

   cat("\nSummary for SIMCA multiple classes classification result\n")
   if (is.null(object$c.ref)) {
      cat("No reference data provided to calculate prediction performance.")
      return()
   }

   fprintf("\nNumber of classes: %d\n", length(nc))
   print(as.matrix.simcamres(object))
   cat("\n")
}

#' Print method for SIMCAM results object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' SIMCAM results (object of class \code{simcamres})
#' @param ...
#' other arguments
#'
#' @export
print.simcamres <- function(x, ...) {
   cat("\nResult for SIMCA multiple classes classification (class simcamres)\n\n")
   print.classres(x, "")
   cat("\n")
}


################################
#  Plotting methods            #
################################


#' Cooman's plot for SIMCAM results
#'
#' @description
#' Shows a Cooman's plot for a pair of SIMCA models
#'
#' @param obj
#' SIMCAM results (object of class \code{simcamres})
#' @param nc
#' vector with two values - classes (SIMCA models) to show the plot for
#' @param main
#' main plot title
#' @param cgroup
#' vector of values to use for color grouping of plot points
#' @param show.plot
#' logical, show plot or just return plot data
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' The plot is similar to \code{\link{plotCooman.simcam}} but shows points only for this result
#' object and does not show critical limits (which are part of a model).
#'
#' @export
plotCooman.simcamres <- function(obj, nc = c(1, 2), main = "Cooman's plot",
   cgroup = obj$c.ref, show.plot = TRUE, ...) {

   attrs <- mda.getattr(obj$c.pred)
   res1 <- obj$pred.res[[nc[1]]]
   res2 <- obj$pred.res[[nc[2]]]

   plot_data <- cbind(
      res1$Q[, res1$ncomp.selected],
      res2$Q[, res2$ncomp.selected]
   )

   attr(plot_data, "exclrows") <- attrs$exclrows
   attr(plot_data, "name") <- "Cooman's plot"
   rownames(plot_data) <- rownames(obj$c.pred)
   colnames(plot_data) <- c(
      paste0("Distance to class ", obj$classnames[nc[1]]),
      paste0("Distance to class ", obj$classnames[nc[2]])
   )

   if (!show.plot) {
      return(plot_data)
   }

   mdaplot(plot_data, type = "p", cgroup = cgroup, ...)
}

#' Prediction plot for SIMCAM results
#'
#' @description
#' Makes a plot with predicted class values for classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' vector with classes to show predictions for.
#' @param main
#' title of the plot
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} or \code{\link{mdaplot}} function
#' can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotPredictions.simcamres <- function(obj, nc = seq_len(obj$nclasses), main = "Predictions", ...) {
   return(plotPredictions.classres(obj, nc = nc, ncomp = 1, main = main, ...))
}

#' Model overview plot for SIMCAM results
#'
#' @description
#' Just shows a prediction plot for SIMCAM results.
#'
#' @param x
#' SIMCAM results (object of class \code{simcamres})
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{simcamres}} function.
#'
#' @export
plot.simcamres <- function(x, ...) {
   plotPredictions(x)
}
