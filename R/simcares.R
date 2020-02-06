#' Results of SIMCA one-class classification
#'
#'  @description
#' \code{simcares} is used to store results for SIMCA one-class classification.
#' @param class.res
#' results of classification (class \code{classres}).
#' @param pca.res
#' results of PCA decomposition of data (class \code{pcares}).
#'
#' @details
#' Class \code{simcares} inherits all properties and methods of class \code{\link{pcares}}, and
#' has additional properties and functions for representing of classification results, inherited
#' from class \code{\link{classres}}.
#'
#' There is no need to create a \code{simcares} object manually, it is created automatically when
#' build a SIMCA model (see \code{\link{simca}}) or apply the model to a new data (see
#' \code{\link{predict.simca}}). The object can be used to show summary and plots for the results.
#'
#' @return
#' Returns an object (list) of class \code{simcares} with the same fields as \code{\link{pcares}}
#' plus extra fields, inherited from \code{\link{classres}}:
#' \item{c.pred}{predicted class values (+1 or -1).}
#' \item{c.ref}{reference (true) class values if provided.}
#'
#' The following fields are available only if reference values were provided.
#' \item{tp}{number of true positives.}
#' \item{fp}{nmber of false positives.}
#' \item{fn}{number of false negatives.}
#' \item{specificity}{specificity of predictions.}
#' \item{sensitivity}{sensitivity of predictions.}
#'
#' @seealso
#' Methods for \code{simcares} objects:
#' \tabular{ll}{
#'  \code{print.simcares} \tab shows information about the object.\cr
#'  \code{summary.simcares} \tab shows statistics for results of classification.\cr
#' }
#'
#' Methods, inherited from \code{\link{classres}} class:
#' \tabular{ll}{
#'  \code{\link{showPredictions.classres}} \tab show table with predicted values.\cr
#'  \code{\link{plotPredictions.classres}} \tab predicted classes plot.\cr
#'  \code{\link{plotSensitivity.classres}} \tab sensitivity plot.\cr
#'  \code{\link{plotSpecificity.classres}} \tab specificity plot.\cr
#'  \code{\link{plotPerformance.classres}} \tab performance plot.\cr
#' }
#'
#' Methods, inherited from \code{\link{ldecomp}} class:
#' \tabular{ll}{
#'  \code{\link{plotResiduals.ldecomp}} \tab makes Q2 vs. T2 residuals plot.\cr
#'  \code{\link{plotScores.ldecomp}} \tab makes scores plot.\cr
#'  \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
#'  \code{\link{plotCumVariance.ldecomp}} \tab makes cumulative explained variance plot.\cr
#' }
#' Check also \code{\link{simca}} and \code{\link{pcares}}.
#'
#' @examples
#' ## make a SIMCA model for Iris setosa class and show results for calibration set
#' library(mdatools)
#'
#' data = iris[, 1:4]
#' class = iris[, 5]
#'
#' # take first 30 objects of setosa as calibration set
#' se = data[1:30, ]
#'
#' # make SIMCA model and apply to test set
#' model = simca(se, 'Se')
#' model = selectCompNum(model, 1)
#'
#' # show infromation and summary
#' print(model$calres)
#' summary(model$calres)
#'
#' # show plots
#' layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
#' plotPredictions(model$calres, show.labels = TRUE)
#' plotResiduals(model$calres, show.labels = TRUE)
#' plotPerformance(model$calres, show.labels = TRUE, legend.position = 'bottomright')
#' layout(1, 1, 1)
#'
#' # show predictions table
#' showPredictions(model$calres)
#' @export
simcares <- function(class.res, pca.res = NULL) {
   res <- c(pca.res, class.res)
   class(res) <- c("simcares", "classres", if (!is.null(pca.res)) c("pcares", "ldecomp"))
   return(res)
}


#' as.matrix method for SIMCA classification results
#'
#' @description
#' Generic \code{as.matrix} function for classification results. Returns matrix with performance
#' values for specific class.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param ncomp
#' model complexity (number of components) to show the parameters for.
#' @param ...
#' other arguments
#'
#' @export
as.matrix.simcares <- function(x, ncomp = NULL, ...) {

   out <- cbind(
      (if (!is.null(x$scores)) as.matrix.ldecomp(x) else matrix(NA, ncol = 2, nrow = x$ncomp)),
      as.matrix.classres(x)
   )

   if (!is.null(ncomp)) {
      out <- out[ncomp, , drop = FALSE]
   }

   out[, 1:2] <- round(out[, 1:2], 2)

   colnames(out) <- c("Expvar", "Cumexpvar", "TP", "FP", "TN", "FN",
      "Spec.", "Sens.", "Accuracy")

   return(out)
}


#' Summary method for SIMCA results object
#'
#' @description
#' Shows performance statistics for the results.
#'
#' @param object
#' SIMCA results (object of class \code{simcares})
#' @param ...
#' other arguments
#'
#' @export
summary.simcares <- function(object, ...) {

   cat("\nSummary for SIMCA one-class classification result\n")
   cat(sprintf("\nClass name: %s\n", object$classname))
   cat(sprintf("Number of selected components: %d\n", object$ncomp.selected))
   cat("\n")

   print(as.matrix.simcares(object, ...))
   cat("\n")
}

#' Print method for SIMCA results object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' SIMCA results (object of class \code{simcares})
#' @param ...
#' other arguments
#'
#' @export
print.simcares <- function(x, ...) {

   cat("Result for SIMCA one-class classification (class simcares)\n")
   cat(sprintf("Method for critical limits: %s\n", x$lim.type))
   print.ldecomp(x, "")
   print.classres(x, "")
   cat("\n")
}
