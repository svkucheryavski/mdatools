#' Results of PCA decomposition
#'
#' @description
#' \code{pcares} is used to store and visualise results for PCA decomposition.
#'
#' @param ...
#' all arguments supported by \code{ldecomp}.
#'
#' @details
#' In fact \code{pcares} is a wrapper for \code{\link{ldecomp}} - general class for storing
#' results for linear decomposition X = TP' + E. So, most of the methods, arguments and
#' returned values are inherited from \code{ldecomp}.
#'
#' There is no need to create a \code{pcares} object manually, it is created automatically when
#' build a PCA model (see \code{\link{pca}}) or apply the model to a new data (see
#' \code{\link{predict.pca}}). The object can be used to show summary and plots for the results.
#'
#' It is assumed that data is a matrix or data frame with I rows and J columns.
#'
#' @return
#' Returns an object (list) of class \code{pcares} and \code{ldecomp} with following fields:
#' \item{scores }{matrix with score values (I x A).}
#' \item{residuals }{matrix with data residuals (I x J).}
#' \item{T2 }{matrix with score distances (I x A).}
#' \item{Q }{matrix with orthogonal distances (I x A).}
#' \item{ncomp.selected }{selected number of components.}
#' \item{expvar }{explained variance for each component.}
#' \item{cumexpvar }{cumulative explained variance.}
#'
#' @seealso
#' Methods for \code{pcares} objects:
#' \tabular{ll}{
#'  \code{print.pcares} \tab shows information about the object.\cr
#'  \code{summary.pcares} \tab shows statistics for the PCA results.\cr
#'  \code{\link{plotResiduals.pcares}} \tab makes Q vs. T2 distance plot.\cr
#' }
#'
#' Methods, inherited from \code{\link{ldecomp}} class:
#' \tabular{ll}{
#'  \code{\link{plotScores.ldecomp}} \tab makes scores plot.\cr
#'  \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
#'  \code{\link{plotCumVariance.ldecomp}} \tab makes cumulative explained variance plot.\cr
#' }
#'
#' Check also \code{\link{pca}} and \code{\link{ldecomp}}.
#'
#' @examples
#' ### Examples for PCA results class
#'
#' library(mdatools)
#'
#' ## 1. Make a model for every odd row of People data
#' ## and apply it to the objects from every even row
#'
#' data(people)
#' x = people[seq(1, 32, 2), ]
#' x.new = people[seq(1, 32, 2), ]
#'
#' model = pca(people, scale = TRUE, info = "Simple PCA model")
#' model = selectCompNum(model, 4)
#'
#' res = predict(model, x.new)
#' summary(res)
#' plot(res)
#'
#' ## 1. Make PCA model for People data with autoscaling
#' ## and full cross-validation and get calibration results
#'
#' data(people)
#' model = pca(people, scale = TRUE, info = "Simple PCA model")
#' model = selectCompNum(model, 4)
#'
#' res = model$calres
#' summary(res)
#' plot(res)
#'
#' ## 2. Show scores plots for the results
#' par(mfrow = c(2, 2))
#' plotScores(res)
#' plotScores(res, cgroup = people[, "Beer"], show.labels = TRUE)
#' plotScores(res, comp = c(1, 3), show.labels = TRUE)
#' plotScores(res, comp = 2, type = "h", show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' ## 3. Show residuals and variance plots for the results
#' par(mfrow = c(2, 2))
#' plotVariance(res, type = "h")
#' plotCumVariance(res, show.labels = TRUE)
#' plotResiduals(res, show.labels = TRUE, cgroup = people[, "Sex"])
#' plotResiduals(res, ncomp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' @export
pcares <- function(...) {
   # Creates an object of pcares class. In fact the class is a wrapper for ldecomp and
   # uses its methods and attributes.

   res <- ldecomp(...)
   class(res) <- c("pcares", "ldecomp")

   return(res)
}

#' Residual distance plot
#'
#' @description
#' Shows a plot with orthogonal (Q, q) vs score (T2, h) distances for data objects. By default the
#' distance values are normalize using corresponding means (q/q0 and h/h0).
#'
#' @param obj
#' object of \code{ldecomp} class.
#' @param ncomp
#' number of components to show the plot for (if NULL, selected by model value will be used).
#' @param norm
#' logical, if TRUE disance values be normalized (u/u0)
#' @param log
#' logical, if TRUE, then log(1 + u) transformation is applied
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param labels
#' what to show as labels if necessary
#' @param main
#' main title for the plot
#' @param show.plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotResiduals.pcares <- function(obj, ncomp = obj$ncomp.selected,
   norm = TRUE, log = FALSE, show.labels = FALSE, labels = "names", main = NULL,
   show.plot = TRUE, ...) {

   # function for transforming distances
   transform <- function(u, u0, norm, log) {
      if (norm) u <- u / u0
      if (log) u <- log(1 + u)
      return(u)
   }

   # function for creating labels depending on transformation
   get_label <- function(lab, norm, log) {
      if (norm) lab <- paste0(lab, "/", lab, "0")
      if (log) lab <- paste0("log(1 + ", lab, ")")
      return(lab)
   }

   # get scale factors
   h0 <- attr(obj$T2, "u0")[[ncomp]]
   q0 <- attr(obj$Q, "u0")[[ncomp]]

   # check that scaling values exist
   if (norm && (is.null(h0) || is.null(q0))) {
      warning("Can not normalize distances as scaling values are absent.")
      norm <- FALSE
   }

   # prepare plot data
   h <- transform(obj$T2[, ncomp], h0, norm, log)
   q <- transform(obj$Q[, ncomp], q0, norm, log)

   # default values for local labels
   lxlab <- get_label("h", norm, log)
   lylab <- get_label("q", norm, log)

   # combine everything to dataset and assign attributes
   plot_data <- mda.cbind(h, q)
   plot_data <- mda.setattr(plot_data, mda.getattr(obj$Q), "row")
   rownames(plot_data) <- rownames(obj$Q)
   colnames(plot_data) <- c(
      paste0("Score distance, ", lxlab),
      paste0("Orthogonal distance, ", lylab)
   )

   # if no plot required - return plot series object
   if (!show.plot) {
      return(plot_data)
   }

   # set up main title for the plot
   if (is.null(main)) {
      main <- if (is.null(ncomp)) "Residuals" else sprintf("Residuals (ncomp = %d)", ncomp)
   }

   # show plot
   return(mdaplot(plot_data, main = main, show.labels = show.labels, labels = labels, ...))
}

#' Plot method for PCA results object
#'
#' @description
#' Show several plots to give an overview about the PCA results
#'
#' @param x
#' PCA results (object of class \code{pcares})
#' @param comp
#' which components to show the scores plot for (can be one value or vector with two values).
#' @param ncomp
#' how many components to use for showing the residual distance plot
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#'
#' @export
plot.pcares <- function(x, comp = c(1, 2), ncomp = x$ncomp.selected, show.labels = T, ...) {
   par(mfrow = c(2, 2))
   plotScores(x, comp = comp, show.labels = show.labels, ...)
   plotResiduals(x, ncomp, show.labels = show.labels, ...)
   plotVariance(x, type = "h", show.labels = show.labels, ...)
   plotCumVariance(x, type = "h", show.labels = show.labels, ...)
   par(mfrow = c(1, 1))
}


#' Summary method for PCA results object
#'
#' @description
#' Shows some statistics (explained variance, eigenvalues) about the results.
#'
#' @param object
#' PCA results (object of class \code{pcares})
#' @param ...
#' other arguments
#'
#' @export
summary.pcares <- function(object, ...) {
   summary.ldecomp(object, "Summary for PCA results", ...)
}


#' Print method for PCA results object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' PCA results (object of class \code{pcares})
#' @param ...
#' other arguments
#'
#' @export
print.pcares <- function(x, ...) {
   print.ldecomp(x, "Results for PCA decomposition (class pcares)", ...)
   cat("\n")
}
