#' Results of PCA decomposition
#' @description
#' \code{pcares} is used to store results for PCA decomposition of data.
#'
#' @param ...
#' other arguments supported by \code{ldecomp}.
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
#' @return
#' Returns an object (list) of class \code{pcares} and \code{ldecomp} with following fields:
#' \item{scores }{matrix with score values (nobj x ncomp).}
#' \item{residuals }{matrix with data residuals (nobj x nvar).}
#' \item{T2 }{matrix with T2 distances (nobj x ncomp).}
#' \item{Q }{matrix with Q residuals (nobj x ncomp).}
#' \item{tnorm }{vector with singular values used for scores normalization.}
#' \item{ncomp.selected }{selected number of components.}
#' \item{expvar }{explained variance for each component.}
#' \item{cumexpvar }{cumulative explained variance.}
#'
#' @seealso
#' Methods for \code{pcares} objects:
#' \tabular{ll}{
#'  \code{print.pcares} \tab shows information about the object.\cr
#'  \code{summary.pcares} \tab shows statistics for the PCA results.\cr
#' }
#'
#' Methods, inherited from \code{\link{ldecomp}} class:
#' \tabular{ll}{
#'  \code{\link{plotScores.ldecomp}} \tab makes scores plot.\cr
#'  \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
#'  \code{\link{plotCumVariance.ldecomp}} \tab makes cumulative explained variance plot.\cr
#'  \code{\link{plotResiduals.ldecomp}} \tab makes Q vs. T2 residuals plot.\cr
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
#' model = pca(people, scale = TRUE, cv = 1, info = 'Simple PCA model')
#' model = selectCompNum(model, 4)
#'
#' res = predict(model, x.new)
#' summary(res)
#' plot(res)
#'
#' ## 1. Make PCA model for People data with autoscaling
#' ## and full cross-validation and get calibration results
#'
#'
#' data(people)
#' model = pca(people, scale = TRUE, cv = 1, info = 'Simple PCA model')
#' model = selectCompNum(model, 4)
#'
#' res = model$calres
#' summary(res)
#' plot(res)
#'
#' ## 2. Show scores plots for the results
#' par(mfrow = c(2, 2))
#' plotScores(res)
#' plotScores(res, cgroup = people[, 'Beer'], show.labels = TRUE)
#' plotScores(res, comp = c(1, 3), show.labels = TRUE)
#' plotScores(res, comp = 2, type = 'h', show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' ## 3. Show residuals and variance plots for the results
#' par(mfrow = c(2, 2))
#' plotVariance(res, type = 'h')
#' plotCumVariance(res, show.labels = TRUE, legend.position = 'bottomright')
#' plotResiduals(res, show.labels = TRUE, cgroup = people[, 'Sex'])
#' plotResiduals(res, ncomp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' @export
pcares <- function(...) {
   # Creates an object of pcares class. In fact the class is a wrapper for ldecomp and
   # uses its methods and attributes.

   res <- ldecomp(...)
   res$categories = NULL
   class(res) = c("pcares", "ldecomp")

   res
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
#' logical, shall values be normalized or not
#' @param log
#' logical, if TRUE, then log(1 + u) transformation is applied
#' @param cgroup
#' color grouping of distance points, if "type" then they will be grouped by type ("normal", "extreme", "outliers")
#' @param show.limits
#' logical, shall extreme and outlier limits be shown or not
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param labels
#' what to show as labels if necessary
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param plot
#' logical, shall plot be created or just plot series object is needed
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotResiduals.pcares <- function(obj, ncomp = obj$ncomp.selected, cgroup = obj$categories,
   norm = TRUE, log = FALSE, show.labels = FALSE, labels = "names", main = NULL,
   show.plot = TRUE, ...) {

   # function for transforming distances
   transform <- function(u, u0, norm, log) {
      if (norm) u <- u/u0
      if (log) u <- log(1 + u)
      return(u)
   }

   # function for creating labels depending on transformation
   get_label <- function(lab, norm, log) {
      if (norm) lab <- paste0(lab, "/", lab, "0")
      if (log) lab <- paste0("log(1 + ", lab, ")")
      return(lab)
   }

   attrs <- mda.getattr(obj$Q)
   h0 <- attr(obj$T2, "u0")
   q0 <- attr(obj$Qlim, "u0")

   # check that scaling values exist
   if (norm && (is.null(h0) || is.null(q0))) {
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
   plot_data <- mda.setattr(plot_data, attrs, "row")
   rownames(plot_data) <- rownames(obj$Q)
   colnames(plot_data) <- c(
      paste0("Orthogonal distance, ", lxlab),
      paste0("Score distance, ", lylab)
   )

   # if no plot required - return plot series object
   if (!show.plot) {
      return(plotseries(plot_data, type = "p", labels = labels, cgroup = cgroup))
   }

   # set up main title for the plot
   if (is.null(main)) {
      main <- if (is.null(ncomp)) "Residuals" else sprintf("Residuals (ncomp = %d)", ncomp)
   }

   # show plot
   mdaplot(plot_data, main = main, cgroup = cgroup, show.labels = show.labels, labels = labels,
      cgroup = cgroup, ...)
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
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#'
#' @export
plot.pcares = function(x, comp = c(1, 2), show.labels = T, ...) {
   par(mfrow = c(2, 2))
   plotScores(x, comp = comp, show.labels = show.labels, ...)
   plotResiduals(x, show.labels = show.labels, ...)
   plotVariance(x, show.labels = show.labels, ...)
   plotCumVariance(x, show.labels = show.labels, ...)
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
   print.ldecomp(x, 'Results for PCA decomposition (class pcares)', ...)
   cat('\n')
}