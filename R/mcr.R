#' General class for Multivariate Curve Resolution model
#'
#' @description
#' \code{mcr} is used to store and visualise general MCR data and results.
#'
#' @param x
#' spectra of mixtures (as matrix or data frame)
#' @param ncomp
#' number of pure components to resolve
#' @param method
#' function for computing spectra of pure components
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values)
#' @param info
#' text with information about the MCR model
#' @param ...
#' other parameters realted to specific method
#'
#' @export
mcr <- function(x, ncomp, method, exclrows = NULL, exclcols = NULL, info = "", ...) {


   model <- list()
   model$ncomp <- ncomp
   model$info <- info
   model$call <- match.call()

   class(model) <- c("mcr")
   model$spectra <- method(model, x, ...)
   model$contr <- predict(model, x)

   return(model)
}

#' MCR predictions
#'
#' @description
#' Applies MCR model to a new set of spectra and returns matrix with contributions.
#'
#' @param object
#' an MCR model (object of class \code{mcr}).
#' @param x
#' spectral values (matrix or data frame).
#' @param ...
#' other arguments.
#'
#' @return
#' Matrix with contributions
#'
#' @export
predict.mcr <- function(object, x, ...) {

   # convert to matrix
   x <- mda.df2mat(x)
   attrs <- attributes(x)

   if (is.null(dim(x))) {
      stop("Test set should be a matrix or a data frame.")
   }

   if (ncol(x) != nrow(object$loadings)) {
      stop("Number and type of data columns should be the same as in calibration dataset.")
   }

   # compute scores and residuals
   x <- prep.autoscale(x, center = object$center, scale = object$scale)
   scores <- x %*% object$loadings

   # set names
   rownames(scores) <- rownames(residuals) <- attrs$dimnames[[1]]
   colnames(scores) <- colnames(object$loadings)

   # set attributes
   scores <- mda.setattr(scores, attrs, "row")
   residuals <- mda.setattr(residuals, attrs)
   attr(scores, "name") <- "Scores"
   attr(scores, "xaxis.name") <- "Components"
   attr(residuals, "name") <- "Residuals"

   if (is.null(attrs$yaxis.name)) {
      attr(scores, "yaxis.name") <- "Objects"
   }

   return(contributions)
}


summary.mcr <- function(obj) {

}

print.mcr <- function(obj) {

}


########################
#  Plotting methods    #
########################

plotSpectra.mcr <- function(obj, comp, ...) {

}

plotContributions.mcr <- function(obj, comp, ...) {

}

plotVariance.mcr <- function(obj, ...) {

}

#' Plot summary for MCR model
#'
#' @param obj
#' \code{mcr} model object
#'
#' @export
plot.mcr <- function(obj) {

}
