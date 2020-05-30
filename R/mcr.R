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
   attrs <- mda.getattr(x)
   St <- object$resspec
   Ct <- x %*% St %*% solve(crossprod(St))
   f <- as.matrix(rowSums(x))
   a <- solve(crossprod(Ct)) %*% t(Ct) %*% f
   A <- if (length(a) == 1) a else diag(as.vector(a))
   Ct <- Ct %*% A

   Ct <- mda.setattr(Ct, attrs, "rows")
   return(Ct)
}


########################
#  Static methods      #
########################

#' Compute explained variance for MCR case
#'
#' @param obj
#' object of class \code{mcr}
#' @param x
#' original spectral data
#'
#' @export
getVariance.mcr <- function(obj, x) {

   cumresvar <- rep(0, obj$ncomp)
   for (i in seq_len(obj$ncomp)) {
      cumresvar[i] <- sum((x - tcrossprod(
         obj$rescont[, seq_len(i), drop = FALSE],
         obj$resspec[, seq_len(i), drop = FALSE]
      ))^2)
   }

   cumexpvar <- 100 - cumresvar / sum(x^2) * 100
   expvar <- c(cumexpvar[1], diff(cumexpvar))

   names(cumexpvar) <- names(expvar) <- colnames(obj$resspec)
   attr(expvar, "name") <- "Variance"
   attr(cumexpvar, "name") <- "Cumulative variance"
   attr(expvar, "xaxis.name") <- attr(cumexpvar, "xaxis.name") <- "Components"
   attr(expvar, "yaxis.name") <- attr(cumexpvar, "yaxis.name") <- "Explained variance, %"
   return(list(cumexpvar = cumexpvar, expvar = expvar))
}

########################
#  Plotting methods    #
########################

#' Show plot with resolved spectra
#'
#' @param obj
#' object of clacc \code{mcr}
#' @param comp
#' vector with number of components to make the plot for
#' @param type
#' type of the plot
#' @param col
#' vector with colors for individual components
#'
#' @export
plotSpectra.mcr <- function(obj, comp = seq_len(obj$ncomp), type = "l",
   col = mdaplot.getColors(obj$ncomp), ...) {
   stopifnot("Parameter 'comp' has wrong value." = min(comp) > 0 && max(comp) <= obj$ncomp)

   mdaplotg(mda.subset(mda.t(obj$resspec), comp), type = type, col = col[comp], ...)
}

#' Show plot with resolved contributions
#'
#' @param obj
#' object of clacc \code{mcr}
#' @param comp
#' vector with number of components to make the plot for
#' @param type
#' type of the plot
#' @param col
#' vector with colors for individual components
#'
#' @export
plotContributions.mcr <- function(obj, comp = seq_len(obj$ncomp), type = "l",
   col = mdaplot.getColors(obj$ncomp), ...) {
   stopifnot("Parameter 'comp' has wrong value." = min(comp) > 0 && max(comp) <= obj$ncomp)

   mdaplotg(mda.subset(mda.t(obj$rescont), comp), type = type, ...)
}

#' Show plot with explained variance
#'
#' @param obj
#' object of clacc \code{mcr}
#' @param type
#' type of the plot
#' @param labels
#' what to use as data labels
#' @param xticks
#' vector with ticks for x-axis
#'
#' @export
plotVariance.mcr <- function(obj, type = "h", labels = "values",
   xticks = seq_len(obj$ncomp), ...) {

   mdaplot(obj$expvar, type = type, labels = labels, ...)
}

#' Show plot with cumulative explained variance
#'
#' @param obj
#' object of clacc \code{mcr}
#' @param type
#' type of the plot
#' @param labels
#' what to use as data labels
#' @param xticks
#' vector with ticks for x-axis
#'
#' @export
plotCumVariance.mcr <- function(obj, type = "b", labels = "values",
   xticks = seq_len(obj$ncomp), ...) {

   mdaplot(obj$cumexpvar, type = type, labels = labels, ...)
}

#' Plot summary for MCR model
#'
#' @param obj
#' \code{mcr} model object
#'
#' @export
plot.mcr <- function(obj) {
}
