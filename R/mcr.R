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
   Ct <- x %*% St %*% solve(crossprod(St))
   f <- as.matrix(rowSums(x))
   a <- solve(crossprod(Ct)) %*% t(Ct) %*% f
   A <- if (length(a) == 1) a else diag(as.vector(a))
   Ct <- Ct %*% A

   Ct <- mda.setattr(Ct, attrs, "rows")
   return(Ct)
}


summary.mcr <- function(obj) {
}

print.mcr <- function(obj) {
}


########################
#  Plotting methods    #
########################

plotSpectra.mcr <- function(obj, comp = seq_len(obj$ncomp), type = "l",
   col = mdaplot.getColors(obj$ncomp), ...) {
   stopifnot("Parameter 'comp' has wrong value." = min(comp) > 0 && max(comp) <= obj$ncomp)

   mdaplotg(mda.subset(obj$purespec, comp), type = type, col = col[comp], ...)
}

plotContributions.mcr <- function(obj, comp = seq_len(obj$ncomp), type = "l",
   col = mdaplot.getColors(obj$ncomp), ...) {
   stopifnot("Parameter 'comp' has wrong value." = min(comp) > 0 && max(comp) <= obj$ncomp)

   mdaplotg(mda.subset(mda.t(obj$pureconc), comp), type = type, ...)
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
