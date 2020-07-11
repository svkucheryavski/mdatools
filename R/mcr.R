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

   variance <- rbind(expvar, cumexpvar)
   colnames(variance) <- colnames(obj$resspec)
   rownames(variance) <- c("Variance", "Cumulative variance")
   attr(variance, "name") <- "Explained variance"
   attr(variance, "xaxis.name") <- "Components"
   attr(variance, "yaxis.name") <- "Explained variance, %"

   return(variance)
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
#' @param ...
#' other parameters suitable for \code{mdaplotg}
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
#' @param ...
#' other parameters suitable for \code{mdaplotg}
#'
#' @export
plotContributions.mcr <- function(obj, comp = seq_len(obj$ncomp), type = "l",
   col = mdaplot.getColors(obj$ncomp), ...) {
   stopifnot("Parameter 'comp' has wrong value." = min(comp) > 0 && max(comp) <= obj$ncomp)

   mdaplotg(mda.subset(mda.t(obj$rescont), comp), type = type, col = col[comp], ...)
}

#' Show plot with explained variance
#'
#' @param obj
#' object of clacc \code{mcr}
#' @param type
#' type of the plot
#' @param labels
#' what to use as data labels
#' @param main
#' title of the plot
#' @param xticks
#' vector with ticks for x-axis
#' @param ...
#' other parameters suitable for \code{mdaplot}
#'
#' @export
plotVariance.mcr <- function(obj, type = "h", labels = "values", main = "Variance",
   xticks = seq_len(obj$ncomp), ...) {

   mdaplot(mda.subset(obj$variance, 1), type = type, labels = labels, xticks = xticks,
      main = main, ...)
}

#' Show plot with cumulative explained variance
#'
#' @param obj
#' object of clacc \code{mcr}
#' @param type
#' type of the plot
#' @param labels
#' what to use as data labels
#' @param main
#' title of the plot
#' @param xticks
#' vector with ticks for x-axis
#' @param ...
#' other parameters suitable for \code{mdaplot}
#'
#' @export
plotCumVariance.mcr <- function(obj, type = "b", labels = "values", main = "Cumulative variance",
   xticks = seq_len(obj$ncomp), ...) {

   mdaplot(mda.subset(obj$variance, 2), type = type, labels = labels, xticks = xticks,
      main = main, ...)
}

#' Plot summary for MCR model
#'
#' @param x
#' \code{mcr} model object
#' @param ...
#' other parameters
#'
#' @export
plot.mcr <- function(x, ...) {
   par(mfrow = c(2, 2))
   plotSpectra(x)
   plotContributions(x)
   plotVariance(x)
   plotCumVariance(x)
}
