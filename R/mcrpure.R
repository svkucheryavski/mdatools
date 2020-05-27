#' Multivariate curve resolution with purity approach
#'
#' @description
#' \code{purity} is used to resolve spectroscopic data into pure spectra and contributions using
#'  purity concept.
#'
#' @param x
#' spectra of mixtures (matrix or data frame).
#' @param ncomp
#' maximum number of components to calculate.
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values)
#' @param savgol
#' list with parameters for Savitzky-Golay preprocessing
#' @param info
#' a short text with model description.
#'
#' @details
#'
#' @return
#' Returns an object of \code{\link{mcrdata}} class.
#'
#' More details and examples can be found in the Bookdown tutorial.
#'
#' @references
#' 1. .
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @seealso
#' Methods for \code{pca} objects:
#' \tabular{ll}{
#'    \code{plot.purity} \tab makes an overview of PCA model with four plots.\cr
#'    \code{summary.purity} \tab shows some statistics for the model.\cr
#'    \code{\link{predict.mcrpure}} \tab applies PCA model to a new data.\cr
#' }
#'
#' Plotting methods for \code{pca} objects:
#' \tabular{ll}{
#'    \code{\link{plotScores.pca}} \tab shows scores plot.\cr
#'    \code{\link{plotLoadings.pca}} \tab shows loadings plot.\cr
#'    \code{\link{plotVariance.pca}} \tab shows explained variance plot.\cr
#'    \code{\link{plotCumVariance.pca}} \tab shows cumulative explained variance plot.\cr
#'    \code{\link{plotResiduals.pca}} \tab shows plot for residual distances (Q vs. T2).\cr
#'    \code{\link{plotBiplot.pca}} \tab shows bi-plot.\cr
#'    \code{\link{plotExtreme.pca}} \tab shows extreme plot.\cr
#'    \code{\link{plotT2DoF}} \tab plot with degrees of freedom for score distance.\cr
#'    \code{\link{plotQDoF}} \tab plot with degrees of freedom for orthogonal distance.\cr
#'    \code{\link{plotDistDoF}} \tab plot with degrees of freedom for both distances.\cr
#' }
#'
#' Most of the methods for plotting data are also available for PCA results (\code{\link{pcares}})
#' objects. Also check \code{\link{pca.mvreplace}}, which replaces missing values in a data matrix
#' with approximated using iterative PCA decomposition.
#'
#' @examples
#' library(mdatools)
#'
#' @export
mcrpure <- function(x, ncomp, purevars = NULL, offset = 0.05, savgol = NULL, exclrows = NULL,
   exclcols = NULL, info = "") {

   stopifnot("Parameter 'offset' should be between 0 and 1." = offset < 1 && offset > 0)
   stopifnot("Provided pure variables have wrong values." =
      is.null(purevars) || (min(purevars) > 0 & max(purevars) <= ncol(x)))


   # exclude columns if "exclcols" is provided
   if (length(exclcols) > 0) {
      x <- mda.exclcols(x, exclcols)
   }

   # exclude rows if "exclrows" is provided
   if (length(exclrows) > 0) {
      x <- mda.exclrows(x, exclrows)
   }

   # get pure variables and unmix data
   model <- getPureVariables(x, ncomp, purevars = purevars, offset = offset)
   #res2 = unmix(spectra, D[, res1$purevars, drop = F], by.spec = by.spec)

   class(model) <- c("mcr", "mcrpure")
   return(model)
}

#' Identifies pure variables
#'
#' @description
#' The method identifies indices of pure variables using the SIMPLISMA
#' algorithm.
#'
#' @param D
#' matrix with the spectra
#' @param ncomp
#' number of pure components
#' @param purevars
#' user provided values gor pure variables (no calculation will be run in this case)
#' @param offset
#' offset (between 0 and 1) for calculation of parameter alpha
#' @param exclcols
#' optinal, indices of columns to be excluded from calculations
#' @param exclrows
#' optinal, indices of columns to be excluded from calculations
#'
#' @return
#' The function returns a list with with following fields:
#' \item{ncomp }{number of pure components.}
#' \item{purvars }{vector with indices for pure variables.}
#' \item{purityspec }{matrix with purity values for each resolved components.}
#' \item{stdspec }{matrix with weighted standard deviation values.}
#' \item{purity }{vector with purity values for resolved components.}
#'
#' @export
getPureVariables <- function(D, ncomp, purevars, offset) {

   attrs <- mda.getattr(D)
   purevals <- NULL
   purityspec <- NULL

   # if indices for pure variables are not provided - compute them
   if (is.null(purevars)) {
      exclrows <- attrs$exclrows
      exclcols <- attrs$exclcols
      D <- prepCalData(D, exclrows, exclcols)

      # get dimensions
      nspec <- nrow(D)
      nvar <- ncol(D)

      # get indices for excluded columns if provided
      colind <- rep(TRUE, nvar)
      colind[exclcols] <- FALSE

      # calculate purity spectrum
      mu <- apply(D, 2, mean)
      sigma <- apply(D, 2, sd) * sqrt((nspec - 1) / nspec)
      alpha <- offset * max(mu)
      purity <- sigma / (mu + alpha)

      # calculate variance-covariance matrix
      Dprime <- sweep(D, 2, (sqrt(mu^2 + (sigma + alpha)^2)), '/')
      R <- crossprod(Dprime) / nrow(Dprime)

      # loop to locate the pure variables
      purityspec <- matrix(0, nrow = ncomp, ncol = nvar)
      purevars <- rep(0, ncomp)
      purevals <- rep(0, ncomp)

      for (i in seq_len(ncomp)) {
         weights <- sapply(
            seq_len(nvar),
            function(j) det(
               R[c(j, purevars[seq_len(i-1)]), c(j, purevars[seq_len(i-1)]), drop = FALSE]
            )
         )

         purityspec[i, colind] <- purity[colind] * weights
         purevars[i] <- which.max(purityspec[i, ])
         purevals[i] <- purityspec[i, purevars[i]]
      }

      dim(purevals) <- c(1, ncomp)
      colnames(purevals) <- rownames(purityspec) <- paste("Comp", seq_len(ncomp))
      purityspec <- mda.setattr(purityspec, attrs, "col")

      attr(purevals, "name") <- "Purity"
      attr(purevals, "xaxis.name") <- "Components"
      attr(purityspec, "name") <- "Purity spectra"
   }

   attr(purevars, "xaxis.name") <- attrs$xaxis.name
   attr(purevars, "xaxis.values") <- purevars
   if (!is.null(attrs$xaxis.values)) {
      attr(purevars, "xaxis.values") <- attrs$xaxis.values[purevars]
   }

   return(
      list(
         ncomp = ncomp,
         purity = purevals,
         purityspec = purityspec,
         purevars = purevars
      )
   )
}

#' Purity values plot
#'
#' @param obj
#' \code{mcrpure} object
#' @param xticks
#' ticks for x axis
#' @param type
#' type of the plot
#' @param ...
#' other parameters suitable for \code{mdaplot}
#'
#' The plot shows largest weighted purity value for each component graphically.
#'
#' @export
plotPurity.mcrpure <- function(obj, xticks = seq_len(obj$ncomp), type = "h",
   labels = "values", ...) {

   if (is.null(obj$purity)) {
      warning("Purity values are not computed when pure variables provided by user.")
      return()
   }

   mdaplot(obj$purity, xticks = xticks, type = type, labels = labels, ...)
}

#' Purity spectra plot
#'
#' @param obj
#' \code{mcrpure} object
#' @param comp
#' vector of components to show the purity spectra for
#' @param type
#' type of the plot
#' @param col
#' colors for the plot (should be a vector with one value for each component in \code{obj})
#' @param show.lines
#' if \code{TRUE} show the selected pure variables as vertical lines
#' @param lines.col
#' color for the selected pure variable lines (by default same as for plots but semitransparent)
#' @param lines.lty
#' line type for the purity lines
#' @param lines.lwd
#' line width for the purity lines
#' @param ...
#' other parameters suitable for \code{mdaplotg}
#'
#' The plot shows weighted purity value of each variable separately for each specified component.
#'
#' @export
plotPuritySpectra.mcrpure <- function(obj, comp = seq_len(obj$ncomp), type = "l",
   col = mdaplot.getColors(obj$ncomp), show.lines = TRUE,
   lines.col = adjustcolor(col, alpha.f = 0.5), lines.lty = 3, lines.lwd = 1, ...) {

   if (is.null(obj$purityspec)) {
      warning("Purity values are not computed when pure variables provided by user.")
      return()
   }

   stopifnot("Parameter 'comp' has wrong value." = min(comp) > 0 && max(comp) <= obj$ncomp)
   mdaplotg(mda.subset(obj$purityspec, comp), type = type, col = col[comp], ...)

   if (show.lines) {
      abline(
         v = attr(obj$purevars, "xaxis.values")[comp],
         col = lines.col[comp],
         lty = lines.lty,
         lwd = lines.lwd
      )
   }
}

#' Print method for PCA model object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' a PCA model (object of class \code{pca})
#' @param ...
#' other arguments
#'
#' @export
print.mcrpure <- function(x, ...) {
}

#' Summary method for PCA model object
#'
#' @description
#' Shows some statistics (explained variance, eigenvalues) for the model.
#'
#' @param object
#' a PCA model (object of class \code{pca})
#' @param ...
#' other arguments
#'
#' @export
summary.mcrpure <- function(object, ...) {
}


