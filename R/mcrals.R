#' Multivariate curve resolution using Alternating Least Squares
#'
#' @description
#' \code{mcralls} allows to resolve spectroscopic data to linear combination of individual spectra
#' and contributions using the alternating least squares (ALS) algorithm with constraints.
#'
#' @param x
#' spectra of mixtures (matrix or data frame).
#' @param ncomp
#' number of components to calculate.
#' @param cont.constraints
#' a list with constraints to be applied to contributions  (see details).
#' @param spec.constraints
#' a list with constraints to be applied to spectra  (see details).
#' @param max.niter
#' maximum number of iterations
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values).
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values).
#' @param info
#' a short text with description of the case (optional).
#'
#' @details
#' The method implements the iterative ALS algorithm, where, at each iteration, spectra and
#' contributions of each chemical component are estimated and then a set of constraints is
#' applied to each. The method is well described in [1, 2].
#'
#' The constraints for contributions and for spectra should be provided as a list with name of the
#' constraint and all necessary parameters. You can see which constraints and parameters are
#' currently supported by running \code{mcrals.constlist()}. See the code examples below or a
#' Bookdown tutorial for more details.
#'
#'
#' @return
#' Returns an object of \code{\link{mcrpure}} class with the following fields:
#' \item{resspec}{matrix with resolved spectra.}
#' \item{rescont}{matrix with resolved contributions.}
#' \item{cont.constraints}{list with contribution constraints provided by user.}
#' \item{spec.constraints}{list with spectra constraints provided by user.}
#' \item{expvar }{vector with explained variance for each component (in percent).}
#' \item{cumexpvar }{vector with cumulative explained variance for each component (in percent).}
#' \item{ncomp}{number of resolved components}
#' \item{max.niter}{maximum number of iterations}
#' \item{info }{information about the model, provided by user when build the model.}
#'
#'
#' More details and examples can be found in the Bookdown tutorial.
#'
#' @references
#' 1.
#'
#' 2.
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @seealso
#' Methods for \code{mcrals} objects:
#' \tabular{ll}{
#'    \code{summary.mcrals} \tab shows some statistics for the case.\cr
#'    \code{\link{unmix.mcrals}} \tab makes unmixing of new set of spectra.\cr
#'    \code{\link{predict.mcr}} \tab computes contributions by projection of new spectra to
#'    the resolved ones.\cr
#' }
#'
#' Plotting methods for \code{mcrals} objects:
#' \tabular{ll}{
#'    \code{\link{plotSpectra.mcr}} \tab shows plot with resolved spectra.\cr
#'    \code{\link{plotContributions.mcr}} \tab shows plot with resolved contributions.\cr
#'    \code{\link{plotVariance.mcr}} \tab shows plot with explained variance.\cr
#'    \code{\link{plotCumVariance.mcr}} \tab shows plot with cumulative explained variance.\cr
#' }
#'
#' @examples
#' library(mdatools)
#'
#' # resolve mixture of carbonhydrates Raman spectra
#'
#' data(carbs)
#' m = mcrpure(carbs$D, ncomp = 3)
#'
#' # plot cumulative and individual explained variance
#' par(mfrow = c(1, 2))
#' plotVariance(m)
#' plotCumVariance(m)
#'
#' # plot resolved spectra (all of them or individually)
#' par(mfrow = c(2, 1))
#' plotSpectra(m)
#' plotSpectra(m, comp = 2:3)
#'
#' # plot resolved contributions (all of them or individually)
#' par(mfrow = c(2, 1))
#' plotContributions(m)
#' plotContributions(m, comp = 2:3)
#'
#' # of course you can do this manually as well, e.g. show original
#' # and resolved spectra
#' par(mfrow = c(1, 1))
#' mdaplotg(
#'    list(
#'       "original" = prep.norm(carbs$D, "area"),
#'       "resolved" = prep.norm(mda.subset(mda.t(m$resspec), 1), "area")
#'    ), col = c("gray", "red"), type = "l"
#' )
#'
#' # in case if you have reference spectra of components you can compare them with
#' # the resolved ones:
#' par(mfrow = c(3, 1))
#' for (i in 1:3) {
#'    mdaplotg(
#'       list(
#'          "pure" = prep.norm(mda.subset(mda.t(carbs$S), 1), "area"),
#'          "resolved" = prep.norm(mda.subset(mda.t(m$resspec), 1), "area")
#'       ), col = c("gray", "red"), type = "l", lwd = c(3, 1)
#'    )
#' }
#'
#' # See bookdown tutorial for more details.
#'
#' @export
mcrals <- function(x, ncomp, cont.constraints = list(), spec.constraints = list(), max.niter = 100,
   exclrows = NULL, exclcols = NULL, verbose = FALSE, info = "") {

   stopifnot("Parameter 'max.niter' should be positive." = max.niter > 0)
   stopifnot("Parameter 'cont.constraints' has incorrect values" =
      check.constraint(cont.constraints))
   stopifnot("Parameter 'spec.constraints' has incorrect values" =
      check.constraint(spec.constraints))


   # get pure variables and unmix data
   x <- prepCalData(x, exclrows, exclcols, min.nrows = 2, min.ncols = 2)
   model <- mcrals.cal(x, ncomp, cont.constraints, spec.constraints, max.niter, verbose)

   # compute explained variance
   model <- c(model, getVariance.mcr(model, x))
   class(model) <- c("mcr", "mcrals")
   model$info <- info

   return(model)
}

#' Print method for mcrpure object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' \code{mcrpure} object
#' @param ...
#' other arguments
#'
#' @export
print.mcrals <- function(x, ...) {
   cat("\nMCRALS unmixing case (class mcrals)\n")

   if (length(x$info) > 1) {
      cat("\nInfo:\n")
      cat(x$info)
   }

   cat("\n\nCall:\n")
   print(x$call)

   cat("\nMajor model fields:\n")
   cat("$ncomp - number of calculated components\n")
   cat("$resspec - matrix with resolved spectra\n")
   cat("$rescons - matrix with resolved concentrations\n")
   cat("$expvar - vector with explained variance\n")
   cat("$cumexpvar - vector with cumulative explained variance\n")
   cat("$info - case info provided by user\n")
}

#' Summary method for mcrals object
#'
#' @description
#' Shows some statistics (explained variance, etc) for the case.
#'
#' @param object
#' \code{mcrals} object
#' @param ...
#' other arguments
#'
#' @export
summary.mcrals <- function(object, ...) {
   cat("\nSummary for MCR ALS case (class mcrals)\n")

   if (length(object$info) > 1) {
      fprintf("\nInfo:\n%s\n", object$info)
   }

   if (length(object$exclrows) > 0) {
      fprintf("Excluded rows: %d\n", length(object$exclrows))
   }

   if (length(object$exclcols) > 0) {
      fprintf("Excluded coumns: %d\n", length(object$exclcols))
   }

   #fprintf("\nOffset: %s\n", object$offset)
   cat("\n")

   data <- cbind(
      round(object$expvar, 2),
      round(object$cumexpvar, 2),
   )
   colnames(data) <- c("Expvar", "Cumexpvar")
   rownames(data) <- colnames(object$resspec)
   show(data)
}


################################
#  Static methods              #
################################


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
#' @param use.deriv
#' a number which tells how to use derivative.
#' @param savgol
#' list with parameters for Savitzky-Golay preprocessing (if \code{use.deriv} is not 0).
#'
#' @return
#' The function returns a list with with following fields:
#' \item{ncomp }{number of pure components.}
#' \item{purvars }{vector with indices for pure variables.}
#' \item{purityspec }{matrix with purity values for each resolved components.}
#' \item{purity }{vector with purity values for resolved components.}
#'
#' @export
mcrals.cal <- function(D, ncomp, cont.constraints, spec.constraints, max.niter, verbose {

   attrs <- mda.getattr(D)
   exclrows <- attrs$exclrows
   exclcols <- attrs$exclcols

   # get dimensions
   nspec <- nrow(D)
   nvar <- ncol(D)
   nvarvis <- nvar - length(exclcols)

   # get indices for excluded columns if provided
   colind <- seq_len(nvar)

   # remove hidden rows and columns
   if (!is.null(exclrows)) D <- D[-exclrows, , drop = FALSE]

   if (!is.null(exclcols)) {
      D <- D[, -exclcols, drop = FALSE]
      colind <- colind[-exclcols]
   }

   # initialize S with random numbers
   S <- matrix(runif(nvarvis * ncomp), nvarvis, ncomp)

   # main loop for ALS
   for (i in seq_len(max.niter)) {
      ## normalize S
      S <- S %*% diag(1 / sqrt(colSums(S^2)), nvarvis, nvarvis)

      ## compute C and apply constraints
      C <- D %*% (S %*% solve(crossprod(S)))
      for (cc in cont.constraints) {
         C <- apply(cc, C)
      }

      ## compute S and apply constraints
      S <- t(D) %*% (C %*% solve(crossprod(C)))
      for (sc in spec.constraints) {
         S <- apply(sc, S)
      }
   }

   # if there were excluded rows or columns, handle this
   Ct <- matrix(0, nobj, ncomp)
   Ct[rowind, ] <- C

   St <- matrix(0, ncomp, nvar)
   St[, colind] <- S

   # add attributes
   Ct <- mda.setattr(Ct, attrs, "row")
   St <- mda.setattr(St, attrs, "col")
   colnames(Ct) <- rownames(St) <- paste("Comp", seq_len(ncomp))

   if (is.null(attr(Ct, "xaxis.name"))) attr(Ct, "xaxis.name") <- "Observations"
   attr(Ct, "name") <- "Resolved contributions"
   attr(St, "name") <- "Resolved spectra"

   # combine the results with model object
   return(list(rescont = Ct, resspec = mda.t(St)))
   return(
      list(
         ncomp = ncomp,
         rescont = rescont,
         resspec = resspec,
         cont.constraints = cont.constraints,
         spec.constraints = spec.constraints,
         max.niter = max.niter,
      )
   )
}


########################
#  Plotting methods    #
########################


