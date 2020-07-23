#' Multivariate curve resolution based on pure variables
#'
#' @description
#' \code{mcrpure} allows to resolve spectroscopic data to linear combination of individual spectra
#' and contributions using the pure variables approach.
#'
#' @param x
#' spectra of mixtures (matrix or data frame).
#' @param ncomp
#' maximum number of components to calculate.
#' @param purevars
#' vector with indices for pure variables (optional, if you want to provide the variables directly).
#' @param offset
#' offset for correcting noise in computing maximum angles (should be value within [0, 1)).
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values).
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values).
#' @param info
#' a short text with description of the case (optional).
#'
#' @details
#' The method estimates purity of each variable and then uses the purest ones to decompose the
#' spectral data into spectra (`resspec`) and contributions (`rescont`) of individual
#' chemical components by ordinary least squares.
#'
#' The pure variabes are identified using stepwise maximum angle calculations and described in
#' detail in [1]. So the purity of a spectral variable (wavelength, wavenumber) is actually an
#' angle (measured in degrees) between the variable and vector of ones for the first component;
#' and between the variable and space formed by previously found pure variables for the other
#' components.
#'
#' @return
#' Returns an object of \code{\link{mcrpure}} class with the following fields:
#' \item{resspec}{matrix with resolved spectra.}
#' \item{rescont}{matrix with resolved contributions.}
#' \item{purevars}{indices of the selected pure variables.}
#' \item{purevals}{purity values for the selected pure variables.}
#' \item{purityspec}{purity spectra (matrix with purity values for each variable and component).}
#' \item{expvar }{vector with explained variance for each component (in percent).}
#' \item{cumexpvar }{vector with cumulative explained variance for each component (in percent).}
#' \item{offset}{offset value used to compute the purity}
#' \item{ncomp}{number of resolved components}
#' \item{info }{information about the model, provided by user when build the model.}
#'
#'
#' More details and examples can be found in the Bookdown tutorial.
#'
#' @references
#' 1. Willem Windig, Neal B. Gallagher, Jeremy M. Shaver, Barry M. Wise. A new approach for
#' interactive self-modeling mixture analysis. Chemometrics and Intelligent Laboratory Systems,
#' 77 (2005) 85–96. DOI: 10.1016/j.chemolab.2004.06.009
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @seealso
#' Methods for \code{mcrpure} objects:
#' \tabular{ll}{
#'    \code{summary.mcrpure} \tab shows some statistics for the case.\cr
#'    \code{\link{unmix.mcrpure}} \tab makes unmixing of new set of spectra.\cr
#'    \code{\link{predict.mcrpure}} \tab computes contributions by projection of new spectra to
#'    the resolved ones.\cr
#' }
#'
#' Plotting methods for \code{mcrpure} objects:
#' \tabular{ll}{
#'    \code{\link{plotPurity.mcrpure}} \tab shows plot with maximum purity of each component.\cr
#'    \code{\link{plotPuritySpectra.mcrpure}} \tab shows plot with purity spectra.\cr
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
#' # examples for purity spectra plot (you can select which components to show)
#' par(mfrow = c(2, 1))
#' plotPuritySpectra(m)
#' plotPuritySpectra(m, comp = 2:3)
#'
#' # you can do it manually and combine e.g. with original spectra
#' par(mfrow = c(1, 1))
#' mdaplotg(
#'    list(
#'       "spectra" = prep.norm(carbs$D, "area"),
#'       "purity" = prep.norm(mda.subset(mda.t(m$resspec), 1), "area")
#'    ), col = c("gray", "red"), type = "l"
#' )
#'
#' # show the maximum purity for each component
#' par(mfrow = c(1, 1))
#' plotPurity(m)
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
mcrpure <- function(x, ncomp, purevars = NULL, offset = 0.05, exclrows = NULL, exclcols = NULL,
   info = "") {

   stopifnot("Parameter 'offset' should be within [0, 1)." = offset < 1 && offset >= 0)
   stopifnot("Provided pure variables have wrong values." =
      is.null(purevars) || (min(purevars) > 0 & max(purevars) <= ncol(x)))

   # get pure variables and unmix data
   x <- prepCalData(x, exclrows, exclcols, min.nrows = 2, min.ncols = 2)
   model <- getPureVariables(x, ncomp, purevars, offset)
   model <- c(model, unmix.mcrpure(model, x))

   # compute explained variance
   model$variance <- getVariance.mcr(model, x)

   # combine everything as an S3 object
   class(model) <- c("mcr", "mcrpure")
   model$call <- match.call()
   model$info <- info

   return(model)
}

#' Unmix spectral data using pure variables estimated before
#'
#' @param obj
#' \code{mcrpure} object
#' @param x
#' matrix with spectra
#'
#' @returns
#' Returns a list with resolved spectra and contributions (matrices).
unmix.mcrpure <- function(obj, x) {
   attrs <- mda.getattr(x)

   # get only values for the purest variables
   Dr <- x[, obj$purevars, drop = FALSE]

   # resolve spectra and concentrations
   St <- tryCatch(
      t(x) %*% Dr %*% solve(crossprod(Dr)),
      error = function(e)
         stop("Unable to resolve the components, perhaps 'ncomp' is too large.", call. = FALSE)
   )

   Ct <- tryCatch(
      x %*% St %*% solve(crossprod(St)),
      error = function(e)
         stop("Unable to resolve the components, perhaps 'ncomp' is too large.", call. = FALSE)
   )

   # scale
   f <- as.matrix(rowSums(x))
   a <- solve(crossprod(Ct)) %*% t(Ct) %*% f
   A <- if (length(a) == 1) a else diag(as.vector(a))
   St <- St %*% A %*% solve(crossprod(A))
   Ct <- Ct %*% A

   # add attributes
   Ct <- mda.setattr(Ct, attrs, "row")
   St <- mda.setattr(t(St), attrs, "col")
   colnames(Ct) <- rownames(St) <- names(obj$purevars)

   if (is.null(attr(Ct, "yaxis.name"))) attr(Ct, "yaxis.name") <- "Observations"
   attr(Ct, "xaxis.name") <- attr(St, "yaxis.name") <- "Arbitrary units"
   attr(Ct, "name") <- "Resolved contributions"
   attr(St, "name") <- "Resolved spectra"

   # combine the results with model object
   return(list(rescont = Ct, resspec = mda.t(St)))
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
predict.mcrpure <- function(object, x, ...) {
   attrs <- mda.getattr(x)
   St <- object$resspec
   Ct <- x %*% St %*% solve(crossprod(St))
   f <- as.matrix(rowSums(x))
   a <- solve(crossprod(Ct)) %*% t(Ct) %*% f
   A <- if (length(a) == 1) a else diag(as.vector(a))
   Ct <- Ct %*% A
   colnames(Ct) <- colnames(object$rescont)
   Ct <- mda.setattr(Ct, attrs, "row")
   return(Ct)
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
print.mcrpure <- function(x, ...) {
   cat("\nMCR Purity unmixing case (class mcrpure)\n")

   if (length(x$info) > 1) {
      cat("\nInfo:\n")
      cat(x$info)
   }

   cat("\n\nCall:\n")
   print(x$call)

   cat("\nMajor model fields:\n")
   cat("$ncomp - number of calculated components\n")
   cat("$resspec - matrix with resolved spectra\n")
   cat("$rescont - matrix with resolved contributions\n")
   cat("$purityspec - matrix with purity spectra\n")
   cat("$purevars - vector with indices of pure variables\n")
   cat("$purevals - purity values for the selected pure variables\n")
   cat("$expvar - vector with explained variance\n")
   cat("$cumexpvar - vector with cumulative explained variance\n")
   cat("$offset - offset value used to compute the purity\n")
   cat("$info - case info provided by user\n")
}

#' Summary method for mcrpure object
#'
#' @description
#' Shows some statistics (explained variance, etc) for the case.
#'
#' @param object
#' \code{mcrpure} object
#' @param ...
#' other arguments
#'
#' @export
summary.mcrpure <- function(object, ...) {
   cat("\nSummary for MCR Purity case (class mcrpure)\n")

   if (length(object$info) > 1) {
      fprintf("\nInfo:\n%s\n", object$info)
   }

   drvstr <- c("no", "only for estimation of pure variables", "for pure variables and unmixing")

   if (length(object$exclrows) > 0) {
      fprintf("Excluded rows: %d\n", length(object$exclrows))
   }

   if (length(object$exclcols) > 0) {
      fprintf("Excluded coumns: %d\n", length(object$exclcols))
   }

   data <- cbind(
      round(t(object$variance), 2),
      object$purevars,
      round(object$purevals, 3)
   )

   colnames(data) <- c("Expvar", "Cumexpvar", "Varindex", "Purity")
   rownames(data) <- colnames(object$purityspec)
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
#'
#' @return
#' The function returns a list with with following fields:
#' \item{ncomp }{number of pure components.}
#' \item{purvars }{vector with indices for pure variables.}
#' \item{purityspec }{matrix with purity values for each resolved components.}
#' \item{purity }{vector with purity values for resolved components.}
#'
#' @export
getPureVariables <- function(D, ncomp, purevars, offset) {

   attrs <- mda.getattr(D)

   # if indices for pure variables are not provided - compute them
   if (is.null(purevars)) purevars <- rep(0, ncomp)

   exclrows <- attrs$exclrows
   exclcols <- attrs$exclcols

   # get dimensions
   nspec <- nrow(D)
   nvar <- ncol(D)
   nspecvis <- nspec - length(exclrows)
   nvarvis <- nvar - length(exclcols)

   # get indices for excluded columns if provided
   colind <- seq_len(nvar)

   # remove hidden rows and columns
   if (!is.null(exclrows)) {
      D <- D[-exclrows, , drop = FALSE]
   }

   if (!is.null(exclcols)) {
      D <- D[, -exclcols, drop = FALSE]
      colind <- colind[-exclcols]
   }

   # calculate statistics and noise correction vector
   mu <- apply(D, 2, mean)
   n <- mu / (mu + offset * max(mu))
   n <- (1 + offset) * n

   # prepare vectores for resultss
   purityspec <- matrix(0, nrow = ncomp, ncol = nvar)
   purevars <- rep(0, ncomp)
   purevals <- rep(0, ncomp)

   # scale columns of the original spectra to unit length
   Dr <- D %*% diag(1 / sqrt(colSums(D^2)))

   # initialize loadings
   for (i in seq_len(ncomp)) {

      # compute loadings (V in the original paper)
      if (i  == 1) {
         V <- svd(matrix(1, nspecvis, 1))$u[, seq_len(i), drop = FALSE]
      } else {
         V <- svd(D[, purevars, drop = FALSE])$u[, seq_len(i - 1), drop = FALSE]
      }

      # project data to the previously selected pure variable space
      s <- crossprod(V, Dr)

      # compute angles and keep them as purity spectrum
      s <- sqrt(colSums(s^2))
      s[s > 1] <- 1
      purityspec[i, colind] <- acos(s) / pi * 180 * n

      # find the largest value in the purity spectrum and save corresponding
      # column index and next pure variable
      if (purevars[i] == 0) purevars[i] <- which.max(purityspec[i, colind])
      purevals[i] <- purityspec[i, purevars[i]]
   }

   purevars <- colind[purevars]
   names(purevars) <- names(purevals) <- rownames(purityspec) <- paste("Comp", seq_len(ncomp))

   purityspec <- mda.setattr(purityspec, attrs, "col")
   attr(purityspec, "name") <- "Purity spectra"
   attr(purityspec, "yaxis.name") <- "Purity"

   attr(purevals, "name") <- "Purity"
   attr(purevals, "xaxis.name") <- "Components"

   attr(purevars, "xaxis.name") <- attrs$xaxis.name
   attr(purevars, "xaxis.values") <- purevars
   if (!is.null(attrs$xaxis.values)) {
      attr(purevars, "xaxis.values") <- attrs$xaxis.values[purevars]
   }

   return(
      list(
         ncomp = ncomp,
         offset = offset,
         purityspec = mda.t(purityspec),
         purevars = purevars,
         purevals = purevals
      )
   )
}


########################
#  Plotting methods    #
########################


#' Purity values plot
#'
#' @param obj
#' \code{mcrpure} object
#' @param xticks
#' ticks for x axis
#' @param type
#' type of the plot
#' @param labels
#' what to use as data labels
#' @param ...
#' other parameters suitable for \code{mdaplot}
#'
#' The plot shows largest weighted purity value for each component graphically.
#'
#' @export
plotPurity.mcrpure <- function(obj, xticks = seq_len(obj$ncomp), type = "h",
   labels = "values", ...) {

   p <- mdaplot(obj$purevals, xticks = xticks, type = type, labels = labels, ...)
   invisible(p)
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
   lines.col = adjustcolor(col, alpha.f = 0.75), lines.lty = 3, lines.lwd = 1, ...) {

   stopifnot("Parameter 'comp' has wrong value." = min(comp) > 0 && max(comp) <= obj$ncomp)
   p <- mdaplotg(mda.subset(mda.t(obj$purityspec), comp), type = type, col = col[comp], ...)

   if (show.lines) {
      abline(
         v = attr(obj$purevars, "xaxis.values")[comp],
         col = lines.col[comp],
         lty = lines.lty,
         lwd = lines.lwd
      )
   }

   invisible(p)
}
