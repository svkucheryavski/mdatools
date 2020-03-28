#' Autoscale values
#'
#' @description
#' Autoscale (mean center and standardize) values in columns of data matrix.
#'
#' @param data
#' a matrix with data values
#' @param center
#' a logical value or vector with numbers for centering
#' @param scale
#' a logical value or vector with numbers for weighting
#' @param max.cov
#' columns that have coefficient of variation (in percent) below or equal to `max.cov` will not
#' be scaled
#'
#' @return
#' data matrix with processed values
#'
#' @description
#'
#' The use of `max.cov` allows to avoid overestimation of inert variables, which vary
#' very little. Note, that the `max.cov` value is already in percent, e.g. if `max.cov = 0.1` it
#' will compare the coefficient of variation of every variable with 0.1% (not 1%). If you do not
#' want to use this option simply keep `max.cov = 0`.
#'
#' @export
prep.autoscale <- function(data, center = TRUE, scale = FALSE, max.cov = 0) {

   attrs <- mda.getattr(data)
   dimnames <- dimnames(data)

   # define values for centering
   if (is.logical(center) && center) center <- apply(data, 2, mean)

   if (is.numeric(center) && length(center) != ncol(data)) {
      stop("Number of values in 'center' should be the same as number of columns in 'daata'")
   }

   # define values for weigting
   if (is.logical(scale) && scale) scale <- apply(data, 2, sd)

   if (is.numeric(scale) && length(scale) != ncol(data)) {
      stop("Number of values in 'scale' should be the same as number of columns in 'daata'")
   }

   # compute coefficient of variation and set scale to 1 if it is below
   # a user defined threshold
   if (is.numeric(scale)) {
      m <- if (is.numeric(center)) center else apply(data, 2, mean)
      cv <- scale / abs(m) * 100
      scale[is.nan(cv) | cv <= max.cov] <- 1
   }

   # make autoscaling and attach preprocessing attributes
   data <- scale(data, center = center, scale = scale)

   data <- mda.setattr(data, attrs)
   attr(data, "scaled:center") <- NULL
   attr(data, "scaled:scale") <- NULL
   attr(data, "prep:center") <- center
   attr(data, "prep:scale") <- scale
   dimnames(data) <- dimnames

   return(data)
}

#' Standard Normal Variate transformation
#'
#' @description
#' Applies Standard Normal Variate (SNV) transformation to the rows of data matrix
#'
#' @param data
#' a matrix with data values
#'
#' @return
#' data matrix with processed values
#'
#' @details
#' SNV is a simple preprocessing to remove scatter effects (baseline offset and slope) from
#' spectral data, e.g. NIR spectra.
#'
#'  @examples
#'
#'  ### Apply SNV to spectra from simdata
#'
#'  library(mdatools)
#'  data(simdata)
#'
#'  spectra = simdata$spectra.c
#'  wavelength = simdata$wavelength
#'
#'  cspectra = prep.snv(spectra)
#'
#'  par(mfrow = c(2, 1))
#'  mdaplot(cbind(wavelength, t(spectra)), type = 'l', main = 'Before SNV')
#'  mdaplot(cbind(wavelength, t(cspectra)), type = 'l', main = 'After SNV')
#'
#' @export
prep.snv <- function(data) {
   attrs <- mda.getattr(data)
   dimnames <- dimnames(data)

   data <- t(scale(t(data), center = TRUE, scale = TRUE))
   data <- mda.setattr(data, attrs)
   dimnames(data) <- dimnames

   return(data)
}

#' Normalization
#'
#' @description
#' Normalizes signals (rows of data matrix) to unit area or unit length
#'
#' @param data
#' a matrix with data values
#' @param type
#' type of normalization \code{"area"} or \code{"length"}
#'
#' @return
#' data matrix with normalized values
#'
#' @export
prep.norm <- function(data, type = "area") {
   attrs <- mda.getattr(data)
   dimnames <- dimnames(data)

   w <- switch(
      type,
      "area" = apply(abs(data), 1, sum),
      "length" = sqrt(apply(data^2, 1, sum))
   )

   if (is.null(w)) {
      stop("Wrong value for argument 'type'.")
   }

   data <- sweep(data, 1, w, "/")
   data <- mda.setattr(data, attrs)
   dimnames(data) <- dimnames

   return(data)
}

#' Savytzky-Golay filter
#'
#' @description
#' Applies Savytzky-Golay filter to the rows of data matrix
#'
#' @param data
#' a matrix with data values
#' @param width
#' width of the filter window
#' @param porder
#' order of polynomial used for smoothing
#' @param dorder
#' order of derivative to take (0 - no derivative)
#'
#' @export
prep.savgol <- function(data, width = 3, porder = 1, dorder = 0) {
   attrs <- mda.getattr(data)
   dimnames <- dimnames(data)

   if (width < 3) {
      stop("Filter width ('width') should be equal at least to 3.")
   }

   if (width %% 2 != 1) {
      stop("Filter width ('width') should be an odd integer number.")
   }

   if (dorder < 0 || dorder > 2) {
      stop("Wrong value for the derivative order (should be 0, 1, or 2).")
   }

   if (porder < 0 || porder > 4) {
      stop("Wrong value for the polynomial degree (should be integer number between 0 and 4).")
   }

   if (porder < dorder) {
      stop("Polynomal degree ('porder') should not be smaller the derivative order ('dorder').")
   }

   nobj <- nrow(data)
   nvar <- ncol(data)
   pdata <- matrix(0, ncol = nvar, nrow = nobj)

   for (i in seq_len(nobj)) {
      d <- data[i, ]
      w <- (width - 1) / 2
      f <- pinv(outer(-w:w, 0:porder, FUN = "^"))
      d <- convolve(d, rev(f[dorder + 1, ]), type = "o")
      pdata[i, ] <- d[(w + 1) : (length(d) - w)]
   }

   pdata <- mda.setattr(pdata, attrs)
   dimnames(pdata) <- dimnames

   return(pdata)
}

#' Multiplicative Scatter Correction transformation
#'
#' @description
#' Applies Multiplicative Scatter Correction (MSC) transformation to data matrix (spectra)
#'
#' @param spectra
#' a matrix with spectra values
#' @param mspectrum
#' mean spectrum (if NULL will be calculated from \code{spectra})
#'
#' @return
#' preprocessed spectra (calculated mean spectrum is assigned as attribut 'mspectrum')
#'
#' @details
#' MSC is used to remove scatter effects (baseline offset and slope) from
#' spectral data, e.g. NIR spectra.
#'
#'  @examples
#'
#'  ### Apply MSC to spectra from simdata
#'
#'  library(mdatools)
#'  data(simdata)
#'
#'  spectra = simdata$spectra.c
#'  cspectra = prep.msc(spectra)
#'
#'  par(mfrow = c(2, 1))
#'  mdaplot(spectra, type = "l", main = "Before MSC")
#'  mdaplot(cspectra, type = "l", main = "After MSC")
#'
#' @export
prep.msc <- function(spectra, mspectrum = NULL) {
   attrs <- mda.getattr(spectra)
   dimnames <- dimnames(spectra)

   if (is.null(mspectrum)) {
      mspectrum <- apply(spectra, 2, mean)
   }

   if (!is.null(mspectrum)) {
      dim(mspectrum) <- NULL
   }

   if (length(mspectrum) != ncol(spectra)) {
      stop("Length of 'mspectrum' should be the same as number of columns in 'spectra'.")
   }

   cspectra <- matrix(0, nrow = nrow(spectra), ncol = ncol(spectra))
   for (i in seq_len(nrow(spectra))) {
      coef <- coef(lm(spectra[i, ] ~ mspectrum))
      cspectra[i, ] <- (spectra[i, ] - coef[1]) / coef[2]
   }

   cspectra <- mda.setattr(cspectra, attrs)
   attr(cspectra, "mspectrum") <- mspectrum
   dimnames(cspectra) <- dimnames

   return(cspectra)
}

#' Pseudo-inverse matrix
#'
#' @description
#' Computes pseudo-inverse matrix using SVD
#'
#' @param data
#' a matrix with data values to compute inverse for
#'
#' @export
pinv <- function(data) {
   # Calculates pseudo-inverse of data matrix
   s <- svd(data)
   s$v %*% diag(1 / s$d) %*% t(s$u)
}
