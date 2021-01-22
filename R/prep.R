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

   f <- function(data, center, scale, max.cov) {
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
      attr(data, "scaled:center") <- NULL
      attr(data, "scaled:scale") <- NULL
      attr(data, "prep:center") <- center
      attr(data, "prep:scale") <- scale

      return(data)
   }

   return(prep.generic(data, f, center = center, scale = scale, max.cov = max.cov))
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

   f <- function(data) t(scale(t(data), center = TRUE, scale = TRUE))
   return(prep.generic(data, f))
}

#' Normalization
#'
#' @description
#' Normalizes signals (rows of data matrix) to unit area, unit length or unit sum
#'
#' @param data
#' a matrix with data values
#' @param type
#' type of normalization \code{"area"}, \code{"length"} or \code{"sum"}.
#'
#' @return
#' data matrix with normalized values
#'
#' @export
prep.norm <- function(data, type = "area") {

   f <- function(data, type) {

      w <- switch(
         type,
         "area" = apply(abs(data), 1, sum),
         "length" = sqrt(apply(data^2, 1, sum)),
         "sum" = apply(data, 1, sum)
      )

      if (is.null(w)) stop("Wrong value for argument 'type'.")
      return(sweep(data, 1, w, "/"))
   }

   return(prep.generic(data, f, type = type))
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

   f <- function(data, width, porder, dorder) {
      stopifnot("Filter width ('width') should be equal at least to 3." = width > 2)
      stopifnot("Filter width ('width') should be an odd integer number." = width %% 2 == 1)
      stopifnot("Wrong value for the derivative order (should be 0, 1, or 2)." = dorder %in% (0:2))
      stopifnot("Wrong value for the polynomial degree (should be integer number between 0 and 4)." = porder %in% (0:4))
      stopifnot("Polynomal degree ('porder') should not be smaller the derivative order ('dorder')." = porder >= dorder)

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

      return(pdata)
   }

   return(prep.generic(data, f, width = width, porder = porder, dorder = dorder))
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

   f <- function(spectra, mspectrum) {
      if (is.null(mspectrum)) {
         mspectrum <- apply(spectra, 2, mean)
      }

      if (!is.null(mspectrum)) {
         dim(mspectrum) <- NULL
      }

      if (length(mspectrum) != ncol(spectra)) {
         stop("Length of 'mspectrum' should be the same as number of columns in 'spectra'.")
      }

      pspectra <- matrix(0, nrow = nrow(spectra), ncol = ncol(spectra))
      for (i in seq_len(nrow(spectra))) {
         coef <- coef(lm(spectra[i, ] ~ mspectrum))
         pspectra[i, ] <- (spectra[i, ] - coef[1]) / coef[2]
      }

      attr(pspectra, "mspectrum") <- mspectrum
      return(pspectra)
   }

   return(prep.generic(spectra, f, mspectrum = mspectrum))
}


#' Kubelka-Munk transformation
#'
#' @description
#' Applies Kubelka-Munk (km) transformation to data matrix (spectra)
#'
#' @param spectra
#' a matrix with spectra values (absolute reflectance values)
#'
#' @return
#' preprocessed spectra.
#'
#' @details
#' Kubelka-Munk is useful preprocessing method for diffuse reflection spectra (e.g. taken for
#' powders or rough surface). It transforms the reflectance spectra R to K/M units as follows:
#' (1 - R)^2 / 2R
#'
#' @export
prep.ref2km <- function(spectra) {
   stopifnot("Can't use Kubelka-Munk transformation as some of the values are zeros or negative." = all(spectra > 0))
   f <- function(x) (1 - x)^2 / (2 * x)

   return(prep.generic(spectra, f))
}

#' Baseline correction using assymetric least squares
#'
#' @param spectra
#' matrix with spectra (rows correspond to individual spectra)
#' @param plambda
#' power of the penalty parameter (e.g. if plambda = 5, lambda = 10^5)
#' @param p
#' assymetry ratio (should be between 0 and 1)
#' @param max.niter
#' maximum number of iterations
#'
#' @details
#' The function implements baseline correction algorithm based on Whittaker smoother. The method
#' was first shown in [1]. The function has two main parameters - power of a penalty parameter
#' (usually varies betwen 2 and 9) and the ratio of assymetry (usually between 0.1 and 0.001). The
#' choice of the parameters depends on how broad the disturbances of the baseline are and how
#' narrow the original spectral peaks are.
#'
#' @examples
#' # take spectra from carbs dataset
#' data(carbs)
#' spectra = mda.t(carbs$S)
#'
#' # apply the correction
#' pspectra = prep.alsbasecorr(spectra, plambda = 3, p = 0.01)
#'
#' # show the original and the corrected spectra individually
#' par(mfrow = c(3, 1))
#' for (i in 1:3) {
#'    mdaplotg(list(
#'       original = mda.subset(spectra, i),
#'       corrected = mda.subset(pspectra, i)
#'    ), type = "l", col = c("black", "red"), lwd = c(2, 1), main = rownames(spectra)[i])
#' }
#'
#' @importFrom Matrix Matrix Diagonal
#'
#' @export
prep.alsbasecorr <- function(spectra, plambda = 5, p = 0.1, max.niter = 10) {
   attrs <- mda.getattr(spectra)
   dimnames <- dimnames(spectra)

   baseline <- function(y) {

      m <- length(y)
      D <- Matrix::Matrix(diff(diag(m), difference = 2), sparse = TRUE)
      LDD <- (10^plambda) * Matrix::crossprod(D)
      w <- Matrix::Matrix(1, nrow = m, ncol = 1, sparse = TRUE)

      for (i in seq_len(max.niter)) {
         W <- Matrix::Diagonal(x = as.numeric(w))
         z <- Matrix::solve(W + LDD, w * y)
         w.old <- w
         w <- p * (y > z) + (1 - p) * (y < z)

         if (sum(abs(w - w.old)) < 10^-12) {
            break
         }
      }

      return(as.numeric(z))
   }

   pspectra <- t(apply(spectra, 1, function(x) x - baseline(x)))
   pspectra <- mda.setattr(pspectra, attrs)
   dimnames(pspectra) <- dimnames

   return(pspectra)
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

#' Generic function for preprocessing
#'
#' @param x
#' data matrix to be preprocessed
#' @param f
#' function for preprocessing
#' @param ...
#' arguments for the function f
#'
#'
prep.generic <- function(x, f, ...) {
   attrs <- mda.getattr(x)
   dimnames <- dimnames(x)
   x.p <- f(x, ...)
   x.p <- mda.setattr(x.p, attrs)
   dimnames(x.p) <- dimnames
   return(x.p)
}