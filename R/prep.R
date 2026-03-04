############################################################
# Main methods for preprocessing.                          #
############################################################


#' Remove spikes from Raman spectra
#'
#' @description
#' Removes spikes from Raman spectra using median absolute deviation computed for signal
#' differences. See [1] for more details
#'
#' @param data
#' a matrix with data values
#' @param width
#' width of the moving median filter
#' @param threshold
#' threshold to compare modified z-score value with to detect the spike.
#'
#' @return
#' data matrix with processed values
#'
#' @references
#' 1. Darren A. Whitaker, Kevin Hayes, A simple algorithm for despiking Raman spectra,
#' Chemometrics and Intelligent Laboratory Systems, 179, 2018, pp. 82-84,
#' 10.1016/j.chemolab.2018.06.009.
#'
#' @export
prep.spikes <- function(data, width = 5, threshold = 6) {

   stopifnot("Parameter 'width' should be at least equal to 3." = width > 2)
   stopifnot("Parameter 'width' should be an odd integer number." = width %% 2 == 1)

   f <- function(data, width, threshold) {

      n <- ncol(data)
      half <- round((width - 1) / 2)

      z <- apply(data, 1, function(x) {
         dx <- diff(x)
         med <- median(dx, na.rm = TRUE)
         mad <- mad(dx, na.rm = TRUE)

         # correct mad value to avoid division by zero
         if (mad < 1e-10) mad <- 1

         # compute z here we do not multiply to constant because mad() does this already
         abs((dx - med) / mad)
      })

      z <- t(rbind(rep(0, ncol(z)), z))
      spikes.map <- z > threshold
      spikes.rows <- which(rowSums(spikes.map) > 0)
      if (length(spikes.rows) > 0) {
         for (rr in spikes.rows) {
            spikes.cols <- which(spikes.map[rr, ] > 0)
            for (cc in spikes.cols) {
               cc.ind <- seq(max(cc - half, 1), min(cc + half, n))
               cc.ind <- cc.ind[spikes.map[rr, cc.ind] == 0]
               if (length(cc.ind) > 0)
                  data[rr, cc] <- mean(data[rr, cc.ind])
            }
         }
      }

      return(data)
   }

   return(prep.generic(data, f, width = width, threshold = threshold))
}


#' Normalization
#'
#' @description
#' Normalizes signals (rows of data matrix).
#'
#' @param data
#' a matrix with data values
#' @param type
#' type of normalization \code{"area"}, \code{"length"}, \code{"sum"}, \code{"snv"}, \code{"is"}, or \code{"pqn"}.
#' @param col.ind
#' indices of columns (can be either integer or logical values) for normalization to internal
#' standard peak.
#' @param ref.spectrum
#' reference spectrum for PQN normalization, if not provided a mean spectrum for data is used
#'
#' @details
#' The \code{"area"}, \code{"length"}, \code{"sum"} types do preprocessing to unit area (sum of
#' absolute values), length or sum of all values in every row of data matrix. Type \code{"snv"}
#' does the Standard Normal Variate normalization. Type
#' \code{"is"} does the normalization to internal standard peak, whose position is defined by
#' parameter `col.ind`. If the position is a single value, the rows are normalized to the height
#' of this peak. If `col.ind` points to several adjacent values, the rows are normalized to the area
#' under the peak - sum of the intensities.
#'
#' The \code{"pqn"} is Probabilistic Quotient Normalization as described in [1]. In this case you also
#' need to provide a reference spectrum (e.g. mean or median of spectra for some reference samples). If
#' reference spectrum is not provided it will be computed as mean of the spectra to be
#' preprocessed (parameter \code{data}).
#'
#' @references
#' 1. F. Dieterle, A. Ross, H. Senn. Probabilistic Quotient Normalization as Robust Method to
#' Account for Dilution of Complex Biological Mixtures. Application in 1 H NMR Metabonomics.
#' Anal. Chem. 2006, 78, 4281–4290.
#'
#' @return
#' data matrix with normalized values
#'
#' @export
prep.norm <- function(data, type = "area", col.ind = NULL, ref.spectrum = NULL) {

   type <- match.arg(type, c("area", "length", "sum", "snv", "is", "pqn"))

   if (type == "is" && is.null(col.ind)) {
      stop("For 'is' normalization you need to provide indices for IS peak.", call. = FALSE)
   }

   if (is.logical(col.ind)) {
      col.ind <- which(col.ind)
   }

   if (!is.null(col.ind) && (min(col.ind) < 1 || max(col.ind) > ncol(data))) {
      stop("Values for 'col.ind' seem to be wrong.", call. = FALSE)
   }

   if (type == "pqn" && is.null(ref.spectrum)) {
      ref.spectrum <- colMeans(data)
   }

   pqn <- function(data, ref.spectrum) {

      if (length(ref.spectrum) != ncol(data)) {
         stop("prep.norm: 'ref.spectrum' should have the same number of values as the number of columns in 'data'.", call. = FALSE)
      }

      # 1. unit area normalization
      ref.spectrum <- as.numeric(ref.spectrum)
      ref.spectrum <- ref.spectrum / sum(abs(ref.spectrum))

      # 2. compute  and return median quotients for each spectrum
      return(apply(sweep(data, 2, ref.spectrum, "/"), 1, median))
   }

   f <- function(data, type, col.ind, ref.spectrum) {

      # preliminary normalize the dataset to unit sum
      if (type == "pqn") {
         data <- prep.norm(data, type = "area")
      }

      if (type == "snv") {
         data <- sweep(data, 1, rowMeans(data), "-")
      }

      w <- switch(
         type,
         "snv" = apply(data, 1, sd),
         "area" = rowSums(abs(data)),
         "length" = sqrt(rowSums(data^2)),
         "sum" = rowSums(data),
         "is" = rowSums(data[, col.ind, drop = FALSE]),
         "pqn" = pqn(data, ref.spectrum)
      )

      if (is.null(w)) stop("Wrong value for argument 'type'.", call. = FALSE)
      return(sweep(data, 1, w, "/"))
   }

   return(prep.generic(data, f, type = type, col.ind = col.ind, ref.spectrum = ref.spectrum))
}


#' Savitzky-Golay filter
#'
#' @description
#' Applies Savitzky-Golay filter to the rows of data matrix
#'
#' @param data
#' a matrix with data values
#' @param width
#' width of the filter window
#' @param porder
#' order of polynomial used for smoothing
#' @param dorder
#' order of derivative to take (0 - no derivative)
#' @param w
#' do not use, required for training of preprocessing model.
#'
#' @details
#' The function implements algorithm described in [1] which handles the edge points correctly and
#' does not require to cut the spectra.
#'
#' @references
#' 1. Peter A. Gorry. General least-squares smoothing and differentiation by the convolution
#' (Savitzky-Golay) method. Anal. Chem. 1990, 62, 6, 570–573, https://doi.org/10.1021/ac00205a007.
#'
#' @export
prep.savgol <- function(data, width = 3, porder = 1, dorder = 0, w = NULL) {

   stopifnot("Filter width ('width') should be equal at least to 3." = width > 2)
   stopifnot("Filter width ('width') should be an odd integer number." = width %% 2 == 1)
   stopifnot("Wrong value for the derivative order (should be 0, 1, or 2)." = dorder %in% (0:2))
   stopifnot("Wrong value for the polynomial degree (should be integer number between 1 and 4)." = porder %in% (1:4))
   stopifnot("Polynomial degree ('porder') should not be smaller than the derivative order ('dorder')." = porder >= dorder)

   f <- function(x) {
      if (is.null(w)) {
         props <- prep.savgol.params(data, width, porder, dorder)
         w <- props$w
      }

      m <- round((ncol(w) - 1)/2)
      nvar <- ncol(x)
      px <- matrix(0.0, nrow(x), ncol(x))

      for (i in seq_len(m)) {
         px[, i] <- apply(x[, seq_len(2 * m + 1), drop = FALSE], 1,
            function(xx) convolve(xx, w[, i], type = "filter")[1])
         px[, nvar - i + 1] <- apply(x[, (nvar - 2 * m):nvar, drop = FALSE], 1,
            function(xx) convolve(xx, w[, width - i + 1], type = "filter")[1])
      }

      px[, (m + 1):(nvar - m)] <- t(apply(x, 1, function(xx) convolve(xx, w[, m + 1], type = "filter")))

      return(px)
   }

   return(prep.generic(data, f))
}


#' Kubelka-Munk transformation
#'
#' @description
#' Applies Kubelka-Munk (km) transformation to data matrix (spectra)
#'
#' @param data
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
prep.ref2km <- function(data) {
   stopifnot("Can't use Kubelka-Munk transformation as some of the values are zeros or negative." = all(data > 0))
   f <- function(x) (1 - x)^2 / (2 * x)

   return(prep.generic(data, f))
}


#' Baseline correction using asymmetric least squares
#'
#' @param data
#' matrix with spectra (rows correspond to individual spectra)
#' @param plambda
#' power of the penalty parameter (e.g. if plambda = 5, lambda = 10^5)
#' @param p
#' asymmetry ratio (should be between 0 and 1)
#' @param max.niter
#' maximum number of iterations
#'
#' @return
#' preprocessed spectra.
#'
#' @details
#' The function implements baseline correction algorithm based on Whittaker smoother. The method
#' was first shown in [1]. The function has two main parameters - power of a penalty parameter
#' (usually varies between 2 and 9) and the ratio of asymmetry (usually between 0.1 and 0.001). The
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
#' @importFrom spam spam diag.spam
#' @importFrom methods as
#'
#' @export
prep.alsbasecorr <- function(data, plambda = 5, p = 0.1, max.niter = 10) {

   f <- function(data, plambda, p, max.niter) {

      n <- nrow(data)
      m <- ncol(data)
      baseline <- matrix(0, n, m)

      # Build the second-difference penalty matrix D
      D <- spam(0, m - 2, m)
      for (i in seq_len(m - 2)) {
         D[i, i]     <- 1
         D[i, i + 1] <- -2
         D[i, i + 2] <- 1
      }

      # Compute LDD = lambda * t(D) %*% D
      lambda <- 10^plambda
      LDD <- lambda * crossprod(D)

      # Precompute initial weight vector
      w.ini <- rep(1, m)

      for (i in seq_len(n)) {
         y <- data[i, ]
         w <- w.ini

         for (j in seq_len(max.niter)) {
            W <- diag.spam(w)
            A <- W + LDD
            b <- w * y

            # Solve (W + LDD) z = w * y
            z <- solve(A, b)

            w.old <- w
            w <- p * (y > z) + (1 - p) * (y < z)

            if (sum(abs(w - w.old)) < 1e-10) break
         }

         baseline[i, ] <- as.numeric(z)
      }

      return(data - baseline)
   }

   return(prep.generic(data, f, plambda = plambda, p = p, max.niter = max.niter))
}


#' Transformation
#'
#' @description
#' Transforms values using any mathematical function (e.g. log).
#'
#' @param data
#' a matrix with data values
#' @param fun
#' reference to a transformation function, e.g. `log` or `function(x) x^2`.
#' @param ...
#' optional parameters for the transformation function
#'
#' @return
#' data matrix with transformed values
#'
#' @examples
#' # generate a matrix with two columns
#' y <- cbind(rnorm(100, 10, 1), rnorm(100, 20, 2))
#'
#' # apply log transformation
#' py1 = prep.transform(y, log)
#'
#' # apply power transformation
#' py2 = prep.transform(y, function(x) x^-1.25)
#'
#' # show distributions
#' par(mfrow = c(2, 3))
#' for (i in 1:2) {
#'    hist(y[, i], main = paste0("Original values, column #", i))
#'    hist(py1[, i], main = paste0("Log-transformed, column #", i))
#'    hist(py2[, i], main = paste0("Power-transformed, column #", i))
#' }
#'
#' @export
prep.transform <- function(data, fun, ...) {
   return(prep.generic(data, fun, ...))
}


#' Variable selection
#'
#' @description
#' Returns dataset with selected variables
#'
#' @param data
#' a matrix with data values
#' @param var.ind
#' indices of variables (columns) to select, can be either numeric or logical
#'
#' @return
#' data matrix with the selected variables (columns)
#'
#' @export
prep.varsel <- function(data, var.ind) {
   if (!is.null(attr(data, "exclcols"))) {
      stop("prep.varsel() can not be used for dataset with excluded (hidden) columns.", call. = FALSE)
   }
   return(mda.subset(data, select = var.ind))
}


#' Centering data columns.
#'
#' @param data
#' a matrix with data values.
#' @param type
#' type of statistic to use for centering ('mean', or 'median').
#' @param center
#' do not use, required for training of preprocessing model.
#'
#' @return
#' preprocessed data matrix
#'
#' @export
prep.center <- function(data, type = "mean", center = NULL) {

   f <- function(data) {

      # define values for centering
      if (is.null(center)) {
         params <- prep.center.params(data, type)
         center <- params$center
      }

      if (length(center) != ncol(data)) {
         stop("Number of values in 'center' should be the same as number of columns in 'data'", call. = FALSE)
      }

      # center data and attach preprocessing attributes
      data <- scale(data, center = center, scale = FALSE)
      attr(data, "scaled:center") <- NULL
      attr(data, "scaled:scale") <- NULL
      attr(data, "prep:center") <- center
      return(data)
   }

   return(prep.generic(data, f))
}


#' Scaling data columns.
#'
#' @param data
#' a matrix with data values.
#' @param type
#' type of statistic to use for scaling ('sd', 'iqr', 'range', 'pareto')
#' @param max.cov
#' columns that have coefficient of variation (in percent) below or equal to `max.cov` will not
#' be scaled.
#' @param scale
#' do not use, required for training of preprocessing model.
#'
#' @return
#' preprocessed data matrix
#'
#' @export
prep.scale <- function(data, type = "sd", max.cov = 0, scale = NULL) {

   f <- function(data) {

      # define values for scaling
      if (is.null(scale)) {
         params <- prep.scale.params(data, type, max.cov, scale)
         scale <- params$scale
      }

      if (length(scale) != ncol(data)) {
         stop("Number of values in 'scale' should be the same as number of columns in 'data'", call. = FALSE)
      }

      # scale data and attach preprocessing attributes
      data <- scale(data, center = FALSE, scale = scale)
      attr(data, "scaled:center") <- NULL
      attr(data, "scaled:scale") <- NULL
      attr(data, "prep:scale") <- scale
      return(data)
   }

   return(prep.generic(data, f))
}


#' Applies Extended Multiplicative Scatter Correction to data rows
#'
#' @param data
#' a matrix with data values.
#' @param degree
#' polynomial degree, if 0 then the result will be the same as for conventional MSC.
#' @param mspectrum
#' optional reference spectrum (if not provided, mean spectrum will be used).
#' @param lnorm
#' do not use, required for training of preprocessing model.
#' @param A
#' do not use, required for training of preprocessing model.
#'
#' @export
prep.emsc <- function(data, degree = 0, mspectrum = NULL, lnorm = NULL, A = NULL) {

   f <- function(data) {

      # define values for EMSC
      if (is.null(A) || is.null(lnorm)) {
         params <- prep.emsc.params(data, degree, mspectrum)
         mspectrum <- params$mspectrum
         A <- params$A
         lnorm <- params$lnorm
      }

      if (length(mspectrum) != ncol(data)) {
         stop("Number of values in 'mspectrum' should be the same as number of columns in 'data'", call. = FALSE)
      }

      # solve A * C = X'
      C <- qr.solve(A, t(data))

      if (ncol(A) > 2) {
         p <- ncol(A)
         B <- t(A[, 3:p, drop = FALSE] %*% C[3:p, , drop = FALSE])
         data <- (data - C[1, ] - B) / C[2, ]
      } else {
         data <- (data - C[1, ]) / C[2, ]
      }

      names(mspectrum) <- NULL
      dim(mspectrum) <- NULL
      attr(data, "mspectrum") <- mspectrum
      return(data)
   }

   return(prep.generic(data, f))
}



############################################################
# Methods for parameters estimation/fitting                #
############################################################

#' Precomputes parameters for normalization
#'
#' @param data
#' a matrix with data values.
#' @param type
#' type of normalization.
#' @param col.ind
#' indices of columns for IS normalization.
#' @param ref.spectrum
#' reference spectrum for PQN normalization.
#'
#' @return
#' list with parameter values
prep.norm.params <- function(data, type = "area", col.ind = NULL, ref.spectrum = NULL) {
   if (type == "pqn" && is.null(ref.spectrum)) {
      ref.spectrum <- colMeans(data)
   }
   return(list(type = type, col.ind = col.ind, ref.spectrum = ref.spectrum))
}


#' Precomputes parameters for centering
#'
#' @param data
#' a matrix with data values.
#' @param type
#' type of statistic to use for centering (\code{'mean'}, or \code{'median'}).
#' @param center
#' vector with precomputed values for centering.
#'
#' @return
#' list with parameter values
prep.center.params <- function(data, type = "mean", center = NULL) {

   type <- match.arg(type, c("mean", "median"))
   if (!is.null(center)) {
      return(list(center = center))
   }

   data <- mda.purgeRows(data)
   f <- c("mean" = mean, "median" = median)
   center <- if (type == "mean") colMeans(data) else apply(data, 2, f[[type]])
   return(list(type = type, center = center))
}


#' Precomputes parameters for scaling
#'
#' @param data
#' a matrix with data values.
#' @param type
#' type of statistic to use for scaling (\code{'sd'}, \code{'iqr'}, \code{'range'}, or \code{'pareto'})
#' @param max.cov
#' columns that have coefficient of variation (in percent) below or equal to `max.cov` will not
#' be scaled.
#' @param scale
#' vector with precomputed values for scaling.
#'
#' @return
#' list with parameter values
prep.scale.params <- function(data, type = "sd", max.cov = 0, scale = NULL) {

   type <- match.arg(type, c("sd", "iqr", "range", "pareto"))
   if (!is.null(scale)) {
      return(list(scale = scale))
   }

   data <- mda.purgeRows(data)
   f <- c("sd" = sd, "iqr" = IQR, "range" = function(x) max(x) - min(x), "pareto" = function(x) sqrt(sd(x)))
   scale <- apply(data, 2, f[[type]])

   if (max.cov > 0) {
      m <- colMeans(data)
      cv <- scale / abs(m) * 100
      scale[is.nan(cv) | cv <= max.cov] <- 1
   }
   return(list(type = type, scale = scale))
}


#' Precomputes parameters for EMSC
#'
#' @param data
#' a matrix with data values.
#' @param degree
#' polynomial degree.
#' @param mspectrum
#' reference spectrum.
#'
#' @return
#' list with parameter values
prep.emsc.params <- function(data, degree = 0, mspectrum = NULL) {

   nx <- ncol(data)
   p <- 2 + degree

   # wavelength/wavenumbers indices from [1, n] to [-1, 1]
   lnorm <- seq(-1, 1, length.out = nx)

   if (is.null(mspectrum)) {
      mspectrum <- colMeans(data)
   }
   dim(mspectrum) <- NULL
   names(mspectrum) <- NULL

   # build design matrix A: nx x p
   A <- matrix(1.0, nx, p)    # first column is just ones for intercept
   A[, 2] <- mspectrum        # second column contains reference spectrum (slope)
   if (degree > 0) {
      for (d in seq_len(degree)) {
         A[, d + 2] <- lnorm^(d)
      }
   }

   return(list(
      mspectrum = mspectrum,
      lnorm = lnorm,
      A = A
   ))
}


#' Precomputes parameters for Savitzky-Golay
#'
#' @param data
#' a matrix with data values.
#' @param width
#' width of the filter window
#' @param porder
#' order of polynomial used for smoothing
#' @param dorder
#' order of derivative to take (0 - no derivative)
#'
#' @return
#' list with parameter values
prep.savgol.params <- function(data, width = 3, porder = 1, dorder = 0) {

   # compute grams polynomials
   gram <- function(i, m, k, s) {
      if (k > 0) {
         return((4 * k - 2) / (k * (2 * m - k + 1)) * (i * gram(i, m, k - 1, s) + s * gram(i, m, k - 1, s - 1)) -
            ((k - 1) * (2 * m + k)) / (k * (2 * m - k + 1)) * gram(i, m, k - 2, s))
      }
      if (k == 0 && s == 0) return(1)
      return(0)
   }

   # compute generalized factorial
   genfact <- function(a, b) {
      f <- 1
      if ((a - b + 1) > a) return(f)

      for (i in (a - b + 1):a) {
         f <- f * i
      }

      return(f)
   }

   # compute weights for convolution depending on position
   weight <- function(i, t, m, n, s) {
      sum <- 0
      for (k in 0:n) {
         sum <- sum + (2 * k + 1) * (genfact(2 * m, k) / genfact(2 * m + k + 1, k + 1)) *
            gram(i, m, k, 0) * gram(t, m, k, s)
      }
      return(sum)
   }

   m <- round((width - 1) / 2)
   w <- outer(-m:m, -m:m, function(x, y) weight(x, y, m, porder, dorder))

   return(list(width = width, porder = porder, dorder = dorder, w = w))
}



############################################################
# Methods for convertion to/from JSON                      #
############################################################

#' Converts preprocessing item from 'prep.varsel' method to JSON elements
#'
#' @param params
#' model parameters precomputed by using prep.fit()
#' @param npred
#' number of predictors in original dataset
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#'
#' @returns list with main elements required for JSON with preprocessing model
#' compatible with preprocessing web-application (ml, mpl, mp, mpl_new, mp_new, info)
#'
prep.varsel.asjson <- function(params, npred, left = NULL, right = NULL) {

   ind <- params$var.ind
   if (is.null(ind)) {
      stop("Parameter 'var.ind' is not specified.", call. = FALSE)
   }

   if (is.logical(ind)) {
      ind <- which(ind)
   }

   if (any(diff(ind) != 1)) {
      stop("Parameter 'var.ind' should be specified as a single interval without 'holes' inside.", call. = FALSE)
   }

   left.local <- if (is.null(left)) min(ind) else min(ind) + left
   right.local <- if (is.null(left)) max(ind) else max(ind) + left
   return(list(
      ml = 4,
      mpl = 2,
      mp = c(left.local - 1, right.local),
      mpl_new = 2,
      mp_new = c(left.local - 1, right.local),
      info = paste0("left: ", left.local - 1, " ", "right: ", right.local)
   ))
}


#' Converts JSON elements to preprocessing item for 'prep.varsel' method
#'
#' @param mp
#' model parameters from JSON (user defined)
#' @param mp_new
#' model parameters from JSON after training
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#' @param left.tot
#' total, accumulated shift on the left part
#'
#' @returns \code{prep} object for the method
#'
prep.varsel.fromjson <- function(mp, mp_new, left = 0, right = 1, left.tot = 0) {
   right.local <- mp[2] - left.tot
   left.local  <- mp[1] - left.tot + 1
   p <- prep("varsel", var.ind = left.local:right.local)
   return(p)
}


#' Converts preprocessing item from 'prep.spikes' method to JSON elements
#'
#' @param params
#' model parameters precomputed by using prep.fit()
#' @param npred
#' number of predictors in original dataset
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#'
#' @returns list with main elements required for JSON with preprocessing model
#' compatible with preprocessing web-application (ml, mpl, mp, mpl_new, mp_new, info)
#'
prep.spikes.asjson <- function(params, npred, left = NULL, right = NULL) {

   width <- params$width
   threshold <- params$threshold

   return(list(
      ml = 5,
      mpl = 2,
      mp = c(width, threshold),
      mpl_new = 2 + 3 * npred + width,
      mp_new = c(width, threshold, rep(0, 3 * npred + width)),
      info = paste0(width, "/", threshold)
   ))
}


#' Converts JSON elements to preprocessing item for 'prep.spikes' method
#'
#' @param mp
#' model parameters from JSON (user defined)
#' @param mp_new
#' model parameters from JSON after training
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#' @param left.tot
#' total, accumulated shift on the left part
#'
#' @returns \code{prep} object for the method
#'
prep.spikes.fromjson <- function(mp, mp_new, left = 0, right = 1, left.tot = 0) {
   width <- mp[1]
   threshold <- mp[2]

   p <- prep("spikes", width = width, threshold = threshold)
   return(p)
}


#' Converts preprocessing item from 'prep.norm' method to JSON elements
#'
#' @param params
#' model parameters precomputed by using prep.fit()
#' @param npred
#' number of predictors in original dataset
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#'
#' @returns list with main elements required for JSON with preprocessing model
#' compatible with preprocessing web-application (ml, mpl, mp, mpl_new, mp_new, info)
#'
prep.norm.asjson <- function(params, npred, left = 0, right = 1) {

   types <- c("snv", "area", "length", "is")
   type <- which(types == params$type) - 1
   if (length(type) == 0) {
      stop("prep.norm.asjson: normalization type '", params$type,
         "' is not supported for JSON export.", call. = FALSE)
   }
   col.ind <- if (!is.null(params$col.ind)) params$col.ind  else 0

   return(list(
      ml = 0,
      mpl = 2,
      mp = c(type, left + col.ind - 1),
      mpl_new = 2,
      mp_new = c(type, left + col.ind - 1),
      info = if (col.ind > 0) paste0("var #", left + col.ind) else params$type
   ))
}


#' Converts JSON elements to preprocessing item for 'prep.norm' method
#'
#' @param mp
#' model parameters from JSON (user defined)
#' @param mp_new
#' model parameters from JSON after training
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#' @param left.tot
#' total, accumulated shift on the left part
#'
#' @returns \code{prep} object for the method
#'
prep.norm.fromjson <- function(mp, mp_new, left = 0, right = 1, left.tot = 0) {
   types <- c("snv", "area", "length", "is")
   type <- types[mp[1] + 1]
   col.ind <- mp[2] - left.tot + 1
   if (type == "is") {
      p <- prep("norm", type = type, col.ind = col.ind)
      p$params <- list(type = type, col.ind = col.ind, ref.spectrum = NULL)
   } else {
      p <- prep("norm", type = type)
      p$params <- list(type = type, col.ind = NULL, ref.spectrum = NULL)
   }
   return(p)
}


#' Converts preprocessing item from 'prep.savgol' method to JSON elements
#'
#' @param params
#' model parameters precomputed by using prep.fit()
#' @param npred
#' number of predictors in original dataset
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#'
#' @returns list with main elements required for JSON with preprocessing model
#' compatible with preprocessing web-application (ml, mpl, mp, mpl_new, mp_new, info)
#'
prep.savgol.asjson <- function(params, npred, left = 0, right = 1) {

   width <- params$width
   porder <- params$porder
   dorder <- params$dorder

   return(list(
      ml = 1,
      mpl = 3,
      mp = c(params$width, params$porder, params$dorder), # esmc, degree, NA
      mpl_new = length(params[["w"]]),
      mp_new = as.numeric(params[["w"]]),
      info = paste0(width, ", ", porder, ", ", dorder)
   ))
}


#' Converts JSON elements to preprocessing item for 'prep.savgol' method
#'
#' @param mp
#' model parameters from JSON (user defined)
#' @param mp_new
#' model parameters from JSON after training
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#' @param left.tot
#' total, accumulated shift on the left part
#'
#' @returns \code{prep} object for the method
#'
prep.savgol.fromjson <- function(mp, mp_new, left = 0, right = 1, left.tot = 0) {

   width <- mp[1]
   porder <- mp[2]
   dorder <- mp[3]
   w <- matrix(mp_new, width, width)
   p <- prep("savgol", width = width, porder = porder, dorder = dorder)
   p$params <- list(width = width, porder = porder, dorder = dorder, w = w)

   return(p)
}


#' Converts preprocessing item from 'prep.emsc' method to JSON elements
#'
#' @param params
#' model parameters precomputed by using prep.fit()
#' @param npred
#' number of predictors in original dataset
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#'
#' @returns list with main elements required for JSON with preprocessing model
#' compatible with preprocessing web-application (ml, mpl, mp, mpl_new, mp_new, info)
#'
prep.emsc.asjson <- function(params, npred, left = 0, right = 1) {
   nterms <- ncol(params$A)
   degree <- nterms - 2

   pad <- NULL
   nvar <- right - left
   vardiff <- npred - nvar
   pad <- rep(0, vardiff * (nterms + 2))

   return(list(
      ml = 3,
      mpl = 3,
      mp = c(0, degree, 0), # esmc, degree, NA
      mpl_new = 2 + npred * nterms + npred + npred + nterms,
      mp_new = c(0, nterms, as.numeric(params$A), as.numeric(params$lnorm), as.numeric(params$mspectrum), rep(0, nterms), pad),
      info = paste0("emsc (", degree, ")")
   ))
}


#' Converts JSON elements to preprocessing item for 'prep.emsc' method
#'
#' @param mp
#' model parameters from JSON (user defined)
#' @param mp_new
#' model parameters from JSON after training
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#' @param left.tot
#' total, accumulated shift on the left part
#'
#' @returns \code{prep} object for the method
#'
prep.emsc.fromjson <- function(mp, mp_new, left = 0, right = 1, left.tot = 0) {

   mpl_new <- length(mp_new)
   nterms <- mp_new[2]
   npred <- round((mpl_new - 2 - nterms) / (2 + nterms))
   degree <- nterms - 2

   nvar <- right - left
   s <- 3
   e <- s + nterms * nvar - 1
   A <- matrix(mp_new[s:e], nvar, nterms)

   s <- e + 1
   e <- s + nvar - 1
   lnorm <- mp_new[s:e]

   s <- e + 1
   e <- s + nvar - 1
   mspectrum <- mp_new[s:e]


   p <- prep("emsc", degree = degree)
   p$params <- list(mspectrum = mspectrum, lnorm = lnorm, A = A)

   return(p)
}


#' Converts preprocessing item from 'prep.center' method to JSON elements
#'
#' @param params
#' model parameters precomputed by using prep.fit()
#' @param npred
#' number of predictors in original dataset
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#'
#' @returns list with main elements required for JSON with preprocessing model
#' compatible with preprocessing web-application (ml, mpl, mp, mpl_new, mp_new, info)
#'
prep.center.asjson <- function(params, npred, left = 0, right = 1) {

   types <- c("mean", "median")
   type <- params[["type"]]
   type.num <- which(types == type)

   pad <- NULL
   nvar <- right - left
   vardiff <- npred - nvar
   pad <- rep(0, vardiff * 2)

   return(list(
      ml = 2,
      mpl = 2,
      mp = c(type.num, 0),
      mpl_new = 2 * npred,
      mp_new = c(as.numeric(params[["center"]]), rep(1, nvar), pad),
      info = paste0(type, "/no")
   ))
}


#' Converts JSON elements to preprocessing item for 'prep.center' method
#'
#' @param mp
#' model parameters from JSON (user defined)
#' @param mp_new
#' model parameters from JSON after training
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#' @param left.tot
#' total, accumulated shift on the left part
#'
#' @returns \code{prep} object for the method
#'
prep.center.fromjson <- function(mp, mp_new, left = 0, right = 1, left.tot = 0) {

   types <- c("mean", "median")
   mpl_new <- length(mp_new)
   npred <- round(mpl_new / 2)
   nvar <- (right - left)
   center <- mp_new[1:nvar]

   type <- types[mp[1]]
   p <- prep("center", type = type)
   p$params <- list(type = type, center = center)

   return(p)
}


#' Converts preprocessing item from 'prep.scale' method to JSON elements
#'
#' @param params
#' model parameters precomputed by using prep.fit()
#' @param npred
#' number of predictors in original dataset
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#'
#' @returns list with main elements required for JSON with preprocessing model
#' compatible with preprocessing web-application (ml, mpl, mp, mpl_new, mp_new, info)
#'
prep.scale.asjson <- function(params, npred, left = 0, right = 1) {

   types.r <- c("sd", "iqr", "range", "pareto")
   type <- params[["type"]]
   type.num <- which(types.r == type)

   pad <- NULL
   nvar <- (right - left)
   vardiff <- npred - nvar
   pad <- rep(0, vardiff * 2)

   return(list(
      ml = 2,
      mpl = 2,
      mp = c(0, type.num), #
      mpl_new = 2 * npred,
      mp_new = c(rep(0, nvar), as.numeric(params[["scale"]]), pad),
      info = paste0("no/", type)
   ))
}


#' Converts JSON elements to preprocessing item for 'prep.scale' method
#'
#' @param mp
#' model parameters from JSON (user defined)
#' @param mp_new
#' model parameters from JSON after training
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#' @param left.tot
#' total, accumulated shift on the left part
#'
#' @returns \code{prep} object for the method
#'
prep.scale.fromjson <- function(mp, mp_new, left = 0, right = 1, left.tot = 0) {
   types.r <- c("sd", "iqr", "range", "pareto")

   mpl_new <- length(mp_new)
   npred <- round(mpl_new / 2)
   nvar <- (right - left)
   scale <- mp_new[(nvar + 1):(2 * nvar)]

   type <- types.r[mp[2]]
   p <- prep("scale", type = type)
   p$params <- list(type = type, scale = scale)

   return(p)
}


#' Converts preprocessing item from 'prep.alsbasecorr' method to JSON elements
#'
#' @param params
#' model parameters precomputed by using prep.fit()
#' @param npred
#' number of predictors in original dataset
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#'
#' @returns list with main elements required for JSON with preprocessing model
#' compatible with preprocessing web-application (ml, mpl, mp, mpl_new, mp_new, info)
#'
prep.alsbasecorr.asjson <- function(params, npred, left = 0, right = 1) {
   return(list(
      ml = 3,
      mpl = 3,
      mp = c(1, params$plambda, params$p), # alsbasecorr, lambda, p
      mpl_new = 3 + npred, # empty space for baseline
      mp_new = c(1, params$plambda, params$p, rep(0, npred)),
      info = paste0("als (", params$plambda, ",", params$p, ")")
   ))
}


#' Converts JSON elements to preprocessing item for 'prep.alsbasecorr' method
#'
#' @param mp
#' model parameters from JSON (user defined)
#' @param mp_new
#' model parameters from JSON after training
#' @param left
#' index of first variable after trimming (if any)
#' @param right
#' index of last variable after trimming (if any)
#' @param left.tot
#' total, accumulated shift on the left part
#'
#' @returns \code{prep} object for the method
#'
prep.alsbasecorr.fromjson <- function(mp, mp_new, left = 0, right = 1, left.tot = 0) {
   p <- prep("alsbasecorr", plambda = mp[2], p = mp[3])
   return(p)
}


######################################################################
# Methods for combining preprocessing methods together,              #
# fitting preprocessing model, and save/load them to/from JSON       #
######################################################################

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
   stopifnot("First argument of preprocessing function should be a matrix." = is.matrix(x))
   attrs <- mda.getattr(x)
   dimnames <- dimnames(x)
   x.p <- f(x, ...)
   x.p <- mda.setattr(x.p, attrs)
   dimnames(x.p) <- dimnames
   return(x.p)
}


#' Shows a list with implemented preprocessing methods
#'
#' @export
getImplementedPrepMethods <- function() {
   list(

      # prep.varsel <- function(data, var.ind)
      "varsel" = list(
         name = "varsel",
         method = prep.varsel,
         params = list(var.ind = NULL),
         params.show = c("var.ind"),
         jmethod = prep.varsel.asjson,
         params.info = list(var.ind = "indices of variables (columns) to select."),
         info = "Select user-defined variables (columns of dataset)."
      ),

      # prep.norm <- function(data, type = "area", col.ind = NULL, ref.spectrum = NULL)
      "norm" = list(
         name = "norm",
         method = prep.norm,
         params = list(type = "area", col.ind = NULL, ref.spectrum = NULL),
         params.show = c("type", "col.ind"),
         pmethod = prep.norm.params,
         jmethod = prep.norm.asjson,
         params.info = list(
            type = "type of normalization ('area', 'sum', 'length', 'is', 'snv', 'pqn').",
            col.ind = "indices of columns (variables) for normalization to internal standard peak.",
            ref.spectrum = "reference spectrum for PQN normalization."
         ),
         info = "Normalization."
      ),

      # prep.savgol <- function(data, width = 3, porder = 1, dorder = 0)
      "savgol" = list(
         name = "savgol",
         method = prep.savgol,
         params = list(width = 3, porder = 1, dorder = 0),
         params.show = c("width", "porder", "dorder"),
         pmethod = prep.savgol.params,
         jmethod = prep.savgol.asjson,
         params.info = list(
            width = "width of the filter.",
            porder = "polynomial order.",
            dorder = "derivative order."
         ),
         info = "Savitzky-Golay filter."
      ),

      # prep.alsbasecorr <- function(spectra, plambda = 5, p = 0.1, max.niter = 10)
      "alsbasecorr" = list(
         name = "alsbasecorr",
         method = prep.alsbasecorr,
         params = list(plambda = 5, p = 0.1, max.niter = 10),
         params.show = c("plambda", "p"),
         jmethod = prep.alsbasecorr.asjson,
         params.info = list(
            plambda = "power of the penalty parameter (e.g. if plambda = 5, lambda = 10^5)",
            p = "asymmetry ratio (should be between 0 and 1)",
            max.niter = "maximum number of iterations"
         ),
         info = "Asymmetric least squares baseline correction."
      ),

      "spikes" = list(
         name = "spikes",
         method = prep.spikes,
         params = list(width = 5, threshold = 6),
         params.show = c("width", "threshold"),
         jmethod = prep.spikes.asjson,
         params.info = list(
            width = "width of the moving median filter",
            threshold = "threshold for spike detection"
         ),
         info = "Removes cosmic spikes."
      ),

      # prep.center <- function(data, type = "mean", center = NULL)
      "center" = list(
         name = "center",
         method = prep.center,
         params = list(type = "mean"),
         params.show = c("type"),
         pmethod = prep.center.params,
         jmethod = prep.center.asjson,
         params.info = list(
            type = "what to use for centering ('mean', or 'median')"
         ),
         info = "Center dataset columns."
      ),

      "scale" = list(
         name = "scale",
         method = prep.scale,
         params = list(type = "sd"),
         params.show = c("type"),
         pmethod = prep.scale.params,
         jmethod = prep.scale.asjson,
         params.info = list(
            type = "what to use for scaling ('sd', 'iqr', 'range', or 'pareto')"
         ),
         info = "Scale dataset columns."
      ),

      # prep.msc <- function(spectra, mspectrum = NULL)
      "emsc" = list(
         name = "emsc",
         method = prep.emsc,
         params = list(degree = 0, mspectrum = NULL),
         params.show = c("degree"),
         pmethod = prep.emsc.params,
         jmethod = prep.emsc.asjson,
         params.info = list(
            degree = "polynomial degree (0 for MSC)",
            mspectrum = "reference spectrum (if NULL mean spectrum will be used)."
         ),
         info = "Multiplicative scatter correction."
      )
   )
}


#' Shows information about all implemented preprocessing methods.
#'
#' @export
prep.list <- function() {
   p <- getImplementedPrepMethods()

   cat("\n\nList of available preprocessing methods:\n")
   lapply(p, function(c) {
      cat("\n\n")
      fprintf(" %s\n", c$info)
      cat(" ---------------\n")
      fprintf(" name: '%s'\n", c$name)

      if (length(c$params.info) == 0) {
         cat(" no parameters required\n")
      } else {
         cat(" parameters:\n")
         for (i in seq_along(c$params.info)) {
            fprintf("  '%s': %s\n", names(c$params.info)[i], c$params.info[[i]])
         }
      }
   })

   invisible()
}


#' Class for preprocessing object/item.
#'
#' @param name
#' short text with name for the preprocessing method.
#' @param ...
#' a list with named parameters for the method (if empty - default parameters will be used).
#'
#' @details
#' Use this class to create a list with a sequence of preprocessing methods to keep them together
#' in right order and with defined parameters. The list/object can be provided as an extra argument
#' to any modelling function (e.g. \code{pca}, \code{pls}, etc), so the optimal model parameters and
#' the optimal preprocessing will be stored together and can be applied to a raw data by using
#' method \code{predict}.
#'
#' For your own preprocessing method you need to create a function, which takes matrix with values
#' (dataset) as the first argument, does something and then return a matrix with the same dimension
#' and same attributes as the result. The method can have any number of optional parameters.
#'
#' See Bookdown tutorial for details.
#'
#' @export
prep <- function(name, ...) {

   params <- list(...)
   pmethod <- NULL
   jmethod <- NULL

   # 1. first check name
   item <- getImplementedPrepMethods()[[name]]
   stopifnot("prep: either name of preprocessing method is wrong or you need to provide a reference to
         function implementing the method if it is user defined." = !is.null(item))

   # 2. check the parameters
   if (length(params) == 0) params <- item$params
   if (length(params) > 0 && !all(names(params) %in% names(item$params))) {
      stop("prep: provided preprocessing parameters have wrong name.", call. = FALSE)
   }

   method <- item$method
   pmethod <- item$pmethod
   jmethod <- item$jmethod

   obj <- list(
      name = name,
      method = method,
      params = params,
      pmethod = pmethod,
      jmethod = jmethod
   )

   class(obj) <- c("prep")
   return(obj)
}


#' Fits preprocessing model
#'
#' @param obj
#' list with preprocessing methods (created using \code{prep} or \code{prep.fit} function).
#' @param x
#' matrix with training set to be used for computing data dependent parameters
#'
#' @return same list but with updated methods parameters computed based on
#' the training set.
#'
#' @export
prep.fit <- function(obj, x) {

   stopifnot("prep.fit: the first argument must be a list with preprocessing methods" =
      is.list(obj) && inherits(obj[[1]], "prep"))
   stopifnot("prep.fit: argument 'x' must be a matrix" =
      !is.null(x) && is.matrix(x))

   # remove excluded rows
   if (!is.null(attr(x, "exclrows"))) x = mda.purgeRows((x))

   npred <- ncol(x)
   out <- list()
   i <- 1
   for (p in obj) {
      if (!is.null(p[["pmethod"]])) {
         p[["params"]] <- do.call(p[["pmethod"]], c(list(data = x), p[["params"]]))
      }
      x <- do.call(p[["method"]], c(list(data = x), p[["params"]]))
      out[[i]] <- p
      i <- i + 1
   }
   out[["_npred"]] <- npred
   class(out) <- c("prepmodel")
   return(out)
}


#' Applies a list with preprocessing methods to a dataset
#'
#' @param obj
#' list with preprocessing methods (created using \code{prep} or \code{prep.fit} function).
#' @param x
#' matrix with dataset
#'
#' @export
prep.apply <- function(obj, x) {
   stopifnot("prep.apply: the first argument must be a list with preprocessing methods" =
      is.list(obj) && inherits(obj[[1]], "prep"))
   stopifnot("prep.apply: argument 'x' must be a matrix" =
      !is.null(x) && is.matrix(x))

   for (p in obj) {
      if (!is.list(p)) next
      x <- do.call(p[["method"]], c(list(data = x), p[["params"]]))
   }
   return(x)
}



#' Converts preprocessing model to JSON elements.
#'
#' @param obj
#' list with preprocessing methods (created using \code{prep.fit} function).
#'
#' @returns stringified JSON.
#'
#' @export
prep.asjson <- function(obj) {

   npred <- obj[["_npred"]]
   if (npred < 1) stop("prep.asjson: preprocessing object does not contain information about number of predictors.", call. = FALSE)

   left <- 0
   right <- npred
   left.tot <- 0

   ml <- NULL
   mp <- NULL
   mpl <- NULL
   mp_new <- NULL
   mpl_new <- NULL
   info <- rep("\'\'", 8)

   n <- 1
   scale.flag <- FALSE
   for (p in obj) {
      if (!is.list(p)) next
      if (is.null(p[["jmethod"]])) stop("prep.asjson: preprocessing list contains method, which can not be converted to JSON.", call. = FALSE)

      out <- do.call(p[["jmethod"]], list(params = p[["params"]], npred = npred, left = left, right = right))

      if (p[["name"]] == "scale" && scale.flag == TRUE) {
         # so we got scaling after centering, in this case we combine them together
         # because in web-app it is a single method

         scale.flag <- FALSE

         # identify number of variables if tails were trimmed
         # npred - total number of variables
         # nvar - number of variables after trimming
         # npad - number of elements to pad (the trimmed elements)
         l <- length(mp)
         mp[l] <- out$mp[2]
         nvar <- right - left
         npad <- (npred - nvar)

         # identify the start and end indices inside mp_new
         # where the scaling vector should be placed (right after centering vector)
         s <- length(mp_new) - npred - npad + 1
         e <- s + nvar - 1

         # save the vector with scaling values to the proper place
         mp_new[s:e] <- out$mp_new[(nvar + 1):(2 * nvar)]

         # identify location of "/" symbol in centering info, e.g. "'mean/none'"
         s1 <- regexec("/", info[n - 1])[[1]]
         s2 <- regexec("/", out[["info"]])[[1]]
         n2 <- nchar(out[["info"]])

         # select all before "/" and combine with scaling type
         info[n - 1] <- paste0(substring(info[n - 1], 1, s1), substring(out[["info"]], s2 + 1, n2), "'")

         # skip the code below
         next
      }

      # if method is centering, we set flag just in case if scaling is next
      # so we can combine them together
      if (p[["name"]] == "center") scale.flag <- TRUE else scale.flag <- FALSE

      ml <- c(ml, out[["ml"]])
      mp <- c(mp, out[["mp"]])
      mpl <- c(mpl, out[["mpl"]])
      mp_new <- c(mp_new, out[["mp_new"]])
      mpl_new <- c(mpl_new, out[["mpl_new"]])
      info[n] <- paste0("\'", out[["info"]], "\'")

      if (out[["ml"]] == 4) {
         # it is trim tail method, modify left and right
         left <- out[["mp"]][1]
         right <- out[["mp"]][2]
         left.tot <- left.tot + left
      }

      n <- n + 1
   }

   m <- paste0(
      "{'ml':{'__type':'Int32Array','data':[",
      paste0(ml, collapse = ","),
      "]}, 'mpl':{'__type':'Int32Array', 'data':[",
      paste0(mpl, collapse = ","),
      "]}, 'mp':{'__type':'Float64Array','data':[",
      paste0(mp, collapse = ","),
      "]}, 'mpl_new':{'__type':'Int32Array','data':[",
      paste0(mpl_new, collapse = ","),
      "]}, 'mp_new':{'__type':'Float64Array','data':[",
      paste0(format(mp_new, digits=14), collapse = ","),
      "]}, 'npred': ", npred, ", 'class': ['prepmodel'],'info': [",
      paste0(info, collapse = ","),
      "]}"
   )

   m <- gsub("\'", "\"", m)
   return(m)
}


#' Converts JSON string to preprocessing model
#'
#' @param str
#' string with JSON
#'
#' @return list with the methods.
prep.fromjson <- function(str) {

   class <- extractStringArray(str, "class")
   if (is.null(class) || length(class) != 1 || !("prepmodel" %in% class)) {
      stop("Selected JSON file does not contain a preprocessing model.", call. = FALSE)
   }

   # list of methods for conversion of every item
   methods <- list(
      prep.norm.fromjson,
      prep.savgol.fromjson,
      prep.center.fromjson,
      list(prep.emsc.fromjson, prep.alsbasecorr.fromjson), # a list of two because in webapp it is one method (baseline)
      prep.varsel.fromjson,
      prep.spikes.fromjson
   )

   # extract values and array from the JSON
   npred <- extractValue(str, "npred")
   ml <- extractArray(str, "ml")
   mp <- extractArray(str, "mp")
   mpl <- extractArray(str, "mpl")
   mp_new <- extractArray(str, "mp_new")
   mpl_new <- extractArray(str, "mpl_new")

   left <- 0
   right <- npred
   left.tot <- 0

   # number of methods and empty list for preprocessing model
   n <- length(ml)
   m <- list()

   # offsets to locate part of "mp" and "mp_new" related to each method
   offset_mp_new <- 1
   offset_mp <- 1

   # we need two counters because some of the JSON methods will end up
   # with two R methods (like "baseline" and "scale")
   im <- 1
   ip <- 1

   # we need to remember that in JSON all indices, like indices for methods
   # start with 0 while in R they start with 1, hence we add 1 to them
   while (ip <= n) {

      # subset local values of "mp" and "mp_new" for the i-th method
      mp_new_i <- mp_new[offset_mp_new:(offset_mp_new + mpl_new[ip] - 1)]
      mp_i <- mp[offset_mp:(offset_mp + mpl[ip] - 1)]

      # if it is "scale" method we need special treatment
      if (ml[ip] == 2) {

         # if it contains both centering and scaling we need
         # to create two separate R methods
         if (mp_i[1] != 0 && mp_i[2] != 0) {
            params <- list(mp = mp_i, mp_new = mp_new_i, left = left, right = right, left.tot = left.tot)
            m[[im]] <- do.call(prep.center.fromjson, params)
            im <- im + 1
            m[[im]] <- do.call(prep.scale.fromjson, params)
            offset_mp_new <- offset_mp_new + mpl_new[ip]
            offset_mp <- offset_mp + mpl[ip]
            im <- im + 1
            ip <- ip + 1
            next
         } else {
            f <- if (mp_i[1] != 0) prep.center.fromjson else prep.scale.fromjson
         }
      } else if (ml[ip] == 3) {
         # this is baseline method, depending on mp[1] value it can be emsc or als
         f <- methods[[ml[ip] + 1]][[mp_i[1] + 1]]
      } else {
         f <- methods[[ml[ip] + 1]]
      }

      # call fromjson method to get all parameters
      m[[im]] <- do.call(f, list(mp = mp_i, mp_new = mp_new_i, left = left, right = right, left.tot = left.tot))

      # move offset to the next method
      offset_mp_new <- offset_mp_new + mpl_new[ip]
      offset_mp <- offset_mp + mpl[ip]

      if (ml[ip] == 4) {
         left.tot <- mp_i[1]
         right <- mp_i[2] - left
         left <-  mp_i[1] - left
      }

      im <- im + 1
      ip <- ip + 1
   }

   m[["_npred"]] <- npred
   class(m) <- c("prepmodel")
   return(m)
}


#' Saves preprocessing model to JSON file which can be loaded to web-application (mda.tools/prep).
#'
#' @param obj
#' list with preprocessing methods (created using \code{prep.fit} function).
#' @param fileName
#' file name (or full path) to JSON file to save the model into.
#'
#' @export
writeJSON.prepmodel <- function(obj, fileName) {
   m <- prep.asjson(obj)
   fileConn <- file(fileName)
   writeLines(m, fileConn)
   close(fileConn)
}


#' Show summary of the preprocessing model.
#'
#' @param object
#' preprocessing model (created by \code{\link{prep.fit}}).
#' @param ...
#' potential further arguments (required for Method/Generic reasons).
#'
#' @return
#' the \code{object} argument (invisibly).
#'
#' @export
summary.prepmodel <- function(object, ...) {
   cat("\nPreprocessing model:\n")
   print(object)
   cat("\n")
   invisible(object)
}


#' Print the information about methods in the preprocessing model.
#'
#' @param x
#' preprocessing model (created by \code{\link{prep.fit}}).
#' @param ...
#' potential further arguments (required for Method/Generic reasons).
#'
#' @return
#' the \code{x} argument (invisibly).
#'
#' @export
print.prepmodel <- function(x, ...) {
   prep.list <- getImplementedPrepMethods()

   par2str <- function(n, p) {
      if (is.null(p)) return("")
      if (is.vector(p) && length(p) > 1) return(paste0(n, " = ", min(p), ":", max(p)))
      return(paste0(n, " = ", p))
   }

   for (p in x) {
      if (!is.list(p)) next
      name <- p[["name"]]
      params <- p[["params"]]
      params.list <- prep.list[[name]][["params.show"]]
      out <- paste0(" - ", name)
      params.out <- sapply(params.list, function(n) par2str(n, params[[n]]))
      if (length(params.out) > 0)
         out <- paste0(out, ": ", paste0(params.out, collapse = ", "))
      cat(out, "\n")
   }
   invisible(x)
}

############################################################
# Service methods                                          #
############################################################

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
   tcrossprod(s$v %*% diag(1 / pmax(s$d, .Machine$double.eps), length(s$d), length(s$d)), s$u)
}



#################################################################################
# Legacy methods - still supported but no development and not part of prep.     #
#################################################################################

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
#' @details
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
      if (is.logical(center) && center) center <- colMeans(data)

      if (is.numeric(center) && length(center) != ncol(data)) {
         stop("Number of values in 'center' should be the same as number of columns in 'data'", call. = FALSE)
      }

      # define values for weighting
      if (is.logical(scale) && scale) scale <- apply(data, 2, sd)

      if (is.numeric(scale) && length(scale) != ncol(data)) {
         stop("Number of values in 'scale' should be the same as number of columns in 'data'", call. = FALSE)
      }

      # compute coefficient of variation and set scale to 1 if it is below
      # a user defined threshold
      if (is.numeric(scale)) {
         m <- if (is.numeric(center)) center else colMeans(data)
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
#' @examples
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
   # .Deprecated("prep.norm")
   return(prep.norm(data, type = "snv"))
}


#' Multiplicative Scatter Correction transformation
#'
#' @description
#' Applies Multiplicative Scatter Correction (MSC) transformation to data matrix (spectra)
#'
#' @param data
#' a matrix with data values (spectra)
#' @param mspectrum
#' mean spectrum (if NULL will be calculated from \code{spectra})
#' @param ...
#' other optional components
#'
#' @return
#' preprocessed spectra (calculated mean spectrum is assigned as attribute 'mspectrum')
#'
#' @details
#' MSC is used to remove scatter effects (baseline offset and slope) from
#' spectral data, e.g. NIR spectra.
#'
#' @examples
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
prep.msc <- function(data, mspectrum = NULL, ...) {
   # .Deprecated("prep.emsc")
   return(prep.emsc(data, degree = 0, mspectrum = mspectrum, ...))
}



#' Applies a list with preprocessing methods to a dataset
#'
#' @param obj
#' list with preprocessing methods (created using \code{prep} function).
#' @param x
#' matrix with dataset
#' @param ...
#' other arguments
#'
#' @export
employ.prep <- function(obj, x, ...) {
   # .Deprecated("prep.apply")
   prep.apply(obj, x)
}
