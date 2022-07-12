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
#' Normalizes signals (rows of data matrix).
#'
#' @param data
#' a matrix with data values
#' @param type
#' type of normalization \code{"area"}, \code{"length"}, \code{"sum"}, \code{"snv"}, \code{"is"}, or \code{"pqn"}.
#' @param col.ind
#' indices of columns (can be either integer or logical valuws) for normalization to internal
#' standard peak.
#' @param ref.spectrum
#' reference spectrum for PQN normalization, if not provided a mean spectrum for data is used
#'
#' @details
#' The \code{"area"}, \code{"length"}, \code{"sum"} types do preprocessing to unit area (sum of
#' absolute values), length or sum of all values in every row of data matrix. Type \code{"snv"}
#' does the Standard Normal Variate normalization, similar to \code{\link{prep.snv}}. Type
#' \code{"is"} does the normalization to internal standard peak, whose position is defined by
#' parameter `col.ind`. If the position is a single value, the rows are normalized to the height
#' of this peak. If `col.ind` points on several adjucent vales, the rows are normalized to the area
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

   if (type == "snv") return(prep.snv(data))

   if (type == "is" && is.null(col.ind)) {
      stop("For 'is' normalization you need to provide indices for IS peak.")
   }

   if (is.logical(col.ind)) {
      col.ind <- which(col.ind)
   }

   if (!is.null(col.ind) && (min(col.ind) < 1 || max(col.ind) > ncol(data))) {
      stop("Values for 'col.ind' seem to be wrong.")
   }

   if (type == "pqn" && is.null(ref.spectrum)) {
      ref.spectrum <- apply(data, 2, mean)
   }

   pqn <- function(data, ref.spectrum) {

      if (length(ref.spectrum) != ncol(data)) {
         stop("prep.norm: 'ref.spectrum' should have the same number of values as the number of columns in 'data'.")
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

      w <- switch(
         type,
         "area" = apply(abs(data), 1, sum),
         "length" = sqrt(apply(data^2, 1, sum)),
         "sum" = apply(data, 1, sum),
         "is" = apply(data[, col.ind, drop = FALSE], 1, sum),
         "pqn" = pqn(data, ref.spectrum)
      )

      if (is.null(w)) stop("Wrong value for argument 'type'.")
      return(sweep(data, 1, w, "/"))
   }

   return(prep.generic(data, f, type = type, col.ind = col.ind, ref.spectrum = ref.spectrum))
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
#' @details
#' The function implements algorithm described in [1] which handles the edge points correctly and
#' does not require to cut the spectra.
#'
#' @references
#' 1. Peter A. Gorry. General least-squares smoothing and differentiation by the convolution
#' (Savitzky-Golay) method. Anal. Chem. 1990, 62, 6, 570–573, https://doi.org/10.1021/ac00205a007.
#'
#' @export
prep.savgol <- function(data, width = 3, porder = 1, dorder = 0) {

   stopifnot("Filter width ('width') should be equal at least to 3." = width > 2)
   stopifnot("Filter width ('width') should be an odd integer number." = width %% 2 == 1)
   stopifnot("Wrong value for the derivative order (should be 0, 1, or 2)." = dorder %in% (0:2))
   stopifnot("Wrong value for the polynomial degree (should be integer number between 0 and 4)." = porder %in% (0:4))
   stopifnot("Polynomal degree ('porder') should not be smaller the derivative order ('dorder')." = porder >= dorder)

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

   f <- function(x, width, porder, dorder) {

      nvar <- ncol(x)
      px <- matrix(0.0, nrow(x), ncol(x))
      m <- round((width - 1) / 2)
      w <- outer(-m:m, -m:m, function(x, y) weight(x, y, m, porder, dorder))

      for (i in seq_len(m)) {
         px[, i] <- apply(x[, seq_len(2 * m + 1), drop = FALSE], 1,
            function(xx) convolve(xx, w[, i], type = "filter")[1])
         px[, nvar - i + 1] <- apply(x[, (nvar - 2 * m):nvar, drop = FALSE], 1,
            function(xx) convolve(xx, w[, width - i + 1], type = "filter")[1])
      }

      px[, (m + 1):(nvar - m)] <- t(apply(x, 1, function(xx) convolve(xx, w[, m + 1], type = "filter")))

      return(px)
   }

   return(prep.generic(data, f, width = width, porder = porder, dorder = dorder))
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
prep.msc <- function(data, mspectrum = NULL) {

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

   return(prep.generic(data, f, mspectrum = mspectrum))
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

#' Baseline correction using asymetric least squares
#'
#' @param data
#' matrix with spectra (rows correspond to individual spectra)
#' @param plambda
#' power of the penalty parameter (e.g. if plambda = 5, lambda = 10^5)
#' @param p
#' assymetry ratio (should be between 0 and 1)
#' @param max.niter
#' maximum number of iterations
#'
#' @return
#' preprocessed spectra.
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
#' @importFrom methods as
#'
#' @export
prep.alsbasecorr <- function(data, plambda = 5, p = 0.1, max.niter = 10) {
   attrs <- mda.getattr(data)
   dimnames <- dimnames(data)

   m <- ncol(data)
   baseline <- matrix(0, nrow(data), ncol(data))
   LDD <- Matrix::Matrix((10^plambda) * crossprod(diff(diag(m), difference = 2)), sparse = TRUE)
   w.ini <- matrix(rep(1, m))

   for (i in seq_len(nrow(data))) {
      y <- data[i, ]
      w <- w.ini
      for (j in seq_len(max.niter)) {
         W <- Matrix::Diagonal(x = as.numeric(w))
         z <- Matrix::solve(as(W + LDD, "dgCMatrix"), w * y, sparse = TRUE)
         w.old <- w
         w <- p * (y > z) + (1 - p) * (y < z)

         if (sum(abs(w - w.old)) < 10^-10) break
      }

      baseline[i, ] <- as.numeric(z)
   }

   pspectra <- data - baseline
   pspectra <- mda.setattr(pspectra, attrs)
   dimnames(pspectra) <- dimnames

   return(pspectra)
}


#' Transformation
#'
#' @description
#' Transforms values from using any mathematical function (e.g. log).
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
#' indices of variables (columns) to select, can bet either numeric or logical
#'
#' @return
#' data matrix with the selected variables (columns)
#'
#' @export
prep.varsel <- function(data, var.ind) {
   if (!is.null(attr(data, "exclcols"))) {
      stop("prep.varsel() can not be used for dataset with excluded (hidden) columns.")
   }
   return(mda.subset(data, select = var.ind))
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


#' Shows a list with implemented preprocessing methods
#'
#' @export
getImplementedPrepMethods <- function() {
   list(
      # prep.snv <- function(data)
      "snv" = list(
         name = "snv",
         method = prep.snv,
         params = list(),
         params.info = list(),
         info = "Standard normal variate normalization."
      ),

      # prep.ref2km <- function(spectra)
      "ref2km" = list(
         name = "ref2km",
         method = prep.ref2km,
         params = list(),
         params.info = list(),
         info = "Transform reflectance spectra to Kubelka-Munk units."
      ),

      # prep.msc <- function(spectra, mspectrum = NULL)
      "msc" = list(
         name = "msc",
         method = prep.msc,
         params = list(mspectrum = NULL),
         params.info = list(mspectrum = "mean spectrum (if NULL will be computed based on data)."),
         info = "Multiplicative scatter correction."
      ),

      # prep.transform <- function(data, fun, ...)
      "transform" = list(
         name = "transform",
         method = prep.transform,
         params = list(fun = log),
         params.info = list(fun = "function to transform the values (e.g. 'log')"),
         info = "Transformation of data values using math functions (log, sqrt, etc.)."
      ),

      # prep.varsel <- function(data, var.ind)
      "varsel" = list(
         name = "varsel",
         method = prep.varsel,
         params = list(var.ind = NULL),
         params.info = list(var.ind = "indices of variables (columns) to select."),
         info = "Select user-defined variables (columns of dataset)."
      ),

      # prep.norm <- function(data, type = "area", col.ind = NULL)
      "norm" = list(
         name = "norm",
         method = prep.norm,
         params = list(type = "area", col.ind = NULL),
         params.info = list(
            type = "type of normalization ('area', 'sum', 'length', 'is', 'snv', 'pqn').",
            col.ind = "indices of columns (variables) for normalization to internal standard peak."
         ),
         info = "Normalization."
      ),

      # prep.savgol <- function(data, width = 3, porder = 1, dorder = 0)
      "savgol" = list(
         name = "savgol",
         method = prep.savgol,
         params = list(width = 3, porder = 1, dorder = 0),
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
         params.info = list(
            plambda = "power of the penalty parameter (e.g. if plambda = 5, lambda = 10^5)",
            p = "assymetry ratio (should be between 0 and 1)",
            max.niter = "maximum number of iterations"
         ),
         info = "Asymmetric least squares baseline correction."
      ),

      # prep.autoscale <- function(data, center = TRUE, scale = FALSE, max.cov = 0)
      "autoscale" = list(
         name = "autoscale",
         method = prep.autoscale,
         params = list(center = TRUE, scale = FALSE, max.cov = 0),
         params.info = list(
            center = "a logical value or vector with numbers for centering.",
            scale = "a logical value or vector with numbers for weighting.",
            max.cov = "columns with coefficient of variation (in percent) below `max.cov` will not be scaled"
         )
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


#' Class for preprocessing object
#'
#' @param name
#' short text with name for the preprocessing method.
#' @param params
#' a list with parameters for the method (if NULL - default parameters will be used).
#' @param method
#' method to call when applying the preprocessing, provide it only for user defined methods.
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
#' and same attributes as the result. The method can have any number of optional  parameters.
#'
#' See Bookdown tutorial for details.
#'
#' @export
prep <- function(name, params = NULL, method = NULL) {

   if (!is.null(params) && !is.list(params)) {
      stop("prep: argument 'params' must be a list with parameter names and values.")
   }

   if (is.null(method)) {
      # assuming it is one of the standard constraints
      # 1. first check name
      item <- getImplementedPrepMethods()[[name]]
      stopifnot("prep: either name of preprocessing method is wrong or you need to provide a reference to
         function implementing the method if it is user defined." = !is.null(item))

      # 2. check the parameters
      if (is.null(params)) params <- item$params
      if (length(params) > 0 && !all(names(params) %in% names(item$params))) {
         stop("prep: provided preprocessing parameters have wrong name.")
      }

      method <- item$method
   } else {
      # user defined method, check that it works
      res <- tryCatch(
         do.call(method, c(list(data = matrix(runif(50, 1, 10), 5, 10)), params)),
         error = function(m) stop("prep: the method you provided raises an error: \n", m),
         warning = function(m) stop("prep: the method you provided raises a warning: \n", m)
      )

      stopifnot("prep: the method you provided does not return matrix with correct dimension." =
         dim(res) == c(5, 10))
   }

   obj <- list(
      name = name,
      method = method,
      params = params
   )

   class(obj) <- c("prep")
   return(obj)
}

#' Applies a list with preprocessing methods to a dataset
#'
#' @param obj
#' list with preprocssing methods (created using \code{prep} function).
#' @param x
#' matrix with dataset
#' @param ...
#' other arguments
#'
#' @export
employ.prep <- function(obj, x, ...) {

   stopifnot("employ.prep: the first argument must be a list with preprocessing methods" =
      is.list(obj) && class(obj[[1]]) == "prep")
   for (p in obj) {
      x <- do.call(p$method, c(list(data = x), p$params))
   }

   return(x)
}
