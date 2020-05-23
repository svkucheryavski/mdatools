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
mcrpure <- function(x, ncomp, exclrows = NULL, exclcols = NULL, offset = offset, savgol = NULL, info = "") {

   x <- prepCalData(x, exclrows, exclcols)

   # get pure variables and unmix data
   model <- getPureVariables(x, ncomp, offset = offset, exclude = exclude)
   res2 = unmix(spectra, D[, res1$purevars, drop = F], by.spec = by.spec)

   class(model) <- c("mcr", "mcrpure")
   return(model)
}


getPureSpectra.mcrpure <- function() {

  # combine all values to the final list
  res = c(res1, res2)
  res$by.spec = by.spec
  res$savgol = savgol
}

#' Identifies pure variables
#'
#' @description
#' The method identifies indices of pure variables using the SIMPLISMA
#' algorithm.
#'
#' @param x
#' matrix with the spectra or their inverted second derivative
#' @param ncomp
#' number of pure components
#' @param offset
#' offset in percent for calculation of parameter alpha
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
getPureVariables <- function(x, ncomp, offset = 5, exclude = NULL) {
  nspec = nrow(D)
  nvar = ncol(D)

  # get indices for excluded columns if provided
  colind = rep(TRUE, nvar)
  if (!is.null(exclude))
  {
    if (max(exclude)> nvar || min(exclude) < 1)
      stop('Wrong values for the parameter "exclude"!')
    colind[exclude] = FALSE
  }


  # calculate purity spectrum
  mu = apply(D, 2, mean)
  sigma = apply(D, 2, sd) * sqrt((nspec - 1) / nspec)
  alpha = offset / 100 * max(mu)
  purity = sigma / (mu + alpha)

  # calculate variance-covariance matrix
  Dprime = sweep(D, 2, (sqrt(mu^2 + (sigma + alpha)^2)), '/')
  R = 1/nrow(Dprime) * (t(Dprime) %*% Dprime)

  # loop to locate the pure variables variables
  purityspec = matrix(0, nrow = ncomp, ncol = nvar)
  stdspec = matrix(0, nrow = ncomp, ncol = nvar)
  purevars = rep(0, ncomp)
  purevals = rep(0, ncomp)
  for (i in 1:ncomp)
  {
    weights = matrix(0, ncol = ncol(D), nrow = 1)
    for (j in 1:ncol(D))
    {
      weights[j] = det(R[c(j, purevars[1:(i-1)]), c(j, purevars[1:(i-1)]), drop = F]);
    }

    purityspec[i, colind] = purity[colind] * weights
    stdspec[i, colind] = sigma * weights
    purevars[i] = which.max(purityspec[i, ])
    purevals[i] = purityspec[i, purevars[i]]
  }

  res = list()
  res$ncomp = ncomp
  res$purity = purevals
  res$purityspec = purityspec
  res$stdspec = stdspec
  res$purevars = purevars

  res
}

#' Unmix spectral data with least squares method
#'
#' @description
#' \code{unmix} decomposes the spectral data to a matrix with concentrations and
#' a matrix with spectra of pure components using set of pure variables.
#'
#' @param D
#' matrix with the spectral data
#' @param Dr
#' a matrix with concentration estimates (e.g. columns of D with pure variables)
#' @param by.spec
#' logical, should the algorithm works with spectra (rows) or with concentrations (columns)
#'
#' @return
#' Returns a list with two matrices, `conc` for concentrations and `spec` for spectra as well
#' as a vector with residual variance values
#'
#' @export
unmix = function(D, Dr, by.spec = T)
{

  # resolve spectra and concentrations
  St = t(D) %*% Dr %*% solve(t(Dr) %*% Dr)
  Ct = D %*% St %*% solve(t(St) %*% St)

  # scale
  #if (savgol[1] > 0)
  #  f = as.matrix(rowSums(deriv))
  #else
  #  f = as.matrix(rowSums(D))
  f = as.matrix(rowSums(D))
  a = solve(t(Ct) %*% Ct) %*% t(Ct) %*% f

  if (length(a) == 1)
    A = a
  else
    A = diag(as.vector(a))

  Ct = Ct %*% A
  St = St %*% A %*% solve(t(A) %*% A)

  # calculate residual variance
  ncomp = ncol(Ct)
  resvar = rep(0, ncomp)
  for (i in 1:ncomp)
  {
    exp = Ct[, 1:i, drop = F] %*% t(St[, 1:i, drop = F])
    res = D - exp
    resvar[i] = sum(res^2) / sum(D^2)
  }

  # add names
  names(resvar) = colnames(Ct) = colnames(St) = paste('C', 1:ncomp, sep = '')
  rownames(Ct) = rownames(D)
  rownames(St) = colnames(D)

  # switch spectra and concentrations if unmixing was in contribution space
  if (by.spec)
  {
    spec = St
    conc = Ct
  }
  else
  {
    spec = Ct
    conc = St
  }

  # add attributes for hyperspectral image if any
  attr(conc, 'width') = attr(D, 'width')
  attr(conc, 'height') = attr(D, 'height')

  # return the results
  res = list(
    conc = conc,
    spec = spec,
    resvar = resvar
  )

  res
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


