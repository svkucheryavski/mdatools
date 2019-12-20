#' Principal Component Analysis
#'
#' @description
#' \code{pca} is used to build and explore a principal component analysis (PCA) model.
#'
#' @param x
#' calibration data (matrix or data frame).
#' @param ncomp
#' maximum number of components to calculate.
#' @param center
#' logical, do mean centering of data or not.
#' @param scale
#' logical, do standardization of data or not.
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values)
#' @param x.test
#' test data (matrix or data frame).
#' @param method
#' method to compute principal components ("svd", "nipals").
#' @param rand
#' vector with parameters for randomized PCA methods (if NULL, conventional PCA is used instead)
#' @param lim.type
#' which method to use for calculation of critical limits for residual distances (see details)
#' @param alpha
#' significance level for extreme limits for T2 and Q disances.
#' @param gamma
#' significance level for outlier limits for T2 and Q distances.
#' @param info
#' a short text with model description.
#'
#' @details
#'
#' Note, that from v. 0.10.0 cross-validation as well as "jm" method for computing residual
#' distance critical limits are not more supported in PCA.
#'
#' If number of components is not specified, a minimum of number of objects - 1 and number of
#' variables in calibration set is used. One can also specified an optimal number of component,
#' once model is calibrated (\code{ncomp.selected}). The optimal number of components is used to
#' build a residuals distance plot, as well as for SIMCA classification.
#'
#' If some of rows of calibration set should be excluded from calculations (e.g. because they are
#' outliers) you can provide row numbers, names, or logical values as parameter \code{exclrows}. In
#' this case they will be completely ignored we model is calibrated. However, score and residuls
#' distances will be computed for these rows as well and then hidden. You can show them
#' on corresponding plots by using parameter \code{show.excluded = TRUE}.
#'
#' It is also possible to exclude selected columns from calculations by provideing parameter
#' \code{exclcols} in form of column numbers, names or logical values. In this case loading matrix
#' will have zeros for these columns. This allows to compute PCA models for selected variables
#' without removing them physically from a dataset.
#'
#' Take into account that if you see other packages to make plots (e.g. ggplot2) you will
#' not be able to distinguish between hidden and normal objects.
#'
#' By default loadings are computed for the original dataset using either SVD or NIPALS algorithm.
#' However, for datasets with large number of rows (e.g. hyperspectral images), there is a
#' possibility to run algorithms based on random permutations [1, 2]. In this case you have
#' to define parameter \code{rand} as a vector with two values: \code{p} - oversampling parameter
#' and \code{k} - number of iterations. Usually \code{rand = c(15, 0)} or  \code{rand = c(5, 1)}
#' are good options, which give quite almost precise solution but much faster.
#'
#' There are several ways to calculate critical limits for orthogonal (Q, q) and score (T2, h)
#' distances. In \code{mdatools} you can specify one of the following methods via parameter
#' \code{lim.type}: \code{"chisq"} - method based on chi-square distribution [4],
#' \code{"ddmoments"} and \code{"ddrobust"} - related to data driven method proposed by Pomerantsev
#' and Rodionova [5]. The \code{"ddmoments"} is based on method of moments for estimation of
#' distribution parameters (also known as "classical" approach) while \code{"ddrobust"} is
#' based in robust estimation.
#'
#' If \code{lim.type="chisq"} is used, only limits for Q-distances are computed based on chi-square
#' distribution. Corresponding limits for T2-distances are computed using Hotelling's T-squared
#' distribution. The methods utilizing the data driven approach calculate limits for both
#' distances bases on chi-square distribution and parameters estimated from the calibration data.
#'
#' The critical limits are calculated for a significance level defined by parameter \code{'alpha'}.
#' You can also specify another parameter, \code{'gamma'}, which is used to calculate acceptance
#' limit for outliers (shown as dashed line on residual distance plot).
#'
#' You can also recalculate the limits for existent model by using different values for alpha and
#' gamme, without recomputing the model itself. In this case use the following code (it is assumed
#' that you current PCA/SIMCA model is stored in variable \code{m}):
#' \code{m = setDistanceLimits(m, lim.type, alpha, gamma)}.
#'
#' In case of PCA the critical limits are just shown on residual plot as lines and can be used for
#' detection of extreme objects (solid line) and outliers (dashed line). When PCA model is used for
#' classification in SIMCA (see \code{\link{simca}}) the limits are also employed for
#' classification of objects.
#'
#' @return
#' Returns an object of \code{pca} class with following fields:
#' \item{ncomp }{number of components included to the model.}
#' \item{ncomp.selected }{selected (optimal) number of components.}
#' \item{loadings }{matrix with loading values (nvar x ncomp).}
#' \item{eigenvals }{vector with eigenvalues for all existent components.}
#' \item{expvar }{vector with explained variance for each component (in percent).}
#' \item{cumexpvar }{vector with cumulative explained variance for each component (in percent).}
#' \item{T2lim }{statistical limit for T2 distance.}
#' \item{Qlim }{statistical limit for Q residuals.}
#' \item{info }{information about the model, provided by user when build the model.}
#' \item{calres }{an object of class \code{\link{pcares}} with PCA results for a calibration data.}
#' \item{testres }{an object of class \code{\link{pcares}} with PCA results for a test data, if it
#' was provided.}
#'
#' More details and examples can be found in the Bookdown tutorial.
#'
#' @references
#' 1. N. Halko, P.G. Martinsson, J.A. Tropp. Finding structure with randomness: probabilistic
#' algorithms for constructing approximate matrix decompositions. SIAM Review, 53 (2010) pp.
#' 217-288.
#' 2. S. Kucheryavskiy, Blessing of randomness against the curse of dimensionality,
#' Journal of Chemometrics, 32 (2018).
#' 3. A.L. Pomerantsev, Acceptance areas for multivariate classification derived by projection
#' methods, Journal of Chemometrics, 22 (2008) pp. 601-609.
#' 4. A.L. Pomerantsev, O.Ye. Rodionova, Concept and role of extreme objects in PCA/SIMCA,
#' Journal of Chemometrics, 28 (2014) pp. 429-438.
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @seealso
#' Methods for \code{pca} objects:
#' \tabular{ll}{
#'    \code{plot.pca} \tab makes an overview of PCA model with four plots.\cr
#'    \code{summary.pca} \tab shows some statistics for the model.\cr
#'    \code{\link{selectCompNum.pca}} \tab set number of optimal components in the model\cr
#'    \code{\link{setDistanceLimits.pca}} \tab set critical limits for residuals\cr
#'    \code{\link{predict.pca}} \tab applies PCA model to a new data.\cr
#'    \code{\link{plotScores.pca}} \tab shows scores plot.\cr
#'    \code{\link{plotLoadings.pca}} \tab shows loadings plot.\cr
#'    \code{\link{plotVariance.pca}} \tab shows explained variance plot.\cr
#'    \code{\link{plotCumVariance.pca}} \tab shows cumulative explained variance plot.\cr
#'    \code{\link{plotResiduals.pca}} \tab shows Q vs. T2 residuals plot.\cr
#' }
#'
#' Most of the methods for plotting data are also available for PCA results (\code{\link{pcares}})
#' objects. Also check \code{\link{pca.mvreplace}}, which replaces missing values in a data matrix
#' with approximated using iterative PCA decomposition.
#'
#' @examples
#' library(mdatools)
#' ### Examples for PCA class
#'
#' ## 1. Make PCA model for People data with autoscaling
#'
#' data(people)
#' model = pca(people, scale = TRUE, info = "Simple PCA model")
#' model = selectCompNum(model, 4)
#' summary(model)
#' plot(model, show.labels = TRUE)
#'
#' ## 2. Show scores and loadings plots for the model
#'
#' par(mfrow = c(2, 2))
#' plotScores(model, comp = c(1, 3), show.labels = TRUE)
#' plotScores(model, comp = 2, type = "h", show.labels = TRUE)
#' plotLoadings(model, comp = c(1, 3), show.labels = TRUE)
#' plotLoadings(model, comp = c(1, 2), type = "h", show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' ## 3. Show residual distance and variance plots for the model
#' par(mfrow = c(2, 2))
#' plotVariance(model, type = "h")
#' plotCumVariance(model, show.labels = TRUE, legend.position = "bottomright")
#' plotResiduals(model, show.labels = TRUE)
#' plotResiduals(model, ncomp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' @export
pca <- function(x, ncomp = min(nrow(x) - 1, ncol(x), 20), center = TRUE, scale = FALSE,
   exclrows = NULL, exclcols = NULL, x.test = NULL, method = "svd", rand = NULL,
   lim.type = "ddmoments", alpha = 0.05, gamma = 0.01, info = "") {

   # exclude columns if "exclcols" is provided
   if (length(exclcols) > 0) {
      x <- mda.exclcols(x, exclcols)
   }

   # exclude rows if "exclrows" is provided
   if (length(exclrows) > 0) {
      x <- mda.exclrows(x, exclrows)
   }

   # calibrate model and set distance limits
   model <- pca.cal(x, ncomp, center = center, scale = scale, method = method, rand = rand)
   model$info <- info
   model$call <- match.call()

   # set distance limits
   model <- setDistanceLimits(model, lim.type = lim.type, alpha = alpha, gamma = gamma)

   # apply the model to calibration set
   model$calres <- predict.pca(model, x)

   # apply model to test set if provided
   if (!is.null(x.test)) {
      model$testres <- predict.pca(model, x.test)
   }

   class(model) <- c("pca")
   return(model)
}

#' Select optimal number of components for PCA model
#'
#' @description
#' Allows user to select optimal number of components for a PCA model
#'
#' @param model
#' PCA model (object of class \code{pca})
#' @param ncomp
#' number of components to select
#' @param ...
#' other parameters if any
#'
#' @return
#' the same model with selected number of components
#'
#' @export
selectCompNum.pca <- function(obj, ncomp, ...) {
   if (ncomp < 1 || ncomp > obj$ncomp) {
      stop("Wrong number of selected components.")
   }

   obj$ncomp.selected <- ncomp
   obj$calres$ncomp.selected <- ncomp

   if (!is.null(obj$testres)) {
      obj$testres$ncomp.selected <- ncomp
   }

   return(obj)
}


#' Compute and set statistical limits for Q and T2 residual distances.
#'
#' @description
#' Computes statisticsl limits for Q and T2 residual distances (based on calibration set)
#' and assign the calculated values as model properties. The method also categorizes objects
#' from calibration and test set accordingly (see \code{categorize.pca} for details).
#'
#' @param obj
#' object with PCA model
#' @param lim.type
#' type of limits ("chisq", "ddmoments", "ddrobust")
#' @param alpha
#' significance level for detection of extreme objects
#' @param gamma
#' significance level for detection of outliers (for data driven approach)
#' @param ...
#' other arguments
#'
#' @details
#'
#' The limits can be accessed as fields of model objects: \code{$Qlim} and \code{$T2lim}. Each
#' is a matrix with four rows and \code{ncomp} columns. First row contains critical limits for
#' extremes, second row - for outliers, third row contains mean value for corresponding distance
#' (or its robust estimate in case of \code{lim.type = "ddrobust"}) and last row contains the
#' degrees of freedom.
#'
#' @return
#' Object models with the two fields updated.
#'
#' @export
setDistanceLimits.pca <- function(obj, lim.type = obj$lim.type, alpha = obj$alpha,
   gamma = obj$gamma, ...) {

   obj$Qlim <- pca.getQLimits(obj, lim.type, alpha = alpha, gamma = gamma)
   obj$T2lim <- pca.getT2Limits(obj, lim.type, alpha = alpha, gamma = gamma)

   obj$alpha <- alpha
   obj$gamma <- gamma
   obj$lim.type <- lim.type

   if (!is.null(obj$calres)) {
      attr(obj$calres$Q, "u0") <- obj$Qlim[, 3]
      attr(obj$calres$T2, "u0") <- obj$T2lim[, 3]
   }

   if (!is.null(obj$testres)) {
      attr(obj$testres$Q, "u0") <- obj$Qlim[, 3]
      attr(obj$testres$T2, "u0") <- obj$T2lim[, 3]
   }

   return(obj)
}

#' Categorize PCA results based on orthogonal and score distances.
#'
#' @description
#' The method compares score and orthogonal distances of PCA results from \code{res} with
#' critical limits computed for the PCA model and categorizes the corresponding objects as
#' "normal", "extreme" or "outlier".
#'
#' @param obj
#' object with PCA model
#' @param res
#' object with PCA results
#' @param ncomp
#' number of components to use for the categorization
#' @param ...
#' other parameters
#'
#' @details
#' The method does not categorize hidden values if any.
#'
#' @return
#' vector (factor) with results of categorization.
#'
#' @export
categorize.pca <- function(obj, res, ncomp = obj$ncomp.selected, ...) {

   create_categories <- function(nobj, extremes_ind, outliers_ind) {
      categories <- rep(1, nobj)
      categories[extremes_ind] <- 2
      categories[outliers_ind] <- 3
      return(factor(categories, levels = 1:3, labels = c("normal", "extreme", "outlier")))
   }

   # get distance values for selected number of components
   h <- res$T2[, ncomp]
   q <- res$Q[, ncomp]

   # remove excluded values if any
   rows_excluded <- attr(res$T2, "exclrows")
   if (length(rows_excluded) > 0) {
      h <- h[-rows_excluded]
      q <- q[-rows_excluded]
   }

   nobj <- length(h)

   # if chisq / hotelling
   if (obj$lim.type == "chisq") {
      outliers_ind <- (h >= obj$T2lim[2, ncomp] | q >= obj$Qlim[2, ncomp])
      extremes_ind <- (h >= obj$T2lim[1, ncomp] | q >= obj$Qlim[1, ncomp])
      return(create_categories(nobj, extremes_ind, outliers_ind))
   }

   # if data driven
   h0 <- obj$T2lim[3, ncomp]
   Nh <- round(obj$T2lim[4, ncomp])
   q0 <- obj$Qlim[3, ncomp]
   Nq <- round(obj$Qlim[4, ncomp])

   f <- Nh * h / h0 + Nq * q / q0
   outliers_ind <- f > (obj$T2lim[2, ncomp] * Nh / h0)
   extremes_ind <- f > (obj$T2lim[1, ncomp] * Nh / h0)
   return(create_categories(nobj, extremes_ind, outliers_ind))
}

#' PCA predictions
#'
#' @description
#' Applies PCA model to a new data.
#'
#' @param object
#' a PCA model (object of class \code{pca}).
#' @param x
#' data values (matrix or data frame).
#' @param ...
#' other arguments.
#'
#' @return
#' PCA results (an object of class \code{pcares})
#'
#' @export
predict.pca <- function(object, x, ...) {
   # convert to matrix
   x <- mda.df2mat(x)
   attrs <- attributes(x)

   if (ncol(x) != nrow(object$loadings)) {
      stop("Number and type of data columns should be the same as in calibration dataset.")
   }

   # compute scores and residuals
   x <- prep.autoscale(x, center = object$center, scale = object$scale)
   scores <- x %*% object$loadings
   residuals <- x - tcrossprod(scores, object$loadings)

   # set names
   rownames(scores) <- rownames(residuals) <- attrs$dimnames[[1]]
   colnames(scores) <- colnames(loadings)
   colnames(residuals) <- attrs$dimnames[[2]]

   # set attributes
   scores <- mda.setattr(scores, attrs, "row")
   residuals <- mda.setattr(residuals, attrs)
   attr(scores, "name") <- "Scores"
   attr(scores, "xaxis.name") <- "Components"
   attr(residuals, "name") <- "Residuals"

   # create and return the results object
   res <- pcares(scores, object$loadings, residuals, object$eigenvals, object$ncomp.selected)
   attr(res$Q, "u0") <- object$Qlim[3, ]
   attr(res$T2, "u0") <- object$T2lim[3, ]

   return(res)
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
print.pca <- function(x, ...) {

   cat("\nPCA model (class pca)\n")

   if (length(x$info) > 1) {
      cat("\nInfo:\n")
      cat(x$info)
   }

   cat("\n\nCall:\n")
   print(x$call)

   cat('\nMajor model fields:\n')
   cat("$loadings - matrix with loadings\n")
   cat("$eigenvals - eigenvalues for components\n")
   cat("$ncomp - number of calculated components\n")
   cat("$ncomp.selected - number of selected components\n")
   cat("$center - values for centering data\n")
   cat("$scale - values for scaling data\n")
   cat("$alpha - significance level for critical limits\n")
   cat("$gamma - significance level for outlier limits\n")
   cat("$Qlim - critical values and parameters for Q distance\n")
   cat("$T2lim - critical values and parameters for T2 distance\n")

   cat("\nFields with results:\n")
   cat("$calres - results (scores, etc) for calibration set\n")

   if (!is.null(x$testres)) {
      cat("$testres - results for test set\n")
   }
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
summary.pca <- function(object, ...) {

   cat("\nSummary for PCA model (class pca)\n")

   if (length(object$info) > 1) {
      fprintf("\nInfo:\n%s\n", object$info)
   }

   if (!is.null(object$rand)) {
      fprintf("\nParameters for randomized algorithm: q = %d, p = %d\n",
         object$rand[1], object$rand[2])
   }

   if (length(object$exclrows) > 0) {
      fprintf("Excluded rows: %d\n", length(object$exclrows))
   }

   if (length(object$exclcols) > 0) {
      fprintf("Excluded coumns: %d\n", length(object$exclcols))
   }

   fprintf("Type of limits: %s\n", object$lim.type)
   fprintf("Alpha: %s\n", object$alpha)
   fprintf("Gamma: %s\n", object$gamma)
   cat("\n")

   data <- cbind(
      round(object$eigenvals, 3),
      round(object$calres$expvar, 2),
      round(object$calres$cumexpvar, 2),
      object$Qlim[4, ],
      object$T2lim[4, ]
   )

   colnames(data) <- c("Eigenvals", "Expvar", "Cumexpvar", "Nq", "Nh")
   rownames(data) <- colnames(object$loadings)
   show(data)
}


################################
#  Static methods              #
################################


#' Replace missing values in data
#'
#' @description
#' \code{pca.mvreplace} is used to replace missing values in a data matrix with
#' approximated by iterative PCA decomposition.
#'
#' @param x
#' a matrix with data, containing missing values.
#' @param center
#' logical, do centering of data values or not.
#' @param scale
#' logical, do standardization of data values or not.
#' @param maxncomp
#' maximum number of components in PCA model.
#' @param expvarlim
#' minimum amount of variance, explained by chosen components (used for selection of optimal number
#' of components in PCA models).
#' @param covlim
#' convergence criterion.
#' @param maxiter
#' maximum number of iterations if convergence criterion is not met.
#'
#' @details
#' The function uses iterative PCA modeling of the data to approximate and impute missing values.
#' The result is most optimal for data sets with low or moderate level of noise and with number of
#' missing values less than 10\% for small dataset and up to 20\% for large data.
#'
#' @return
#' Returns the same matrix \code{x} where missing values are replaced with approximated.
#'
#' @references
#' Philip R.C. Nelson, Paul A. Taylor, John F. MacGregor. Missing data methods in PCA and PLS:
#' Score calculations with incomplete observations. Chemometrics and Intelligent Laboratory
#' Systems, 35 (1), 1996.
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @examples
#' library(mdatools)
#'
#' ## A very simple example of imputing missing values in a data with no noise
#'
#' # generate a matrix with values
#' s = 1:6
#' odata = cbind(s, 2*s, 4*s)
#'
#' # make a matrix with missing values
#' mdata = odata
#' mdata[5, 2] = mdata[2, 3] = NA
#'
#' # replace missing values with approximated
#' rdata = pca.mvreplace(mdata, scale = TRUE)
#'
#' # show all matrices together
#' show(cbind(odata, mdata, round(rdata, 2)))
#'
#' @export
pca.mvreplace <- function(x, center = T, scale = F, maxncomp = 7,
                         expvarlim = 0.95, covlim = 10^-6, maxiter = 100) {
   mvidx <- is.na(x)

   # check if any column has more than 20% values
   nmv_cols <- colSums(mvidx) / nrow(x)
   if (any(nmv_cols > 0.2)) {
      stop("Some of columns have more than 20% missing values.")
   }

   # make initial estimates with mean values for each column
   col_means <- matrix(apply(x, 2, mean, na.rm = TRUE), nrow(x), ncol(x))
   x[mvidx] <- col_means[midx]

   # autoscale
   x <- scale(x, center = center, scale = scale)

   if (scale) {
      gsd <- attr(x.rep, 'scaled:scale')
   }

   if (center) {
      gmean <- attr(x.rep, 'scaled:center');
   }

   n <- 1
   scoresp <- 0
   scores <- 1
   cond <- 1
   maxncomp <- min(maxncomp, nrow(x) - 1, ncol(x))

   while (cond > covlim && n < maxiter) {

      # recenter data on every iteration
      x <- scale(x, center = TRUE, scale = FALSE)
      lmean <- attr(x, 'scaled:center')

      # compute PCA decomposition annd find number of components
      res <- pca.svd(x.rep, maxncomp)
      expvar <- cumsum(res$eigenvals/sum(res$eigenvals))
      ncomp <- min(which(expvar >= expvarlim), maxncomp)

      # correct number of components for border cases
      if (ncomp == length(expvar)) ncomp <- ncomp - 1
      if (ncomp == 0) ncomp <- 1

      # get and trancate scores and loadings and reestimate the values
      scoresp <- scores
      loadings <- res$loadings[, 1:ncomp]
      scores <- x %*% loadings
      x_new <- tcrossprod(scores, loadings)

      # remove centering
      x_new <- sweep(x.new, 2L, lmean, '+', check.margin = F)

      # replace missing values by the calculated
      x[mvidx] <- x_new[mvidx]

      if (n > 2) {
         # calculate difference between scores for convergence
         ncompcond <- min(ncol(scores), ncol(scoresp))
         cond <- sum((scores[, 1:ncompcond] - scoresp[, 1:ncompcond])^2)
      }
   }

   # rescale the data back and return
   if (scale) {
      x <- sweep(x.rep, 2L, gsd, '*', check.margin = FALSE)
   }

   if (center) {
      x <- sweep(x.rep, 2L, gmean, '+', check.margin = FALSE)
   }

   return(x)
}

#' Low-dimensional approximation of data matrix X
#'
#' @param X
#' data matrix
#' @param k
#' rank of X (number of components)
#' @param rand
#' a vector with two values - number of iterations (q) and oversmapling parameter (p)
#' @param dist
#' distribution for generating random numbers, 'unif' or 'norm'
#'
#' @import stats
getB <- function(X, k = NULL, rand = NULL, dist = "unif") {

   if(is.null(rand)) {
      return(X)
   }

   nrows <- nrow(X)
   ncols <- ncol(X)

   q <- rand[1]
   p <- rand[2]
   k <- if (is.null(k)) ncols else 2 * k

   l <- k + p
   Y <- if (dist == 'unif')
            X %*% matrix(runif(ncols * l, -1, 1), ncols, l)
         else
            X %*% matrix(rnorm(ncols * l), ncols, l)

   Q <- qr.Q(qr(Y))
   if (q > 0) {
      for (i in seq_len(q)) {
         Y <- crossprod(X, Q)
         Q <- qr.Q(qr(Y))
         Y <- X %*% Q
         Q <- qr.Q(qr(Y))
      }
   }

   return(crossprod(Q, X))
}


#' Singular Values Decomposition based PCA algorithm
#'
#' @description
#' Computes principal component space using Singular Values Decomposition
#'
#' @param x
#' a matrix with data values (preprocessed)
#' @param ncomp
#' number of components to calculate
#'
#' @return
#' a list with scores, loadings and eigenvalues for the components
#'
pca.svd <- function(x, ncomp = min(ncol(x), nrow(x) - 1)) {

   s <- svd(x, nu = ncomp, nv = ncomp)
   return(
      list(
         loadings = s$v,
         scores = s$u %*% diag(s$d)[1:ncomp, 1:ncomp],
         eigenvals = s$d[1:ncomp]^2 / (nrow(x) - 1),
         ncomp = ncomp
      )
   )
}

#' NIPALS based PCA algorithm
#'
#' @description
#' Calculates principal component space using non-linear iterative partial least squares algorithm
#' (NIPALS)
#'
#' @param x
#' a matrix with data values (preprocessed)
#' @param ncomp
#' number of components to calculate
#' @param tol
#' tolerance (if difference in eigenvalues is smaller - convergence achieved)
#'
#' @return
#' a list with scores, loadings and eigenvalues for the components
#'
#' @references
#' Geladi, Paul; Kowalski, Bruce (1986), "Partial Least Squares
#' Regression:A Tutorial", Analytica Chimica Acta 185: 1-17
#'
pca.nipals <- function(x, ncomp = min(ncol(x), nrow(x) - 1), tol = 10^-10) {
   nobj <- nrow(x)
   nvar <- ncol(x)

   scores <- matrix(0, nrow = nobj, ncol = ncomp)
   loadings <- matrix(0, nrow = nvar, ncol = ncomp)

   E <- x
   for (i in seq_len(ncomp)) {
      ind <- which.max(apply(E, 2, sd))
      t <- E[, ind, drop = F]
      tau <- th <- 99999999
      while (th > tol * tau) {
         p <- crossprod(E, t) / as.vector(crossprod(t))
         p <- p / as.vector(crossprod(p)) ^ 0.5
         t <- (E %*% p)/as.vector(crossprod(p))

         th <- abs(tau - as.vector(crossprod(t)))
         tau <- as.vector(crossprod(t))
      }

      E <- E - tcrossprod(t, p)
      scores[, i] <- t
      loadings[, i] <- p
   }

   return(
      list(
         loadings = loadings,
         scores = scores,
         eigenvals = colSums(scores^2) / (nobj - 1)
      )
   )
}

#' Runs one of the selected PCA methods
#'
#' @param x
#' data matrix
#' @param ncomp
#' number of components
#' @param method
#' name of PCA methods ('svd', 'nipals')
#' @param rand
#' parameters for randomized algorithm (if not NULL)
#'
#' @export
pca.run <- function(x, ncomp, method, rand = NULL) {

   # define which function to use depending on method name
   f <- switch(
      method,
      "svd" = pca.svd,
      "nipals" = pca.nipals,
      stop('Wrong value for PCA method!')
   )

   # use proper function to compute scores, loadings and eigenvalues
   res <- f(getB(x, k = ncomp, rand = rand), ncomp)

   # recompute scores and eigenvalues if randomized algorithm was used
   if (!is.null(rand)) {
      res$scores <- x %*% res$loadings
      res$eigenvals <- colSums(res$scores^2)/(nrow(x) - 1)
   }

   # compute and add residuals
   res$residuals <- x - tcrossprod(res$scores, res$loadings)

   return(res)
}

#' PCA model calibration
#'
#' @description
#' Calibrates (builds) a PCA model for given data and parameters
#'
#' @param x
#' matrix with data values
#' @param ncomp
#' number of principal components to calculate
#' @param center
#' logical, do mean centering or not
#' @param scale
#' logical, do standardization or not
#' @param method
#' algorithm for compiting PC space (only 'svd' and 'nipals' are supported so far)
#' @param rand
#' vector with parameters for randomized PCA methods (if NULL, conventional PCA is used instead)
#'
#' @return
#' an object with calibrated PCA model
#'
pca.cal <- function(x, ncomp, center, scale, method, rand = NULL) {

   # check if data has missing values
   if (any(is.na(x))) {
      stop("Data has missing values, try to fix this using pca.mvreplace.")
   }

   # prepare empty list for model object
   model <- list()

   # save data attributes
   attrs <- mda.getattr(x)

   # convert data to a matrix
   x <- mda.df2mat(x)
   nrows_full <- nrow(x)
   ncols_full <- ncol(x)

   # correct maximum number of components
   ncols <- ncols_full - length(attrs$exclcols)
   nrows <- nrows_full - length(attrs$exclrows)

   # make sure that ncomp is correct
   ncomp <- min(ncomp, ncols, nrows - 1)

   # prepare data for model calibration and cross-validation
   x_cal <- x

   # remove excluded rows
   if (length(attrs$exclrows) > 0) {
      x_cal = x_cal[-attrs$exclrows, , drop = F]
   }

   # autoscale and save the mean and std values for predictions
   x_cal <- prep.autoscale(x_cal, center = center, scale = scale)
   model$center <- attr(x_cal, "prep:center")
   model$scale <- attr(x_cal, "prep:scale")

   # remove excluded columns
   if (length(attrs$exclcols) > 0) {
      x_cal <- x_cal[, -attrs$exclcols, drop = F]
   }

   # compute loadings, scores and eigenvalues for data without excluded elements
   res <- pca.run(x_cal, ncomp, method, rand)

   # correct loadings for missing columns in x
   # corresponding rows in loadings will be set to 0 and excluded
   if (length(attrs$exclcols) > 0) {
      loadings <- matrix(0, nrow = ncols_full, ncol = ncomp)
      loadings[-attrs$exclcols, ] <- res$loadings
      loadings <- mda.exclrows(loadings, attrs$exclcols)
   } else {
      loadings <- res$loadings
   }

   if (is.null(dim(loadings))) {
      loadings <- matrix(loadings, ncol = ncomp)
   }

   # set names and attributes for the loadings
   rownames(loadings) <- colnames(x)
   colnames(loadings) <- paste("Comp", 1:ncol(loadings))
   attr(loadings, "name") <- "Loadings"
   attr(loadings, "xaxis.name") <- "Components"
   attr(loadings, "yaxis.name") <- attrs$xaxis.name
   attr(loadings, "yaxis.values") <- attrs$xaxis.values
   model$loadings <- loadings

   # organize eigen values
   model$eigenvals <- res$eigenvals
   names(model$eigenvals) <- colnames(loadings)

   # compute parameters for distance critical limits
   model$limParams <- pca.getLimParams(res)

   # finalize model
   model$method <- method
   model$rand <- rand

   # setup other fields and return the model
   model$ncomp <- ncomp
   model$ncomp.selected <- model$ncomp

   # save excluded columns and rows
   model$exclcols <- attrs$exclcols
   model$exclrows <- attrs$exclrows
   class(model) <- "pca"

   return(model)
}

#' Compute parameters for critical limits based on calibration results
#'
#' @param res
#' results returned by \code{pca.run()}
#'
pca.getLimParams <- function(res) {

   dist <- ldecomp.getDistances(res$scores, res$loadings, res$residuals, res$eigenvals)
   return(
      list(
         "T2" = list(
            "moments" = ddmoments.param(dist$T2),
            "robust" = ddrobust.param(dist$T2),
            "nobj" = nrow(res$scores)
         ),
         "Q" = list(
            "moments" = ddmoments.param(dist$Q),
            "robust" = ddrobust.param(dist$Q),
            "nobj" = nrow(res$scores)
         )
      )
   )
}

#' Compute critical limits for orthogonal distances (Q)
#'
#' @param model
#' PCA model
#' @param lim.type
#' which method to use for calculation of critical limits for residuals (see details)
#' @param alpha
#' significance level for extreme limits.
#' @param gamma
#' significance level for outlier limits.
#'
pca.getQLimits <- function(model, lim.type, alpha, gamma) {

   params <- model$limParams
   pQ <- if (regexpr("robust", lim.type) > 0) params$Q$robust else params$Q$moments
   pT2 <- if (regexpr("robust", lim.type) > 0) params$T2$robust else params$T2$moments

   DoF <- round(pQ$Nu)
   DoF[DoF < 1] <- 1
   DoF[DoF > 250] <- 250

   lim <- switch(lim.type,
      "chisq" = chisq.crit(pQ, alpha, gamma),
      "ddmoments" = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pQ$u0),
      "ddrobust"  = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pQ$u0),
      stop("Wrong value for 'lim.type' parameter.")
   )

   lim <- rbind(lim, pQ$u0, DoF)
   colnames(lim) <- colnames(model$loadings)
   rownames(lim) <- c("Extremes limits", "Outliers limits", "Mean", "DoF")
   attr(lim, "name") <- "Critical limits for orthogonal distances (Q)"
   attr(lim, "alpha") <- alpha
   attr(lim, "gamma") <- gamma
   attr(lim, "lim.type") <- lim.type

   return(lim)
}

#' Compute critical limits for score distances (T2)
#'
#' @param model
#' PCA model
#' @param lim.type
#' which method to use for calculation ("chisq", "ddmoments", "ddrobust")
#' @param alpha
#' significance level for extreme limits.
#' @param gamma
#' significance level for outlier limits.
#'
pca.getT2Limits <- function(model, lim.type, alpha, gamma) {

   params <- model$limParams
   pQ <- if (regexpr("robust", lim.type) > 0) params$Q$robust else params$Q$moments
   pT2 <- if (regexpr("robust", lim.type) > 0) params$T2$robust else params$T2$moments

   DoF <- round(pT2$Nu)
   DoF[DoF < 1] <- 1
   DoF[DoF > 250] <- 250

   lim <- switch(lim.type,
      "chisq" = hotelling.crit(pT2$nobj, 1:model$ncomp, alpha, gamma),
      "ddmoments" = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pT2$u0),
      "ddrobust"  = scale(dd.crit(pQ, pT2, alpha, gamma), center = FALSE, scale = DoF / pT2$u0),
      stop("Wrong value for 'lim.type' parameter.")
   )

   lim <- rbind(lim, pT2$u0, DoF)
   colnames(lim) <- colnames(model$loadings)
   rownames(lim) <- c("Extremes limits", "Outliers limits", "Mean", "DoF")
   attr(lim, "name") <- "Critical limits for score distances (T2)"
   attr(lim, "alpha") <- alpha
   attr(lim, "gamma") <- gamma
   attr(lim, "lim.type") <- lim.type

   return(lim)
}


################################
#  Plotting methods            #
################################


#' Explained variance plot for PCA
#'
#' @description
#' Shows a plot with explained variance or cumulative explained variance for components.
#'
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param type
#' type of the plot ("b", "l", "h")
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pca}} function.
#'
#' @export
plotVariance.pca <- function(obj, type = "b", main = "Variance", xlab = "Components",
   ylab = "Explained variance, %", show.legend = TRUE, res = getModelRes(obj),
   variance = "expvar", ...) {

   if (length(res) == 1) {
      return(
         plotVariance(res[[1]], type = type, main = main, xlab = xlab, ylab = ylab,
            variance = variance, ...)
      )
   }

   nres <- length(res)
   plot_data <- matrix(0, nres, obj$ncomp)
   for (i in seq_len(nres)) {
      plot_data[i, ] <- plotVariance(res[[i]], type = type, variance = variance, show.plot = FALSE)
   }

   colnames(plot_data) <- colnames(obj$loadings)
   rownames(plot_data) <- names(res)
   mdaplotg(plot_data, main = main, xlab = xlab, xticks = 1:obj$ncomp, ylab = ylab,
            show.legend = show.legend, type = type, ...)
}

#' Cumulative explained variance plot for PCA
#'
#' @description
#' Shows a plot with cumulative explained variance for components.
#'
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pca}} function.
#'
#' @export
plotCumVariance.pca <- function(obj, type = "b", main = "Cumulative variance", xlab = "Components",
   ylab = "Explained variance, %", show.legend = TRUE, res = getModelRes(obj), ...) {

   plotVariance(obj, type = type, main = main, xlab = xlab, show.legend = show.legend,
      res = res, variance = "cumexpvar")
}


#' Scores plot for PCA
#'
#' @description
#' Shows a scores plot for selected components.
#'
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param comp
#' a value or vector with several values - number of components to show the plot for
#' @param type
#' type of the plot ('b', 'l', 'h')
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param cgroup
#' a vector with numeric values or a factor used for color grouping of plot points.
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pca}} function.
#'
#' @export
plotScores.pca <- function(obj, comp = c(1, 2), type = "p", main = "Scores", show.labels = F,
   show.legend = TRUE, legend.position = "topright", show.axes = TRUE,
   res = getModelRes(obj), ...) {

   if (length(res) == 1) {
      return(
         plotScores(res[[1]], comp = comp, type = type, main = main,
            show.labels = show.labels, ...)
      )
   }

   plot_data <- list()
   for (i in seq_len(length(res))) {
      plot_data[[names(res[i])]] <-
         plotScores(res[[i]], comp = comp, type = type, show.plot = FALSE)
   }

   if (type != "p") {
      stop("You have several result objects in model, only scatter plot is available for scores.")
   }


   # set up values for showing axes lines
   show.lines <- FALSE
   if (show.axes) {
      show.lines <- if (length(comp) == 2 && type == "p") c(0, 0) else c(NA, 0)
   }

   mdaplotg(data = plot_data, main = main, show.labels = show.labels, show.legend = show.legend,
      legend.position = legend.position, show.lines = show.lines, ...)
}


#' Residuals plot for PCA
#'
#' @description
#' Shows a plot with Q residuals vs. Hotelling T2 values for selected number of components.
#'
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param ncomp
#' how many components to use (if NULL - user selected optimal value will be used)
#' @param norm
#' logical, show normalized Q vs T2 (\code{norm = T}) values or original ones (\code{norm = F})
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param show.limits
#' logical, show or not lines with statistical limits for the residuals
#' @param xlim
#' limits for x-axis
#' @param ylim
#' limits for y-axis
#' @param lim.col
#' vector with two values - line color for extreme and outlier borders
#' @param lim.lwd
#' vector with two values - line width for extreme and outlier borders
#' @param lim.lty
#' vector with two values - line type for extreme and outlier borders
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pca}} function.
#'
#' @export
plotResiduals.pca = function(obj, ncomp = NULL, norm = F, main = NULL, xlab = NULL, ylab = NULL,
                             show.labels = F,  show.legend = T, show.limits = T,
                             xlim = NULL, ylim = NULL,
                             lim.col = c('#333333', '#333333'),
                             lim.lwd = c(1, 1), lim.lty = c(2, 3),
                             ...) {
   if (is.null(main)) {
      if (is.null(ncomp))
         main = 'Residuals'
      else
         main = sprintf('Residuals (ncomp = %d)', ncomp)
   }

   if (is.null(ncomp))
      ncomp = obj$ncomp.selected

   if (ncomp > obj$ncomp || ncomp < 1)
      stop('Wrong number of components!')

   data = list()
   data$cal = mda.cbind(
      mda.subset(obj$calres$T2, select = ncomp),
      mda.subset(obj$calres$Q, select = ncomp)
   )

   # set values for normalization of residuals if necessary
   if (norm) {
      T2.mean = obj$T2lim[3, ncomp]
      Q.mean = obj$Qlim[3, ncomp]
      if (is.null(xlab))
         xlab = expression(paste('Hotelling ', T^2, ' distance (norm)'))
      if (is.null(ylab))
         ylab = 'Squared residual distance, Q (norm)'
   } else {
      T2.mean = 1
      Q.mean = 1
      if (is.null(xlab))
         xlab = expression(paste('Hotelling ', T^2, ' distance'))
      if (is.null(ylab))
         ylab = 'Squared residual distance, Q'
   }

   data$cal[, 1] = data$cal[, 1] / T2.mean
   data$cal[, 2] = data$cal[, 2] / Q.mean
   x.max = max(data$cal[, 1])
   y.max = max(data$cal[, 2])

   if (!is.null(obj$cvres)) {
      data$cv = mda.cbind(
         obj$cvres$T2[, ncomp]/T2.mean,
         obj$cvres$Q[ , ncomp]/Q.mean
      )
      x.max = max(x.max, data$cv[, 1])
      y.max = max(y.max, data$cv[, 2])
   }

   if (!is.null(obj$testres)) {
      data$test = mda.cbind(
         obj$testres$T2[, ncomp]/T2.mean,
         obj$testres$Q[ , ncomp]/Q.mean
      )
      x.max = max(x.max, data$test[, 1])
      y.max = max(y.max, data$test[, 2])
   }

   if (show.limits == T) {
      # get residual limits, correct if necessary and recalculate axes maximum limit
      lim = cbind(obj$T2lim[1:2, ncomp], obj$Qlim[1:2, ncomp])
      if (substr(obj$lim.type, 1, 2) != 'dd') {
         lim[, 1] = lim[, 1] / T2.mean
         lim[, 2] = lim[, 2] / Q.mean
         x.max = max(x.max, lim[, 1])
         y.max = max(y.max, lim[, 2])
      } else {
         lim[, 1] = lim[, 1] * T2.mean / Q.mean
         lim[, 2] = lim[, 2] / Q.mean
         x.max = 1.5 * x.max
         y.max = 1.5 * y.max
      }
   }

   # use computed max values for axes limits if user did not specify anything
   if (is.null(xlim))
      xlim = c(0, 1.2 * x.max)
   if (is.null(ylim))
      ylim = c(0, 1.2 * y.max)

   # show plot
   if (length(data) == 1) {
      mdaplot(data[[1]], main = main, xlab = xlab, ylab = ylab,
              show.labels = show.labels, xlim = xlim, ylim = ylim, ...)
   } else {
      mdaplotg(data, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
               show.labels = show.labels, show.legend = show.legend, ...)
   }

   # show limits
   if (show.limits) {
      ldecomp.plotLimits(lim, obj$lim.type, lim.col, lim.lwd, lim.lty)
   }
}

#' Loadings plot for PCA
#'
#' @description
#' Shows a loadings plot for selected components.
#'
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param comp
#' a value or vector with several values - number of components to show the plot for
#' @param type
#' type of the plot ('b', 'l', 'h')
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#'
#' @details
#' See examples in help for \code{\link{pca}} function.
#'
#' @export
plotLoadings.pca = function(obj, comp = c(1, 2), type = NULL, main = 'Loadings',
                            xlab = NULL, ylab = NULL, show.labels = NULL, show.legend = TRUE,
                            show.axes = TRUE, ...) {

   if (max(comp) > obj$ncomp || min(comp) < 1)
      stop('Wrong number of components!')

   if (is.null(type)) {
      if (length(comp) == 2)
         type = 'p'
      else
         type = 'l'
   }

   data = mda.subset(obj$loadings, select = comp)

   if (type == 'p') {
      if (show.axes == TRUE) {
         if (length(comp) > 1)
            show.lines = c(0, 0)
         else
            show.lines = c(NA, 0)
      } else {
         show.lines = FALSE
      }

      if (is.null(show.labels))
         show.labels = TRUE

      colnames(data) =
         paste('Comp ', comp, ' (', round(obj$calres$expvar[comp], 2) , '%)', sep = '')
      mdaplot(data, type = type, show.labels = show.labels, show.lines = show.lines,
              main = main, ylab = ylab, xlab = xlab, ...)
   } else {
      if (is.null(show.legend))
         show.legend = TRUE

      if (is.null(show.labels))
         show.labels = FALSE

      if (is.null(ylab))
         ylab = ''

      if (show.axes == TRUE && type != 'h')
         show.lines = c(NA, 0)
      else
         show.lines = FALSE

      mdaplotg(mda.t(data), show.legend = show.legend, type = type, show.labels = show.labels,
               show.lines = show.lines, main = main, ylab = ylab, xlab = xlab, ...)
   }
}


#' PCA biplot
#'
#' @description
#' Shows a biplot for selected components.
#'
#' @param obj
#' a PCA model (object of class \code{pca})
#' @param comp
#' a value or vector with several values - number of components to show the plot for
#' @param pch
#' a vector with two values - markers for scores and loadings
#' @param col
#' a vector with two colors for scores and loadings
#' @param main
#' main title for the plot
#' @param lty
#' line type for loadings
#' @param lwd
#' line width for loadings
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.axes
#' logical, show or not a axes lines crossing origin (0,0)
#' @param show.excluded
#' logical, show or hide rows marked as excluded (attribute `exclrows`)
#' @param lab.col
#' a vector with two colors for scores and loadings labels
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' @export
plotBiplot.pca = function(obj, comp = c(1, 2), pch = c(16, NA), col = mdaplot.getColors(2),
                          main = 'Biplot', lty = 1, lwd = 1, show.labels = FALSE, show.axes = TRUE,
                          show.excluded = FALSE, lab.col = c('#90A0D0', '#D09090'), ...) {

   if (length(comp) != 2)
      stop('Biplot can be made only for two principal components!')

   if (show.axes == TRUE)
      show.lines = c(0, 0)
   else
      show.lines = FALSE

   loadings = mda.subset(obj$loadings, select = comp)
   scores = mda.subset(obj$calres$scores, select = comp)
   attrs = mda.getattr(scores)

   loadsScaleFactor = sqrt(max(rowSums(loadings^2)))

   scoresScaleFactor = max(abs(scores))
   scores = (scores / scoresScaleFactor) * loadsScaleFactor
   scores = mda.setattr(scores, attrs)


   colnames(scores) =
      paste('Comp ', comp, ' (', round(obj$calres$expvar[comp], 2) , '%)', sep = '')

   mdaplotg(list(scores = scores, loadings = loadings), type = 'p', pch = pch,
            show.legend = FALSE, show.labels = show.labels, lab.col = lab.col,
            main = main, colmap = col, show.lines = show.lines, show.excluded = show.excluded, ...)

   if (show.excluded == TRUE && length(attr(loadings, 'exclrows')) > 0)
      loadings = loadings[-attr(loadings, 'exclrows'), , drop = F]

   segments(0, 0, loadings[, 1], loadings[, 2], col = col[2], lty = lty, lwd = lwd)
}

#' Model overview plot for PCA
#'
#' @description
#' Shows a set of plots (scores, loadings, residuals and explained variance) for PCA model.
#'
#' @param x
#' a PCA model (object of class \code{pca})
#' @param comp
#' vector with two values - number of components to show the scores and loadings plots for
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.legend
#' logical, show or not a legend on the plot
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{pca}} function.
#'
#' @export
plot.pca = function(x, comp = c(1, 2), show.labels = FALSE, show.legend = TRUE, ...) {
   obj = x

   par(mfrow = c(2, 2))
   plotScores(obj, comp = comp, show.labels = show.labels, show.legend = show.legend)
   plotLoadings(obj, comp = comp, show.labels = show.labels, show.legend = show.legend)
   plotResiduals(obj, ncomp = obj$ncomp.selected,  show.labels = show.labels,
                 show.legend = show.legend, show.limits = T)
   plotCumVariance(obj, show.legend = show.legend)
   par(mfrow = c(1, 1))
}
