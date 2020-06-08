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
#' @param spec.ini
#' a matrix with initial estimation of the pure components spectra.
#' @param cont.solver
#' which function to use as a solver for resolving of pure components contributions (see detials).
#' @param spec.solver
#' which function to use as a solver for resolving of pure components spectra (see detials).
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values).
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values).
#' @param verbose
#' logical, if TRUE information about every iteration will be shown.
#' @param max.niter
#' maximum number of iterations.
#' @param tol
#' tolerance, when explained variance change is smaller than this value, iterations stop.
#' @param info
#' a short text with description of the case (optional).
#'
#' @details
#' The method implements the iterative ALS algorithm, where, at each iteration, spectra and
#' contributions of each chemical component are estimated and then a set of constraints is
#' applied to each. The method is well described in [1, 2].
#'
#' The method assumes that the spectra (D) is a linear combination of pure components spectra (S)
#' and pure component concentrations (C):
#'
#' D = CS' + E
#'
#' So the task is to get C and S by knowing D. In order to do that you need to provide:
#'
#' 1. Constraints for spectra and contributions. The constraints should be provided as a list
#' with name of the constraint and all necessary parameters. You can see which constraints and
#' parameters are currently supported by running \code{constraintList()}. See the code examples
#' below or a Bookdown tutorial for more details.
#'
#' 2. Initial estimation of the pure components spectra, S. By default method uses a matrix with
#' random numbers but you can provide a better guess (for example by running \code{\link{mcrpure}})
#' as a first step.
#'
#' 3. Which solver to use for resolving spectra and concentrations. There are two built in solvers:
#' \code{mcrals.nnls} (default) and \code{mcrals.ols}. The first implements non-negative least
#' squares method which gives non-negative (thus physically meaningful) solutions. The second is
#' ordinary least squares and if you want to get non-negative spectra and/or contributions in this
#' case you need to provide a non-negativity constraint.
#'
#' The algorithm iteratively resolves C and S and checks how well CS' is to D. The iterations stop
#' either when number exceeds value in \code{max.niter} or when improvements (difference between
#' explained variance on current and previous steps) is smaller than \code{tol} value.
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
#' 1. J. Jaumot, R. Gargallo, A. de Juan, and R. Tauler, "A graphical user-friendly interface for
#' MCR-ALS: a new tool for multivariate curve resolution in MATLAB", Chemometrics and Intelligent #' Laboratory Systems 76, 101-110 (2005).
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
#'
#' # define constraints for contributions
#' cc <- list(
#'    constraint("non-negativity")
#' )
#'
#' # define constraints for spectra
#' cs <- list(
#'    constraint("non-negativity"),
#'    constraint("norm", params = list(type = "area"))
#' )
#'
#' # because by default initial approximation is made by using random numbers
#' # we need to seed the generator in order to get reproducable results
#' set.seed(6)
#'
#' # run ALS
#' m <- mcrals(carbs$D, ncomp = 3, cont.constraints = cc, spec.constraints = cs)
#' summary(m)
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
mcrals <- function(x, ncomp,
   cont.constraints = list(),
   spec.constraints = list(),
   spec.ini = matrix(runif(ncol(x) * ncomp), ncol(x), ncomp),
   cont.solver = mcrals.nnls,
   spec.solver = mcrals.nnls,
   exclrows = NULL, exclcols = NULL, verbose = FALSE,
   max.niter = 100, tol = 10^-6, info = "") {

   stopifnot("Parameter 'max.niter' should be positive." = max.niter > 0)

   # get pure variables and unmix data
   x <- prepCalData(x, exclrows, exclcols, min.nrows = 2, min.ncols = 2)
   model <- mcrals.cal(
      x, ncomp,
      cont.constraints = cont.constraints,
      spec.constraints = spec.constraints,
      spec.ini = spec.ini,
      cont.solver = cont.solver,
      spec.solver = spec.solver,
      max.niter = max.niter,
      tol = tol,
      verbose = verbose
   )

   # compute explained variance
   model <- c(model, getVariance.mcr(model, x))

   # reorder spectra and contributions according to the explained variance
   ind <- order(model$expvar, decreasing = TRUE)
   model$rescont <- mda.subset(model$rescont, select = ind)
   model$resspec <- mda.subset(model$resspec, select = ind)
   model$expvar <- model$expvar[ind, drop = FALSE]
   colnames(model$rescont) <- colnames(model$resspec) <- names(model$expvar) <-
      paste("Comp", seq_len(ncomp))

   # add class name, call and info and return
   class(model) <- c("mcr", "mcrals")
   model$call <- match.call()
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

   if (length(object$spec.constraints) > 0) {
      cat("\nConstraints for spectra:\n")
      cat(paste(" - ", sapply(object$spec.constraints, function(x) x$name), collapse = "\n"))
      cat("\n")
   }

   if (length(object$cont.constraints) > 0) {
      cat("\nConstraints for contributions:\n")
      cat(paste(" - ", sapply(object$cont.constraints, function(x) x$name), collapse = "\n"))
      cat("\n")
   }

   cat("\n")
   data <- cbind(
      round(object$expvar, 2),
      round(object$cumexpvar, 2)
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
#' @param cont.constraints
#' a list with constraints to be applied to contributions  (see details).
#' @param spec.constraints
#' a list with constraints to be applied to spectra  (see details).
#' @param spec.ini
#' a matrix with initial estimation of the pure components spectra.
#' @param cont.solver
#' which function to use as a solver for resolving of pure components contributions (see detials).
#' @param spec.solver
#' which function to use as a solver for resolving of pure components spectra (see detials).
#' @param max.niter
#' maximum number of iterations.
#' @param tol
#' tolerance, when explained variance change is smaller than this value, iterations stop.
#' @param verbose
#' logical, if TRUE information about every iteration will be shown.
#'
#' @return
#' The function returns a list with with following fields:
#' \item{ncomp }{number of pure components.}
#' \item{resspec}{matrix with resolved spectra.}
#' \item{rescont}{matrix with resolved contributions.}
#' \item{cont.constraints}{list with contribution constraints provided by user.}
#' \item{spec.constraints}{list with spectra constraints provided by user.}
#' \item{max.niter}{maximum number of iterations}
#'
#' @export
mcrals.cal <- function(D, ncomp, cont.constraints, spec.constraints, spec.ini, cont.solver,
   spec.solver, max.niter, tol, verbose) {

   attrs <- mda.getattr(D)
   exclrows <- attrs$exclrows
   exclcols <- attrs$exclcols

   # get dimensions
   nobj <- nrow(D)
   nvar <- ncol(D)

   # get indices for excluded columns if provided
   colind <- seq_len(nvar)
   rowind <- seq_len(nobj)

   # remove hidden rows if any
   if (!is.null(exclrows)) {
      D <- D[-exclrows, , drop = FALSE]
      rowind <- rowind[-exclrows]
   }

   # remove hidden columns if any
   if (!is.null(exclcols)) {
      D <- D[, -exclcols, drop = FALSE]
      spec.ini <- spec.ini[-exclcols, , drop = FALSE]
      colind <- colind[-exclcols]
   }

   # main loop for ALS
   totvar <- sum(D^2)
   var <- 0

   if (verbose) {
      cat(sprintf("\nStarting iterations (max.niter = %d):\n", max.niter))
   }

   S <- spec.ini
   for (i in seq_len(max.niter)) {

      ## compute C and apply constraints
      C <- cont.solver(D, S)
      for (cc in cont.constraints) {
         C <- employ(cc, C)
      }

      ## compute S and apply constraints
      S <- spec.solver(D, C)
      for (sc in spec.constraints) {
         S <- employ(sc, S)
      }

      var_old <- var
      var <- 1 - sum((D - tcrossprod(C, S))^2) / totvar
      if ( (var - var_old) < tol) {
         if (verbose) cat("No more improvements.\n")
         break
      }

      if (verbose) {
         cat(sprintf("Iteration %4d, R2 = %7.4f\n", i, var * 100))
      }
   }


   # if there were excluded rows or columns, handle this
   Ct <- matrix(0, nobj, ncomp)
   Ct[rowind, ] <- C

   St <- matrix(0, ncomp, nvar)
   St[, colind] <- t(S)

   # add attributes
   Ct <- mda.setattr(Ct, attrs, "row")
   St <- mda.setattr(St, attrs, "col")

   if (is.null(attr(Ct, "xaxis.name"))) attr(Ct, "xaxis.name") <- "Observations"
   attr(Ct, "name") <- "Resolved contributions"
   attr(St, "name") <- "Resolved spectra"

   # combine the results with model object
   return(
      list(
         ncomp = ncomp,
         rescont = Ct,
         resspec = mda.t(St),
         cont.constraints = cont.constraints,
         spec.constraints = spec.constraints,
         max.niter = max.niter
      )
   )
}

#' Ordinary least squares
#'
#' @param D
#' a matrix
#' @param B
#' a matrix
#'
#' @details
#' Computes OLS solution for D = AB' (or D' = AB'), where D, A are known
#'
#' @export
mcrals.ols <- function(D, A) {
   if (ncol(D) != nrow(A)) D <- t(D)
   return(D %*% (A %*% solve(crossprod(A))))
}


#' Non-negative least squares
#'
#' @param D
#' a matrix
#' @param A
#' a matrix
#' @param tol
#' tolerance parameter for algorithm convergence
#'
#' @details
#' Computes NNLS solution for B: D = AB' subject to B >= 0. Implements
#' the active-set based algorithm proposed by Lawson and Hanson [1].
#'
#' @references
#' 1.  Lawson, Charles L.; Hanson, Richard J. (1995). Solving Least Squares Problems. SIAM.
#'
#' @export
mcrals.nnls <- function(D, A,
   tol = 10 * .Machine$double.eps * as.numeric(sqrt(crossprod(A[, 1]))) * nrow(A)) {

   if (nrow(D) != nrow(A)) D <- t(D)

   nvars <- ncol(D)
   ncomp <- ncol(A)
   itmax <- 30 * ncomp
   B <- matrix(0, nvars, ncomp)

   for (v in seq_len(nvars)) {

      d <- D[, v, drop = FALSE]

      # initialize vector of zeros and infinites
      nz <- rep(0, ncomp)
      wz <- nz

      # vector of non-active columns (positive set)
      p <- rep(FALSE, ncomp)

      # vector of active columns (zero or active set)
      z <- rep(TRUE, ncomp)

      # initialize vector with current solution
      b <- nz

      # compute initial residuals and values for Lagrange multipliers
      E <- d - A %*% b
      w <- t(A) %*% E

      iter <- 0

      # outer loop
      while (any(z) && any(w[z] > tol)) {

         # set intermediate solution for b to zeros
         b.hat <- nz

         # create vector with Lagrange multipliers for active set (-Inf for positive set)
         wz[p] <- -Inf
         wz[z] <- w[z]

         # find index with larges multiplier
         t <- which.max(wz)

         # move variable corrsponding to the index from active set to positive set
         p[t] <- TRUE
         z[t] <- FALSE

         # compute intermediate solution using only positive set
         b.hat[p] <- mcrals.ols(d, A[, p, drop = FALSE])

         # inner loop to remove elements from the positive set
         while (any(b.hat[p] <= 0)) {
            iter <- iter + 1

            if (iter > itmax) {
               # too many iterations, exitint with current solution
               return(b.hat)
            }

            # find which values in current solution are negative but they are in positive lest
            q <- (b.hat <= 0) & p

            # adjust solution to make them non-negative
            alpha <- min(b[q]/(b[q] - b.hat[q]))
            b <- b + alpha * (b.hat - b)

            # reset the active set and positive set according to the current solution
            z <- ((abs(b) < tol) & p) | z
            p <- !z

            # compute intermedate solution again
            b.hat <- nz
            b.hat[p] <- mcrals.ols(d, A[, p, drop = FALSE])
         }

         # set current solution to intermediate solution and recompute the multipliers
         b <- b.hat
         E <- d - A %*% b
         w <- crossprod(A, E)
      }

      B[v, ] <- b
   }

   return(B)
}
