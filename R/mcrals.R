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
#' @param cont.forced
#' a matrix which allows to force some of the concentration values (see details).
#' @param spec.forced
#' a matrix which allows to force some of the spectra values (see details).
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
#' Parameters \code{cont.force} and \code{spec.force} allows you to force some parts of the
#' contributions or the spectra to be equal to particular pre-defined values. In this case you need
#' to provide the parameters (or just one of them) in form of a matrix. For example \code{cont.force}
#' should have as many rows as many you have in the original spectral data \code{x} and as many
#' columns as many pure components you want to resolve. Feel all values of this matrix with
#' \code{NA} and the values you want to force with real numbers. For example if you know that in
#' the first measurement concentration of 2 and 3 components was zero, set the corresponding
#' values of \code{cont.force} to zero. See also the last case in the examples section.
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
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @seealso
#' Methods for \code{mcrals} objects:
#' \tabular{ll}{
#'    \code{summary.mcrals} \tab shows some statistics for the case.\cr
#'    \code{\link{predict.mcrals}} \tab computes contributions by projection of new spectra to
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
#'
#' \donttest{
#'
#' library(mdatools)
#'
#' # resolve mixture of carbonhydrates Raman spectra
#'
#' data(carbs)
#'
#' # define constraints for contributions
#' cc <- list(
#'    constraint("nonneg")
#' )
#'
#' # define constraints for spectra
#' cs <- list(
#'    constraint("nonneg"),
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
#' # This example shows how to force some of the contribution values
#' # First of all we combine the matrix with mixtures and the pure spectra, so the pure
#' # spectra are on top of the combined matrix
#' Dplus <- mda.rbind(mda.t(carbs$S), carbs$D)
#'
#' # since we know that concentration of C2 and C3 is zero in the first row (it is a pure
#' # spectrum of first component), we can force them to be zero in the optimization procedure.
#' # Similarly we can do this for second and third rows.
#'
#' cont.forced <- matrix(NA, nrow(Dplus), 3)
#' cont.forced[1, ] <- c(NA, 0, 0)
#' cont.forced[2, ] <- c(0, NA, 0)
#' cont.forced[3, ] <- c(0, 0, NA)
#'
#' m <- mcrals(Dplus, 3, cont.forced = cont.forced, cont.constraints = cc, spec.constraints = cs)
#' plot(m)
#'
#' # See bookdown tutorial for more details.
#'
#' }
#'
#' @export
mcrals <- function(x, ncomp,
   cont.constraints = list(),
   spec.constraints = list(),
   spec.ini = matrix(runif(ncol(x) * ncomp), ncol(x), ncomp),
   cont.forced = matrix(NA, nrow(x), ncomp),
   spec.forced = matrix(NA, ncol(x), ncomp),
   cont.solver = mcrals.nnls,
   spec.solver = mcrals.nnls,
   exclrows = NULL, exclcols = NULL, verbose = FALSE,
   max.niter = 100, tol = 10^-6, info = "") {

   stopifnot("Parameter 'max.niter' should be positive." = max.niter > 0)

   if (any(spec.ini < 0)) {
      warning("Initial estimation of pure spectra has negative numbers, it can cause problems.")
   }

   # get pure variables and unmix data
   x <- prepCalData(x, exclrows, exclcols, min.nrows = 2, min.ncols = 2)

   # check dimensions of cont.forced and spec.forced
   stopifnot("Number of rows in 'cont.forced' should be the same as number of rows in 'x'." = nrow(cont.forced) == nrow(x))
   stopifnot("Number of columns in 'cont.forced' should be the same as 'ncomp'." = ncol(cont.forced) == ncomp)
   stopifnot("Number of rows in 'spec.forced' should be the same as number of columns in 'x'." = nrow(spec.forced) == ncol(x))
   stopifnot("Number of columns in 'spec.forced' should be the same as 'ncomp'." = ncol(cont.forced) == ncomp)

   model <- mcrals.cal(
      x, ncomp,
      cont.constraints = cont.constraints,
      spec.constraints = spec.constraints,
      spec.ini = spec.ini,
      cont.forced = cont.forced,
      spec.forced = spec.forced,
      cont.solver = cont.solver,
      spec.solver = spec.solver,
      max.niter = max.niter,
      tol = tol,
      verbose = verbose
   )

   # compute explained variance
   model$variance <- getVariance.mcr(model, x)

   # add class name, call and info and return
   class(model) <- c("mcr", "mcrals")
   model$call <- match.call()
   model$info <- info

   return(model)
}

#' MCR ALS predictions
#'
#' @description
#' Applies MCR-ALS model to a new set of spectra and returns matrix with contributions.
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
predict.mcrals <- function(object, x, ...) {
   attrs <- mda.getattr(x)
   Ct <- object$cont.solver(x, object$resspec)
   Ct <- mda.setattr(Ct, attrs, "rows")
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

   cat("\nConstraints for spectra:\n")
   if (length(object$spec.constraints) > 0) {
      cat(paste(" - ", sapply(object$spec.constraints, function(x) x$name), collapse = "\n"))
      cat("\n")
   } else {
      cat(" - none\n\n")
   }

   cat("\nConstraints for contributions:\n")
   if (length(object$cont.constraints) > 0) {
      cat(paste(" - ", sapply(object$cont.constraints, function(x) x$name), collapse = "\n"))
      cat("\n")
   } else {
      cat(" - none\n\n")
   }

   cat("\n")
   out <- round(t(object$variance[1:2, ]), 2)
   rownames(out) <- colnames(object$variance)
   colnames(out) <- c("Expvar", "Cumexpvar")
   show(out)
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
#' @param cont.forced
#' a matrix which allows to force some of the concentration values (see details).
#' @param spec.forced
#' a matrix which allows to force some of the spectra values (see details).
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
mcrals.cal <- function(D, ncomp, cont.constraints, spec.constraints, spec.ini,
   cont.forced, spec.forced, cont.solver, spec.solver, max.niter, tol, verbose) {

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

   cont.forced.ind <- !is.na(cont.forced)
   spec.forced.ind <- !is.na(spec.forced)

   St <- spec.ini
   for (i in seq_len(max.niter)) {

      ## resolve contributions
      Ct <- tryCatch(
         cont.solver(D, St),
         error = function(e) {
            #print(e)
            stop("Unable to resolve the components, perhaps 'ncomp' is too large or initial estimates for spectra are not good enough.", call. = FALSE)
         }
      )

      ## force the user defined values
      if (any(cont.forced.ind)) {
         Ct[cont.forced.ind] <- cont.forced[cont.forced.ind]
      }

      ## apply constraints to the resolved contributions
      for (cc in cont.constraints) {
         Ct <- employ.constraint(cc, x = Ct, d = D)
      }

      ## resolve spectra
      St <- tryCatch(
         spec.solver(D, Ct),
         error = function(e) {
            #print(e)
            stop("Unable to resolve the components, perhaps 'ncomp' is too large or initial estimates for spectra are not good enough.", call. = FALSE)
         }
      )

      ## force the user defined values
      if (any(spec.forced.ind)) {
         St[spec.forced.ind] <- spec.forced[spec.forced.ind]
      }

      ## apply constraints to the resolved spectra
      for (sc in spec.constraints) {
         St <- employ.constraint(sc, x = St, d = D)
      }

      if (verbose) {
         cat(sprintf("Iteration %4d, R2 = %.6f\n", i, var))
      }

      var_old <- var
      var <- tryCatch(
         1 - sum((D - tcrossprod(cont.solver(D, St), St))^2) / totvar,
         error = function(e) {
            print(e)
            stop("Unable to resolve the components, perhaps 'ncomp' is too large.\n
               or initial estimates for spectra are not good enough.", call. = FALSE)
         }
      )

      if ( (var - var_old) < tol) {
         if (verbose) cat("No more improvements.\n")
         break
      }

   }

   # predict final Ct values based on last version of the St
   Ct <- tryCatch(
      cont.solver(D, St),
      error = function(e) {
         print(e)
         stop("Unable to resolve the components, perhaps 'ncomp' is too large.\n
            or initial estimates for spectra are not good enough.", call. = FALSE)
      }
   )

   # if there were excluded rows or columns, handle this
   cont <- matrix(0, nobj, ncomp)
   cont[rowind, ] <- Ct

   spec <- matrix(0, ncomp, nvar)
   spec[, colind] <- t(St)

   # add names for components
   rownames(spec) <- colnames(cont) <- paste("Comp", seq_len(ncomp))

   # add attributes
   cont <- mda.setattr(cont, attrs, "row")
   spec <- mda.setattr(spec, attrs, "col")

   if (is.null(attr(cont, "yaxis.name"))) attr(cont, "yaxis.name") <- "Observations"
   attr(cont, "name") <- "Resolved contributions"
   attr(spec, "name") <- "Resolved spectra"

   # combine the results with model object
   return(
      list(
         ncomp = ncomp,
         rescont = cont,
         resspec = mda.t(spec),
         cont.constraints = cont.constraints,
         spec.constraints = spec.constraints,
         cont.solver = cont.solver,
         spec.colver = spec.solver,
         max.niter = max.niter
      )
   )
}

#' Ordinary least squares
#'
#' @param D
#' a matrix
#' @param A
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

   nvar <- ncol(D)
   ncomp <- ncol(A)
   itmax <- 30 * ncomp
   B <- matrix(0, nvar, ncomp)

   for (v in seq_len(nvar)) {

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


#' Fast combinatorial non-negative least squares
#'
#' @param D
#' a matrix
#' @param A
#' a matrix
#' @param tol
#' tolerance parameter for algorithm convergence
#'
#' @details
#' Computes Fast combinatorial NNLS solution for B: D = AB' subject to B >= 0. Implements
#' the method described in [1].
#'
#' @references
#' 1. Van Benthem, M.H. and Keenan, M.R. (2004), Fast algorithm for the solution of large scale
#' non-negativity-constrained least squares problems. J. Chemometrics, 18: 441-450.
#' doi:10.1002/cem.889
#'
#' @export
mcrals.fcnnls <- function(D, A,
   tol = 10 * .Machine$double.eps * as.numeric(sqrt(crossprod(A[, 1]))) * nrow(A)) {

   ind2sub <- function(size, ind) {
      return(list(
         r = ((ind - 1) %% size[1]) + 1,
         c = floor((ind - 1) / size[1]) + 1
      ))
   }

   sub2ind <- function(size, r, c) {
      return((c - 1) * size[1] + r)
   }

   cssls <- function(AtD, AtA, P.set = NULL) {

      if (is.null(P.set) || length(P.set) == 0 || all(P.set)) {
         return(solve(AtA, AtD))
      }

      B <- matrix(0, nrow(AtD), ncol(AtD))
      ncomp <- nrow(P.set)
      nvar <- ncol(P.set)

      codedPset <- (2^((ncomp - 1):0)) %*% P.set

      sortedEset <- order(codedPset)
      sortedPset <- codedPset[sortedEset]

      breaks <- diff(sortedPset)
      breakIdx <- c(0, which(breaks > 0), nvar)

      for (k in seq_len(length(breakIdx) - 1)) {
         cols2solve <- sortedEset[(breakIdx[k] + 1) : breakIdx[k+1]]
         vars <- P.set[, sortedEset[breakIdx[k] + 1]]
         B[vars, cols2solve] <- solve(
            AtA[vars, vars, drop = FALSE],
            AtD[vars, cols2solve, drop = FALSE]
         )
      }

      return(B)
   }


   if (nrow(D) != nrow(A)) D <- t(D)

   nvar <- ncol(D)
   ncomp <- ncol(A)
   nobj <- nrow(A)

   W <- matrix(0, ncomp, nvar)
   iter <- 0
   maxiter <- 3 * ncomp

   # precompute parts of pseudoinverse
   AtA <- crossprod(A)
   AtD <- crossprod(A, D)

   # Obtain the initial feasible solution and corresponding passive set
   B <- cssls(AtD, AtA)
   P.set <- B > 0
   B[!P.set] <- 0
   B.hat <- B

   # F.set contains indices of variables where there is at least one
   # negative coefficient (active columns)
   F.set <- which(!apply(P.set, 2, all))

   # active set algorithm for NNLS
   while (length(F.set) != 0) {
      # get solution for the variables from active columns
      B[, F.set] <- cssls(AtD[, F.set, drop = FALSE], AtA, P.set[, F.set, drop = FALSE])

      # find variables from active columns where at least one coefficient is negative
      # and keep them in H.set (so it is a part of F.set which was not improved)
      H.set <- F.set[which(apply(B[, F.set, drop = FALSE] < 0, 2, any))]

      # make infeasible solutions feasible (standard NNLS inner loop)
      # get number of currently active columns
      nHset <- length(H.set)
      if (nHset > 0) {

         while (length(H.set) != 0 && (iter < maxiter)) {
            iter <- iter + 1
            alpha <- matrix(Inf, ncomp, nHset)

            # find rows and columns of all negative coefficients in the passive set
            nvInd <- which(P.set[, H.set, drop = FALSE] & (B[, H.set, drop = FALSE] < 0),
               arr.ind = TRUE)
            nvR <- nvInd[, 1]
            nvC <- nvInd[, 2]

            # convert rows and columns into indices for every negative coefficient in
            # the passive set - indices based on the size of H.set and alpha
            hIdx <- sub2ind(dim(alpha), nvR, nvC)

            # get similar indices but based on the size of B
            bIdx <- sub2ind(dim(B), nvR, H.set[nvC])

            # compute alpha for the negative coefficients
            alpha[hIdx] <- B.hat[bIdx] / (B.hat[bIdx] - B[bIdx])

            # find smallest values for alpha in each column
            alphaMin <- apply(alpha, 2, min)

            # find row indices which contain the smallest values
            alphaMinIdx <- apply(alpha, 2, which.min)

            # replace alpha values with the smallest ones for each column
            alpha <- matrix(alphaMin, ncomp, nHset, byrow = TRUE)

            # adjust coefficients from the problematic columns
            B.hat[, H.set] <- B.hat[, H.set] - alpha * (B.hat[, H.set] - B[, H.set])

            # get indices of the coefficients with smallest alpha in each column
            # and set them to 0 plus remove them from the passive set
            idx2zero <- sub2ind(dim(B.hat), alphaMinIdx, H.set)
            B.hat[idx2zero] <- 0
            P.set[idx2zero] <- FALSE

            B[, H.set] <- cssls(AtD[, H.set, drop = FALSE], AtA, P.set[, H.set, drop = FALSE])
            H.set <- which(apply(B < 0, 2, any))
            nHset <- length(H.set)
         }

      }

      # Check solutions for optimality
      W[, F.set] <- AtD[, F.set, drop = FALSE] - AtA %*% B[, F.set, drop = FALSE]
      J.set <- apply(((!P.set[, F.set, drop = FALSE]) * W[, F.set, drop = FALSE]) <= 0, 2, all)
      F.set <- F.set[!J.set]

      # For non-optimal solutions, add the appropriate variable to Pset
      if (length(F.set) > 0) {
         mxidx <- apply((!P.set[, F.set, drop = FALSE]) * W[, F.set, drop = FALSE], 2, which.max)
         P.set[sub2ind(c(ncomp, nvar), mxidx, F.set)] <- TRUE
         B.hat[, F.set] <- B[, F.set]
      }

      if (iter == maxiter) {
         warning("Maximum number of iterations is reached, quit with current solution.")
         return(t(B))
      }
   }

   return(t(B))
}