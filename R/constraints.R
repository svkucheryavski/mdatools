###################
# Static methods  #
###################

#' Shows a list with implemented constraints
#'
#' @export
getImplementedConstraints <- function() {
   list(
      "nonneg" = list(
         name = "nonneg",
         method = constraintNonNegativity,
         params = list(),
         params.info = list(),
         info = "Non-negativity (sets negative values to zero)"
      ),
      "unimod" = list(
         name = "unimod",
         method = constraintUnimod,
         params = list(tol = 0),
         params.info = list(tol = "tolerance (between 0 and 1)"),
         info = "Unimodality (forces contribution or spectral profile to have a single maximum)"
      ),
      "closure" = list(
         name = "closure",
         method = constraintClosure,
         params = list(sum = 1),
         params.info = list(sum = "value, the data rows should sum up to"),
         info = "Closure (forces contributions or spectral profiles sum up to constant value)"
      ),
      "norm" = list(
         name = "norm",
         method = constraintNorm,
         params = list(type = "length"),
         params.info = list(type = "type of normalization: 'length', 'area' or 'sum'"),
         info = "Normalization (normalize spectra or contributions)"
      ),
      "angle" = list(
         name = "angle",
         method = constraintAngle,
         params = list(weight = 0.05),
         params.info = list(weight = "how much of mean will be added: between 0 and 1"),
         info = "Angle (increases contrast among resolved spectra or contributions)"
      )
   )
}

#' Shows information about all implemented constraints
#'
#' @export
constraints.list <- function() {
   constraints <- getImplementedConstraints()

   cat("\nList of constraints available for mcrals():\n")
   lapply(constraints, function(c) {
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

#' Method for closure constraint
#'
#' @description
#' Force rows of data sum up to given value
#'
#' @param x
#' data matrix (spectra or contributions)
#' @param d
#' matrix with the original spectral values
#' @param sum
#' which value the specra or contributions should sum up to
#'
#' @export
constraintClosure <- function(x, d, sum = 1) {
   stopifnot("Parameter 'sum' should be positive number." = sum > 0 )

   # if all values in a row are zeros set sums to one so they will not be scaled
   rsums <- rowSums(x)
   rsums[rsums == 0] <- 1

   # scale the values so for evert row they sum up to 1
   s <- diag(sum / rsums, nrow(x), nrow(x))
   return(s %*% x)
}

#' Method for unimodality constraint
#'
#' @description
#' forces column of matrix to have one maximum each
#'
#' @param x
#' data matrix (spectra or contributions)
#' @param d
#' matrix with the original spectral values
#' @param tol
#' tolerance (value between 0 and 1) to take make method stable to small fluctuations
#'
#' @export
constraintUnimod <- function(x, d, tol = 0) {

   f <- function(y, max, indseq, step) {

      for (i in indseq) {
         if (y[i] <= max) {
            max <- y[i]
         } else if (y[i] > max * (1 + tol)) {
            y[i] <- y[i + step]
            max <- y[i]
         }
      }

      return(y)
   }

   # find maximum for each component
   peak.ind <- apply(x, 2, which.max)
   nvar <- nrow(x)

   # process each component separately
   for (a in seq_len(ncol(x))) {

      if (peak.ind[a] != 1) {
          # flatten peaks to the left of maximum
          left_part <- (peak.ind[a] - 1):1
          x[, a] <- f(x[, a], max = x[peak.ind[a], a], indseq = left_part, step = +1)
      }

      if (peak.ind[a] != nvar) {
          # flatten peaks to the right of maximum
          right_part <- (peak.ind[a] + 1):nvar
          x[, a] <- f(x[, a], max = x[peak.ind[a], a], indseq = right_part, step = -1)
      }
   }

   return(x)
}

#' Method for non-negativity constraint
#'
#' @description
#' Set all negative values in the matrix to 0
#'
#' @param x
#' data matrix (spectra or contributions)
#' @param d
#' matrix with the original spectral values
#'
#' @export
constraintNonNegativity <- function(x, d) {
   x[x < 0] <- 0
   return(x)
}

#' Method for normalization constraint
#'
#' @description
#' Normalize rows of matrix to unit length or area
#'
#' @param x
#' data matrix (spectra or contributions)
#' @param d
#' matrix with the original spectral values
#' @param type
#' type of normalization ("area", "length" or "sum")
#'
#' @export
constraintNorm <- function(x, d, type = "length") {
   types <- c("area", "length", "sum")
   stopifnot("Parameter 'type' should be either 'area', 'length' or 'sum'." = type %in% types )
   return(t(prep.norm(t(x), type)))
}

#' Method for angle constraint
#'
#' @description
#' Adds a small portion of mean to contributions or spectra to increase contrast
#'
#' @param x
#' data matrix (spectra or contributions)
#' @param d
#' matrix with the original spectral values
#' @param weight
#' how many percent of mean to add (between 0 and 1)
#'
#' @export
constraintAngle <- function(x, d, weight = 0.05) {
   stopifnot("Parameter 'weight' should be between 0 and 1." = weight >= 0 && weight <= 1 )

   # compute mean values for rows or columns of d
   m <- apply(d, ifelse(nrow(x) == ncol(d), 2, 1), mean)

   # scale mean values to unit length
   m <- m / sqrt(sum(m^2))

   # scale data to unit length
   x <- t(prep.norm(t(x), "length"))

   return((1 - weight) * x + matrix(m * weight, nrow(x), ncol(x)))
}

#' Class for MCR-ALS constraint
#'
#' @param name
#' short text with name for the constraint
#' @param params
#' a list with parameters for the constraint method (if NULL - default parameters will be used)
#' @param method
#' method to call when applying the constraint, provide it only for user defined constraints
#'
#' @details
#' Use this class to create constraints and add them to a list for MCR-ALS curve resuliton (see
#' \code{\link{mcrals}}). Either provide name and parameters to one of the existing constraint
#' implementations or make your own. See the list of implemented constraints by running
#' \code{constraints()}
#'
#' For your own constraint you need to create a method, which takes matrix with values (either
#' spectra or contributions being resolved) as the first argument, does something and then return
#' a matrix with the same dimension as the result. The method can have any number of optional
#' parameters.
#'
#' See help for \code{\link{mcrals}} or Bookdown tutorial for details.
#'
#' @export
constraint <- function(name, params = NULL, method = NULL) {

   if (is.null(method)) {
      # assuming it is one of the standard constraints
      # 1. first check name
      item <- getImplementedConstraints()[[name]]
      stopifnot("Either name of constraint is wrong or you need to provide a method if the
         constraint is user defined." = !is.null(item))

      # 2. check the parameters
      if (is.null(params)) params <- item$params
      if (length(params) > 0 && !(names(params) %in% names(item$params))) {
         stop("Provided constraint parameters have wrong name.")
      }

      method <- item$method
   } else {
      # user defined constraint, check that it works
      res <- tryCatch(
         do.call(method, c(matrix(runif(25, 5, 10)), params)),
         error = function(m) stop("The method you provided raises an error: \n", m),
         warning = function(m) stop("The method you provided raises a warning: \n", m)
      )

      stopifnot("The method you provided does not return matrix with correct dimension." =
         dim(res) == c(5, 10))
   }

   obj <- list(
      name = name,
      method = method,
      params = params
   )

   class(obj) <- c("constraint")
   return(obj)
}

#' Applies constraint to a dataset
#'
#' @param obj
#' object with constraint
#' @param x
#' matrix with pure spectra or contributions
#' @param d
#' matrix with original spectral values
#' @param ...
#' other arguments
#'
#' @export
employ.constraint <- function(obj, x, d, ...) {
   return(do.call(obj$method, c(list(x = x, d = d), obj$params)))
}

