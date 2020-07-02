###################
# Static methods  #
###################

#' Shows a list with implemented constraints
#'
#' @export
getImplementedConstraints <- function() {
   list(
      "non-negativity" = list(
         name = "non-negativity",
         method = constraintNonNegativity,
         params = list(),
         params.info = list(),
         info = "Non-negativity (sets negative values to zero)"
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

#' Method for non-negativity constraint
#'
#' @description
#' Set all negative values in the matrix to 0
#'
#' @param x
#' data matrix (spectra or contributions)
#'
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
#' @param type
#' type of normalization ("area", "length" or "sum")
#'
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
      stopifnot("Either name of constraint is wrong you you need to provide a method if the
         consttaint is user defined." = !is.null(item))

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
#'
#' @export
employ.constraint <- function(obj, x, d) {
   return(do.call(obj$method, c(list(x = x, d = d), obj$params)))
}

