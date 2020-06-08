###################
# Static methods  #
###################
getImplementedConstraints <- function() {
   list(
      "non-negativity" = list(
         name = "non-negativity",
         method = constraintNonNegativity,
         params = list(),
         info = "Non-negativity constraint."
      ),
      "norm" = list(
         name = "norm",
         method = constraintNorm,
         params = list(type = "area"),
         info = "Normalization constraint (type can be either 'area' or 'length')."
      )
   )
}


#' Method for non-negativity constraint
#'
#' @description
#' Set all negative values in the matrix to 0
#'
#' @param x
#' data matrix (spectra or contributions)
#'
constraintNonNegativity <- function(x) {
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
#' type of normalization ("area" or "length")
#'
constraintNorm <- function(x, type = "area") {
   types <- c("area", "length")
   stopifnot("Parameter 'type' should be either 'area' or 'length'." = type %in% types )
   norma <- if(type == "length") sqrt(colSums(x^2)) else colSums(abs(x))
   return(x %*% diag(1/norma, nrow = ncol(x), ncol = ncol(x)))
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
#' matrix with data values
#'
#' @export
employ.constraint <- function(obj, x) {
   return(do.call(obj$method, c(list(x = x), obj$params)))
}

