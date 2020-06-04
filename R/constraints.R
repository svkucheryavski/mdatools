###################
# Static methods  #
###################
getImplementedConstraints <- function() {
   list(
      "non-negativity" = list(
         method = constraintsNN,
         params = list(),
         info = "Non-negativity constraint."
      )
   )
}


#' Method for non-negativity constraint
#'
#' @param x
#' data matrix (spectra or contributions)
constraintsNN <- function(x) {
   x[x < 0] = 0
   return(x)
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

      method <- item$method
   } else {
      # user defined constraint, check that it works
      res <- tryCatch(
         do.call(method, c(matrix(runif(25, 5, 5)), params)),
         error = function(m) stop("The method you provided raises an error: \n", m),
         warning = function(m) stop("The method you provided raises a warning: \n", m)
      )
   }

   obj <- list(
      name = name,
      method = method,
      params = params
   )

   class(obj) <- c("constraint")
}

apply.constraint <- function(obj, x) {
   return(do.call(obj$method, c(x, obj$params)))
}

