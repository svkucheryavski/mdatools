#' Define parameters based on 'cv' value
#'
#' @param cv
#' settings for cross-validation provided by user
#' @param nobj
#' number of objects in calibration set
#'
crossval.getParams <- function(cv, nobj) {

   nrep <- 1

   # random
   if (is.numeric(cv)) {
      return(
         list(
            type = "rand",
            nrep = 1,
            nseg = if (cv == 1) nobj else cv
         )
      )
   }

   # leave one out
   type <- cv[[1]]
   if (type == "loo") {
      return(
         list(
            type = "loo",
            nrep = 1,
            nseg = nobj
         )
      )
   }

   # venetian blinds
   nseg <- cv[[2]]
   if (type == "ven") {
      return(
         list(
            type = "ven",
            nrep = nrep,
            nseg = nseg
         )
      )
   }

   nrep <- if (length(cv) == 3) cv[[3]] else 1
   return(
      list(
         type = type,
         nrep = nrep,
         nseg = nseg
      )
   )
}

#' Generate sequence of indices for cross-validation
#'
#' @description
#' Generates and returns sequence of object indices for each segment in random segmented
#' cross-validation
#'
#' @param nobj
#' number of objects in a dataset
#' @param cv
#' cross-validation settings, can be a number or a list. If cv is a number, it will be
#' used as a number of segments for random cross-validation (if cv = 1, full cross-validation
#' will be preformed), if it is a list, the following syntax can be used:
#' cv = list('rand', nseg, nrep) for random repeated cross-validation with nseg segments and nrep
#' repetitions or cv = list('ven', nseg) for systematic splits to nseg segments ('venetian blinds').
#' @param resp
#' vector with response values to use in case of venetian blinds
#'
#' @return
#' matrix with object indices for each segment
#'
#' @export
crossval <- function(cv = 1, nobj = NULL, resp = NULL) {

   # get cross-validation parameters
   if (is.null(nobj)) nobj <- length(resp)

   # if user already provided matrix with values - return it
   if (is.numeric(cv) && length(cv) == nobj) return(as.matrix(cv))

   p <- crossval.getParams(cv = cv, nobj = nobj)
   if (!(p$type %in% c("rand", "ven", "loo"))) {
      stop("Wrong name for cross-validation method.")
   }

   # check number of repetitions
   if (p$nrep < 1 || p$nrep > 100) {
      stop("Wrong value for cv repetitions (should be between 1 and 100).")
   }

   # check number of segments
   if (p$nseg < 2 || p$nseg > nobj) {
      stop("Wrong value for number of segments (should be between 2 and number of objects).")
   }

   if (p$type == "loo") {
      return(matrix(seq_len(nobj), ncol = 1))
   }

   if (p$type == "rand") {
      return(sapply(seq_len(p$nrep), function(i) rep(seq_len(p$nseg), length.out = nobj)[sample(nobj)]))
   }

   if (p$type == "ven") {

      ind <- if (is.null(resp)) seq_len(nobj) else order(order(resp))
      return(matrix(rep(seq_len(p$nseg), length.out = nobj)[ind], ncol = 1))
   }

   stop("Something went wrong.")
}

#' String with description of cross-validation method
#'
#' @param cv
#' a list with cross-validation settings
#'
#' @return
#' a string with the description text
#'
crossval.str <- function(cv) {

   if (length(cv) == 0) return("none")

   if (is.numeric(cv)) {

      if (length(cv) > 1) {
         return (sprintf("user defined with %d segments", length(unique((cv)))))
      }

      return(
         if (cv == 1) "full (leave one out)"
         else sprintf("random with %d segments", cv)
      )
   }


   type <- cv[[1]]
   if (type == "loo") {
      return("full (leave one out)")
   }

   if (type == "ven") {
      return(sprintf("venetian blinds with %.0f segments", cv[[2]]))
   }

   return(
      sprintf("random with %.0f segments%s",
         cv[[2]], if (length(cv) == 3) paste(" and", cv[[3]], "repetitions") else "")
   )
}
