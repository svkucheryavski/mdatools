#' Results of PCA decomposition
#'
#' @description
#' \code{pcares} is used to store and visualise results for PCA decomposition.
#'
#' @param ...
#' all arguments supported by \code{ldecomp}.
#'
#' @details
#' In fact \code{pcares} is a wrapper for \code{\link{ldecomp}} - general class for storing
#' results for linear decomposition X = TP' + E. So, most of the methods, arguments and
#' returned values are inherited from \code{ldecomp}.
#'
#' There is no need to create a \code{pcares} object manually, it is created automatically when
#' build a PCA model (see \code{\link{pca}}) or apply the model to a new data (see
#' \code{\link{predict.pca}}). The object can be used to show summary and plots for the results.
#'
#' It is assumed that data is a matrix or data frame with I rows and J columns.
#'
#' @return
#' Returns an object (list) of class \code{pcares} and \code{ldecomp} with following fields:
#' \item{scores }{matrix with score values (I x A).}
#' \item{residuals }{matrix with data residuals (I x J).}
#' \item{T2 }{matrix with score distances (I x A).}
#' \item{Q }{matrix with orthogonal distances (I x A).}
#' \item{ncomp.selected }{selected number of components.}
#' \item{expvar }{explained variance for each component.}
#' \item{cumexpvar }{cumulative explained variance.}
#'
#' @seealso
#' Methods for \code{pcares} objects:
#' \tabular{ll}{
#'  \code{print.pcares} \tab shows information about the object.\cr
#'  \code{summary.pcares} \tab shows statistics for the PCA results.\cr
#' }
#'
#' Methods, inherited from \code{\link{ldecomp}} class:
#' \tabular{ll}{
#'  \code{\link{plotScores.ldecomp}} \tab makes scores plot.\cr
#'  \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
#'  \code{\link{plotCumVariance.ldecomp}} \tab makes cumulative explained variance plot.\cr
#'  \code{\link{plotDistances.ldecomp}} \tab makes Q vs. T2 distance plot.\cr
#' }
#'
#' Check also \code{\link{pca}} and \code{\link{ldecomp}}.
#'
#' @examples
#' ### Examples for PCA results class
#'
#' library(mdatools)
#'
#' ## 1. Make a model for every odd row of People data
#' ## and apply it to the objects from every even row
#'
#' data(people)
#' x = people[seq(1, 32, 2), ]
#' x.new = people[seq(2, 32, 2), ]
#'
#' model = pca(people, scale = TRUE, info = "Simple PCA model")
#' model = selectCompNum(model, 4)
#'
#' res = predict(model, x.new)
#' summary(res)
#' plot(res)
#'
#' ## 2. Make PCA model for People data with autoscaling
#' ## and full cross-validation and get calibration results
#'
#' data(people)
#' model = pca(people, scale = TRUE, info = "Simple PCA model")
#' model = selectCompNum(model, 4)
#'
#' res = model$calres
#' summary(res)
#' plot(res)
#'
#' ## 3. Show scores plots for the results
#' par(mfrow = c(2, 2))
#' plotScores(res)
#' plotScores(res, cgroup = people[, "Beer"], show.labels = TRUE)
#' plotScores(res, comp = c(1, 3), show.labels = TRUE)
#' plotScores(res, comp = 2, type = "h", show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' ## 4. Show residuals and variance plots for the results
#' par(mfrow = c(2, 2))
#' plotVariance(res, type = "h")
#' plotCumVariance(res, show.labels = TRUE)
#' plotResiduals(res, show.labels = TRUE, cgroup = people[, "Sex"])
#' plotResiduals(res, ncomp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' @export
pcares <- function(...) {
   # Creates an object of pcares class. In fact the class is a wrapper for ldecomp and
   # uses its methods and attributes.

   obj <- ldecomp(...)
   class(obj) <- c("pcares", "ldecomp")

   return(obj)
}

#' Plot method for PCA results object
#'
#' @description
#' Show several plots to give an overview about the PCA results
#'
#' @param x
#' PCA results (object of class \code{pcares})
#' @param comp
#' which components to show the scores plot for (can be one value or vector with two values).
#' @param ncomp
#' how many components to use for showing the residual distance plot
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#'
#' @export
plot.pcares <- function(x, comp = c(1, 2), ncomp = x$ncomp.selected, show.labels = TRUE, ...) {
   op <- par(mfrow = c(2, 2))
   on.exit(par(op))
   plotScores(x, comp = comp, show.labels = show.labels, ...)
   plotResiduals(x, ncomp, show.labels = show.labels, ...)
   plotVariance(x, type = "h", show.labels = show.labels, ...)
   plotCumVariance(x, type = "h", show.labels = show.labels, ...)
}

#' Summary method for PCA results object
#'
#' @description
#' Shows some statistics (explained variance, eigenvalues) about the results.
#'
#' @param object
#' PCA results (object of class \code{pcares})
#' @param ...
#' other arguments
#'
#' @export
summary.pcares <- function(object, ...) {
   summary.ldecomp(object, "Summary for PCA results", ...)
   invisible(object)
}

#' Print method for PCA results object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' PCA results (object of class \code{pcares})
#' @param ...
#' other arguments
#'
#' @export
print.pcares <- function(x, ...) {
   print.ldecomp(x, "Results for PCA decomposition (class pcares)", ...)
   cat("\n")
   invisible(x)
}


#' Save PCA results to CSV file
#'
#' @description
#' Saves PCA results to CSV file with structure identical to
#' the one generated by web-application at https://mda.tools/pca
#'
#' @param res
#' PCA results (object of class \code{pcares})
#' @param fileName
#' name (or full path) to CSV file to be created.
#' @param name
#' short name of the result object (e.g. \code{"cal"}, \code{"test"}. etc.).
#' @param sep
#' values separator (either \code{","} or \code{";"}).
#' @param dataFile
#' optional, name of the data file used to create the results.
#' @param classes
#' vector with class names for every object in the dataset (if any)
#' @param model
#' optional, if PCA model is provided, loadings and vectors for centering/scaling will be added
#' @param ...
#' other optional parameters
#'
#' @export
writeCSV.pcares <- function(res, fileName, name, sep = ",", dataFile = "", classes = NULL, model = NULL, ...) {

   # for future use
   outliers <- res$exclrows
   exclvars <- res$exclcols
   rowLabels <- rownames(res$scores)

   # number of rows with/without outliers
   nAll <- nrow(res$scores)
   nOut <- length(outliers)
   n <- nAll - nOut

   # number of components
   nComp <- res$ncomp
   nStat <- if (is.null(model)) 4 else 5 # var, expvar, resvar, resexpvar [eigenvals]

   # number of empty columns to add in addition to row name and first value
   nFiller <- max(nComp, nStat)
   filler <- paste1(rep('', nFiller), sep)
   fillerStat <- if (nStat < nComp) paste1(rep('', nComp - nStat), sep) else NULL
   fillerComp <- if (nComp < nStat) paste1(rep('', nStat - nComp), sep) else NULL

   # vector with component names/labels and values separator
   compNames <- paste1("PC", seq_len(nComp), sep)
   compNames <- substr(compNames, 1, nchar(compNames) - 1)

   # function for adding component based outcomes for every row
   addChunk <- function (name, values, out) {

      if (!is.null(name)) {
         out <- c(out, paste1(name, ':', sep, '', sep, filler))
         out <- c(out, paste1('object', sep, 'group', sep, 'outlier', sep, compNames, fillerComp))
      }


      for (i in seq_len(nAll)) {
         groupname <- if (is.null(classes)) '' else classes[i]

         # outlier
         if (nOut > 0 && i %in% outliers) {
            str <- paste1(rep('-', nComp), sep)
            out <- c(out, paste1(rowLabels[i], sep, groupname, sep, 'yes', sep, substr(str, 1, nchar(str) - 1), fillerComp))
         } else {
            str <- paste1(values[i, ], sep)
            out <- c(out, paste1(rowLabels[i], sep, groupname, sep, '', sep, substr(str, 1, nchar(str) - 1), fillerComp))
         }
      }

      return(out)
   }

   # general information
   center <- if (is.logical(res$center) && res$center == FALSE) "false" else "true"
   scale <- if (is.logical(res$scale) && res$scale == FALSE) "false" else "true"

   out <- NULL
   out <- c(out, paste1('Number of components:', sep, nComp, sep, '', filler))
   out <- c(out, paste1('Center:', sep, center, sep, '', filler))
   out <- c(out, paste1('Scale:', sep, scale, sep, '', filler))
   out <- c(out, ' ')

   # variances
   out <- c(out, paste1('Statistics:', sep, '', sep, '', filler))

   if (!is.null(model)) {
      out <- c(out, paste1('', sep, '', sep, 'PC', sep, 'Expvar', sep, 'Cumexpvar', sep, 'Resvar', sep, 'Cumresvar', sep, 'Eigenvalues', fillerStat))

      for (a in seq_len(nComp)) {
         out <- c(out, paste1('', sep, '', sep, a, sep,
               res$expvar[a] / 100, sep,
               res$cumexpvar[a] / 100, sep,
               1 - res$expvar[a] / 100, sep,
               1 - res$cumexpvar[a] / 100, sep,
               model$eigenvals[a],
               fillerStat
            )
         )
      }
   } else {
      out <- c(out, paste1('', sep, '', sep, 'PC', sep, 'Expvar', sep, 'Cumexpvar', sep, 'Resvar', sep, 'Cumresvar', fillerStat))
      for (a in seq_len(nComp)) {
         out <- c(out, paste1('', sep, a, sep,
               res$expvar[a] / 100, sep,
               res$cumexpvar[a] / 100, sep,
               1 - res$expvar[a] / 100, sep,
               1 - res$cumexpvar[a] / 100,
               fillerStat
            )
         )
      }
   }
   out <- c(out, ' ')

   if (!is.null(model)) {
      varlabels <- rownames(model$loadings)
      nCols <- nrow(model$loadings)
      center <- if(is.vector(model$center) && length(model$center) == nCols) model$center else rep(0, nCols)
      scale <- if(is.vector(model$scale) && length(model$scale) == nCols) model$scale else rep(1, nCols)

      out <- c(out, ' ')
      out <- c(out, paste1('Loading, center and scale values:', sep, '', sep, '', filler))
      out <- c(out, paste1('variable', sep, 'center', sep, 'scale', sep, compNames, fillerComp))

      for (i in seq_len(nCols)) {
         loadstr <- paste1(model$loadings[i, ], sep)
         out <- c(out, paste1(
            varlabels[i], sep,
            center[i], sep,
            scale[i], sep,
            substr(loadstr, 1, nchar(loadstr) - 1),
            fillerComp)
         )
      }
      out <- c(out, ' ')
   }


   # scores and distances
   out <- c(out, ' ')
   out <- addChunk('Scores', res$scores, out)
   out <- c(out, ' ')
   out <- addChunk('Score distances, h', res$T2, out)
   out <- c(out, ' ')
   out <- addChunk('Orthogonal distances, q', res$Q, out)


   # if decimal separator is not ".", replace all "." with the correct separator
   if (sep == ';') {
      out <- gsub("\\.", ",", out)
   }

   # add header
   out <- c(paste1('Data filename:', sep, dataFile, sep, '', filler), out)
   out <- c(paste1(paste1('PCA results (',  name, ')'), sep, 'https://mda.tools/pca/', sep, '', filler), ' ', out)

   writeLines(out, fileName)
}