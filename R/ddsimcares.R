#' Results of DD-SIMCA one-class classification
#'
#'  @description
#' \code{ddsimcares} is used to store results for DD-SIMCA one-class classification. Do not create
#' this object manually, it will be created automatically by applying DD-SIMCA model.
#'
#' @param pcares
#' results of PCA decomposition of data (class \code{pcares}).
#' @param outcomes
#' outcomes of DD-SIMCA classification procedure.
#' @param indices
#' list with the object indices (members, strangers, unknown).
#' @param numbers
#' list with the object numbers in each subset (members, strangers, unknown).
#' @param alpha
#' significance level used for making the predictions.
#' @param classname
#' short text (up to 20 symbols) with class name.
#' @param c.ref
#' optional, vector with reference classes.
#'
#' @details
#' Class \code{ddsimcares} inherits all properties and methods of class \code{\link{pcares}}, and
#' has additional properties and functions for representing of classification results and other
#' DD-SIMCA outcomes.
#'
#' Do not create a \code{ddsimcares} object manually, it is created automatically when
#' a DD-SIMCA model is developed (see \code{\link{ddsimca}}) or when the model is applied to a new
#' data (see \code{\link{predict.ddsimca}}). The object can be used to show summary and plots for
#' the results.
#'
#' @return
#' Returns an object (list) of class \code{ddsimcares} with the same fields as \code{\link{pcares}}
#' plus additional field \code{simcares} which is a list with all DD-SIMCA outcomes and related
#' properties:
#'
#' @seealso
#' Methods specific for \code{ddsimcares} objects:
#' \tabular{ll}{
#'  \code{print.ddsimcares} \tab shows information about the object.\cr
#'  \code{summary.ddsimcares} \tab shows statistics for results of classification.\cr
#'  \code{as.data.frame.ddsimcares} \tab shows converts DD-SIMCA results into data frame.\cr
#'  \code{as.matrix.ddsimcares} \tab shows converts summary of DD-SIMCA results into matrix.\cr
#'  \code{writeCSV.ddsimcares} \tab saves DD-SIMCA results into a CSV file.\cr
#'  \code{plotAcceptance.ddsimcares} \tab shows acceptance plot (q/q0 vs h/h0) with decision and outlier boundaries.\cr
#'  \code{plotExtremes.ddsimcares} \tab shows extremes plot.\cr
#'  \code{plotAliens.ddsimcares} \tab shows aliens plot.\cr
#'  \code{plotDistance.ddsimcares} \tab shows plot with individual distances (q, h or f).\cr
#' }
#'
#' Methods, inherited from \code{\link{ldecomp}} class:
#' \tabular{ll}{
#'  \code{\link{plotScores.ldecomp}} \tab makes scores plot.\cr
#'  \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
#'  \code{\link{plotCumVariance.ldecomp}} \tab makes cumulative explained variance plot.\cr
#' }
#' Check also \code{\link{ddsimca}} and \code{\link{pcares}}.
#'
#' @examples
#' ## make a DD-SIMCA model for Iris setosa class and show results for calibration set
#' library(mdatools)
#'
#' data = iris[, 1:4]
#' class = iris[, 5]
#'
#' # take every second of first 50 objects (setosa) as calibration set
#' se = data[1:2:50, ]
#'
#' # take the rest as test set
#' ind.test = c(2:2:50, 51:150)
#' x.test = data[ind.test, ]
#' c.test = class[ind.test]
#'
#' # make DD-SIMCA model and set optimal number of components to 1
#' model = ddsimca(se, 'setosa', scale = TRUE)
#' model = selectCompNum(model, 1)
#'
#' # apply model to test set
#' r = predict(model, x.test, c.test)
#' print(r)
#'
#' # show summary
#' summary(r)
#'
#' # show plots
#' par(mdrow = c(2, 2))
#' plotAcceptance(r)
#' plotFoMs(r)
#' plotExtremes(r)
#' plotSelectivity(r)
#'
#' @export
ddsimcares <- function(pcares, outcomes, classname = classname, indices, numbers, alpha, c.ref = NULL) {
   res <- pcares
   res[["simca"]] <- list(
      alpha = alpha,
      outcomes = outcomes,
      indices = indices,
      numbers = numbers,
      classname = classname,
      c.ref = c.ref
   )
   class(res) <- c("ddsimcares", "pcares", "ldecomp")
   return(res)
}


#' Creates a data frame from DD-SIMCA classification results.
#'
#' @description
#' The data frame contains object label, class label (if provided), role, decision, as well
#' ass valued for distances (h/h0, q/q0, f) for every object from dataset used to
#' create the results object.
#'
#' @param x
#' classification results (object of class \code{ddsimcamres}).
#' @param ncomp
#' model complexity (number of components) to compute the classiciation results for.
#' @param limType
#' limit type to use ('classic' or 'robust')
#' @param ...
#' other arguments
#'
#' @returns a data frame.
#'
#' @export
as.data.frame.ddsimcares <- function(x, ncomp = x$ncomp.selected, limType = "classic", ...) {

  if (ncomp < 1 || ncomp > ncol(x$scores)) {
      stop("Wrong value for 'ncomp' parameter.")
   }

   limType <- processLimType(limType)
   res <- x$simca
   v <- res$outcomes[[limType]]$values[[ncomp]]


   out <- data.frame(
      'decision' = ifelse(v$decisions, "in", "out"),
      'h/h0' = v$h,
      'q/q0' = v$q,
      'f' = v$f,
      check.names = FALSE
   )
   if (!is.null(res$c.ref)) {
      out.extra <- data.frame(
         'class' = res$c.ref,
         'role' = v$roles
      )
      out <- cbind(out.extra, out)
   }

   rownames(out) <- rownames(x$scores)
   return(out)
}


#' Creates a matrix from DD-SIMCA classification results.
#'
#' @description
#' The matrix contains results characteristics computed for different number
#' of components, including number of objects accepted/rejected by the model,
#' explained and cumulative explained variance, as well as (if reference classes
#' are provided) main figures of merits (TP, FP, TN, FN, sensitivity, specificity,
#' efficiency and selectivity).
#'
#' If dataset uuse to create the results contained only target class members or only
#' objects from alternative classes, the FoMs will be reduced accordingly.
#'
#' @param x
#' classification results (object of class \code{ddsimcamres}).
#' @param limType
#' limit type to use ('classic' or 'robust')
#' @param ...
#' other arguments
#'
#' @returns a matrix.
#'
#' @export
as.matrix.ddsimcares <- function(x, limType = "classic", ...) {

   limType <- processLimType(limType)
   res <- x$simca$outcomes[[limType]]
   ncomp <- ncol(x$scores)

   labels <- c(
      "nin" = "In",
      "nout" = "Out",
      "sns" = "Sensitivity",
      "spc" = "Specificity",
      "sel" = "Selectivity",
      "eff" = "Efficiency",
      "acc" = "Accuracy",
      "TP" = "TP",
      "FN" = "FN",
      "TN" = "TN",
      "FP" = "FP",
      "Expvar" = "Expvar",
      "Cumexpvar" = "Cumexpvar"
   )

   names <- c("nin", "nout")

   if (x$simca$numbers$members > 0) {
      names <- c(names, c("TP", "FN", "sns"))
   }

   if (x$simca$numbers$strangers > 0) {
      names <- c(names, c("TN", "FP", "spc", "sel"))
   }

   if (x$simca$numbers$members > 0 && x$simca$numbers$strangers > 0) {
      names <- c(names, c("acc", "eff"))
   }

   out <- do.call(cbind, res[names])
   colnames(out) <- labels[colnames(out)]
   rownames(out) <- colnames(x$scores)
   return(cbind(as.matrix.ldecomp(x), out))
}



################################
#  Static and local methods    #
################################

#' Make correction to limit types names.
#'
#' @param limType
#' the limit type provided by user
#'
processLimType <- function(limType) {
   if (is.null(limType)) return("moments")
   if (limType == "classic") return("moments")
   if (!(limType == "moments" || limType == "robust") ) {
      stop("Wrong value for parameter 'limType'.")
   }
   return(limType)
}


#' Save DD-SIMCA results to CSV file
#'
#' @description
#' Saves DD-SIMCA results to CSV file with structure identical to
#' the one generated by web-application at https://mda.tools/ddsimca
#'
#' @param res
#' DD-SIMCA results (object of class \code{ddsimcares})
#' @param fileName
#' name (or full path) to CSV file to be created.
#' @param ncomp
#' model complexity (number of components) to compute the classiciation results for.
#' @param limType
#' limit type to use ('classic' or 'robust')
#' @param sep
#' values separator (either \code{","} or \code{";"}).
#' @param name
#' short name of the result object (e.g. \code{"cal"}, \code{"test"}. etc.).
#' @param dataFile
#' optional, name of the data file used to create the results.
#' @param src
#' optional, source of the results
#' @param ...
#' other optional parameters
#'
#' @export
writeCSV.ddsimcares <- function(res, fileName, name = "cal", sep = ",", dataFile = "",
   ncomp = res$ncomp.selected, limType = "classic", src = "mdatools for R", ...) {

   paste1 <- function(...) {
      paste0(..., collapse = "")
   }

   limType <- processLimType(limType)

   # reference classes
   rowlabels <- rownames(res$scores)
   classes <- res$simca$c.ref
   hasClasses <- !is.null(classes)

   # for future use
   outliers <- res$exclrows
   exclvars <- res$exclcols
   rowLabels <- rownames(res$scores)

   # number of rows with/without outliers
   nAll <- nrow(res$scores)
   nOut <- length(outliers)
   n <- nAll - nOut

   # number of components
   filler <- paste1(rep('', 3 + hasClasses), sep)


   # DD-SIMCA outcomes
   v <- res$simca$outcomes[[limType]]
   va <- v$values[[ncomp]]

   # general information
   center <- if (is.logical(res$center) && res$center == FALSE) "false" else "true"
   scale <- if (is.logical(res$scale) && res$scale == FALSE) "false" else "true"

   out <- NULL
   out <- c(out, paste1('Target class:', sep, res$simca$classname, sep, '', filler))
   out <- c(out, paste1('Reference classes:', sep, (if (hasClasses) "provided" else "not provided"), sep, '', filler))
   out <- c(out, paste1('Number of components:', sep, ncomp, sep, '', filler))
   out <- c(out, paste1('Alpha:', sep, res$simca$alpha, sep, '', filler))
   out <- c(out, paste1('Limits estimator:', sep, if (limType == "moments") "classic" else limType, sep, '', filler))
   out <- c(out, ' ')

   out <- c(out, paste1('Nh:', sep, va$Nh, sep, '', filler))
   out <- c(out, paste1('Nq:', sep, va$Nq, sep, '', filler))
   out <- c(out, paste1('Nf:', sep, va$Nf, sep, '', filler))
   out <- c(out, ' ')

   out <- c(out, paste1('F crit (extremes):', sep, va$fce, sep, '', filler))
   out <- c(out, paste1('F crit (outliers):', sep, va$fco, sep, '', filler))
   out <- c(out, ' ')


   if (length(exclvars) > 0) {
      out <- c(out, paste1('Number of excluded variables:', sep, length(exclvars), sep, '', filler))
      out <- c(out, paste1('Indices of excluded variables:', sep, arr2int(exclvars), sep, '', filler))
      out <- c(out, ' ')
   }

   out <- c(out, paste1('Objects:', sep, nAll, sep, '', filler))
   out <- c(out, paste1('In:', sep, sum(va$decisions & !res$simca$indices$excluded), sep, '', filler))
   out <- c(out, paste1('Out:', sep, sum(!va$decisions & !res$simca$indices$excluded), sep, '', filler))
   out <- c(out, paste1('Excluded:', sep, sum(res$simca$indices$excluded), sep, '', filler))

   # statistics for members if any
   if (res$simca$numbers$members > 0) {
      out <- c(out, ' ')
      out <- c(out, paste1('Members:', sep, res$simca$numbers$members, sep, '', filler))
      out <- c(out, paste1('Regular:', sep, sum(va$roles == "regular"), sep, '', filler))
      out <- c(out, paste1('Extreme:', sep, sum(va$roles == "extreme"), sep, '', filler))
      out <- c(out, paste1('Outliers:', sep, sum(va$roles == "outlier"), sep, '', filler))
      out <- c(out, paste1('Sensitivity:', sep, v$sns[ncomp], sep, '', filler))
   }

   # statistics for strangers if any
   if (res$simca$numbers$strangers > 0) {
      out <- c(out, ' ')
      out <- c(out, paste1('Strangers:', sep, res$simca$numbers$strangers, sep, '', filler))
      out <- c(out, paste1('In:', sep, va$FP, sep, '', filler))
      out <- c(out, paste1('Aliens:', sep, sum(va$roles == "alien"), sep, '', filler))
      out <- c(out, paste1('External:', sep, sum(va$roles == "external"), sep, '', filler))
      out <- c(out, paste1('Specificity:', sep, v$spc[ncomp], sep, '', filler))
      out <- c(out, paste1('Selectivity:', sep, v$sel[ncomp], sep, '', filler))

      out <- c(out, ' ')
      out <- c(out, paste1('Beta:', sep, va$beta, sep, '', filler))
      out <- c(out, paste1('Non-centrality (s):', sep, va$s, sep, '', filler))
      out <- c(out, paste1('Scaling (f0):', sep, va$f0, sep, '', filler))

      if (res$simca$numbers$members > 0) {

         PPV <- v$TP[ncomp] / (v$TP[ncomp] + v$FP[ncomp])
         TPR <- v$sns[ncomp];
         FPR <- 1 - v$spc[ncomp];
         F1 <- 2 * PPV * TPR / (PPV + TPR)
         PLR <- TPR / FPR;
         NLR <- (1 - TPR) / (1 - FPR);

         out <- c(out, ' ')
         out <- c(out, paste1('Efficiency:', sep, v$eff[ncomp], sep, '', filler))
         out <- c(out, paste1('Accuracy:', sep, v$acc[ncomp], sep, '', filler))
         out <- c(out, paste1('Precision:', sep, PPV, sep, '', filler))
         out <- c(out, paste1('F1-score:', sep, F1, sep, '', filler))
         out <- c(out, paste1('PLR:', sep, PLR, sep, '', filler))
         out <- c(out, paste1('NLR:', sep, NLR, sep, '', filler))
      }
   }

   # outcomes for every object
   out <- c(out, ' ')
   out <- c(out, paste1('', sep, 'object', sep, 'class', sep, 'role', sep, 'decision', sep, 'h/h0', sep, 'q/q0', sep, 'f'))

   odf <- as.data.frame.ddsimcares(res, limType = limType)
   decision <- ifelse(va$decisions, "in", "out")

   for (i in seq_len(nAll)) {

      groupname <- if (is.null(classes)) '' else classes[i]

      # outlier
      if (nOut > 0 && i %in% outliers) {
         out <- c(out, paste1('', sep, rowlabels[i], sep, '-', sep, 'removed', sep, '-', sep, '-', sep, '-', sep, '-'))
      } else {
         out <- c(out, paste1('', sep, rowlabels[i], sep, classes[i], sep, va$roles[i], sep, decision[i], sep, va$h[i], sep, va$q[i], sep, va$f[i]))
      }
   }


   # if decimal separator is not ".", replace all "." with the corect separator
   if (sep == ';') {
      out = gsub("\\.", ",", out);
   }

   # add header
   out <- c(paste1('Data filename:', sep, dataFile, sep, '', filler), out)
   out <- c (paste1(paste1('DD-SIMCA results (',  name, ')'), sep, src, sep, '', filler), ' ', out);

   writeLines(out, fileName)
}


#################################
# Plotting methods              #
#################################


#' Selectivity area plot (similar to ROC curve).
#'
#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares}).
#' @param ncomp
#' number of components in the model to show the plot for.
#' @param limType
#' limit type to show the plot for ('classic' or 'robust').
#' @param col
#' color of the plot elements.
#'
#' @description
#' The method generates a sequence of values for significance level (alpha) and then
#' computes corresponding selectivity (1 - beta) for each value from the sequence using model parameters
#' for given number of components and limit type. Then the values are shown as line plot using
#' coordinates (1 - alpha) vs beta and area under the curve is shown in the plot title.
#'
#' Current values of alpha and beta are shown on the plot as a marker.
#'
#' @export
plotSelectivityArea.ddsimcares <- function(obj,
   ncomp = obj$ncomp.selected,
   limType = "classic",
   col = "#2679b2",
   xlab = bquote("1 - selectivity, "~beta),
   ylab = bquote("Sensitivity, 1 - "~alpha),
   ...
   ) {

   nstrangers <- obj$simca$numbers$strangers
   if (nstrangers < 1) {
      stop("This result object does not contain information about non-target class members.")
   }

   # get the limits
   limType <- processLimType(limType)
   v <- obj$simca$outcomes[[limType]]$values[[ncomp]]

   da <- 0.0001
   alpha <- seq(0, 1, by = da)
   norm1 <- 1 / (v$k + v$s)
   norm2 <- 1 / (v$Sz)
   fcrit <- qchisq(1 - alpha, v$k)
   z <- fcrit / v$f0

   beta <- pnorm(((z * norm1)^v$hz - v$Mz)*norm2)
   area <- sum(beta) * da
   pd <- matrix(c(0, 1, 0, 1), ncol = 2, nrow = 2)
   attr(pd, "name") <- sprintf("Selecitivity (A = %d, AUC = %.4f)", ncomp, 1 - area)

   mdaplot(pd, type = "p", cex = 0, xlab = xlab, ylab = ylab, ...)
   lines(beta, 1 - alpha, col = col)
   points(v$beta, 1 - obj$simca$alpha, pch = 16, col = col)
}

#' Sensitivity plot.
#'
#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares}).
#' @param show.ci
#' logical, show or not the expected sensitivity and 95% confidence interval for it.
#' @param ci.col
#' color for lines showing expected sensitivity and confidence interval.
#' @param ci.lty
#' type of lines (lower, expected, upper) for the expected sensitivity and confidence interval.
#' @param ...
#' any parameters suitable for the \code{\link{plotFoM.ddsimcares}} method.
#'
#' The method is a wrapper for more general \code{\link{plotFoM.ddsimcares}} method with a
#' possibility to show confidence interval for the expected sensitivity.
#'
#' @export
plotSensitivity.ddsimcares <- function(obj, show.ci = FALSE, ci.col = "#a0a0a0", ci.lty = c(3, 2, 3), ...) {

   nmembers <- obj$simca$numbers$members
   if (nmembers < 1) {
      stop("This result object does not contain information about class members.")
   }

   p <- plotFoM(obj, fom = "sens", ...)

   if (show.ci) {
      ci.lo <- qbinom(0.025, nmembers, 1 - obj$simca$alpha) / nmembers
      ci.up <- qbinom(0.975, nmembers, 1 - obj$simca$alpha) / nmembers
      abline(h = c(ci.lo, 1 - obj$simca$alpha, ci.up), col = ci.col, lty = ci.lty)
   }

   return(invisible(p))
}


#' Figure of merit plot.
#'
#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares})
#' @param fom
#' name of the figure of merit to show the plot for ('sens', 'spec', 'sel', 'eff', 'acc').
#' @param limType
#' limit type to show the plot for ('classic' or 'robust').
#' @param type
#' type of plot (can be 'b', 'l', or 'h').
#' @param ylim
#' limits for y-axis.
#' @param pch
#' marker number or symbol for the plot points (if 'b' type is selected).
#' @param show.plot
#' logical, if FALSE returns only values and does not show the plot.
#' @param ...
#' any parameters for \code{\link{mdaplot}} method can be provided.
#'
#' Shows the selected figure of merit values vs. number of components. Pay attention, that
#' some of the FoMs may not be available, for example if result object is made only for
#' members of the target class, only sensitivity ('sens') is available.
#'
#' @export
plotFoM.ddsimcares <- function(obj, fom = "sens", limType = "classic", type = "b",
   ylim = c(0, 1.1), pch = 1, show.plot = TRUE, ...) {

   labels <- c(
      "sens" = "Sensitivity",
      "spec" = "Specificity",
      "sel" = "Selectivity",
      "eff" = "Efficiency",
      "acc" = "Accuracy"
   )

   stat <- as.matrix.ddsimcares(obj, limType = limType)

   if (!(fom %in% names(labels)) ||  !(labels[fom] %in% colnames(stat))) {
      stop("Wrong value for 'fom' argument.")
   }

   pd <- matrix(stat[, labels[fom]], nrow = 1)
   attr(pd, "xaxis.values") <- seq_len(obj$ncomp)
   attr(pd, "xaxis.name") <- "Number of components, A"
   attr(pd, "name") <- labels[fom]
   rownames(pd) <- labels[fom]

   if (!show.plot) {
      return(pd)
   }

   return(invisible(mdaplot(pd, type = type, ylim = ylim, pch = pch, ...)))
}


#' Figures of merit plot (multiplt FoMs).
#'
#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares})
#' @param foms
#' vector with FoM names to show on the plot (any combination of 'sens', 'spec', 'sel', 'eff', 'acc').
#' @param limType
#' limit type to show the plot for ('classic' or 'robust').
#' @param type
#' type of plot (can be 'b', 'l', or 'h').
#' @param pch
#' marker number or symbol for the plot points (if 'b' type is selected), can be single or individual for each FoM.
#' @param ylim
#' limits for y-axis.
#' @param legend.position
#' position of legend on the plot.
#' @param main
#' main title.
#' @param show.plot
#' logical, if FALSE returns only values and does not show the plot.
#' @param ...
#' any parameters for \code{\link{mdaplotg}} method can be provided.
#'
#' Shows several figures of merit values vs. number of components as a group plot (see
#' details in \code{\link{mdaplotg}} method desciption.
#'
#' @export
plotFoMs.ddsimcares <- function(obj, foms = NULL, limType = "classic", type = "b", pch = 1,
   ylim = c(0, 1.1), legend.position = "top", main = "Figures of merit", ...) {

   nm <- obj$simca$numbers$members
   ns <- obj$simca$numbers$strangers

   if (is.null(foms)) {
      if (nm > 0 && ns > 0) {
         foms <- c("sens", "spec", "eff")
      } else if (nm > 0) {
         foms <- c("sens")
      } else {
         foms <- c("spec", "sel")
      }
   }

   labels <- c(
      "sens" = "Sensitivity",
      "spec" = "Specificity",
      "sel" = "Selectivity",
      "eff" = "Efficiency",
      "acc" = "Accuracy"
   )

   stat <- as.matrix.ddsimcares(obj, limType = limType)


   pds <- list()
   for (fom in foms) {
      if (!(fom %in% names(labels)) ||  !(labels[fom] %in% colnames(stat))) {
         stop("Wrong value for 'fom' argument.")
      }
      pd <- matrix(stat[, labels[fom]], nrow = 1)
      attr(pd, "xaxis.values") <- seq_len(obj$ncomp)
      attr(pd, "xaxis.name") <- "Number of components, A"
      attr(pd, "name") <- labels[fom]
      pds[[fom]] <- pd
   }

   return(invisible(mdaplotg(pds, type = type, ylim = ylim, pch = pch,
      legend.position = legend.position, main = main, ...)))
}

#' Extremes plot (shortcut to \code{\link{plotExtremes.ddsimcares}}).
#'
#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares})
#' @param ...
#' any parameter acceptable by \code{\link{plotExtremes.ddsimcares}}.
#'
#' @export
plotExtreme.ddsimcares <- function(obj, ...) {
   return(plotExtremes.ddsimcares(obj, ...))
}

#' Extremes plot.
#'
#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares})
#' @param ncomp
#' number of components in the model to show the plot for.
#' @param limType
#' limit type to show the plot for ('classic' or 'robust').
#' @param main
#' main title.
#' @param xlab
#' label for x-axis.
#' @param ylab
#' label for y-axis.
#' @param pch
#' marker number or symbol for the plot points.
#' @param col
#' color of the plot points.
#' @param ellipse.col
#' color of the confidence ellipse elements.
#' @param show.plot
#' logical, if FALSE returns only values and does not show the plot.
#'
#' @description
#' The method generates a sequence of values for significance level (alpha) and then
#' computes how many extremes the model produces by applying it to the target class members. After
#' that the computed number is plotted against the expected number.
#'
#' The plot also shows a confidence ellipse for the expected values.
#'
#' @export
plotExtremes.ddsimcares <- function(obj, ncomp = obj$ncomp.selected, limType = "classic",
   main = sprintf("Extremes (A = %d)", ncomp),
   xlab = "Number of extremes (expected)",
   ylab = "Number of extremes (observed)",
   pch = 1, col = "#2679b2",
   ellipse.col = "#eeeeee",
   show.plot = TRUE, ...) {


   res <- obj$simca
   if (is.null(res$indices) || is.null(res$indices$members)) {
      stop("No information about target class objects is found.")
   }

   ind <- if (!is.null(res$excluded)) res$indices$members && (!res$excluded) else res$indices$members
   if (sum(ind) < 3) {
      stop("Not enough target class members to make the plot (at least 3 is required).")
   }

   # get the limits
   limType <- processLimType(limType)
   v <- res$outcomes[[limType]]$values[[ncomp]]

   # get the DD-SIMCA outcomes for particular limit type and number of components
   f <- v$f[ind]
   p <- 1 - pchisq(f, v$Nf)

   nobj <- length(f)
   expected <- seq_len(nobj)
   alpha <- expected / nobj

   observed <- vapply(alpha, function(x) sum(p < x), 0)
   pd <- cbind(expected, observed)

   colnames(pd) <- c("Number of extremes (expected)", "Number of extremes (observed)")

   if (!show.plot) {
      return(pd)
   }

   # show axes, grid and diagonal

   if (is.null(xlab))
   if (is.null(ylab)) ylab <- "Number of extremes (observed)"

   plot(0, type = "n", xlim = c(0, nobj), ylim = c(0, nobj), main = main, xlab =  xlab, ylab = ylab)
   grid()
   lines(c(0, expected), c(0, expected), type = "l", col = ellipse.col)

   # compute and show the tolerance ellipse
   Nm <- qbinom(0.025, nobj, alpha)
   Np <- qbinom(0.975, nobj, alpha)
   segments(expected, Nm, expected, Np, col = ellipse.col)
   lines(c(0, expected), c(0, Nm), col = ellipse.col)
   lines(c(0, expected), c(0, Np), col = ellipse.col)

   # show the points
   points(expected, observed, pch = pch, col = col, ...)
}



#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares})
#' @param ncomp
#' number of components in the model to show the plot for.
#' @param limType
#' limit type to show the plot for ('classic' or 'robust')
#' @export
plotAliens.ddsimcares <- function(obj, ncomp = obj$ncomp.selected, limType = "classic",
   main = sprintf("Aliens (A = %d)", ncomp),
   xlab = "Number of aliens (expected)",
   ylab = "Number of aliens (observed)",
   pch = 1, col = "#2679b2",
   ellipse.col = "#eeeeee",
   show.plot = TRUE, ...) {


   res <- obj$simca
   if (is.null(res$indices) || is.null(res$indices$strangers)) {
      stop("No information about non-target class objects is found.")
   }

   ind <- if (!is.null(res$excluded)) res$indices$strangers && (!res$excluded) else res$indices$strangers
   if (sum(ind) < 3) {
      stop("Not enough non-target class members to make the plot (at least 3 is required).")
   }

   # get the limits
   limType <- processLimType(limType)
   v <- res$outcomes[[limType]]$values[[ncomp]]

   # get the DD-SIMCA outcomes for particular limit type and number of components
   f <- v$f[ind]

   nobj <- length(f)
   expected <- seq_len(nobj)
   beta <- expected / nobj

   zb <- qnorm(beta);
   eCrit <- v$m * (v$Sz * zb + v$Mz)^(1.0 / v$hz)
   observed <- vapply(eCrit, function(x) sum(f < x), 0)
   pd <- cbind(expected, observed)

   colnames(pd) <- c("Number of aliens (expected)", "Number of aliens (observed)")

   if (!show.plot) {
      return(pd)
   }

   # show axes, grid and diagonal

   if (is.null(xlab))
   if (is.null(ylab)) ylab <- "Number of extremes (observed)"

   plot(0, type = "n", xlim = c(0, nobj), ylim = c(0, nobj), main = main, xlab =  xlab, ylab = ylab)
   grid()
   lines(c(0, expected), c(0, expected), type = "l", col = ellipse.col)

   # compute and show the tolerance ellipse
   i <- seq_len(nobj)
   D <- 2 * sqrt(i * (1 - beta))
   Nm <- qbinom(0.025, nobj, beta)
   Np <- qbinom(0.975, nobj, beta)
   segments(expected, Nm, expected, Np, col = ellipse.col)
   lines(c(0, expected), c(0, Nm), col = ellipse.col)
   lines(c(0, expected), c(0, Np), col = ellipse.col)

   # show the points
   points(expected, observed, pch = pch, col = col, ...)
}


#' Acceptance plot for DD-SIMCA results object.
#'
#' @description
#' Shows a plot with h/h0 and q/q0 distance values for given number of components and type of critical
#' limits.
#'
#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares})
#' @param ncomp
#' number of components in the model to show the plot for.
#' @param limType
#' limit type to show the plot for ('classic' or 'robust')
#' @param show
#' name of subset to show (\code{"members"}, \code{"strangers"}, \code{"all"}, or \code{"auto"}).
#' @param log
#' logical, apply log transformation or not.
#' @param show.legend
#' logica, show or not colorbar legend on the plot.
#' @param pch
#' marker type
#' @param colors
#' named vector with colors for every role
#' @param lim.lty
#' vector with two values - types of the lines showing critical limits (extreme, outlier).
#' @param lim.col
#' vector with two values - colors of the lines showing critical limits (extreme, outlier).
#' @param lim.lwd
#' vector with two values - thickness of the lines showing critical limits (extreme, outlier).
#' @param ylim
#' limits for y-axis, if NULL will be estimated automatically.
#' @param xlim
#' limits for x-axis, if NULL will be estimated automatically.
#' @param show.excluded
#' logical, show or not the excluded rows (objects).
#' @param ...
#' other parameters compatible with \code{\link{mdaplot}} method.
#'
#' @export
plotAcceptance.ddsimcares <- function(obj, ncomp = obj$ncomp.selected, limType = "classic",
      show = "auto",
      log = FALSE,
      show.legend = TRUE,
      pch = 1,
      colors = c("regular" = "#2679B2", "extreme" = "#F2B825", "outlier" = "#D22C2F", "alien" = "#2679B2", "external" = "#D22C2F", "unknown" = "#3c6784", "excluded" = "#00ff00"),
      lim.lty = c(2, 3),
      lim.col = c("darkgray", "darkgray"),
      lim.lwd = c(1, 1),  lab.cex = 0.75, lab.col = "#808080",
      ylim = NULL, xlim = NULL, res.name = NULL,
      show.excluded = FALSE, ...) {

   if (ncomp < 1 || ncomp > ncol(obj$scores)) {
      stop("Wrong value for 'ncomp' parameter.")
   }

   limType <- processLimType(limType)

   # get the DD-SIMCA outcomes for particular limit type and number of components
   res <- obj$simca
   v <- res$outcomes[[limType]]$values[[ncomp]]

   # save attributes
   attrs <- mda.getattr(obj$scores)
   nobj <- nrow(obj$scores)

   # check what to show
   if (show == "auto") {
      if (res$numbers$members > 0) {
         show <- "members"
      } else if (res$numbers$strangers > 0) {
         show <- "strangers"
      } else {
         show <- "all"
      }
   }

   # get object indices, and color grouping settings depending on the "show" value
   if (show == "members") {
      stopifnot("Result object does not contain any class members." = res$numbers$members > 0)
      ind <- res$indices$members
      levels <- c("regular", "extreme", "outlier")
      groups <- v$roles
      colmap <- colors
   } else if (show == "strangers") {
      stopifnot("Result object does not contain any members of non-target classes." = res$numbers$strangers > 0)
      ind <- res$indices$strangers
      levels <- c("alien", "external")
      groups <- v$roles
      colmap <- colors
   } else if (show == "all") {
      ind <- rep(TRUE, nobj)
      groups <- if (!is.null(res$c.ref)) res$c.ref else rep("unknown", nobj)
      levels <- unique(groups)
      colmap <- mdaplot.getColors(length(levels))
      names(colmap) <- levels
      colmap["excluded"] <- "#909090"
   } else {
      stop("Wrong value for parameter 'show' (can be: 'strangers', 'members', 'all' and 'auto').")
   }

   # distance values (already scalued) and roles
   x <- v$h
   y <- v$q

   # compute coordinates of the lines for critical limits
   slope <- - v$Nh / v$Nq
   xle <- seq(0, v$fce / v$Nh, length.out = 200)
   xlo <- seq(0, v$fco / v$Nh, length.out = 200)
   yle <- v$fce / v$Nq + slope * xle
   ylo <- v$fco / v$Nq + slope * xlo

   # apply log transformation if necessary
   if (log == TRUE) {
      x <- log(1 + x)
      y <- log(1 + y)
      xle <- log(1 + xle)
      xlo <- log(1 + xlo)
      yle <- log(1 + yle)
      ylo <- log(1 + ylo)
   }


   # show main plot
   pd <- cbind(x, y)
   rownames(pd) <- rownames(obj$scores)

   # remove excluded if not necessary
   if (!show.excluded && any(v$roles == "excluded")) {
      indExcluded <- v$roles == "excluded"
      pd <- subset(pd, !indExcluded)
      ind <- ind[!indExcluded]
      groups <- groups[!indExcluded]
   } else if (any(v$roles == "excluded")) {
      groups[v$roles == "excluded"] <- "excluded"
   }

   if (sum(ind) != length(ind)) {
      pd <- subset(pd, ind)
      groups <- groups[ind]
   }

   if (any(groups == "excluded")) {
      levels <- c(levels, "excluded")
   }

   colnames(pd) <- if (log) c("Score distance, ln(1 + h/h0)", "Orthogonal distance, ln(1 + q/q0)")
      else c("Score distance, h/h0", "Orthogonal distance, q/q0")

   res.str <- if(is.null(res.name)) "" else paste0(res.name, ", ")
   attr(pd, "name") <- sprintf("Acceptance (%sA = %d)", res.str, ncomp)

   colmap <- colmap[levels]

   # re-assign values for color group categories, so they include the number of
   # objects in each group
   if (!is.null(groups)) {
      sums <- sapply(levels, function(x) sum(groups == x))
      cgroup <- factor(paste0(groups, ": ", sums[groups], ""), levels = paste0(names(sums), ": ", sums, ""))
   } else {
      cgroup <- NULL
   }

   # amend the limits
   if (is.null(xlim)) xlim <- c(0, max(c(pd[, 1], xlo)))
   if (is.null(ylim)) ylim <- c(0, max(c(pd[, 2], ylo)) * 1.2)

   if (length(colmap) == 1) {
      p <- mdaplot(pd, ylim = ylim, xlim = xlim, pch = pch, ...)
   } else {
      p <- mdaplot(pd, ylim = ylim, xlim = xlim, colmap = colmap, pch = pch, cgroup = cgroup, lab.cex = lab.cex, lab.col = lab.col,...)
   }

   # show the critical limits
   if (!is.null(lim.lty) && !is.na(lim.lty[1])) {
      lines(xle, yle, col = lim.col[1], lty = lim.lty[1], lwd = lim.lwd[1])
   }

   if (!is.null(lim.lty) && !is.na(lim.lty[2])) {
      lines(xlo, ylo, col = lim.col[2], lty = lim.lty[2], lwd = lim.lwd[2])
   }

   return(p)
}


#' Show with distance values (score, orthogonal or full) vs object indices for DD-SIMCA results.
#'
#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares})
#' @param ncomp
#' number of components to show the plot for.
#' @param limType
#' limit type to show the plot for ('classic' or 'robust')
#' @param distance
#' which distance to show (\code{"h"} for score, \code{"q"} for orthogonal or \code{"f"} for full distance).
#' @param log
#' logical, apply log transformation or not.
#' @param show.legend
#' logica, show or not colorbar legend on the plot.
#' @param pch
#' marker type
#' @param lim.lty
#' vector with two values - types of the lines showing critical limits (extreme, outlier).
#' @param lim.col
#' vector with two values - colors of the lines showing critical limits (extreme, outlier).
#' @param lim.lwd
#' vector with two values - thickness of the lines showing critical limits (extreme, outlier).
#' @param ylim
#' limits for y-axis, if NULL will be estimated automatically.
#' @param xlim
#' limits for x-axis, if NULL will be estimated automatically.
#' @param ...
#' other parameters compatible with \code{\link{mdaplot}} method.
#'
#' @export
plotDistances.ddsimcares <- function(obj, ncomp = obj$ncomp.selected,
      limType = "classic", distance = "h",
      log = FALSE,
      show.legend = TRUE,
      pch = 1,
      lim.lty = c(2, 3),
      lim.col = c("darkgray", "darkgray"),
      lim.lwd = c(1, 1),
      ylim = NULL, xlim = NULL,
      show.plot = TRUE,
      show.excluded = FALSE,
      ...)
   {

   if (ncomp < 1 || ncomp > ncol(obj$scores)) {
      stop("Wrong value for 'ncomp' parameter.")
   }

   if (!(distance %in% c("h", "q", "f"))) {
      stop("Wrong value for parameter 'distance' (can be 'h', 'q', or 'f').")
   }

   limType <- processLimType(limType)
   res <- obj$simca
   v <- res$outcomes[[limType]]$values[[ncomp]]

   groups <- res$c.ref
   levels <- unique(groups)
   colmap <- "default"

   if (!is.null(groups) && length(levels) > 1) {
      sums <- sapply(levels, function(x) sum(groups == x))
      cgroup <- factor(paste0(groups, " (", sums[groups], ")"), levels = paste0(names(sums), " (", sums, ")"))
   } else {
      cgroup <- NULL
   }


   y <- v[[distance]]
   x <- attr(obj$scores, "yaxis.values")
   yle <- v$fce / v$Nf
   ylo <- v$fco / v$Nf

   if (distance != "f") {
      yle <- 0
      ylo <- 0
   } else {
      y <- y / v$Nf
   }

   if (log == TRUE) {
      y <- log(1 + y)
      yle <- log(1 + yle)
      ylo <- log(1 + ylo)
   }

   exclrows <- attr(obj$scores, "exclrows")
   objnames <- rownames(obj$scores)
   if (!show.excluded && length(exclrows) > 0) {
      y <- y[-exclrows]
      objnames <- objnames[-exclrows]

      if (!is.null(cgroup)) cgroup <- cgroup[-exclrows]
      if (!is.null(x)) x <- x[-exclrows]
      exclrows <- NULL
   }

   if (is.null(x)) {
      x <- seq_along(y)
   }

   if (is.null(ylim)) {
      ylim <- c(0, max(c(y, ylo)) * 1.2)
   }

   if (is.null(xlim)) {
      xlim <- c(0, max(x))
   }

   pd <- cbind(x, y)
   attr(pd, "name") <- sprintf("Distances (A = %d)", ncomp)

   n1 <- c("h" = "Score", "q" = "Orthogonal", "f" = "Full")
   n2 <- c("h" = "h/h0", "q" = "q/q0", "f" = "f/f0")


   colnames(pd) <- c(
      if (is.null(attr(obj$scores, "yaxis.name"))) "Objects" else attr(obj$scores, "yaxis.name"),
      if (log) sprintf("%s distance, log(1 + %s)", n1[distance], n2[distance]) else sprintf("%s distance, %s", n1[distance], n2[distance])
   )
   rownames(pd) <- objnames

   attr(pd, "exclrows") <- exclrows
   if (!show.plot) {
      attr(pd, "ylim") <- ylim
      return(pd)
   }

   mdaplot(pd, ylim = ylim, xlim = xlim, colmap = colmap, pch = pch,
      cgroup = cgroup, show.excluded = show.excluded, ...)


   if (yle > 0 && !is.null(lim.lty) && !is.na(lim.lty[1])) {
      abline(h = yle, col = lim.col[1], lty = lim.lty[1], lwd = lim.lwd[1])
   }

   if (ylo > 0 && !is.null(lim.lty) && !is.na(lim.lty[2])) {
      abline(h = ylo, col = lim.col[2], lty = lim.lty[2], lwd = lim.lwd[2])
   }
}


#' Summary method for DD-SIMCA results.
#'
#' @description
#' Shows performance statistics for the DD-SIMCA results.
#'
#' @param object
#' DD-SIMCA results (object of class \code{ddsimcares})
#' @param ...
#' other arguments
#'
#' @export
summary.ddsimcares <- function(object, ncomp = object$ncomp.selected, limType = "classic", ...) {

   fprintf("\nDD-SIMCA results\n\n")

   if (!is.null(object$info) && nchar(object$info)) {
      fprintf("Info: %s\n", object$info)
   }

   fprintf("Target class: %s\n", object$simca$classname)
   fprintf("Number of components: %d\n", object$ncomp)
   fprintf("Number of selected components: %d\n", ncomp)
   fprintf("Type of limits: %s\n", limType)
   fprintf("Alpha: %.3f\n", object$simca$alpha)

   cat("\n")
   fprintf("Number of objects: %d\n", nrow(object$scores))
   fprintf(" - target class members: %d\n", object$simca$numbers$members)
   fprintf(" - from non-target classes: %d\n", object$simca$numbers$strangers)
   fprintf(" - from unknown classses: %d\n", object$simca$numbers$unknown)


   cat("\n")

   stat <- as.matrix.ddsimcares(object, limType = limType)
   rownames(stat)[ncomp] <- paste0(rownames(stat)[ncomp], " (s)")
   print(stat, 4)
   cat("\n")}


#' Plot method for DD-SIMCA results.
#'
#' @param obj
#' DD-SIMCA results (object of class \code{ddsimcares})
#'
#' @export
plot.ddsimcares <- function(obj, ...) {

   if (!is.null(obj$simca$c.ref)) {
      op <- par(mfrow = c(1, 2))
      on.exit(par(op))
      plotAcceptance.ddsimcares(obj)
      plotFoMs.ddsimcares(obj)
   } else {
      plotAcceptance.ddsimcares(obj)
   }
}

#' Print method for DD-SIMCA results
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' SIMCA results (object of class \code{simcares})
#' @param ...
#' other arguments
#'
#' @export
print.ddsimcares <- function(x, ...) {

   cat("Result for DD-SIMCA one-class classification (class simcares)\n\n")

   cat("\nFields inherited from 'pcares'/'ldecomp' classes:\n")
   cat(" $scores - matrix with score values\n")
   cat(" $T2 - matrix with T2 distances\n")
   cat(" $Q - matrix with Q residuals\n")
   cat(" $ncomp.selected - selected number of components\n")
   cat(" $expvar - explained variance for each component\n")
   cat(" $cumexpvar - cumulative explained variance\n")

   cat("\nFields from 'ddsimcares' class:\n")
   cat(" $simca$alpha - significance level for critical limits\n")
   cat(" $simca$outcomes - list with all classification outcomes\n")
   cat(" $simca$indices - list with object indices (members, strangers, unknown).\n")
   cat(" $simca$numbers - list with object numbers (members, strangers, unknown).\n")
   cat(" $simca$classname - name of target class.\n")
   cat(" $simca$c.ref - vector with reference class labels (if provided).\n")

   cat("\n")
}
