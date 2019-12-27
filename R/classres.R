#' Results of classification
#'
#' @description
#' \code{classres} is used to store results classification for one or multiple classes.
#'
#' @param c.pred
#' matrix with predicted values (+1 or -1) for each class.
#' @param c.ref
#' matrix with reference values for each class.
#' @param p.pred
#' matrix with probability values for each class.
#' @param ncomp.selected
#' vector with selected number of components for each class.
#'
#' @details
#' There is no need to create a \code{classres} object manually, it is created automatically when
#' build a classification model (e.g. using \code{\link{simca}} or \code{\link{plsda}}) or apply
#' the model to new data. For any classification method from \code{mdatools}, a class using to
#' represent results of classification (e.g. \code{\link{simcares}}) inherits fields and methods of
#' \code{classres}.
#'
#' @return
#' \item{c.pred}{predicted class values (+1 or -1).}
#' \item{p.pred}{predicted class probabilities.}
#' \item{c.ref}{reference (true) class values if provided.}
#'
#' The following fields are available only if reference values were provided.
#' \item{tp}{number of true positives.}
#' \item{tn}{number of true negatives.}
#' \item{fp}{nmber of false positives.}
#' \item{fn}{number of false negatives.}
#' \item{specificity}{specificity of predictions.}
#' \item{sensitivity}{sensitivity of predictions.}
#' \item{misclassified}{ratio of misclassified objects.}
#'
#' @seealso
#' Methods \code{classres} class:
#' \tabular{ll}{
#'  \code{\link{showPredictions.classres}} \tab shows table with predicted values.\cr
#'  \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
#'  \code{\link{plotSensitivity.classres}} \tab makes sn plot.\cr
#'  \code{\link{plotSpecificity.classres}} \tab makes specificity plot.\cr
#'  \code{\link{plotMisclassified.classres}} \tab makes ms ratio plot.\cr
#'  \code{\link{plotPerformance.classres}} \tab makes plot with misclassified ratio, specificity
#'  and sensitivity values.\cr
#' }
#'
#' @export
classres <- function(c.pred, c.ref = NULL, p.pred = NULL, ncomp.selected = 1) {

   if (length(dim(c.pred)) != 3) {
      stop("Wrong number of dimensions for 'c.pred' array (should be 3-way array).")
   }

   obj <- list()
   obj$c.ref <- if(!is.null(c.ref)) as.factor(c.ref)
   obj$c.pred <- c.pred
   obj$p.pred <- p.pred
   obj$nclasses <- dim(c.pred)[3]
   obj$ncomp <- dim(c.pred)[2]
   obj$classnames <- dimnames(c.pred)[[3]]

   # check that ncomp.selected is correct
   if (is.null(ncomp.selected)) ncomp.selected <- obj$ncomp
   if (ncomp.selected < 1 || ncomp.selected > obj$ncomp) {
      stop("Wrong value for 'ncomp.selected' parameer.")
   }

   obj$ncomp.selected <- ncomp.selected

   if (!is.null(c.ref)) {
      obj <- c(obj, classres.getClassificationPerformance(c.ref, c.pred))
   }


   obj$call <- match.call()
   class(obj) <- "classres"

   return(obj)
}

#' Calculation of  classification performance parameters
#'
#' @description
#' Calculates and returns performance parameters for classification result (e.g. number of false
#' negatives, false positives, sn, specificity, etc.).
#'
#' @param c.ref
#' reference class values for objects (vector with numeric or text values)
#' @param c.pred
#' predicted class values for objects (array nobj x ncomponents x nclasses)
#'
#' @return
#' Returns a list with following fields:
#' \tabular{ll}{
#'    \code{$fn} \tab number of false negatives (nclasses x ncomponents) \cr
#'    \code{$fp} \tab number of false positives (nclasses x ncomponents) \cr
#'    \code{$tp} \tab number of true positives (nclasses x ncomponents) \cr
#'    \code{$sensitivity} \tab sn values (nclasses x ncomponents) \cr
#'    \code{$specificity} \tab specificity values (nclasses x ncomponents) \cr
#'    \code{$specificity} \tab ms ratio values (nclasses x ncomponents) \cr
#' }
#'
#' @details
#' The function is called automatically when a classification result with reference values is
#' created, for example when applying a \code{plsda} or \code{simca} models.
#'
classres.getClassificationPerformance <- function(c.ref, c.pred) {

   if (is.null(c.ref) || is.null(c.pred)) {
      stop("Both reference and predicted class values are required.")
   }

   if (length(c.ref) != dim(c.pred)[1]) {
      stop("Number of objects in reference and predicted results should be the same.")
   }

   # remove excluded rows for correct calculation of performance
   dim(c.ref) <- NULL
   attrs <- mda.getattr(c.pred)
   if (length(attrs$exclrows) > 0) {
      c.pred <- c.pred[-attrs$exclrows, , , drop = F]
      c.ref <- c.ref[-attrs$exclrows]
   }

   nobj <- dim(c.pred)[1]
   ncomp <- dim(c.pred)[2]
   nclasses <- dim(c.pred)[3]

   tp <- matrix(0, nrow = nclasses, ncol = ncomp)
   fp <- matrix(0, nrow = nclasses, ncol = ncomp)
   fn <- matrix(0, nrow = nclasses, ncol = ncomp)
   tn <- matrix(0, nrow = nclasses, ncol = ncomp)

   # compute main performance indicators
   classnames <- dimnames(c.pred)[[3]]
   for (i in seq_len(nclasses)) {
      fn[i, ] <- colSums((c.ref == classnames[i]) & (c.pred[, , i, drop = F] == -1))
      fp[i, ] <- colSums((c.ref != classnames[i]) & (c.pred[, , i, drop = F] == 1))
      tp[i, ] <- colSums((c.ref == classnames[i]) & (c.pred[, , i, drop = F] == 1))
      tn[i, ] <- colSums((c.ref != classnames[i]) & (c.pred[, , i, drop = F] == -1))
   }

   # compute main statistics
   sn <- tp / (tp + fn)
   sp <- tn / (tn + fp)
   ms <- (fp + fn) / (tp + tn + fp + fn)

   # add row with summary for all classes
   sn <- rbind(sn, colSums(tp) / colSums(tp + fn))
   sp <- rbind(sp, colSums(tn) / colSums(tn + fp))
   ms <- rbind(ms, colSums(fp + fn) / colSums(tp + tn + fp + fn))

   # add names
   row_names <- dimnames(c.pred)[[3]]
   col_names <- dimnames(c.pred)[[2]]
   rownames(fn) <- rownames(fp) <- rownames(tp) <- rownames(tn) <- row_names
   colnames(fn) <- colnames(fp) <- colnames(tp) <- colnames(sn) <- colnames(sp) <- col_names
   rownames(sn) <- rownames(sp) <- rownames(ms) <- c(row_names, "Total")

   return(
      list(
         "fn" = fn,
         "fp" = fp,
         "tp" = tp,
         "tn" = tn,
         "sensitivity" = sn,
         "specificity" = sp,
         "misclassified" = ms
      )
   )
}

#' Confusion matrix for classification results
#'
#' @details
#' Returns confusion matrix for classification results represented by the object.
#'
#' @param obj
#' classification results (object of class \code{simcares}, \code{simcamres}, etc)
#' @param ncomp
#' number of components to make the matrix for (NULL - use selected for a model).
#' @param ...
#' other arguments
#'
#' @description
#' The columns of the matrix correspond to classification results, rows - to the real classes. In
#' case of soft classification with multiple classes (e.g. SIMCAM) sum of values for every row
#' will not correspond to the total number of class members as the same object can be classified
#' as a member of several classes or non of them.
#'
#' @export
getConfusionMatrix.classres <- function(obj, ncomp = obj$ncomp.selected, ...) {

   if (is.null(obj$c.ref)) {
      stop("Reference classes are not available!")
   }

   attrs <- mda.getattr(obj$c.pred)
   c.pred <- obj$c.pred[, ncomp, , drop = FALSE]
   c.ref <- obj$c.ref

   # remove excluded rows
   if (length(attrs$exclrows) > 0) {
      c.pred <- c.pred[-attrs$exclrows, , ,drop = FALSE]
      c.ref <- c.ref[-attrs$exclrows]
   }

   # get class names and numbers
   ref.classes <- levels(c.ref)
   ref.nclasses <- length(ref.classes)

   # compute the confusion matrix
   out <- matrix(0, nrow = ref.nclasses, ncol = obj$nclasses + 1)
   none <- rep(TRUE, length(c.ref))
   for (i in seq_len(obj$nclasses)) {
      ind <- c.pred[, , i] > 0
      out[, i] <- table(c.ref[ind], exclude = FALSE)
      none[ind] <- FALSE
   }

   # find ones that were not classified as member of any class and names
   out[, obj$nclasses + 1] <- table(c.ref[none], exclude = FALSE)
   rownames(out) <- ref.classes
   colnames(out) <- c(obj$classnames, "None")
   return(out)
}

#' Show predicted class values
#'
#' @description
#' Shows a table with predicted class values for classification result.
#'
#' @param obj
#' object with classification results (e.g. \code{plsdares} or \code{simcamres}).
#' @param ncomp
#' number of components to show the predictions for (NULL - use selected for a model).
#' @param ...
#' other parameters
#'
#' @details
#' The function prints a matrix where every column is a class and every row is an data object.
#' The matrix has either -1 (does not belong to the class) or +1 (belongs to the class) values.
#'
#' @export
showPredictions.classres <- function(obj, ncomp = obj$ncomp.selected, ...) {

   if (obj$nclasses == 1) {
      pred = obj$c.pred[ , ncomp[1], , drop = F]
   } else {
      if (length(ncomp) == 1)
         ncomp = matrix(ncomp, nrow = 1, ncol = obj$nclasses)

      pred = NULL
      for (i in 1:obj$nclasses) {
         pred = cbind(pred, obj$c.pred[, ncomp[i], i, drop = F])
      }
   }

   dim(pred) <- c(nrow(pred), obj$nclasses)
   dimnames(pred) <- list(dimnames(obj$c.pred)[[1]], dimnames(obj$c.pred)[[3]])

   print(pred)
}

#' Plot for class belonging probability
#'
#' @description
#' Makes a plot with class belonging probabilities for each object of the classification results.
#' Works only with classification methods, which compute this probability (e.g. SIMCA).
#'
#' @param obj
#' classification results (e.g. object of class \code{simcamres}).
#' @param ncomp
#' number of components to use the probabilities for.
#' @param nc
#' if there are several classes, which class to make the plot for.
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ylim
#' vector with limits for y-axis
#' @param show.lines
#' shows a horizontal line at p = 0.5
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @export
plotProbabilities.classres = function(obj, ncomp = obj$ncomp.selected, nc = 1, type = 'h',
                                      xlab = 'Objects', ylab = 'Probability', main = NULL,
                                      ylim = c(0, 1.1), show.lines = c(NA, 0.5), ...) {
   if (is.null(obj$p.pred)) {
      stop('No probability values are available!')
   }

   if (nc > obj$nclasses || nc < 1) {
      stop('Wrong value for argument "nc"!')
   }

   classname = dimnames(obj$p.pred)[[3]][nc]

   if (is.null(main)) {
      main = sprintf('Class probabilities, %s (ncomp = %d)', classname, ncomp);
   }

   mdaplot(obj$p.pred[, ncomp, nc], show.lines = show.lines, type = type,
           ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
}


#' Sensitivity plot for classification results
#'
#' @description
#' Makes a plot with sn values vs. model complexity (e.g. number of components) for
#' classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotSensitivity.classres <- function(obj, nc = NULL, ...) {
   return(plotPerformance(obj, nc = nc, param = "sensitivity", ...))
}

#' Specificity plot for classification results
#'
#' @description
#' Makes a plot with specificity values vs. model complexity (e.g. number of components) for
#' classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotSpecificity.classres <- function(obj, nc = NULL, ...) {
   return(plotPerformance(obj, nc = nc, param = "specificity", ...))
}

#' Misclassified ratio plot for classification results
#'
#' @description
#' Makes a plot with ms ratio values vs. model complexity (e.g. number of components) for
#' classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotMisclassified.classres = function(obj, nc = NULL, ...) {
   return(plotPerformance(obj, nc = nc, param = "misclassified", ...))
}


#' Performance plot for classification results
#'
#' @description
#' Makes a plot with classification performance parameters vs. model complexity (e.g. number of
#' components) for classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ncomp
#' number of components to make the plot for
#' @param param
#' which performance parameter to make the plot for (\code{'sn'}, \code{'specificity'},
#' \code{'ms'}, \code{'all'}).
#' @param type
#' type of the plot
#' @param legend
#' vector with legend items
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ylim
#' vector with two values - limits for y axis
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplot}} function can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotPerformance.classres = function(obj, nc = NULL, ncomp = obj$ncomp.selected,
      param = 'all', type = 'h', legend = NULL, main = NULL, xlab = 'Components',
      ylab = '', ylim = c(0, 1.1), ...) {

   if (is.null(nc)) {
      nc =  obj$nclasses + 1

      if (obj$nclasses > 1)
         classname = '(all classes)'
      else
         classname = ''
   } else {
      if (nc > obj$nclasses || nc < 1)
         stop('Wrong value for argument "nc"!')

      classname = sprintf('(%s)', dimnames(obj$c.pred)[[3]][nc])
   }

   if (param == 'all') {
      if (is.null(main))
         main = sprintf('Prediction performance %s', classname);

      data = list()
      if (!any(is.na(obj$sn[nc, ])))
          data$sn = obj$sn[nc, ]
      if (!any(is.na(obj$specificity[nc, ])))
         data$specificity = obj$specificity[nc, ]
      if (!any(is.na(obj$ms[nc, ])))
         data$ms = obj$ms[nc, ]

      mdaplotg(data, type = type, legend = legend, main = main, xticks = 1:obj$ncomp,
               xlab = xlab, ylim = ylim, ylab = ylab, ...)
   } else {
      if (is.null(main))
         main = sprintf('%s%s %s', toupper(substring(param, 1, 1)), substring(param, 2), classname)

      data = obj[[param]][nc, , drop = F]
      if (any(is.na(data)))
         stop('This performance parameter has NA values!')

      mdaplot(data, type = type, main = main, xticks = 1:obj$ncomp, xlab = xlab, ylab = ylab,
              ylim = ylim, ...)
   }
}

#' Prediction plot for classification results
#'
#' @description
#' Makes a plot with predicted class values for classification results.
#'
#' @param obj
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ncomp
#' which number of components to make the plot for (one value, if NULL - model selected number will
#' be used).
#' This parameter shal not be used for multiclass models or results as predictions in this case are
#'  only
#' for optimal number of components
#' @param type
#' type of the plot
#' @param main
#' main title for the plot
#' @param ylab
#' label for y axis
#' @param ...
#' most of the graphical parameters from \code{\link{mdaplotg}} or \code{\link{mdaplot}} function
#' can be used.
#'
#' @details
#' See examples in description of \code{\link{plsdares}}, \code{\link{simcamres}}, etc.
#'
#' @export
plotPredictions.classres <- function(obj, nc = NULL, ncomp = NULL, type = 'p',
                                    main = NULL, ylab = '', ...) {

   if (is.null(obj))
      stop('No classification results were provided!')

   # get classnames
   classnames = dimnames(obj$c.pred)[[3]]

   # if vector with class number was not provided show plot for all of them
   if (is.null(nc))
      nc = 1:obj$nclasses

   # if vector with names were provided instead of numbers - convert to numbers
   if (is.character(nc))
      nc = which(classnames %in% nc)

   # take unique values and check correctness
   nc = unique(nc)
   if (length(nc) == 0 || min(nc)< 1 || max(nc) > dim(obj$c.pred)[3])
      stop('Incorrect class numbers or names!')
   nclasses = length(nc)

   # set main title
   if (is.null(main)) {
      main = 'Predictions'
      if (!is.null(ncomp) && length(ncomp) == 1)
         main = sprintf('%s (ncomp = %d)', main, ncomp)
   }

   ncomp = getSelectedComponents.classres(obj, ncomp)

   if (max(ncomp) > dim(obj$c.pred)[2])
      stop('Wrong value for ncomp parameter!')

   # extract predicted values for particular component
   attrs = mda.getattr(obj$c.pred)
   c.pred = as.matrix(obj$c.pred[, ncomp, nc])

   if (nclasses > 1) {
      # prepare matrix for the results
      pdata = matrix(0, nrow = nrow(c.pred), ncol = nclasses + 1)
      # find objects that were not classified for any of the class and set rows = 1 in first column
      pdata[, 1] = apply(c.pred, 1, function(x)(all(x == -1)))
      # set values for the other
      for (i in 1:nclasses)
         pdata[, i + 1] = (c.pred[, i] == 1) * (i + 1)
      # unfold the matrix
      dim(pdata) = NULL
      # add evector with indices of objects
      pdata = cbind(rep(1:nrow(c.pred), nclasses + 1), pdata)
      # remove all rows with zeros in the second column
      pdata = pdata[pdata[, 2] != 0, ]
      # assign proper rownames
      rownames(pdata) = rownames(c.pred)[pdata[, 1]]
      # assign other attributes
      pdata = mda.setattr(pdata, attrs)
      # exclude rows
      attr(pdata, 'exclrows') = NULL
      pdata = mda.exclrows(pdata, pdata[, 1] %in% attrs$exclrows)
      row.ind = pdata[, 1]
   } else {
      pdata = cbind((1:nrow(c.pred)), (c.pred == 1) + 1)
      pdata = mda.setattr(pdata, attrs)
   }

   colnames(pdata) = c(ifelse(is.null(attrs$yaxis.name), 'Objects', attrs$yaxis.name), 'Classes')

   if (is.null(obj$c.ref)) {
      mdaplot(pdata, type = 'p', main = main, ylab = ylab, yticks = c(1:(nclasses + 1)),
              yticklabels = c('None', classnames[nc]), ylim = c(0.8, nclasses + 1.2), ...)
   } else {
      pdata.g = as.factor(obj$c.ref[pdata[, 1]])
      mdaplotg(pdata, type = 'p', main = main, ylab = ylab, yticks = c(1:(nclasses + 1)),
               groupby = pdata.g,
              yticklabels = c('None', classnames[nc]), ylim = c(0.8, nclasses + 1.2), ...)
   }
}

#' Plot function for classification results
#'
#' @description
#' Generic plot function for classification results. Shows predicted class values.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' other arguments
#'
#' @export
plot.classres = function(x, nc = NULL, ...){
   plotPredictions.classres(x, nc = nc, ...)
}

#' as.matrix method for classification results
#'
#' @description
#' Generic \code{as.matrix} function for classification results. Returns matrix with performance
#' values for specific class.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param ncomp
#' model complexity (number of components) to show the parameters for.
#' @param nc
#' if there are several classes, which class to show the parameters for.
#' @param ...
#' other arguments
#'
#' @export
as.matrix.classres = function(x, ncomp = NULL, nc = 1, ...) {
   obj = x

   if (is.null(obj$c.ref))
      return()

   if (is.null(ncomp))
      res = cbind(obj$tp[nc, ], obj$fp[nc, ], obj$tn[nc, ], obj$fn[nc, ],
               round(obj$specificity[nc, ], 3), round(obj$sn[nc, ], 3))
   else
      res = cbind(obj$tp[nc, ncomp], obj$fp[nc, ncomp], obj$tn[nc, ncomp], obj$fn[nc, ncomp],
               round(obj$specificity[nc, ncomp], 3), round(obj$sn[nc, ncomp], 3))

   colnames(res) = c('TP', 'FP', 'TN', 'FN', 'Spec', 'Sens')

   if (any(is.na(obj$specificity)))
      res = res[, -5]

   res
}

#' Print information about classification result object
#'
#' @description
#' Generic \code{print} function for classification results. Prints information about major fields
#' of the object.
#'
#' @param x
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param str
#' User specified text (e.g. to be used for particular method, like PLS-DA, etc).
#' @param ...
#' other arguments
#'
#' @export
print.classres = function(x, str = NULL, ...)
{
   if (is.null(str))
      str = 'Classification results (class classres)\nMajor fields:'

   if (nchar(str) > 0)
      cat(sprintf('\n%s\n', str))

   cat('$c.pred - predicted class values\n')
   if (!is.null(x$c.ref))
   {
      cat('$c.ref - reference (true) class values\n')
      cat('$tp - number of true positives\n')
      cat('$fp - number of false positives\n')
      cat('$fn - number of false negatives\n')
      cat('$specificity - specificity of predictions\n')
      cat('$sn - sn of predictions\n')
      cat('$ms - misclassification ratio for predictions\n')
   }
}

#' Summary statistics about classification result object
#'
#' @description
#' Generic \code{summary} function for classification results. Prints performance values for the
#' results.
#'
#' @param object
#' classification results (object of class \code{plsdares}, \code{simcamres}, etc.).
#' @param ncomp
#' which number of components to make the plot for (can be one value for all classes or vector with
#' separate values for each, if NULL - model selected number will be used).
#' @param nc
#' if there are several classes, which class to make the plot for (NULL - summary for all classes).
#' @param ...
#' other arguments
#'
#' @export
summary.classres = function(object, ncomp = NULL, nc = NULL, ...) {
   cat('\nClassiciation results (class classres) summary\n')
   if (!is.null(object$c.ref)) {

      if (is.null(nc))
         nc = 1:dim(object$c.pred)[3]

      if (!is.null(ncomp)) {
         cat(sprintf('\nNumber of selected components: %d', ncomp))
      } else {
         if (is.null(object$ncomp.selected)) {
            ncomp = 1
         } else {
            ncomp = object$ncomp.selected
            cat(sprintf('\nNumber of selected components: %d', ncomp))
         }
      }

      cat(sprintf('\nNumber of classes: %d\n', ncol(object$c.ref)))

      for (i in nc) {
         cat(sprintf('\nClass "%s":\n', colnames(object$c.ref)[i]))
         res = as.matrix.classres(object, nc = i)
         print(res)
      }
   } else {
      cat('No reference data provided to calculate prediction performance.')
   }
}