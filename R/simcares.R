#'  Results of SIMCA one-class classification
#'  
#'  @description 
#' \code{simcares} is used to store results for SIMCA one-class classification.
#' @param pres
#' results of PCA decomposition of data (class \code{pcares}).
#' @param cres
#' results of classification (class \code{classres}).
#'
#' @details 
#' Class \code{simcares} inherits all properties and methods of class \code{\link{pcares}}, and 
#' has additional properties and functions for representing of classification results, inherited 
#' from class \code{\link{classres}}.
#'   
#' There is no need to create a \code{simcares} object manually, it is created automatically when 
#' build a SIMCA model (see \code{\link{simca}}) or apply the model to a new data (see 
#' \code{\link{predict.simca}}). The object can be used to show summary and plots for the results.
#' 
#' @return
#' Returns an object (list) of class \code{simcares} with the same fields as \code{\link{pcares}} 
#' plus extra fields, inherited from \code{\link{classres}}:
#' \item{c.pred}{predicted class values (+1 or -1).}
#' \item{c.ref}{reference (true) class values if provided.}
#' 
#' The following fields are available only if reference values were provided.
#' \item{tp}{number of true positives.}
#' \item{fp}{nmber of false positives.}
#' \item{fn}{number of false negatives.}
#' \item{specificity}{specificity of predictions.}
#' \item{sensitivity}{sensitivity of predictions.}
#'
#' @seealso 
#' Methods for \code{simcares} objects:
#' \tabular{ll}{
#'  \code{print.simcares} \tab shows information about the object.\cr
#'  \code{summary.simcares} \tab shows statistics for results of classification.\cr
#' }
#' 
#' Methods, inherited from \code{\link{classres}} class:
#' \tabular{ll}{
#'  \code{\link{showPredictions.classres}} \tab show table with predicted values.\cr
#'  \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
#'  \code{\link{plotSensitivity.classres}} \tab makes plot with sensitivity vs. components values.\cr
#'  \code{\link{plotSpecificity.classres}} \tab makes plot with specificity vs. components values.\cr
#'  \code{\link{plotPerformance.classres}} \tab makes plot with both specificity and sensitivity 
#'  values.\cr
#' }
#' 
#' Methods, inherited from \code{\link{ldecomp}} class:   
#' \tabular{ll}{
#'  \code{\link{plotResiduals.ldecomp}} \tab makes Q2 vs. T2 residuals plot.\cr
#'  \code{\link{plotScores.ldecomp}} \tab makes scores plot.\cr
#'  \code{\link{plotVariance.ldecomp}} \tab makes explained variance plot.\cr
#'  \code{\link{plotCumVariance.ldecomp}} \tab makes cumulative explained variance plot.\cr
#' }
#' Check also \code{\link{simca}} and \code{\link{pcares}}.    
#' 
#' @examples
#' ## make a SIMCA model for Iris setosa class and show results for calibration set
#' library(mdatools)
#' 
#' data = iris[, 1:4]
#' class = iris[, 5]
#' 
#' # take first 30 objects of setosa as calibration set 
#' se = data[1:30, ]
#' 
#' # make SIMCA model and apply to test set
#' model = simca(se, 'Se')
#' model = selectCompNum(model, 1)
#' 
#' # show infromation and summary
#' print(model$calres)
#' summary(model$calres)
#' 
#' # show plots
#' layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
#' plotPredictions(model$calres, show.labels = TRUE)
#' plotResiduals(model$calres, show.labels = TRUE)
#' plotPerformance(model$calres, show.labels = TRUE, legend.position = 'bottomright')
#' layout(1, 1, 1)
#' 
#' # show predictions table
#' showPredictions(model$calres)
#' @export
simcares = function(pres, cres) {
   res = c(pres, cres)
   res$classname = dimnames(cres$c.pred)[[3]][1]
   class(res) = c('simcares', 'classres', 'pcares', 'ldecomp')   
   
   res
}   

#' Residuals plot for SIMCA results
#' 
#' @description
#' Shows a plot with Q vs. T2 residuals for SIMCA results
#' 
#' @param obj
#' SIMCA results (object of class \code{simcares})
#' @param ncomp
#' which principal components to show the plot for
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param norm
#' logical, show normalized Q vs T2 (\code{norm = T}) values or original ones (\code{norm = F})
#' @param show.limits
#' logical, show or not lines with statistical limits for the residuals
#' @param legend
#' vector with legend items
#' @param lim.col
#' vector with two values - line color for extreme and outlier borders 
#' @param lim.lwd
#' vector with two values - line width for extreme and outlier borders 
#' @param lim.lty
#' vector with two values - line type for extreme and outlier borders 
#' @param ...
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcares}} function.
#' 
#' @export
plotResiduals.simcares = function(obj, ncomp = NULL, main = NULL, xlab = NULL, ylab = NULL,
                                  norm = F, show.limits = T, legend = NULL, 
                                  lim.col = c('#c0a0a0', '#906060'), 
                                  lim.lwd = c(1, 1), lim.lty = c(2, 3), ...) {
   # Shows residuals plot (T2 vs Q) 
   #
   # Arguments:
   #  obj: SIMCA results (an object of class simcares)
   #  ncomp: number of components to make the plot for
   #  show.limits: logical, show or not statistical limits for the residual values
   #  ...: standard arguments for plots

   # set main title
   if (is.null(main)) {
      if (is.null(ncomp))
         main = 'Residuals'
      else
         main = sprintf('Residuals (ncomp = %d)', ncomp)
   }

   # get selected components and ncomp is NULL
   ncomp = getSelectedComponents(obj, ncomp)

   if (is.null(obj$c.ref)) {
      # if no reference values or all objects are from the same class 
      # show standard plot for PCA
      plotResiduals.ldecomp(obj, ncomp, main = main, xlab = xlab, ylab = ylab, norm = norm,
                            show.limits = show.limits, lim.col = lim.col, lim.lty = lim.lty,
                            lim.lwd = lim.lwd, ...)
   } else {   
      # if objects include members and non-members show plot with
      # color differentiation and legend
      
      c.ref = obj$c.ref     
      if (!is.character(c.ref)) {   
         c.ref[obj$c.ref == 1] = obj$classname
         c.ref[obj$c.ref != 1] = 'Others'
      }
      
      if (sum(c.ref == obj$classname) == length(c.ref)) {
         plotResiduals.pcares(obj, ncomp, main = main, xlab = xlab, ylab = ylab, norm = norm,
                               show.limits = show.limits, lim.col = lim.col, lim.lty = lim.lty,
                               lim.lwd = lim.lwd, ...)
      } else {
         classes = unique(c.ref)
         nclasses = length(classes)

         pdata = list()
         legend.str = NULL
         
         # set values for normalization of residuals if necessary
         if (norm) {
            T2.mean = obj$T2lim[3, ncomp]
            Q.mean = obj$Qlim[3, ncomp]
            if (is.null(xlab))
               xlab = expression(paste('Hotelling ', T^2, ' distance (norm)'))
            if (is.null(ylab))
               ylab = 'Squared residual distance, Q (norm)'      
         } else {
            T2.mean = 1
            Q.mean = 1
            if (is.null(xlab))
               xlab = expression(paste('Hotelling ', T^2, ' distance'))
            if (is.null(ylab))
               ylab = 'Squared residual distance, Q'      
         }
         
         for (i in 1:nclasses) {   
            idx = c.ref == classes[i]
            data = cbind(obj$T2[idx, ncomp]/T2.mean, obj$Q[idx, ncomp]/Q.mean)
            rownames(data) = rownames(obj$c.ref[idx])
            legend.str = c(legend.str, classes[i])
            pdata[[i]] = data
         }

         if (is.null(legend))
            legend = legend.str
         
         mdaplotg(pdata, type = 'p', xlab = xlab, ylab = ylab, main = main, legend = legend, ...)
         
         if (show.limits) {
            # get residual limits, correct if necessary and recalculate axes maximum limit
            lim = cbind(obj$T2lim[1:2, ncomp], obj$Qlim[1:2, ncomp])
            if (substr(obj$lim.type, 1, 2) != 'dd') {
               lim[, 1] = lim[, 1] / T2.mean
               lim[, 2] = lim[, 2] / Q.mean
            } else {
               lim[, 1] = lim[, 1] * T2.mean / Q.mean
               lim[, 2] = lim[, 2] / Q.mean
            }
            ldecomp.plotLimits(lim, obj$lim.type, lim.col, lim.lwd, lim.lty)   
         }
      }
   }   
}


#' Summary method for SIMCA results object
#' 
#' @description
#' Shows performance statistics for the results.
#' 
#' @param object
#' SIMCA results (object of class \code{simcares})
#' @param ...
#' other arguments
#' 
#' @export
summary.simcares = function(object, ...) {
   # Show summary for simcares object
   
   obj = object
   
   cat('\nSummary for SIMCA one-class classification result\n')
   cat(sprintf('\nClass name: %s\n', obj$classname))
   cat(sprintf('Number of selected components: %d\n', obj$ncomp.selected))
   cat('\n')
   pcares = as.matrix.ldecomp(obj)
   calres = as.matrix.classres(obj)    
   print(cbind(round(pcares, 2), calres))
}

#' Print method for SIMCA results object
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
print.simcares = function(x, ...) {   
   # Print information on simcares object
   obj = x
   
   cat('Result for SIMCA one-class classification (class simcares)\n')
   cat(sprintf('Method for critical limits: %s\n', obj$lim.type))
   print.ldecomp(obj, '')
   print.classres(obj, '')
   cat('\n')
}

