#' Results of PCA decomposition
#' @description 
#' \code{pcares} is used to store results for PCA decomposition of data.
#'
#' @param ...
#' other arguments supported by \code{ldecomp}.
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
#' @return 
#' Returns an object (list) of class \code{pcares} and \code{ldecomp} with following fields:
#' \item{scores }{matrix with score values (nobj x ncomp).}
#' \item{residuals }{matrix with data residuals (nobj x nvar).}
#' \item{T2 }{matrix with T2 distances (nobj x ncomp).}
#' \item{Q }{matrix with Q residuals (nobj x ncomp).}
#' \item{tnorm }{vector with singular values used for scores normalization.}
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
#'  \code{\link{plotResiduals.ldecomp}} \tab makes Q vs. T2 residuals plot.\cr
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
#' x.new = people[seq(1, 32, 2), ]
#' 
#' model = pca(people, scale = TRUE, cv = 1, info = 'Simple PCA model')
#' model = selectCompNum(model, 4)
#' 
#' res = predict(model, x.new)
#' summary(res)
#' plot(res)
#' 
#' ## 1. Make PCA model for People data with autoscaling
#' ## and full cross-validation and get calibration results
#' 
#' 
#' data(people)
#' model = pca(people, scale = TRUE, cv = 1, info = 'Simple PCA model')
#' model = selectCompNum(model, 4)
#' 
#' res = model$calres
#' summary(res)
#' plot(res)
#' 
#' ## 2. Show scores plots for the results
#' par(mfrow = c(2, 2))
#' plotScores(res)
#' plotScores(res, cgroup = people[, 'Beer'], show.labels = TRUE)
#' plotScores(res, comp = c(1, 3), show.labels = TRUE)
#' plotScores(res, comp = 2, type = 'h', show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' ## 3. Show residuals and variance plots for the results
#' par(mfrow = c(2, 2))
#' plotVariance(res, type = 'h')
#' plotCumVariance(res, show.labels = TRUE, legend.position = 'bottomright')
#' plotResiduals(res, show.labels = TRUE, cgroup = people[, 'Sex'])
#' plotResiduals(res, ncomp = 2, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#'
#' @export 
pcares = function(...) {
   # Creates an object of pcares class. In fact the class is a wrapper for ldecomp and
   # uses its methods and attributes.
   
   res = ldecomp(...)
   class(res) = c('pcares', 'ldecomp')   
   
   res
}   


#' Residuals plot for PCA results
#' 
#' @description
#' Shows a plot with T2 vs Q values for data objects.
#' 
#' @param obj
#' object of \code{ldecomp} class.
#' @param ncomp
#' what number of components to show the plot for (if NULL, model selected value will be used).
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param show.limits
#' logical, show or not lines for statistical limits of the residuals
#' @param norm
#' logical, show normalized Q vs T2 (\code{norm = T}) values or original ones (\code{norm = F})
#' @param xlim
#' limits for x-axis
#' @param ylim
#' limits for y-axis
#' @param lim.col
#' vector with two values - line color for extreme and outlier borders 
#' @param lim.lwd
#' vector with two values - line width for extreme and outlier borders 
#' @param lim.lty
#' vector with two values - line type for extreme and outlier borders 
#' @param ...
#' most of graphical parameters from \code{\link{mdaplot}} function can be used.
#' 
#' @export
plotResiduals.pcares = function(obj, ncomp = NULL, main = NULL, xlab = NULL, ylab = NULL, 
                                 show.labels = F, show.limits = T, norm = F, 
                                 xlim = NULL, ylim = NULL, 
                                 lim.col = c('#333333', '#333333'), 
                                 lim.lwd = c(1, 1), lim.lty = c(2, 3), ...) {
   if (is.null(main)) {
      if (is.null(ncomp))
         main = 'Residuals'
      else
         main = sprintf('Residuals (ncomp = %d)', ncomp)
   }
   
   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   
   data = mda.cbind(
      mda.subset(obj$T2, select = ncomp), 
      mda.subset(obj$Q, select = ncomp)
   )
   
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
   
   data[, 1] = data[, 1] / T2.mean
   data[, 2] = data[, 2] / Q.mean
   x.max = max(data[, 1])
   y.max = max(data[, 2])
   
   if (show.limits == T) {
      # get residual limits, correct if necessary and recalculate axes maximum limit
      lim = cbind(obj$T2lim[1:2, ncomp], obj$Qlim[1:2, ncomp])
      if (substr(obj$lim.type, 1, 2) != 'dd') {
         lim[, 1] = lim[, 1] / T2.mean
         lim[, 2] = lim[, 2] / Q.mean
         x.max = max(x.max, lim[, 1])
         y.max = max(y.max, lim[, 2])
      } else {
         lim[, 1] = lim[, 1] * T2.mean / Q.mean
         lim[, 2] = lim[, 2] / Q.mean
         x.max = 1.5 * x.max
         y.max = 1.5 * y.max
      }
   }
   
   # use computed max values for axes limits if user did not specify anything
   if (is.null(xlim))
      xlim = c(0, 1.2 * x.max)
   if (is.null(ylim))
      ylim = c(0, 1.2 * y.max)
   
   # show plot
   mdaplot(data, main = main, xlab = xlab, ylab = ylab, show.labels = show.labels, 
           xlim = xlim, ylim = ylim, ...)
   
   # show limits
   if (show.limits) {
      ldecomp.plotLimits(lim, obj$lim.type, lim.col, lim.lwd, lim.lty)   
   }   
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
#' @param show.labels
#' logical, show or not labels for the plot objects
#' @param ...
#' other arguments
#'
#' @export 
plot.pcares = function(x, comp = c(1, 2), show.labels = T, ...) {   
   par(mfrow = c(2, 2))
   plotScores(x, comp = comp, show.labels = show.labels, ...)
   plotResiduals(x, show.labels = show.labels, ...)
   plotVariance(x, show.labels = show.labels, ...)
   plotCumVariance(x, show.labels = show.labels, ...)
   par(mfrow = c(1, 1))
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
summary.pcares = function(object, ...) {
   summary.ldecomp(object, 'Summary for PCA results', ...)
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
print.pcares = function(x, ...) {   
   print.ldecomp(x, 'Results for PCA decomposition (class pcares)', ...)
   cat('\n')
}