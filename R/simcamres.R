#' Results of SIMCA multiclass classification
#' 
#' @description
#' \code{simcamres} is used to store results for SIMCA multiclass classification.
#' 
#' @param cres
#' results of classification (class \code{classres}).
#' @param T2
#' matrix with T2 values for each object and class.
#' @param Q
#' matrix with Q values for each object and class.
#' @param T2lim
#' vector with T2 statistical limits for each class.
#' @param Qlim
#' vector with Q statistical limits for each class.
#' 
#' @details 
#' Class \code{simcamres} inherits all properties and methods of class \code{\link{classres}}, plus 
#' store values necessary to visualise prediction decisions (e.g. Cooman's plot or Residuals plot).
#' 
#' In cotrast to \code{simcares} here only values for optimal (selected) number of components in 
#' each individual SIMCA models are presented.
#' 
#' There is no need to create a \code{simcamres} object manually, it is created automatically when 
#' make a SIMCAM model (see \code{\link{simcam}}) or apply the model to a new data (see 
#' \code{\link{predict.simcam}}). The object can be used to show summary and plots for the results.
#'
#' @return
#' Returns an object (list) of class \code{simcamres} with the same fields as \code{\link{classres}} 
#' plus extra fields for Q and T2 values and limits:
#' 
#' \item{c.pred}{predicted class values.}
#' \item{c.ref}{reference (true) class values if provided.}
#' \item{T2}{matrix with T2 values for each object and class.}
#' \item{Q}{matrix with Q values for each object and class.}
#' \item{T2lim}{vector with T2 statistical limits for each class.}
#' \item{Qlim}{vector with Q statistical limits for each class.}
#' 
#' The following fields are available only if reference values were provided.
#' \item{tp}{number of true positives.}
#' \item{fp}{nmber of false positives.}
#' \item{fn}{number of false negatives.}
#' \item{specificity}{specificity of predictions.}
#' \item{sensitivity}{sensitivity of predictions.}
#'
#' @seealso 
#' Methods for \code{simcamres} objects:
#' \tabular{ll}{
#'  \code{print.simcamres} \tab shows information about the object.\cr
#'  \code{summary.simcamres} \tab shows statistics for results of classification.\cr
#'  \code{\link{plotResiduals.simcamres}} \tab makes Q vs. T2 residuals plot.\cr
#'  \code{\link{plotCooman.simcamres}} \tab makes Cooman's plot.\cr
#' }
#' 
#' Methods, inherited from \code{\link{classres}} class:
#' \tabular{ll}{
#'  \code{\link{showPredictions.classres}} \tab show table with predicted values.\cr
#'  \code{\link{plotPredictions.classres}} \tab makes plot with predicted values.\cr
#' }
#' 
#' Check also \code{\link{simcam}}.    
#'
#' @examples
#' ## make a multiclass SIMCA model for Iris data and apply to test set
#' library(mdatools)
#' 
#' # split data 
#' caldata = iris[seq(1, nrow(iris), 2), 1:4]
#' se = caldata[1:25, ]
#' ve = caldata[26:50, ]
#' vi = caldata[51:75, ]
#' 
#' testdata = iris[seq(2, nrow(iris), 2), 1:4]
#' testdata.cref = iris[seq(2, nrow(iris), 2), 5]
#' 
#' # create individual models
#' semodel = simca(se, classname = 'setosa')
#' semodel = selectCompNum(semodel, 1)
#' 
#' vimodel = simca(vi, classname = 'virginica')
#' vimodel = selectCompNum(vimodel, 1)
#' 
#' vemodel = simca(ve, classname = 'versicolor')
#' vemodel = selectCompNum(vemodel, 1)
#' 
#' # combine models into SIMCAM object, show statistics 
#' model = simcam(list(semodel, vimodel, vemodel), info = 'Iris data')
#' res = predict(model, testdata, testdata.cref)
#' summary(res)
#' 
#' # show predicted values
#' showPredictions(res)
#' 
#' # plot predictions
#' par(mfrow = c(2, 2))
#' plotPredictions(res)
#' plotPredictions(res, nc = 1)
#' plotPredictions(res, nc = c(1, 2))
#' plotPredictions(res, show.labels = TRUE)
#' par(mfrow = c(1, 1))
#' 
#' # show residuals and Cooman's plot
#' 
#' par(mfrow = c(2, 2))
#' plotCooman(res)
#' plotCooman(res, nc = c(1, 3))
#' plotResiduals(res)
#' plotResiduals(res, nc = 3)
#' par(mfrow = c(1, 1))
#' 
#' @export
simcamres = function(cres, T2, Q, T2lim, Qlim)
{
  res = cres
   res$T2 = T2
   res$Q = Q
   res$T2lim = T2lim
   res$Qlim = Qlim
   res$classnames = dimnames(cres$c.pred)[[3]]
   class(res) = c('simcamres', 'classres')   
   
   res
}   

#' Residuals plot for SIMCAM results
#' 
#' @description
#' Shows a plot with Q vs. T2 residuals for SIMCAM results
#' 
#' @param obj
#' SIMCAM results (object of class \code{simcamres})
#' @param nc
#' which class (SIMCA model) to show the plot for
#' @param show.limits
#' logical, show or not lines with statistical limits for the residuals
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param legend
#' vector with legend items
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcamres}} function.
#' 
#' @export
plotResiduals.simcamres = function(obj, nc = 1, show.limits = T, type = 'p', main = NULL, 
                                  xlab = 'T2', ylab = 'Squared residual distance (Q)', legend = NULL, ...)
{
   # set main title
   if (is.null(main))
      main = sprintf('Residuals (%s)', obj$classnames[nc])

   if (!is.null(obj$c.ref))
   {   
      classes = unique(obj$c.ref)
      data = list()
      for (i in 1:length(classes))
         data[[i]] = cbind(obj$T2[obj$c.ref == classes[i], nc], obj$Q[obj$c.ref == classes[i], nc])
      
      if (is.null(legend))
         legend = classes
      
      if (show.limits == T)
         show.lines = c(obj$T2lim[nc], obj$Qlim[nc])
      else
         show.lines = F
      
      mdaplotg(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines,
               legend = legend, ...)
   }
   else
   {

      data = cbind(obj$T2[, nc], obj$Q[, nc])
      if (show.limits == T)
         show.lines = c(obj$T2lim[nc], obj$Qlim[nc])
      else
         show.lines = F
      
      mdaplot(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines, ...)
   }
}

#' Cooman's plot for SIMCAM results
#' 
#' @description
#' Shows a Cooman's plot for a pair of SIMCA models
#' 
#' @param obj
#' SIMCAM results (object of class \code{simcamres})
#' @param nc
#' vector with two values - classes (SIMCA models) to show the plot for
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param show.limits
#' logical, show or not lines with statistical limits for the residuals
#' @param legend
#' vector with legend items 
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcamres}} function.
#' 
#' @export
plotCooman.simcamres = function(obj, nc = c(1, 2), type = 'p', main = "Cooman's plot", xlab = NULL, 
                                ylab = NULL, show.limits = T, legend = NULL, ...)
{
   # set labels for axes
   if (is.null(xlab))
      xlab = sprintf('Distance to class %s', obj$classnames[nc[1]])

   if (is.null(ylab))
      ylab = sprintf('Distance to class %s', obj$classnames[nc[2]])
   
   if (!is.null(obj$c.ref))
   {   
      classes = unique(obj$c.ref)
      data = list()
      for (i in 1:length(classes))
         data[[i]] = cbind(sqrt(obj$Q[obj$c.ref == classes[i], nc[1]]), 
                           sqrt(obj$Q[obj$c.ref == classes[i], nc[2]])
                           )
      
      if (is.null(legend))
         legend = classes
      
      if (show.limits == T)
         show.lines = c(obj$Qlim[nc[1]], obj$Qlim[nc[2]])
      else
         show.lines = F
      
      mdaplotg(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines,
               legend = legend, ...)
   }
   else
   {
      data = cbind(sqrt(obj$Q[, nc[1]]), sqrt(obj$Q[, nc[2]]))
      
      if (show.limits == T)
         show.lines = c(obj$Qlim[nc[1]], obj$Qlim[nc[2]])
      else
         show.lines = F
      
      mdaplot(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines, ...)
   }
}

#' Model overview plot for SIMCAM results
#' 
#' @description
#' Just shows a prediction plot for SIMCAM results.
#' 
#' @param x
#' SIMCAM results (object of class \code{simcamres})
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{simcamres}} function.
#' 
#' @export
plot.simcamres = function(x, ...)
{
   plotPredictions(x)
}

#' Summary method for SIMCAM results object
#' 
#' @description
#' Shows performance statistics for the results.
#' 
#' @param object
#' SIMCAM results (object of class \code{simcamres})
#' @param ...
#' other arguments
#' 
#' @export
summary.simcamres = function(object, ...)
{
   obj = object
   
   cat('\nSummary for SIMCA multiple classes classification result\n')
   if (!is.null(obj$c.ref))
   {   
      classres = NULL
      for (i in 1:obj$nclasses)
         classres = rbind(classres, as.matrix.classres(obj, nc = i))
      rownames(classres) = paste(obj$classnames, ' (', obj$ncomp.selected, ' comp)', sep = '')   
      print(classres)
   }
   else
   {
      cat('\nReference values are not provided.\n')
   }
}  

#' Print method for SIMCAM results object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' SIMCAM results (object of class \code{simcamres})
#' @param ...
#' other arguments
#'
#' @export 
print.simcamres = function(x, ...)
{   
   obj = x
   
   cat('Result for SIMCA multiple classes classification (class simcamres)\n')
   print.classres(obj, '')
   cat('\n')
}
