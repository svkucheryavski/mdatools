#' SIMCA multiclass classification
#' 
#' @description 
#' \code{simcam} is used to combine several one-class SIMCA models for multiclass classification. 
#'
#' @param models
#' list with SIMCA models (\code{simca} objects).
#' @param info
#' text with information about the the object.
#' 
#' @details 
#' Besides the possibility for multiclass classification, SIMCAM also provides tools for
#' investigation of relationship among individual models (classes), such as discrimination power of 
#' variables, Cooman's plot, model distance, etc.
#' 
#' When create \code{simcam} object, the calibration data from all individual SIMCA models is 
#' extracted and combined for making predictions and calculate performance of the multi-class model. 
#' The results are stored in \code{$calres} field of the model object.
#' 
#' @return 
#' Returns an object of \code{simcam} class with following fields:
#' \item{models }{a list with provided SIMCA models.} 
#' \item{dispower }{an array with discrimination power of variables for each pair of individual 
#' models.} 
#' \item{moddist }{a matrix with distance between each each pair of individual models.} 
#' \item{classnames }{vector with names of individual classes.} 
#' \item{nclasses }{number of classes in the object.} 
#' \item{info }{information provided by user when create the object.} 
#' \item{calres }{an object of class \code{\link{simcamres}} with classification results for a 
#' calibration data.} 
#'
#' @seealso 
#' Methods for \code{simca} objects:
#' \tabular{ll}{
#'  \code{print.simcam} \tab shows information about the object.\cr
#'  \code{summary.simcam} \tab shows summary statistics for the models.\cr
#'  \code{plot.simcam} \tab makes an overview of SIMCAM model with two plots.\cr
#'  \code{\link{predict.simcam}} \tab applies SIMCAM model to a new data.\cr
#'  \code{\link{plotModelDistance.simcam}} \tab shows plot with distance between individual 
#'  models.\cr
#'  \code{\link{plotDiscriminationPower.simcam}} \tab shows plot with discrimination power.\cr
#'  \code{\link{plotModellingPower.simcam}} \tab shows plot with modelling power for individual 
#'  model.\cr
#'  \code{\link{plotCooman.simcam}} \tab shows Cooman's plot for calibration data.\cr
#'  \code{\link{plotResiduals.simcam}} \tab shows plot with Q vs. T2 residuals for calibration 
#'  data.\cr
#' }
#' 
#' Methods, inherited from \code{classmodel} class:
#' \tabular{ll}{
#'  \code{\link{plotPredictions.classmodel}} \tab shows plot with predicted values.\cr
#'  \code{\link{plotSensitivity.classmodel}} \tab shows sensitivity plot.\cr
#'  \code{\link{plotSpecificity.classmodel}} \tab shows specificity plot.\cr
#'  \code{\link{plotMisclassified.classmodel}} \tab shows misclassified ratio plot.\cr
#' }
#' 
#' Since SIMCAM objects and results are calculated only for optimal number of components, there is no
#' sense to show such plots like sensitivity or specificity vs. number of components. However they
#' are available as for any other classification model.
#'
#' @examples 
#' ## make a multiclass SIMCA model for Iris data
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
#' # combine models into SIMCAM objects, show statistics and plots
#' model = simcam(list(semodel, vimodel, vemodel), info = 'Iris data')
#' summary(model)
#' plot(model)
#' 
#' # show predictions and residuals for calibration data
#' par(mfrow = c(2, 2))
#' plotPredictions(model)
#' plotCooman(model, nc = c(1, 2))
#' plotResiduals(model, nc = 1)
#' plotResiduals(model, nc = 2)
#' par(mfrow = c(1, 1))
#' 
#' # show different plots for the model
#' par(mfrow = c(2, 2))
#' plotModelDistance(model, nc = 1)
#' plotDiscriminationPower(model, nc = c(1, 2))
#' plotModellingPower(model, nc = 1)
#' plotModellingPower(model, nc = 2)
#' par(mfrow = c(1, 1))
#' 
#' # apply the SIMCAM model to test set and show statistics and plots
#' res = predict(model, testdata, testdata.cref)
#' summary(res)
#' plotPredictions(res)
#'
#' @export     
simcam = function(models, info = '')
{
   nclasses = length(models)
   classnames = unlist(lapply(models, function(x){x$classname}))
   
   model = list()
   model$models = models
   model$classnames = classnames
   model$nclasses = nclasses
   model$info = info
   
   # calculate statistics
   stat = simcam.getPerformanceStatistics(model)
   model$dispower = stat$dispower
   model$moddist = stat$moddist
   
   model$call = match.call()   
   class(model) = c("simcam", "classmodel")
   
   # make predictions for calibration set
   caldata = getCalibrationData(model)
   model$calres = predict(model, caldata$x, caldata$c.ref)
   
   model
}

#' SIMCA multiple classes predictions
#' 
#' @description
#' Applies SIMCAM model (SIMCA for multiple classes) to a new data set
#' 
#' @param object
#' a SIMCAM model (object of class \code{simcam})
#' @param x
#' a matrix with x values (predictors)
#' @param c.ref
#' a vector with reference class names (same as class names in models)
#' @param cv
#' logical, are predictions for cross-validation or not
#' @param ...
#' other arguments
#' 
#' @return
#' SIMCAM results (an object of class \code{simcamres})
#'
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
#' @export
predict.simcam = function(object, x, c.ref = NULL, cv = F, ...) {
   attrs = mda.getattr(x)
   
   c.pred = array(0, dim = c(nrow(x), 1, object$nclasses))
   pred.res = list()
   for (i in 1:object$nclasses) {  
      pred.res[[i]] = predict(object$models[[i]], x, c.ref)
      c.pred[, , i] = pred.res[[i]]$c.pred[, object$models[[i]]$ncomp.selected, ]
   }
   dimnames(c.pred) = list(rownames(x), paste('Comp'), object$classnames)
   c.pred = mda.setattr(c.pred, attrs, 'row')
   attr(c.pred, 'name') = 'SIMCAM predictions'
   cres = classres(c.pred, c.ref)
   res = simcamres(cres, pred.res)
              
   res
}

#' Get calibration data
#' 
#' @description
#' Get data, used for calibration of the SIMCAM model.
#' 
#' @param obj
#' SIMCAM model (object of class \code{simcam})
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
#' @export
getCalibrationData.simcam = function(obj, ...) {
   x = getCalibrationData(obj$models[[1]])
   c.ref = matrix(obj$models[[1]]$classname, nrow = nrow(x), ncol = 1)            
   for (i in 2:obj$nclasses) {
      classdata = getCalibrationData(obj$models[[i]])
      x = mda.rbind(x, classdata)
      c.ref = rbind(c.ref, matrix(obj$models[[i]]$classname, nrow = nrow(classdata), ncol = 1))            
   }
   
   res = list(x = x, c.ref = c.ref)
   
   res
}

#' Performance statistics for SIMCAM model
#' 
#' @description
#' Calculates discrimination power and distance between models for SIMCAM model.
#' 
#' @param model
#' SIMCAM model (object of class \code{simcam})
#' 
simcam.getPerformanceStatistics = function(model) {
   nvar = nrow(model$models[[1]]$loadings)
   nc = length(model$models)
   
   dispower = array(0, dim = c(nc, nc, nvar))
   moddist = array(0, dim = c(nc, nc))
   
   # loop through all combinations of classes
   for (nc1 in 1:nc) {   
      for (nc2 in 1:nc) {     
         m1 = model$models[[nc1]]
         d1 = getCalibrationData(m1)
         m2 = model$models[[nc2]]
         d2 = getCalibrationData(m2)
         
         # apply model 1 to data 2 and vice versa
         m12 = predict.pca(m1, d2)
         m21 = predict.pca(m2, d1)   
         
         # calculate residuals for projections 
         if (m1$ncomp.selected < m1$ncomp) {   
            res1 = 
               m1$calres$scores[, (m1$ncomp.selected + 1):m1$ncomp] %*% 
               t(m1$loadings[, (m1$ncomp.selected + 1):m1$ncomp]) + 
               m1$calres$residuals
            res12 = 
               m12$scores[, (m1$ncomp.selected + 1):m1$ncomp] %*% 
               t(m1$loadings[, (m1$ncomp.selected + 1):m1$ncomp]) + 
               m12$residuals
         } else {
            res1 = m1$calres$residuals
            res12 = m12$residuals            
         }   
         
         if (m2$ncomp.selected < m2$ncomp) {   
            res2 = 
               m2$calres$scores[, (m2$ncomp.selected + 1):m2$ncomp] %*% 
               t(m2$loadings[, (m2$ncomp.selected + 1):m2$ncomp]) + 
               m2$calres$residuals
            
            res21 = 
               m21$scores[, (m2$ncomp.selected + 1):m2$ncomp] %*% 
               t(m2$loadings[, (m2$ncomp.selected + 1):m2$ncomp]) + 
               m21$residuals
         } else {
            res2 = m2$calres$residuals
            res21 = m21$residuals               
         }   
         
         # calculate standard deviations for the residuals
         s1 = colSums(res1^2)/(nrow(d1) - m1$ncomp.selected - 1)
         s2 = colSums(res2^2)/(nrow(d2) - m2$ncomp.selected - 1)
         s12 = colSums(res12^2)/(nrow(d2))
         s21 = colSums(res21^2)/(nrow(d1))
         
         # calculate model distance and discrimination power
         dispower[nc1, nc2, ] = sqrt((s12 + s21)/(s1 + s2))
         moddist[nc1, nc2] = sqrt(sum(s12 + s21)/sum(s1 + s2))
      }
   }
   
   dimnames(dispower) = list(model$classnames, model$classnames, rownames(model$models[[1]]$loadings))
   dimnames(moddist) = list(model$classnames, model$classnames)
  
   attrs = mda.getattr(model$models[[1]]$loadings)
   attr(dispower, 'name') = 'Discrimination power'
   attr(dispower, 'xaxis.name') = attrs$yaxis.name
   
   stat = list(
      dispower = dispower,
      moddist = moddist
   )
   
   stat
}

#' Modelling distance plot for SIMCAM model
#' 
#' @description
#' Shows a plot with distance from data objects to a SIMCA model
#' 
#' @param obj
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' for which class (SIMCA model) to show the plot for
#' @param type
#' type of the plot
#' @param main
#' main plot title
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
#' @export
plotModelDistance.simcam = function(obj, nc = 1, type = 'h', main = NULL, xlab = 'Models', ylab = '', ...) {
   if (is.null(main))
      main = sprintf('Model distance (%s)', obj$classnames[nc])
   
   mdaplot(mda.t(obj$moddist[, nc, drop = F]), type = type, main = main, xlab = xlab, ylab = ylab, ...)
}

#' Discrimination power plot for SIMCAM model
#' 
#' @description
#' Shows a plot with discrimination power of predictors for a pair of SIMCA models
#' 
#' @param obj
#' a SIMCAM model (object of class \code{simcam})
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
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
#' @export
plotDiscriminationPower.simcam = function(obj, nc = c(1, 2), type = 'h', main = NULL, 
                                          xlab = 'Variables', ylab = '', ...) {
   if (is.null(main))
      main = sprintf('Discrimination power (%s vs %s)', obj$classnames[nc[1]], obj$classnames[nc[2]])
   
   nvar = dim(obj$dispower)[[3]]
   
   if (is.null(type))
   {   
      if (nvar < 20)
         type = 'h'
      else
         type = 'l'
   }   
   
   attrs = mda.getattr(obj$dispower)
   data = obj$dispower[nc[1], nc[2], , drop = F]
   varnames = dimnames(obj$dispower)[[3]]
   dim(data) = c(dim(data)[3], 1)
   rownames(data) = varnames
   data = mda.setattr(data, attrs)
   
   mdaplot(data, type = type, main = main, xlab = xlab, ylab = ylab, ...)
}   

#' Cooman's plot for SIMCAM model
#' 
#' @description
#' Shows a Cooman's plot for a pair of SIMCA models
#' 
#' @param obj
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' vector with two values - classes (SIMCA models) to show the plot for
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
#' @export
plotCooman.simcam = function(obj, nc = c(1, 2), ...)
{
   plotCooman(obj$calres, nc = nc, ...)
}

#' Residuals plot for SIMCAM model
#' 
#' @description
#' Shows a plot with residuals for SIMCAM model
#' 
#' @param obj
#' a SIMCAM model (object of class \code{simcam})
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
#' @export
plotResiduals.simcam = function(obj, ...)
{
   plotResiduals(obj$calres, ...)
}

#' Modelling power plot for SIMCAM model
#' 
#' @description
#' Shows a plot with modelling power values for each predictor of selected SIMCA model
#' 
#' @param obj
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' which classe (SIMCA model) to show the plot for
#' @param main
#' main plot title
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
#' @export
plotModellingPower.simcam = function(obj, nc = 1, main = NULL, ...)
{
   if (is.null(main))
      main = sprintf('Modelling power (%s)', obj$classnames[nc])
   
   plotModellingPower.simca(obj$models[[nc]], main = main, ...)
}

#' Model overview plot for SIMCAM
#' 
#' @description
#' Shows a set of plots for SIMCAM model.
#' 
#' @param x
#' a SIMCAM model (object of class \code{simcam})
#' @param nc
#' vector with two values - classes (SIMCA models) to show the plot for
#' @param ...
#' other arguments
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
#' @export
plot.simcam = function(x, nc = c(1, 2), ...)
{
   obj = x
   
   par(mfrow = c(2, 1))
   plotDiscriminationPower(obj, nc)
   plotModelDistance(obj, nc[1])
   par(mfrow = c(1, 1))
}  

#' Summary method for SIMCAM model object
#' 
#' @description
#' Shows performance statistics for the model.
#' 
#' @param object
#' a SIMCAM model (object of class \code{simcam})
#' @param ...
#' other arguments
#' 
#' @export
summary.simcam = function(object, ...)
{
   obj = object
   
   cat('\nSIMCA multiple classes classification (class simcam)\n')
   cat(sprintf('Nmber of classes: %d\n', length(obj$models)))
   
   if (!is.null(obj$info))
      cat(sprintf('Info: %s\n', obj$info))
   
   for (i in 1:length(obj$models))
      summary(obj$models[[i]])
}  

#' Print method for SIMCAM model object
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a SIMCAM model (object of class \code{simcam})
#' @param ...
#' other arguments
#'
#' @export 
print.simcam = function(x, ...)
{
   obj = x
   
   cat('\nSIMCA multiple classes classification (class simcam)\n')
   
   cat('\nCall:\n')
   print(obj$call)
   
   cat('\nMajor fields:\n')   
   cat('$models - list wth individual SIMCA models for each class\n')
   cat('$classnames - vector with names of classes\n')
   cat('$moddist - matrix with distance between the models\n')
   cat('$dispower - matrix with discrimination power values\n')
   cat('$info - information about the object\n')
}  

