## class and methods for SIMCA classification with multiple classes (SIMCA models) ##

simcam = function(models, info = '')
{
   # Create a SIMCA multiple classes classification obj
   #
   # Arguments:
   #   models: a list with SIMCA one-class models (simca)
   #   info: information about the model
   #
   # Returns:
   #   model: a multiple classes SIMCA model (object of simcam class)

   nclasses = length(models)
   
   classnames = models[[1]]$classname
   for (i in 2:nclasses)
   {   
      classnames = c(classnames, models[[i]]$classname)
   }
   
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
#' a vector with reference class values
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
predict.simcam = function(object, x, c.ref = NULL, cv = F, ...)
{
   x = as.matrix(x)
   nobj = nrow(x)

   c.pred = array(0, dim = c(nobj, 1, object$nclasses))
   Q2 = array(0, dim = c(nobj, object$nclasses))
   T2 = array(0, dim = c(nobj, object$nclasses))
   Q2lim = array(0, dim = c(1, object$nclasses))
   T2lim = array(0, dim = c(1, object$nclasses))

   ncomp.selected = matrix(0, nrow = 1, ncol = object$nclasses)
   
   for (i in 1:object$nclasses)
   {  
      if (!is.null(c.ref))
      {   
         if (is.numeric(c.ref))
            res = predict(object$models[[i]], x, c.ref == i)
         else
            res = predict(object$models[[i]], x, c.ref == object$models[[i]]$classname)
      }
      else
         res = predict(object$models[[i]], x)

      ncomp.selected[i] = object$models[[i]]$ncomp.selected
      c.pred[, , i] = res$c.pred[, ncomp.selected[i], ]
      Q2[, i] = res$Q2[, ncomp.selected[i], drop = F]
      T2[, i] = res$T2[, ncomp.selected[i], drop = F]
      Q2lim[i] = res$Q2lim[ncomp.selected[i]]
      T2lim[i] = res$T2lim[ncomp.selected[i]]
   }
   
   dimnames(c.pred) = list(rownames(x), paste('Comp', ncomp.selected[[i]]), object$classnames)
   cres = classres(c.pred, c.ref, ncomp.selected = ncomp.selected)
   res = simcamres(cres, T2, Q2, T2lim, Q2lim)
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
getCalibrationData.simcam = function(obj, ...)
{
   x = NULL
   c.ref = NULL
   for (i in 1:obj$nclasses)
   {
      classdata = getCalibrationData(obj$models[[i]])
      x = rbind(x, classdata)
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
simcam.getPerformanceStatistics = function(model)
{
   nvar = nrow(model$models[[1]]$loadings)
   nc = length(model$models)
   
   dispower = array(0, dim = c(nc, nc, nvar))
   moddist = array(0, dim = c(nc, nc))

   # loop through all combinations of classes
   for (nc1 in 1:nc)
   {   
      for (nc2 in 1:nc)
      {     
         m1 = model$models[[nc1]]
         d1 = getCalibrationData(m1)
         m2 = model$models[[nc2]]
         d2 = getCalibrationData(m2)
         
         # apply model 1 to data 2 and vice versa
         m12 = predict.pca(m1, d2)
         m21 = predict.pca(m2, d1)   

         # calculate residuals for projections 
         if (m1$ncomp.selected < m1$ncomp)
         {   
            res1 = 
               m1$calres$scores[, (m1$ncomp.selected + 1):m1$ncomp] %*% 
               t(m1$loadings[, (m1$ncomp.selected + 1):m1$ncomp]) + 
               m1$calres$residuals
            res12 = 
               m12$scores[, (m1$ncomp.selected + 1):m1$ncomp] %*% 
               t(m1$loadings[, (m1$ncomp.selected + 1):m1$ncomp]) + 
               m12$residuals
         }
         else
         {
            res1 = m1$calres$residuals
            res12 = m12$residuals            
         }   

         if (m2$ncomp.selected < m2$ncomp)
         {   
            res2 = 
               m2$calres$scores[, (m2$ncomp.selected + 1):m2$ncomp] %*% 
               t(m2$loadings[, (m2$ncomp.selected + 1):m2$ncomp]) + 
               m2$calres$residuals
            
            
            res21 = 
               m21$scores[, (m2$ncomp.selected + 1):m2$ncomp] %*% 
               t(m2$loadings[, (m2$ncomp.selected + 1):m2$ncomp]) + 
               m21$residuals
         }
         else
         {
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
   
   stat = list(
      dispower = dispower,
      moddist = moddist
      )
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
#' @param xticklabels
#' labels for x axis ticks
#' @param ...
#' other plot parameters (see \code{mdaplotg} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcam}} function.
#' 
plotModelDistance.simcam = function(obj, nc = 1, type = 'h', main = NULL, xlab = 'Models', ylab = '', 
                                    xticklabels = NULL, ...)
{
   if (is.null(main))
      main = sprintf('Model distance (%s)', obj$classnames[nc])

   if (is.null(xticklabels) && length(obj$moddist[, nc]) < 8)
      xticklabels = obj$classnames

   data = cbind(1:length(obj$models), obj$moddist[, nc, drop = F])
   mdaplot(data, type = type, main = main, xlab = xlab, ylab = ylab, xticklabels = xticklabels, ...)
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
plotDiscriminationPower.simcam = function(obj, nc = c(1, 2), type = 'h', main = NULL, 
                                          xlab = 'Variables', ylab = '', ...)
{
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
   
   data = obj$dispower[nc[1], nc[2], , drop = F]
   varnames = dimnames(obj$dispower)[[3]]
   dim(data) = c(dim(data)[3], 1)
   data = cbind(1:nvar, data)
   rownames(data) = varnames
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
#' @method summary simcam
#' @S3method summary simcam
#'
#' @description
#' Shows performance statistics for the model.
#' 
#' @param object
#' a SIMCAM model (object of class \code{simcam})
#' @param ...
#' other arguments
#' 
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
#' @method print simcam
#' @S3method print simcam
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a SIMCAM model (object of class \code{simcam})
#' @param ...
#' other arguments
#' 
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

