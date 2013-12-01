## class and methods for SIMCA classification with multiple classes (SIMCA models) ##

simcam = function(models, info = '')
{
   # Create a SIMCA multiple classes classification model
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
   class(model) = c("simcam")
   
   # make predictions for calibration set
   caldata = getCalibrationData(model)
   model$calres = predict(model, caldata$data, caldata$c.ref)
   
   model
}

predict.simcam = function(model, data, c.ref = NULL, cv = F)
{
   # Apply the SIMCAM model to a new data set
   #
   # Arguments:
   #   model: a SIMCAM model (object of class simcam)
   #   data: a matrix with new data values
   #   c.ref: reference class values for the data set (optional)
   #
   # Returns:
   #   res: results of SIMCAM classification (object of simcamres class)
   
   data = as.matrix(data)
   nobj = nrow(data)

   c.pred = array(0, dim = c(nobj, 1, model$nclasses))
   Q2 = array(0, dim = c(nobj, model$nclasses))
   T2 = array(0, dim = c(nobj, model$nclasses))
   Q2lim = array(0, dim = c(1, model$nclasses))
   T2lim = array(0, dim = c(1, model$nclasses))

   ncomp.selected = matrix(0, nrow = 1, ncol = model$nclasses)
   
   for (i in 1:model$nclasses)
   {  
      if (!is.null(c.ref))
      {   
         if (is.numeric(c.ref))
            res = predict(model$models[[i]], data, c.ref == i)
         else
            res = predict(model$models[[i]], data, c.ref == model$models[[i]]$classname)
      }
      else
         res = predict(model$models[[i]], data)

      ncomp.selected[i] = model$models[[i]]$ncomp.selected
      c.pred[, , i] = res$c.pred[, ncomp.selected[i], ]
      Q2[, i] = res$Q2[, ncomp.selected[i], drop = F]
      T2[, i] = res$T2[, ncomp.selected[i], drop = F]
      Q2lim[i] = res$Q2lim[ncomp.selected[i]]
      T2lim[i] = res$T2lim[ncomp.selected[i]]
   }
   
   dimnames(c.pred) = list(rownames(data), paste('Comp', ncomp.selected[[i]]), model$classnames)
   cres = classres(c.pred, c.ref, ncomp.selected = ncomp.selected)
   res = simcamres(cres, T2, Q2, T2lim, Q2lim)
   res
}

getCalibrationData.simcam = function(model)
{
   data = NULL
   c.ref = NULL
   for (i in 1:model$nclasses)
   {
      classdata = getCalibrationData(model$models[[i]])
      data = rbind(data, classdata)
      c.ref = rbind(c.ref, matrix(model$models[[i]]$classname, nrow = nrow(classdata), ncol = 1))            
   }

   res = list(data = data,
              c.ref = c.ref
              )

   res
}

simcam.getPerformanceStatistics = function(model)
{
   # calculates discrimination power and distance between models
   #
   # Arguments:
   #  model: SIMCAM model (object of class simcam)
   
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

plotModelDistance.simcam = function(model, nc = 1, type = 'h', main = NULL, xlab = 'Models', ylab = '', 
                                    xticklabels = NULL, ...)
{
   # makes a plot with distance from one model to the others
   #
   # Argumens:
   #  model: a SIMCAM model (object of class simcam)
   #  nc: number of class to show the distance from
   #  ...: standard plot parameters


   if (is.null(main))
      main = sprintf('Model distance (%s)', model$classnames[nc])

   if (is.null(xticklabels) && length(model$moddist[, nc]) < 8)
      xticklabels = model$classnames

   data = cbind(1:length(model$models), model$moddist[, nc, drop = F])
   mdaplot(data, type = type, main = main, xlab = xlab, ylab = ylab, xticklabels = xticklabels, ...)
}

plotDiscriminationPower.simcam = function(model, nc = c(1, 2), type = 'h', xlab = 'Variables', 
                                          ylab = '', main = NULL, ...)
{
   # makes a plot with discrimination power of variables
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   nc: vector with two valies - models to show discrimination power for
   #   ...: stadard plot parameters
   
   if (is.null(main))
      main = sprintf('Discrimination power (%s vs %s)', model$classnames[nc[1]], model$classnames[nc[2]])

   nvar = dim(model$dispower)[[3]]
   
   if (is.null(type))
   {   
      if (nvar < 20)
         type = 'h'
      else
         type = 'l'
   }   
   
   data = model$dispower[nc[1], nc[2], , drop = F]
   varnames = dimnames(model$dispower)[[3]]
   dim(data) = c(dim(data)[3], 1)
   data = cbind(1:nvar, data)
   rownames(data) = varnames
   mdaplot(data, type = type, main = main, xlab = xlab, ylab = ylab, ...)
}   

plotCooman.simcam = function(model, ...)
{
   plotCooman(model$calres, ...)
}

plotResiduals.simcam = function(model, ...)
{
   plotResiduals(model$calres, ...)
}

plotModellingPower.simcam = function(model, nc = 1, main = NULL, ...)
{
   if (is.null(main))
      main = sprintf('Modelling power (%s)', model$classnames[nc])
   
   plotModellingPower(model$models[[nc]], main = main, ...)
}

plotPredictions.simcam = function(model, ...)
{
   plotPredictions(model$calres, ...)
}

plot.simcam = function(model, nc = c(1, 2), ...)
{
   # makes a plot overview of  a SIMCAM model
   #
   # Arguments:
   #  model: a SIMCAM model (object of class simcam)
   #  nc: a vector with two values - number of models to show the plot for

   par(mfrow = c(2, 1))
   plotDiscriminationPower(model, nc)
   plotModelDistance(model, nc[1])
   par(mfrow = c(1, 1))
}  

print.simcam = function(model, ...)
{
   cat('\nSIMCA multiple classes classification (class simcam)\n')
   
   cat('\nCall:\n')
   print(model$call)
   
   cat('\nMajor fields:\n')   
   cat('$models - list wth individual SIMCA models for each class\n')
   cat('$classnames - vector with names of classes\n')
   cat('$moddist - matrix with distance between the models\n')
   cat('$dispower - matrix with discrimination power values\n')
   cat('$info - information about the object\n')
}  

summary.simcam = function(model, ...)
{
   cat('\nSIMCA multiple classes classification (class simcam)\n')
   cat(sprintf('Nmber of classes: %d\n', length(model$models)))
   
   for (i in 1:length(model$models))
      summary(model$models[[i]])
}  
