## class and methods for SIMCA multi class classification results ##

simcamres = function(cres, T2, Q2, T2lim, Q2lim, ...)
{
   # Creates an object of simcamres class. 
   #
   # Arguments:
   #  cres: an object of classres class (results for classification)
   #  T2: T2 values for the objects and classes (selected component only)
   #  Q2: Q2 values for the objects and classes (selected component only)
   #  T2lim: T2 limits for the classes (selected component only)
   #  Q2lim: Q2 limits for the classes (selected component only)

   res = cres
   res$T2 = T2
   res$Q2 = Q2
   res$T2lim = T2lim
   res$Q2lim = Q2lim
   res$classnames = dimnames(cres$c.pred)[[3]]
   class(res) = c('simcamres', 'classres')   
   
   res
}   

plotResiduals.simcamres = function(obj, nc = 1, show.limits = T, type = 'p', main = NULL, 
                                  xlab = 'T2', ylab = 'Q2', legend = NULL, ...)
{
   # Shows residuals plot (T2 vs Q) for a selected model 
   #
   # Arguments:
   #  obj: SIMCA results (an object of class simcares)
   #  nc: which model to make the plot for
   #  show.limits: logical, show or not statistical limits for the residual values
   #  ...: standard arguments for plots

   # set main title
   if (is.null(main))
      main = sprintf('Residuals (%s)', obj$classnames[nc])

   if (!is.null(obj$c.ref))
   {   
      classes = unique(obj$c.ref)
      data = list()
      for (i in 1:length(classes))
         data[[i]] = cbind(obj$T2[obj$c.ref == classes[i], nc], obj$Q2[obj$c.ref == classes[i], nc])
      
      if (is.null(legend))
         legend = classes
      
      if (show.limits == T)
         show.lines = c(obj$T2lim[nc], obj$Q2lim[nc])
      else
         show.lines = F
      
      mdaplotg(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines,
               legend = legend, ...)
   }
   else
   {

      data = cbind(obj$T2[, nc], obj$Q2[, nc])
      if (show.limits == T)
         show.lines = c(obj$T2lim[nc], obj$Q2lim[nc])
      else
         show.lines = F
      
      mdaplot(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines, ...)
   }
}

plotCooman.simcamres = function(obj, nc = c(1, 2), type = 'p', main = "Cooman's plot", xlab = NULL, 
                                ylab = NULL, show.limits = T, legend = NULL, ...)
{

   # Shows Cooman's plot (distance to model 1 vs to model 2)
   #
   # Arguments:
   #  obj: SIMCA results (an object of class simcares)
   #  nc:  vector with two values - models to show the plot for
   #  show.limits: logical, show or not statistical limits for the residual values
   #  ...: standard arguments for plots

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
         data[[i]] = cbind(sqrt(obj$Q2[obj$c.ref == classes[i], nc[1]]), 
                           sqrt(obj$Q2[obj$c.ref == classes[i], nc[2]])
                           )
      
      if (is.null(legend))
         legend = classes
      
      if (show.limits == T)
         show.lines = c(obj$Q2lim[nc[1]], obj$Q2lim[nc[2]])
      else
         show.lines = F
      
      mdaplotg(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines,
               legend = legend, ...)
   }
   else
   {
      data = cbind(sqrt(obj$Q2[, nc[1]]), sqrt(obj$Q2[, nc[2]]))
      
      if (show.limits == T)
         show.lines = c(obj$Q2lim[nc[1]], obj$Q2lim[nc[2]])
      else
         show.lines = F
      
      mdaplot(data, type = type, xlab = xlab, ylab = ylab, main = main, show.lines = show.lines, ...)
   }
}

plot.simcamres = function(obj, ...)
{
   plotPredictions(obj)
}

print.simcamres = function(obj)
{   
   # Print information on simcares object
   
   cat('Result for SIMCA multiple classes classification (class simcamres)\n')
   print.classres(obj, '')
   cat('\n')
}

summary.simcamres = function(obj)
{
   # Show summary for simcamres object
   
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
