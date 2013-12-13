## class and methods for classification models ##

plotPredictions.classmodel = function(model, res = NULL, nc = NULL, ncomp = NULL, main = NULL, ...)
{
   # makes a plot with predictions
   #
   # Arguments:
   #   model: a SIMCA model (object of class simca)
   #   which: which results to show (not used so far)
   #   ncomp: number of components to show the predictions for (default - selected)
   #   type: plot type
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   legend: legend strings
   #
   
   if (is.null(res))
   {
      if (!is.null(model$testres))
         res = 'testres'
      else if (!is.null(model$cvres))
         res = 'cvres'
      else
         res = 'calres'
   }

   resnames = c('calres', 'cvres', 'testres')
   resstrings = c('calbiration', 'cross-validation', 'test set')
   resobj = model[[res]]

   if (!(res %in% resnames))
      stop('Wrong value for argument "res" (use "calres", "cvres" or "testres")!')
   else if (is.null(resobj))
      stop('The result does not exist!')
   else
   {
      if (is.null(main))
      {
         if (is.null(ncomp))
            main = sprintf('Predictions for %s results', resstrings[resnames == res])
         else
            main = sprintf('Predictions for %s results (ncomp = %d)', resstrings[resnames == res], ncomp)
      }
      plotPredictions.classres(resobj, nc = nc, ncomp = ncomp, main = main, ...)
   }
}

plotSpecificity.simca = function(model, ...)
{
   # makes a plot with specificity values vs. number of components
   #
   # Arguments:
   #   model: a classification model (object of class plsda, simca or similar)
   #
   
   plotPerformance(model, param = 'specificity', ...)
}

plotSensitivity.classmodel = function(model, ...)
{
   # makes a plot with sensitivity values vs. number of components 
   #
   # Arguments:
   #   model: a classification model (object of class plsda, simca or similar)
   #

   plotPerformance(model, param = 'sensitivity', ...)
}


plotMisclassified.classmodel = function(model, ...)
{
   # makes a plot with misclassified samples vs. number of components 
   #
   # Arguments:
   #   model: a classification model (object of class plsda, simca or similar)
   #

   plotPerformance(model, param = 'misclassified', ...)
}


plotPerformance.classmodel = function(model, nc = NULL, param = 'specificity', type = 'h', legend = NULL, 
                                 main = NULL, xlab = 'Components', ylab = '', 
                                 ylim = c(0, 1.15), ...)
{
   # makes a plot with classification performance parameters vs. number of components 
   #
   # Arguments:
   #   model: a classification model (object of class plsda, simca or similar)
   #   param: which parameter to make the plot for ('specificity', 'sensitivity', 'misclassified')
   #   type: plot type
   #   legend: legend strings
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   ylim: limits for y axis
   #
  
   if (is.null(nc))
   {
      nc =  model$nclasses + 1
      
      if (model$nclasses > 1)
         classname = '(all classes)'
      else
         classname = ''
   }
   else
   {
      if (nc > model$nclasses || nc < 1)
         stop('Wrong value for argument "nc"!')

      classname = sprintf('(%s)', model$classnames[nc])
   }

   data = cbind(1:model$ncomp, model$calres[[param]][nc, ])
   labels = matrix(mdaplot.formatValues(model$calres[[param]][nc, ]), ncol = 1)
   legend_str = 'cal'
   
   if (!is.null(model$cvres))
   {
      data = cbind(data, model$cvres[[param]][nc, ])   
      labels = cbind(labels, mdaplot.formatValues(model$cvres[[param]][nc, ]))
      legend_str = c(legend_str, 'cv')
   }   
   
   if (!is.null(model$testres))
   {
      data = cbind(data, model$testres[[param]][nc, ])   
      labels = cbind(labels, mdaplot.formatValues(model$testres[[param]][nc, ]))
      legend_str = c(legend_str, 'test')
   }
  
   if (is.null(main))
      main = sprintf('%s%s %s', toupper(substring(param, 1, 1)), substring(param, 2, ), toString(classname))

   if (!is.null(legend))
      legend_str = legend
   
   mdaplotg(data, type = type, main = main, xlab = xlab, ylab = ylab, legend = legend_str,
            ylim = ylim, labels = labels, ...)
}

