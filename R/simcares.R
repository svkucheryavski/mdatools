## class and methods for SIMCA classification results ##

simcares = function(pres, cres)
{
   # Creates an object of simcares class. 
   #
   # Arguments:
   #  pres: an object of pcares class (results for PCA decomposition)
   #  cres: an object of classres class (results for classification)
   
   res = c(pres, cres)
   res$classname = dimnames(cres$c.pred)[[3]][1]
   class(res) = c('simcares', 'classres', 'pcares', 'ldecomp')   
   
   res
}   

#' Residuals plot for SIMCA results
#' 
#' @description
#' Shows a plot with Q2 vs. T2 residuals for SIMCA results
#' 
#' @param obj
#' SIMCA results (object of class \code{simcares})
#' @param ncomp
#' which principal components to show the plot for
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
#' other plot parameters (see \code{mdaplot} for details)
#' 
#' @details
#' See examples in help for \code{\link{simcares}} function.
#' 
plotResiduals.simcares = function(obj, ncomp = NULL, show.limits = T, type = 'p', main = NULL, 
                                  xlab = 'T2', ylab = 'Q2', legend = NULL, ...)
{
   # Shows residuals plot (T2 vs Q2) 
   #
   # Arguments:
   #  obj: SIMCA results (an object of class simcares)
   #  ncomp: number of components to make the plot for
   #  show.limits: logical, show or not statistical limits for the residual values
   #  ...: standard arguments for plots

   # set main title
   if (is.null(main))
   {
      if (is.null(ncomp))
         main = 'Residuals'
      else
         main = sprintf('Residuals (ncomp = %d)', ncomp)
   }

   # get selected components and ncomp is NULL
   ncomp = getSelectedComponents(obj, ncomp)

   if (is.null(obj$c.ref)) 
   {
      # if no reference values or all objects are from the same class 
      # show standard plot for PCA
      plotResiduals.ldecomp(obj, ncomp, main = main, ...)
   }   
   else
   {   
      # if objects include members and non-members show plot with
      # color differentiation and legend
      
      c.ref = obj$c.ref     
      if (!is.character(c.ref))
      {   
         c.ref[obj$c.ref == 1] = obj$classname
         c.ref[obj$c.ref != 1] = 'Others'
      }
      
      if (sum(c.ref == obj$classname) == length(c.ref))
      {
         plotResiduals.ldecomp(obj, ncomp, main = main, ...)
      }
      else
      {
         classes = unique(c.ref)
         nclasses = length(classes)

         pdata = list()
         legend.str = NULL

         for (i in 1:nclasses)
         {   
            idx = c.ref == classes[i]
            data = cbind(obj$T2[idx, ncomp, drop = F], obj$Q2[idx, ncomp, drop = F])
            colnames(data) = c('T2', 'Q2')
            rownames(data) = rownames(obj$c.ref[idx])
           
            legend.str = c(legend.str, classes[i])
            pdata[[i]] = data
         }

         if (is.null(legend))
            legend = legend.str
         
            
         if (show.limits == T)
            show.lines = c(obj$T2lim[1, ncomp], obj$Q2lim[1, ncomp])
         else
            show.lines = F
         
         mdaplotg(pdata, type = type, xlab = xlab, ylab = ylab, main = main,
                  legend = legend, show.lines = show.lines, ...)
      }
   }   
}


#' Summary method for SIMCA results object
#' 
#' @method summary simcares
#' @S3method summary simcares
#'
#' @description
#' Shows performance statistics for the results.
#' 
#' @param object
#' SIMCA results (object of class \code{simcares})
#' @param ...
#' other arguments
#' 
summary.simcares = function(object, ...)
{
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
#' @method print simcares
#' @S3method print simcares
#' 
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' SIMCA results (object of class \code{simcares})
#' @param ...
#' other arguments
#' 
print.simcares = function(x, ...)
{   
   # Print information on simcares object
   obj = x
   
   cat('Result for SIMCA one-class classification (class simcares)\n')
   print.ldecomp(obj, '')
   print.classres(obj, '')
   cat('\n')
}

