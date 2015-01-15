randtest = function(x, y, ncomp = 15, center = T, scale = F, nperm = 1000, 
                    sig.level = 0.05, silent = TRUE)
{   
   # Randomization test for PLS regression
   # 
   # Carries out randomization/permutation test for a PLS regression model 
   # with given data and parameters
   #
   # Arguments
   # ---------
   # x - a matrix with x values (predictors)
   # y - a vector with y values (responses)
   # ncomp - maximum number of components to calculate
   # center - logical, do mean centering or not
   # scale - logical, do standardization or not
   #
   # Return
   # ------
   # res - an object with test results
   #
      
   x = as.matrix(x)
   y = as.matrix(y)
   nobj = nrow(x)
   
   x = prep.autoscale(as.matrix(x), center = center, scale = scale)
   y = prep.autoscale(as.matrix(y), center = center, scale = scale)
   
   stat = matrix(0, ncol = ncomp, nrow = 1)
   alpha = matrix(0, ncol = ncomp, nrow = 1)
   statperm = matrix(0, ncol = ncomp, nrow = nperm)
   corrperm = matrix(0, ncol = ncomp, nrow = nperm)
   
   for (icomp in 1:ncomp)
   {
      if ( !silent )
      {
         cat(sprintf('Permutations for component #%d...\n', icomp))
      }
      
      if (icomp > 1)
      {
         x = x - m$xscores %*% t(m$xloadings)
         y = y - m$xscores %*% t(m$yloadings)         
      }   
      
      m = pls.simpls(x, y, 1)      
      stat[icomp] = (t(m$xscores) %*% y) / nobj
            
      for (iperm in 1:nperm)
      {
         yp = y[sample(1:nobj)]
         mp = pls.simpls(x, yp, 1)      
         statperm[iperm, icomp] = (t(mp$xscores) %*% yp) / nobj      
         corrperm[iperm, icomp] = cor(y, yp)      
      }
      
      alpha[icomp] = sum(statperm[, icomp] > stat[icomp])/nperm
   }   

   ncomp.selected = max(which(alpha <= sig.level))
   colnames(alpha) = colnames(stat) = paste('Comp', 1:ncomp)
   colnames(statperm) = colnames(corrperm) = paste('Comp', 1:ncomp)   
   rownames(statperm) = rownames(corrperm) = 1:nperm
   rownames(alpha) = 'Alpha'
   rownames(stat) = 'Statistic'
   
   res = list(
      nperm = nperm,
      stat = stat,
      alpha = alpha,
      statperm = statperm,
      corrperm = corrperm,
      ncomp.selected = ncomp.selected
      )

   res$call = match.call()
   class(res) = "randtest"   
   
   res
}   
   

#' Histogram plot for randomization test results
#' 
#' @description
#' Makes a histogram for statistic values distribution for particular component, also
#' show critical value as a vertical line.
#' 
#' @param obj
#' results of randomization test (object of class `randtest`)
#' @param comp
#' number of component to make the plot for
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other optional arguments
#' 
#' @details
#' See examples in help for \code{\link{randtest}} function.
#' 
plotHist.randtest = function(obj, comp = NULL, main = NULL, xlab = 'Test statistic', ylab = 'Frequency', ...)    
{
   if (is.null(comp))
      comp = obj$ncomp.selected

   if (is.null(main))
      main = sprintf('Distribution for permutations (ncomp = %d)', comp)
   
   h = hist(obj$statperm[, comp], plot = F)
   
   dx = h$mids[2] - h$mids[1]
   sx = h$mids[1]

   stat = (obj$stat[comp] - sx) / dx
   
   xticks = 1:length(h$mids)
   data = cbind(xticks, h$counts)
   lim = mdaplot.getAxesLim(data, show.lines = c(stat, NA), xticks = xticks)
   mdaplot.plotAxes(xticks, h$mids, NULL, NULL, lim, main = main, xlab = xlab, ylab = ylab)
   mdaplot(data, type = 'h', show.lines = c(stat, NA), show.axes = F)
}

#' Correlation plot for randomization test results
#' 
#' @description
#' Makes a plot with statistic values vs. coefficient of determination between permuted 
#' and reference y-values.
#' 
#' @param obj
#' results of randomization test (object of class `randtest`)
#' @param comp
#' number of component to make the plot for
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other optional arguments
#' 
#' @details
#' See examples in help for \code{\link{randtest}} function.
#' 
plotCorr.randtest = function(obj, comp = NULL, main = NULL, xlab = expression(r^2), ylab = 'Test statistic', ...)
{
   if (is.null(comp))
      comp = obj$ncomp.selected

   if (is.null(main))
      main = sprintf('Permutations (ncomp = %d)', comp)
   
   data = list(cbind(obj$corrperm[, comp]^2, obj$statperm[, comp]), cbind(1, obj$stat[, comp]))   
   fitdata = rbind(apply(data[[1]], 2, mean), data[[2]])
   
   mdaplotg(data, type = 'p', main = main, xlab = xlab, ylab = ylab, ...)
   mdaplot.showRegressionLine(fitdata, col = rgb(0.6, 0.6, 0.6), lty = 2, lwd = 0.75)
}

#' Plot for randomization test results
#' 
#' @method plot randtest
#' @S3method plot randtest
#' 
#' @description
#' Makes a bar plot with alpha values for each component.
#' 
#' @param x
#' results of randomization test (object of class `randtest`)
#' @param main
#' main title for the plot
#' @param xlab
#' label for x axis
#' @param ylab
#' label for y axis
#' @param ...
#' other optional arguments
#' 
#' @details
#' See examples in help for \code{\link{randtest}} function.
#' 
plot.randtest = function(x, main = 'Alpha', xlab = 'Components', ylab = '', ...)
{
   obj = x
   mdaplot(t(rbind(1:length(obj$alpha), obj$alpha)), show.lines = c(NA, 0.05), type = 'h', 
           main = main, xlab = xlab, ylab = ylab, ...)      
}

#' Summary method for randtest object
#' 
#' @method summary randtest
#' @S3method summary randtest
#'
#' @description
#' Shows summary for randomization test results.
#' 
#' @param object
#' randomization test results (object of class \code{randtest})
#' @param ...
#' other arguments
#' 
summary.randtest = function(object, ...)
{
   obj = object
   data = rbind(obj$alpha, obj$stat)
   cat('Summary for permutation test results\n')
   cat(sprintf('Number of permutations: %d\n', obj$nperm))
   cat(sprintf('Suggested number of components: %d\n', obj$ncomp.selected))
   cat('\nStatistics and alpha values:\n')
   show(data)
}  

#' Print method for randtest object
#' 
#' @method print randtest
#' @S3method print randtest
#'
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a randomization test results (object of class \code{randtest})
#' @param ...
#' other arguments
#' 
print.randtest = function(x, ...)
{   
   obj = x
   
   cat('\nRandomization test results (class randtest)\n')
   cat('\nCall:\n')
   print(obj$call)
   cat('\nMajor fields:\n')
   cat('$nperm - number of permutations\n')
   cat('$ncomp.selected - number of selected components (suggested)\n')   
   cat('$alpha - vector with alpha values calculated for each component.\n')
   cat('$stat - vector with statistic values calculated for each component.\n')
   cat('$statperm - matrix with statistic values for each permutation.\n')
   cat('$corrperm - matrix with correlation between predicted and reference y-vales for each permutation.\n')
   cat('\nTry summary(obj) and plot(obj) to see the test results.\n')   
}

