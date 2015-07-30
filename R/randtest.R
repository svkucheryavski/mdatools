#' Randomization test for PLS regression
#' 
#' @description
#' \code{randtest} is used to carry out randomization/permutation test for a PLS regression model 
#' 
#' @param x
#' matrix with predictors.
#' @param y
#' vector or one-column matrix with response.
#' @param ncomp
#' maximum number of components to test.
#' @param center
#' logical, center or not predictors and response values.
#' @param scale
#' logical, scale (standardize) or not predictors and response values.
#' @param nperm
#' number of permutations.
#' @param sig.level
#' significance level.
#' @param silent
#' logical, show or not test progress.
#' 
#' @return 
#' Returns an object of \code{randtest} class with following fields:
#' \item{nperm }{number of permutations used for the test.} 
#' \item{stat }{statistic values calculated for each component.} 
#' \item{alpha }{alpha values calculated for each component.} 
#' \item{statperm }{matrix with statistic values for each permutation.} 
#' \item{corrperm }{matrix with correlation between predicted and reference y-vales for each 
#' permutation.} 
#' \item{ncomp.selected }{suggested number of components.} 
#' 
#' @details    
#' The class implements a method for selection of optimal number of components in PLS1 regression 
#' based on the randomization test [1]. The basic idea is that for each component from 1 to
#' \code{ncomp} a statistic T, which is a covariance between t-score (X score, derived from a PLS 
#' model) and the reference Y values, is calculated. By repeating this for randomly permuted 
#' Y-values a distribution of the statistic is obtained. A parameter \code{alpha} is computed to 
#' show how often the statistic T, calculated for permuted Y-values, is the same or higher than 
#' the same statistic, calculated for original data without permutations.
#' 
#' If a component is important, then the covariance for unpermuted data should be larger than the 
#' covariance for permuted data and therefore the value for \code{alpha} will be quie small (there 
#' is still a small chance to get similar covariance). This makes \code{alpha} very similar to 
#' p-value in a statistical test.
#' 
#' The \code{randtest} procedure calculates alpha for each component, the values can be observed 
#' using \code{summary} or \code{plot} functions. There are also several function, allowing e.g. 
#' to show distribution of statistics and the critical value for each component.
#'
#' @references 
#' S. Wiklund et al. Journal of Chemometrics 21 (2007) 427-439.
#'
#' @seealso 
#' Methods for \code{randtest} objects:
#' \tabular{ll}{
#'  \code{print.randtest} \tab prints information about a \code{randtest} object.\cr
#'  \code{\link{summary.randtest}} \tab shows summary statistics for the test.\cr
#'  \code{\link{plot.randtest}} \tab shows bar plot for alpha values.\cr
#'  \code{\link{plotHist.randtest}} \tab shows distribution of statistic plot.\cr
#'  \code{\link{plotCorr.randtest}} \tab shows determination coefficient plot.\cr
#' }
#' 
#' @examples
#' ### Examples of using the test
#' 
#' ## Get the spectral data from Simdata set and apply SNV transformation
#' 
#' data(simdata)
#' 
#' y = simdata$conc.c[, 3]
#' x = simdata$spectra.c
#' x = prep.snv(x)
#' 
#' ## Run the test and show summary 
#' ## (normally use higher nperm values > 1000)
#' r = randtest(x, y, ncomp = 4, nperm = 200, silent = FALSE)
#' summary(r)
#' 
#' ## Show plots
#' 
#' par( mfrow = c(3, 2))
#' plot(r)
#' plotHist(r, comp = 3)
#' plotHist(r, comp = 4)
#' plotCorr(r, 3)
#' plotCorr(r, 4)
#' par( mfrow = c(1, 1))
#' 
#' @export
randtest = function(x, y, ncomp = 15, center = T, scale = F, nperm = 1000, 
                    sig.level = 0.05, silent = TRUE)
{   
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
#' @export
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
#' @export
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
#' @export
plot.randtest = function(x, main = 'Alpha', xlab = 'Components', ylab = '', ...)
{
   obj = x
   mdaplot(t(rbind(1:length(obj$alpha), obj$alpha)), show.lines = c(NA, 0.05), type = 'h', 
           main = main, xlab = xlab, ylab = ylab, ...)      
}

#' Summary method for randtest object
#' 
#' @description
#' Shows summary for randomization test results.
#' 
#' @param object
#' randomization test results (object of class \code{randtest})
#' @param ...
#' other arguments
#'
#' @export 
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
#' @description
#' Prints information about the object structure
#' 
#' @param x
#' a randomization test results (object of class \code{randtest})
#' @param ...
#' other arguments
#'
#' @export 
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

