## class and methods for iPLS ##

#' Variable selection with interval PLS
#' 
#' @description 
#' Applies iPLS alrogithm to find variable intervals most important for
#' prediction
#' 
#' @param x
#' a matrix with predictor values
#' @param y
#' a vector with response values
#' @param glob.ncomp
#' maximum number of components for a global PLS model
#' @param center
#' logical, center or not the data values
#' @param scale
#' logical, standardize or not the data values
#' @param cv
#' number of segments for cross-validation (1 - full CV)
#' @param exclcols
#' columns of x to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param int.ncomp
#' maximum number of components for interval PLS models
#' @param int.num
#' number of intervals
#' @param int.width
#' width of intervals
#' @param int.limits
#' a two column matrix with manual intervals specification
#' @param int.niter
#' maximum number of iterations (if NULL it will be the same as number of intervals)
#' @param ncomp.selcrit
#' criterion for selecting optimal number of components ('min' for minimum of RMSECV)
#' @param method
#' iPLS method (\code{'forward'} or \code{'backward'})
#' @param silent
#' logical, show or not information about selection process
#'  
#' @return 
#' object of 'ipls' class with several fields, including:
#'    \item{var.selected}{a vector with indices of selected variables}
#'    \item{int.selected}{a vector with indices of selected intervals }
#'    \item{int.num}{total number of intervals}
#'    \item{int.width}{width of the intervals}
#'    \item{int.limits}{a matrix with limits for each interval}
#'    \item{int.stat}{a data frame with statistics for the selection algorithm}
#'    \item{glob.stat}{a data frame with statistics for the first step (individual intervals)}
#'    \item{gm}{global PLS model with all variables included}
#'    \item{om}{optimized PLS model with selected variables}
#'      
#' @details 
#' The algorithm splits the predictors into several intervals and tries to find a combination 
#' of the intervals, which gives best prediction performance. There are two selection methods:
#' "forward" when the intervals are successively included, and "backward" when the intervals
#' are successively excluded from a model. On the first step the algorithm finds the best
#' (forward) or the worst (backward) individual interval. Then it tests the others to find the 
#' one which gives the best model in a combination with the already selected/excluded one. The 
#' procedure continues until the maximum number of iteration is reached. 
#'  
#' There are several ways to specify the intervals. First of all either number of intervals 
#' (\code{int.num}) or width of the intervals (\code{int.width}) can be provided. Alternatively
#' one can specify the limits (first and last variable number) of the intervals manually 
#' with \code{int.limits}.
#' 
#' @references 
#' [1] Lars Noergaard at al.  Interval partial least-squares regression (iPLS): a
#' comparative chemometric study with an example from near-infrared spectroscopy.
#' Appl.Spec. 2000; 54: 413-419
#'  
#' @examples 
#' library(mdatools)
#'
#' ## forward selection for simdata
#'  
#' data(simdata)
#' Xc = simdata$spectra.c
#' yc = simdata$conc.c[, 3, drop = FALSE]
#'
#' # run iPLS and show results  
#' im = ipls(Xc, yc, int.ncomp = 5, int.num = 10, cv = 4, method = "forward")
#' summary(im)
#' plot(im)
#'  
#' # show "developing" of RMSECV during the algorithm execution
#' plotRMSE(im)
#'  
#' # plot predictions before and after selection
#' par(mfrow = c(1, 2))
#' plotPredictions(im$gm)
#' plotPredictions(im$om)
#'  
#' # show selected intervals on spectral plot
#' ind = im$var.selected
#' mspectrum = apply(Xc, 2, mean)
#' plot(simdata$wavelength, mspectrum, type = 'l', col = 'lightblue')
#' points(simdata$wavelength[ind], mspectrum[ind], pch = 16, col = 'blue')
#'  
#' @export
ipls = function(x, y, glob.ncomp = 10, center = T, scale = F, cv = 10,
                exclcols = NULL, exclrows = NULL,  int.ncomp = glob.ncomp, int.num = NULL, 
                int.width = NULL, int.limits = NULL, int.niter = NULL, ncomp.selcrit = 'min', 
                method = 'forward', silent = F) {
   
   # process names and values for xaxis
   xaxis.name = attr(x, 'xaxis.name')
   if (is.null(xaxis.name))
      xaxis.name = 'Variables'

   xaxis.values = attr(x, 'xaxis.values')
   if (is.null(xaxis.values))
      xaxis.values = 1:ncol(x)

   x = as.matrix(x)
   
   if (is.null(dim(y))) {
      y = matrix(y, ncol = 1)       
      rownames(y) = rownames(x)
   }    
  
   # remove excluded columns and rows   
   if (length(exclcols) > 0) {
      x = mda.exclcols(x, exclcols)
      exclcols = attr(x, 'exclcols')
      x = x[, -exclcols, drop = F]
      xaxis.values = xaxis.values[-exclcols]
      attr(x, 'exclcols') = NULL
   }
   
   if (length(exclrows) > 0) {
      x = mda.exclrows(x, exclrows)
      exclrows = attr(x, 'exclrows')
      x = x[-exclrows, , drop = F]
      y = y[-exclrows, , drop = F]
      attr(x, 'exclrows') = NULL
      attr(y, 'exclrows') = NULL
   }
   
   attr(x, 'xaxis.values') = xaxis.values
   attr(x, 'xaxis.name') = xaxis.name
   
   if (length(dim(y)) != 2) {
      stop('Response variable (y) should be a matrix or a sequence!')
   }   
   
   if (ncol(y) > 1) {
      warning('iPLS can work with one y-variable at time, selecting first column.')
      y = y[, 1, drop = F]
   }
   
   if (is.null(colnames(y)))
      colnames(y) = paste('y', 1:ncol(y), sep = '')
  
   nx = ncol(x)
   
   # check and define width and number of intervals
   if (is.null(int.limits)) {   
      if (is.null(int.width)) {   
         if (is.null(int.num))
            stop('Specify either interval width or number of intervals!')
   
         if (int.num < 2 || int.num > nx)
            stop('Wrong value for number of intervals')
      
         int.width = nx / int.num
      } else {
         if (int.width < 1 || int.width > nx)
            stop('Wrong value for interval width!')

         int.num = round(nx / int.width)      
      }   
   
      # generate interval limits similar to the way from [1]
      if (int.num == nx) {
         int1 = 1:nx
         int2 = 1:nx
      } else {   
         varrem = nx %% int.num
         nvarrem = trunc(nx/int.num)
         if (varrem == 0) {
            int1 = seq(1, nx, by = nvarrem)
         } else {   
            int1 = c(
               seq(1, (varrem - 1) * (nvarrem + 1) + 1, by = nvarrem + 1),
               seq((varrem - 1) * (nvarrem + 1) + 2 + nvarrem, nx, by = nvarrem)
            )
         }
         int2 = c(int1[2:int.num] - 1, nx);
      }
      int.limits = cbind(int1, int2);
   } else {
      if (is.null(dim(int.limits)) || ncol(int.limits) != 2 || nrow(int.limits) < 2)
         stop('Interval limits shall be provided as a matrix with two columns!')
      
      df = int.limits[, 2] - int.limits[, 1]    
      if (min(int.limits) < 1 || max(int.limits) > nx || any(df < 0))
         stop('Wrong values for interval limits!')
      
      int.num = nrow(int.limits)
      int.width = mean(df)
   }   
   rownames(int.limits) = 1:int.num
   colnames(int.limits) = c('Left', 'Right')

   ## uncomment code below is for checking how intervals are generated
   #show(int.num)
   #show(int.width)
   #show(int1)
   #show(int2)
   #show(cbind(int.limits, int.limits[, 2] - int.limits[, 1] + 1))
   #stop()
   
  # define number of iterations
   if (is.null(int.niter)) {
      if (int.num < 30)
         int.niter = int.num
      else
         int.niter = 30
   }   
   
   # correct maximum number of components for interval models
   if (!is.numeric(cv))
      nseg = cv[[2]]
   else
      nseg = cv
   
   if (nseg == 1)
      nobj.cv = 1
   else
      nobj.cv = ceiling(nrow(x)/nseg)  
   
   int.ncomp = min(round(int.width), nrow(x) - 1 - nobj.cv, int.ncomp)   
   
   # build an object
   obj = list(
      cv = cv,
      glob.ncomp = glob.ncomp,
      int.ncomp = int.ncomp,
      xmean = apply(x, 2, "mean"),
      int.width = int.width,
      int.num = int.num,
      int.niter = int.niter,
      int.limits = int.limits,
      center = center,
      scale = scale,
      method = method,
      ncomp.selcrit = ncomp.selcrit,
      silent = silent,
      xaxis.name = xaxis.name,
      xaxis.values = xaxis.values
   )
   
   # make a global model
   obj$gm = pls(x, y, ncomp = obj$glob.ncomp, center = obj$center, scale = obj$scale, 
                cv = obj$cv, light = TRUE, ncomp.selcrit = obj$ncomp.selcrit)
   
   if (method == 'forward')
      obj = ipls.forward(x, y, obj)
   else if (method == 'backward')
      obj = ipls.backward(x, y, obj)
   else
      stop('Wrong value for parameter "method"!')
   
   if (!is.null(obj$int.stat)) {
      rownames(obj$int.stat) = 1:nrow(obj$int.stat)
      obj$int.stat$R2 = round(obj$int.stat$R2, 3)
   }   
   
   if (!is.null(obj$glob.stat)) {
      rownames(obj$glob.stat) = 1:nrow(obj$glob.stat)
      obj$glob.stat$R2 = round(obj$glob.stat$R2, 3)
   }   
   
   obj$var.selected = NULL
   for (i in obj$int.selected)
      obj$var.selected = c(obj$var.selected, seq(obj$int.limits[i, 1], obj$int.limits[i, 2]))

   obj$om = pls(x[, obj$var.selected], y, ncomp = obj$int.ncomp, center = obj$center, 
                scale = obj$scale, cv = obj$cv, light = TRUE, ncomp.selcrit = obj$ncomp.selcrit)
   
   obj$call = match.call()   
   class(obj) = "ipls"
   
   obj
}

#' Runs the forward iPLS algorithm
#' 
#' @param x
#' a matrix with predictor values
#' @param y
#' a vector with response values
#' @param obj
#' object with initial settings for iPLS algorithm
#'  
ipls.forward = function(x, y, obj) {
   # define vectors with status, selected and non-selected intervals
   int.nonselected = 1:obj$int.num
   int.selected = NULL
   
   glob.stat = data.frame(
      'n' = 0,
      'start' = 1,
      'end' = ncol(x),
      'nComp' = obj$gm$ncomp.selected,
      'RMSE' = obj$gm$cvres$rmse[1, obj$gm$ncomp.selected],
      'R2' = obj$gm$cvres$r2[1, obj$gm$ncomp.selected]
   )
   
   int.stat = data.frame(
      'n' = 0,
      'start' = 1,
      'end' = ncol(x),
      'selected' = F,
      'nComp' = obj$gm$ncomp.selected,
      'RMSE' = obj$gm$cvres$rmse[1, obj$gm$ncomp.selected],
      'R2' = obj$gm$cvres$r2[1, obj$gm$ncomp.selected]
   )
 
   if (!obj$silent)
      cat(sprintf('Model with all intervals: RMSECV = %f, nLV = %d\n', 
                  int.stat$RMSE[1], int.stat$nComp[1]))
   
   # do loop for max number of intervals
   selind = NULL 
   rmse = Inf
   for (i in 1:obj$int.niter) {
      if (!obj$silent)
         cat(sprintf('Iteration %3d/%3d... ', i, obj$int.niter))
      
      sel = NULL
      for (l in int.nonselected) {
         # combine already selected intervals with the current
         ind = obj$int.limits[l, 1]:obj$int.limits[l, 2]
         Xc = x[, c(selind, ind), drop = F]
         
         # build a model
         m = pls(Xc, y, ncomp = obj$int.ncomp, center = obj$center, scale = obj$scale, 
                 cv = obj$cv, light = TRUE, ncomp.selcrit = obj$ncomp.selcrit)
         
         # if first round, build a data frame with statistics for each interval
         if (i == 1) { 
            glob.stat = rbind(glob.stat, 
                              data.frame(
                                 'n' = l,
                                 'start' = obj$int.limits[l, 1],
                                 'end' = obj$int.limits[l, 2],
                                 'nComp' = m$ncomp.selected,
                                 'RMSE' = m$cvres$rmse[1, m$ncomp.selected],
                                 'R2' = m$cvres$r2[1, m$ncomp.selected]
                              )
            )      
         }
         
         ## monitoring of the algorithm
         #cat(sprintf('%4d %4d %4d %8.3f %8.3f %4d %4d %4d %4d\n',
         #            i, l, m$ncomp.selected, rmse, m$cvres$rmse[1, m$ncomp.selected], 
         #            dim(Xc)[1], dim(Xc)[2],
         #            obj$int.limits[l, 1],obj$int.limits[l, 2]))
         
         # else check if rmse has been improved
         if (rmse > m$cvres$rmse[1, m$ncomp.selected]) {
            ncomp = m$ncomp.selected
            rmse = m$cvres$rmse[1, m$ncomp.selected]
            r2 = m$cvres$r2[1, m$ncomp.selected]
            sel = l
         }
      }
      
      
      if (!is.null(sel)) {
         selind = c(selind, obj$int.limits[sel, 1]:obj$int.limits[sel, 2])
         int.nonselected = int.nonselected[int.nonselected != sel]
         int.selected = c(int.selected, sel)
         
         int.stat = rbind(int.stat, data.frame(
            'n' = sel,
            'start' = obj$int.limits[sel, 1],
            'end' = obj$int.limits[sel, 2],
            'selected' = T,
            'nComp' = ncomp,
            'RMSE' = rmse,
            'R2' = r2
         )) 
         
         if (!obj$silent)
            cat(sprintf('selected interval %3d (RMSECV = %f)\n', sel, rmse))
      } else {   
         # no improvements, quit the outer loop
         if (!obj$silent)
            cat('no improvements, stop.\n')
         break
      }
   }   
   
   # return the selection results
   obj$glob.stat = glob.stat
   obj$int.stat = int.stat
   obj$int.selected = int.selected
   
   obj
}

#' Runs the backward iPLS algorithm
#' 
#' @param x
#' a matrix with predictor values
#' @param y
#' a vector with response values
#' @param obj
#' object with initial settings for iPLS algorithm
#'  
ipls.backward = function(x, y, obj) {
   # define vectors with status, selected and non-selected intervals
   int.selected = 1:obj$int.num
   int.nonselected = NULL
  
   glob.stat = data.frame(
      'n' = 0,
      'start' = 1,
      'end' = ncol(x),
      'nComp' = obj$gm$ncomp.selected,
      'RMSE' = obj$gm$cvres$rmse[1, obj$gm$ncomp.selected],
      'R2' = obj$gm$cvres$r2[1, obj$gm$ncomp.selected]
   )
   
   int.stat = data.frame(
      'n' = 0,
      'start' = 1,
      'end' = ncol(x),
      'selected' = F,
      'nComp' = obj$gm$ncomp.selected,
      'RMSE' = obj$gm$cvres$rmse[1, obj$gm$ncomp.selected],
      'R2' = obj$gm$cvres$r2[1, obj$gm$ncomp.selected]
   )
   
   if (!obj$silent)
      cat(sprintf('Model with all intervals: RMSECV = %f, nLV = %d\n', 
                  int.stat$RMSE[1], int.stat$nComp[1]))
   
   # do loop for max number of intervals
   unselind = NULL 
   rmse = Inf
   for (i in 1:obj$int.niter) {
      if (length(int.selected) == 1)
         break
     
      if (!obj$silent)
         cat(sprintf('Iteration %3d/%3d... ', i, obj$int.niter))
      
      # do loop to select an interval
      unsel = NULL
      for (l in int.selected) {
         # combine already selected intervals with the current
         ind = obj$int.limits[l, 1]:obj$int.limits[l, 2]
         Xc = x[, -c(unselind, ind), drop = F]
         
         # build a model
         m = pls(Xc, y, ncomp = obj$int.ncomp, center = obj$center, scale = obj$scale, 
                 cv = obj$cv, light = TRUE, ncomp.selcrit = obj$ncomp.selcrit)
         
         # if first round, build a data frame with statistics for each interval
         if (i == 1) {
            glob.stat = rbind(glob.stat, 
                             data.frame(
                                'n' = l,
                                'start' = obj$int.limits[l, 1],
                                'end' = obj$int.limits[l, 2],
                                'nComp' = m$ncomp.selected,
                                'RMSE' = m$cvres$rmse[1, m$ncomp.selected],
                                'R2' = m$cvres$r2[1, m$ncomp.selected]
                             )
            )      
         }
         
         ## monitoring of the algorithm
         #cat(sprintf('%4d %4d %4d %8.3f %8.3f %4d %4d %4d %4d\n',
         #            i, l, m$ncomp.selected, rmse, m$cvres$rmse[1, m$ncomp.selected], 
         #            dim(Xc)[1], dim(Xc)[2],
         #            obj$int.limits[l, 1],obj$int.limits[l, 2]))
        
         # if last two intervals are left keep them both
         if (length(int.selected) == 2) {
            int.stat = rbind(int.stat, data.frame(
               'n' = l,
               'start' = obj$int.limits[l, 1],
               'end' = obj$int.limits[l, 2],
               'selected' = F,
               'nComp' = m$ncomp.selected,
               'RMSE' = m$cvres$rmse[1, m$ncomp.selected],
               'R2' = m$cvres$r2[1, m$ncomp.selected]
            ))
            unsel = NULL           
         }  else if (rmse > m$cvres$rmse[1, m$ncomp.selected]) {
            # else check if rmse has been improved
            ncomp = m$ncomp.selected
            rmse = m$cvres$rmse[1, m$ncomp.selected]
            r2 = m$cvres$r2[1, m$ncomp.selected]
            unsel = l
         }
      }
      
      if (!is.null(unsel)) {
         unselind = c(unselind, obj$int.limits[unsel, 1]:obj$int.limits[unsel, 2])
         int.selected = int.selected[int.selected != unsel]
         int.nonselected = c(int.nonselected, unsel)
         
         int.stat = rbind(int.stat, data.frame(
            'n' = unsel,
            'start' = obj$int.limits[unsel, 1],
            'end' = obj$int.limits[unsel, 2],
            'selected' = F,
            'nComp' = ncomp,
            'RMSE' = rmse,
            'R2' = r2
         ))
         
         if (!obj$silent)
            cat(sprintf('excluded interval %3d (RMSECV = %f)\n', unsel, rmse))
         
      } else {   
         # no improvements, quit the outer loop
         if (!obj$silent)
            cat('no improvements, stop.\n')
         break
      }
   }   
   
   # return the selection results
   obj$glob.stat = glob.stat
   obj$int.stat = int.stat
   obj$int.selected = int.selected
   
   obj
}

#' iPLS performance plot
#'
#' @description 
#' Shows PLS performance for each selected or excluded intervals at the
#' first iteration 
#' 
#' @param obj 
#' iPLS results (object of class ipls)
#' @param glob.ncomp
#' number of components for global PLS model with all intervals
#' @param main 
#' main title for the plot
#' @param xlab 
#' label for x-axis
#' @param ylab 
#' label for y-axis
#' @param xlim 
#' limits for x-axis
#' @param ylim 
#' limits for y-axis
#' @param ... 
#' other arguments
#'
#' @details
#' The plot shows intervals as bars, which height corresponds to RMSECV obtained when particular
#' interval was selected (forward) or excluded (backward) from a model at the first iteration. 
#' The intervals found optimal after backward/forward iPLS selection are shown with green color 
#' while the other intervals are gray.
#' 
#' See examples in help for \code{\link{ipls}} function.
#'
#'  @seealso 
#'  \code{\link{summary.ipls}}, \code{\link{plotRMSE.ipls}}
#'    
plotSelection.ipls = function(obj, glob.ncomp = NULL, main = 'iPLS results', 
                         xlab = obj$xaxis.name, ylab = 'RMSECV', xlim = NULL, ylim = NULL, ...) {
   int = obj$int.limits
   xlabels = obj$xaxis.values
   if (!is.null(xlabels)) {
      if (!is.numeric(xlabels))
         stop('Parameter "xlabels" should be a vector with numbers!')
      
      if (length(xlabels) != (max(int) - min(int) + 1)) {
         stop('Wrong values for "xlabels" parameter!')   
      } else {
         int[, 1] = xlabels[int[, 1]]
         int[, 2] = xlabels[int[, 2]]
      }   
   } else {
       xlabels = 1:max(int)
   }   
   
   bwd = (int[, 2] - int[, 1] + 1)
   mids = (int[, 2] + int[, 1])/2
   rmse = obj$glob.stat$RMSE[2:nrow(obj$glob.stat)]
   ncomp = obj$glob.stat$nComp[2:nrow(obj$glob.stat)]
   
   if (is.null(xlim))
      xlim = c(min(int), max(int))
   if (is.null(ylim))
      ylim = c(0, max(rmse) * 1.1)
  
   # rescale mean X values to fit the plot
   xmean = (obj$xmean - min(obj$xmean)) / (max(obj$xmean) - min(obj$xmean)) * (ylim[2] - ylim[1])

   # make plot
   plot(0, 0, type = 'n', main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
   
   # gray and green bars
   bars(mids, rmse, col = rgb(0.9, 0.9, 0.9), bwd = bwd, border = rgb(0.8, 0.8, 0.8)) 
   bars(mids[obj$int.selected], rmse[obj$int.selected], col = rgb(0.5, 1.0, 0.6), 
      bwd = bwd[obj$int.selected], border = rgb(0.4, 0.9, 0.5))

   # mean signal
   lines(xlabels, xmean, col = rgb(1.0, 0.7, 0.7), lwd = 2)

   # number of components for each interval
   text(mids,  matrix(0.05 * ylim[2], ncol = length(mids)), ncomp, 
        col = rgb(0.4, 0.4, 0.4), cex = 0.85) 

   # global model   
   if (is.null(glob.ncomp))
      glob.ncomp = obj$gm$ncomp.selected
   else if (glob.ncomp < 1 || glob.ncomp > obj$gm$ncomp)
      stop('Wrong value for number of components!')

   dx = (xlim[2] - xlim[1])/50
   abline(h = obj$gm$cvres$rmse[1, glob.ncomp], lty = 2, col = rgb(0.5, 0.5, 0.5))
   text(xlim[2] + dx, obj$gm$cvres$rmse[1, glob.ncomp], glob.ncomp, cex = 0.85, 
        col = rgb(0.3, 0.3, 0.3), font = 2, pos = 3)
}

#' RMSE development plot
#'
#' @description 
#' Shows how RMSE develops for each iteration of iPLS selection 
#' algorithm
#' 
#' @param obj 
#' iPLS results (object of class ipls)
#' @param glob.ncomp
#' number of components for global PLS model with all intervals
#' @param main 
#' main title for the plot
#' @param xlab 
#' label for x-axis
#' @param ylab 
#' label for y-axis
#' @param xlim 
#' limits for x-axis
#' @param ylim 
#' limits for y-axis
#' @param ... 
#' other arguments
#'
#' @details
#' The plot shows RMSE values obtained at each iteration of the iPLS algorithm as bars. The first 
#' bar correspond to the global model with all variables included, second - to the model obtained 
#' at the first iteration and so on. Number at the bottom of each bar corresponds to the interval 
#' included or excluded at the particular iteration. 
#' 
#' @seealso 
#' \code{\link{summary.ipls}}, \code{\link{plotSelection.ipls}}
#' 
#' @export      
plotRMSE.ipls = function(obj, glob.ncomp = NULL, main = 'RMSE development', xlab = 'Iterations', 
                    ylab = 'RMSECV', xlim = NULL, ylim = NULL, ...) {

   if (is.null(glob.ncomp))
      glob.ncomp = obj$gm$ncomp.selected
   else if (glob.ncomp < 1 || glob.ncomp > obj$gm$ncomp)
      stop('Wrong value for number of components!')
   
   rmse = obj$int.stat$RMSE
   n = obj$int.stat$n
   mids = 0:(length(n) - 1)

   if (is.null(xlim))
      xlim = c(min(mids) - 0.5, max(mids) + 0.5)
   if (is.null(ylim))
      ylim = c(0, max(rmse) * 1.1)

   # make plot
   plot(0, 0, type = 'n', main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)

   # gray and green bars
   bars(mids, rmse, col = rgb(0.5, 1.0, 0.6), bwd = 1, border = rgb(0.4, 0.9, 0.5)) 
   bars(mids[1], rmse[1], col = rgb(0.98, 0.98, 0.98), bwd = 1, border = rgb(0.85, 0.85, 0.85)) 

   # interval numbers 
   text(mids[-1],  ylim[1] + (ylim[2] - ylim[1])/25, n[-1], col = rgb(0.4, 0.4, 0.4), cex = 0.80) 
}

#' Overview plot for iPLS results
#' 
#' @description
#' Shows a plot for iPLS results.
#' 
#' @param x
#' a  (object of class \code{pca})
#' @param ...
#' other arguments
#' 
#' @details
#' See details for \code{\link{plotSelection.ipls}}.
#' 
#' @export
plot.ipls = function(x, ...) {  
   plotSelection(x, ...)
}

#' Print method for iPLS
#' 
#' @description
#' Prints information about the iPLS object structure
#' 
#' @param x
#' a iPLS (object of class \code{ipls})
#' @param ...
#' other arguments
#'
#' @export 
print.ipls = function(x, ...) {
   obj = x
   
   cat('\niPLS results (class ipls)\n')
   cat('\nCall:\n')
   print(obj$call)
   cat('\nMajor fields:\n')
   cat('$var.selected - vector with selected variables\n')
   cat('$int.selected - vector with selected intervals\n')
   cat('$int.num - number of intervals\n')
   cat('$int.width - width of the intervals\n')
   cat('$int.limits - limits for the intervals\n')
   cat('$int.stat - table with statistics for the interval selection results\n')
   cat('$glob.stat - table with statistics for the first iteration of the algorithm\n')
   cat('\nTry summary(obj) and plot(obj) to see details.\n')   
}

#' Summary for iPLS results
#' 
#' @description
#' Shows statistics and algorithm parameters for iPLS results.
#' 
#' @param object
#' a iPLS (object of class \code{ipls})
#' @param glob.ncomp
#' number of components for global PLS model with all intervals
#' @param ...
#' other arguments
#' 
#' @details 
#' The method shows information on the algorithm parameters as well as a table with selected or
#' excluded interval. The table has the following columns: 'step' showing on which iteration
#' an interval was selected or excluded, 'start and 'end' show variable indices for the interval,
#' 'nComp' is a number of components used in a model, 'RMSE' is RMSECV for the model and 'R2' is 
#' coefficient of determination for the same model.
#'
#' @export 
summary.ipls = function(object, glob.ncomp = NULL, ...) {
   # statistics for global model   
   if (is.null(glob.ncomp))
      glob.ncomp = object$gm$ncomp.selected
   else if (glob.ncomp < 1 || glob.ncomp > object$gm$ncomp)
      stop('Wrong value for number of components!')
   
   glob.rmse = object$gm$cvres$rmse[1, glob.ncomp]
   opt.ncomp = object$om$ncomp.selected
   opt.rmse = object$om$cvres$rmse[1, opt.ncomp]
   
   cat('\niPLS variable selection results\n')
   cat(sprintf('  Method: %s\n', object$method))
   cat(sprintf('  Validation: %s\n', crossval.str(object$cv)))
   cat(sprintf('  Number of intervals: %d\n', object$int.num))
   cat(sprintf('  Number of selected intervals: %d\n', length(object$int.selected)))
   cat(sprintf('  RMSECV for global model: %f (%d LVs)\n', glob.rmse, glob.ncomp))
   cat(sprintf('  RMSECV for optimized model: %f (%d LVs)\n', opt.rmse, opt.ncomp))
   cat('\nSummary for selection procedure:\n')
   show(object$int.stat)
}
