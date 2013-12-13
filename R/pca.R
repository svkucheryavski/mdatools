## class and methods for Principal Component Analysis based methods ##

pca = function(x, ncomp = 20, center = T, scale = F, cv = NULL, x.test = NULL, 
               alpha = 0.05, method = 'svd', info = '', ...)
{
   # Calibrate and validate a PCA model.
   #
   # Arguments:
   #   x: a matrix with data values
   #   ncomp: maximum number of components to calculate
   #   center: logical, mean center or not data values
   #   scale: logical, standardize or not data values
   #   cv: number of segments for random cross-validation (1 - for full CV)
   #   x.test: a matrix with data values for test set validation
   #   alpha: a significance level for Q2 residuals
   #   method: method to estimate principal component space (only SVD is supported so far)
   #   info: a short text with information about the model
   #
   # Returns:
   #   model: a PCA model (object of pca class)
      
   x = as.matrix(x)

   # check if data has missing values
   if (sum(is.na(x)) > 0)
   {
      warning('Data has missing values, will try to fix using pca.mvreplace.')
      x = pca.mvreplace(x, center = center, scale = scale)
   }   
   
   # correct maximum number of components
   ncomp = min(ncomp, ncol(x), nrow(x) - 1)

   # calibrate model  
   model = pca.cal(x, ncomp, center = center, scale = scale, method = method)
   model$ncomp = ncomp
   model$ncomp.selected = model$ncomp
   model$info = info
   model$alpha = alpha
   
   # apply model to calibration set
   model$calres = predict.pca(model, x)
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = pca.crossval(model, x, cv, center = center, scale = scale)
   
   # apply model to test set if provided
   if (!is.null(x.test))
      model$testres = predict.pca(model, x.test)
   
   # calculate and assign limit values for T2 and Q2 residuals
   lim = ldecomp.getResLimits(model$eigenvals, nrow(x), model$ncomp, model$alpha)
   model$T2lim = lim$T2lim
   model$Q2lim = lim$Q2lim
   
   model$call = match.call()   
   class(model) = "pca"
   
   model
}

getCalibrationData.pca = function(model)
{
   x = model$calres$scores %*% t(model$loadings) + model$calres$residuals
   
   if (is.numeric(attr(x, 'prep:scale')))
      x = sweep(x, 2L, attr(x, 'prep:scale'), '*', check.margin = F)
   
   if (is.numeric(attr(x, 'prep:center')))
      x = sweep(x, 2L, attr(x, 'prep:center'), '+', check.margin = F)
   
}

selectCompNum.pca = function(model, ncomp)
{
   # Sets user defined number of optimal components.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   ncomp: number of components to set as optimal
   #
   # Returns:
   #   model: the same model with selected number of optimal components 
   
   if (ncomp < 1 || ncomp > model$ncomp)
      stop('Wrong number of selected components!')
   
   model$ncomp.selected = ncomp   
   
   model$calres$ncomp.selected = ncomp
   
   if (!is.null(model$testres))
      model$testres$ncomp.selected = ncomp

   if (!is.null(model$cvres))
      model$cvres$ncomp.selected = ncomp

   model
}

pca.mvreplace = function(x, center = T, scale = F, maxncomp = 7,
                         expvarlim = 0.95, covlim = 10^-6, maxiter = 100)
{
   # Makes missing values replacement with PCA based estimations.
   #
   # Arguments:
   #   x: a matrix with data values
   #   center: logical, mean center data or not
   #   scale: logical, standardize data or not
   #   maxncomp: maximum number of components to calculate for PCA estimation
   #   expvarlim: limit for explained variance to chose optimal components
   #   covlim: limit for covergence 
   #   maxiter: maximum number of iterations if covergence is not reached
   #
   # Returns:
   #   x.rep: data matrix with replaced missing values 
   
   x.rep = x
   mvidx = is.na(x.rep)

   # calculate number of missing values for every variable
   # and make initial estimates with mean values
   for (i in 1:ncol(x))
   {
      mv = is.na(x[, i])
      
      if (sum(mv)/length(x[, i]) > 0.2)
         stop(sprintf('To many missing values in column #%d', i))
      
      x.rep[mv, i] = mean(x[, i], na.rm = T)        
   }  
   
   # autoscale 
   x.rep = scale(x.rep, center = center, scale = scale)
   
   if (scale == T)
      gsd = attr(x.rep, 'scaled:scale')
   
   if (center == T)
      gmean = attr(x.rep, 'scaled:center');         
   
   x = x.rep
   
   n = 1
   scoresp = 0
   scores = 1
   cond = 1
   while (cond > covlim && n < maxiter)
   {    
      n = n + 1
      
      # rescale data on every iteration
      x.rep = scale(x.rep, center = T, scale = F)
      lmean = attr(x.rep, 'scaled:center')
      
      res = pca.svd(x.rep, maxncomp)
      
      expvar = cumsum(res$eigenvals/sum(res$eigenvals))
      ncomp = min(which(expvar >= expvarlim), maxncomp)
            
      if (ncomp == 0)
         ncomp = 1
      if (ncomp == length(expvar))
         ncomp = ncomp - 1
      
      # get and trancate scores and loadings and reestimate the values
      scoresp = scores
      loadings = res$loadings[, 1:ncomp]      
      scores = x.rep %*% loadings
      x.new = scores %*% t(loadings)   
      
      # remove centering
      x.new = sweep(x.new, 2L, lmean, '+', check.margin = F)

      x.rep = x
      x.rep[mvidx] = x.new[mvidx]
      
      if (n > 2)
      {
         # calculate difference between scores for convergence 
         ncompcond = min(ncol(scores), ncol(scoresp))
         cond = sum((scores[, 1:ncompcond] - scoresp[, 1:ncompcond])^2)
      }      
   }   
   
   # rescale the data back and return
   if (scale == T)
      x.rep = sweep(x.rep, 2L, gsd, '*', check.margin = F)
   
   if (center == T)
      x.rep = sweep(x.rep, 2L, gmean, '+', check.margin = F)

   x.rep
}

pca.cal = function(x, ncomp, center = T, scale = F, method = 'svd', simca = F)
{
   # Calibrates a PCA model.
   #
   # Arguments:
   #   x: a matrix with data values  
   #   ncomp: number of principal components to calculate
   #   center: logical, mean center the data values or not
   #   scale: logical, standardize the data values or not
   #   method: which method to use for computing principal component space
   #   simca: logical, is model build for SIMCA or not
   #
   # Returns:
   #   model: a calibrated PCA model 
   
   
   x = prep.autoscale(x, center = center, scale = scale)
   model = pca.svd(x, ncomp)
   
   model$tnorm = sqrt(colSums(model$scores ^ 2)/(nrow(model$scores) - 1));   
   
   rownames(model$loadings) = colnames(x)
   colnames(model$loadings) = paste('Comp', 1:ncol(model$loadings))
   model$center = attr(x, 'prep:center')
   model$scale = attr(x, 'prep:scale')
   
   if (simca)
      model$modpower = pca.getModellingPower(model, x)
   
   model
}  

pca.svd = function(x, ncomp = NULL)
{
   # Singulare Value Decomposition based PCA algorithm.
   #
   # Arguments:
   #   x: a matrix with data values (preprocessed)  
   #   ncomp: number of components to calculate
   #
   # Returns:
   #   res: a list with scores, loadings and eigenvalues of the components 
   
   if (is.null(ncomp)) 
      ncomp = min(ncol(x), nrow(x) - 1)
   else
      ncomp = min(ncomp, ncol(x), nrow(x) - 1)
   
   s = svd(x)
   loadings = s$v[, 1:ncomp]
      
   res = list(
      loadings = loadings,
      scores = x %*% loadings,
      eigenvals = (s$d^2)/(nrow(x) - 1)
   )
}

pca.nipals = function(x, ncomp)
{
   # NIPALS based PCA algorithm.
   #
   # Arguments:
   #   x: a matrix with data values (preprocessed)  
   #   ncomp: number of components to calculate
   #
   # Returns:
   #   res: a list with scores, loadings and eigenvalues of the components 
   
   
   nobj = nrow(x)
   nvar = ncol(x)   
   ncomp = min(ncomp, nobj - 1, nvar)
   
   scores = matrix(0, nrow = nobj, ncol = ncomp)
   loadings = matrix(0, nrow = nvar, ncol = ncomp)
   eigenvals = rep(0, ncomp);
   
   E = x
   for (i in 1:ncomp)
   {      
      ind = which.max(apply(E, 2, sd))
      t = E[, ind, drop = F]
      tau = 99999
      th = 9999

      while (th > 0.000001)
      {      
         p = (t(E) %*% t) / as.vector((t(t) %*% t))
         p = p / as.vector(t(p) %*% p) ^ 0.5
         t = (E %*% p)/as.vector(t(p) %*% p)
         th = abs(tau - as.vector(t(t) %*% t))
         tau = as.vector(t(t) %*% t)
      }
            
      E = E - t %*% t(p)
      scores[, i] = t
      loadings[, i] = p
      eigenvals[i] = tau / (nobj - 1)
   }

   res = list(
      loadings = loadings,
      scores = scores,
      eigenvals = eigenvals
   )   
}

pca.pp = function(x, ncomp)
{
   # Projection Pursuite based PCA algorithm.
   #
   # Arguments:
   #   x: a matrix with data values (preprocessed)  
   #   ncomp: number of components to calculate
}

pca.crossval = function(model, x, cv, center = T, scale = F)
{
   # Cross-validates a PCA model
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   x: a matrix with data values
   #   cv: number of segments for cross-validation (1 for full CV)
   #   center: logical, mean center data values or not
   #   scale: logical, standardize data values or not
   #
   # Returns:
   #   res: results of cross-validation (object of class pcares) 
      
   ncomp = model$ncomp   
   nobj = nrow(x)
   nvar = ncol(x)
   
   # get matrix with indices for cv segments
   idx = crossval(nobj, cv)
   
   seglen = ncol(idx);
   
   Q2 = matrix(0, ncol = ncomp, nrow = nobj)   
   T2 = matrix(0, ncol = ncomp, nrow = nobj)   
   
   # loop over segments
   for (i in 1:nrow(idx))
   {
      ind = na.exclude(idx[i,])
      
      if (length(ind) > 0)
      {   
         x.cal = x[-ind, , drop = F]
         x.val = x[ind, , drop = F]
         
         m = pca.cal(x.cal, ncomp, center, scale)               
         res = predict.pca(m, x.val, cv = T)
         Q2[ind, ] = res$Q2
         T2[ind, ] = res$T2
      }
   }  
   
   rownames(Q2) = rownames(T2) = rownames(x)
   colnames(Q2) = colnames(T2) = colnames(model$scores)

   # in CV results there are no scores only residuals and variances
   res = pcares(NULL, NULL, NULL, model$calres$totvar, model$tnorm, model$ncomp.selected,
                T2, Q2)
   res$Q2lim = model$Q2lim
   res$T2lim = model$T2lim
   
   res
}  

predict.pca = function(model, x, cv = F)
{
   # Applies PCA model to a data.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   x: a matrix with data values
   #   cv: logical, will prediction be used for cross-validation or not
   #
   # Returns:
   #   res: list with PCA results or (for CV) residual distances for each object 
      
   x = prep.autoscale(x, model$center, model$scale)
   scores = x %*% model$loadings
   residuals = x - scores %*% t(model$loadings)
   
   if (cv == F)
   {   
      totvar = sum(x^2)
      res = pcares(scores, model$loadings, residuals, totvar, model$tnorm, model$ncomp.selected)
      res$Q2lim = model$Q2lim
      res$T2lim = model$T2lim
   }   
   else
   {
      res = ldecomp.getDistances(scores, model$loadings, residuals, model$tnorm)   
   }

   res
}  


plotVariance.pca = function(model, type = 'b', variance = 'expvar', 
                            main = 'Variance', xlab = 'Components', 
                            ylab = 'Explained variance, %',
                            show.legend = T, show.labels = F, ...)
{
   # Makes explained variance plot.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   type: type of the plot ('l', 'b', 'h')
   #   variance: which variance to show the plot for ('expvar', 'cumexpvar')
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   main: main title for the plot
   #   show.legend: logical, show or not legend on the plot
   #   show.labels: logical, show or not labels for the data objects
   #   legend.position: position of legend on the plot   
   #   ...: other possible graphical parameters (see also mdaplotg for details)
   
   data = cbind(1:length(model$calres[[variance]]), model$calres[[variance]])
   labels = mdaplot.formatValues(model$calres[[variance]])
   legend  = 'cal'
   
   if (!is.null(model$cvres))
   {
      data = cbind(data, model$cvres[[variance]])
      labels = cbind(labels, mdaplot.formatValues(model$cvres[[variance]]))
      legend = c(legend, 'cv')
   }      

   if (!is.null(model$testres))
   {
      data = cbind(data, model$testres[[variance]])
      labels = cbind(labels, mdaplot.formatValues(model$testres[[variance]]))
      legend = c(legend, 'test')
   }      
      
   if (show.legend == F)
      legend = NULL
   
   if (show.labels == F)
      labels = NULL
   
   mdaplotg(data, main = main, xlab = xlab, ylab = ylab, 
            labels = labels, legend = legend, type = type, ...)   
}

plotCumVariance.pca = function(model, xlab = 'Components', ylab = 'Explained variance, %', 
                               main = 'Cumulative variance', ...)
{
   # Makes a cumulative explained variance plot
   #
   # Arguments:
   #   model: a PCA model (object of class pca)
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   main: main title for the plot
   #   ...: other graphical parameters (see also mdaplotg for details)
   
   plotVariance.pca(model, variance = 'cumexpvar', xlab = xlab, ylab = ylab, main = main, ...)   
}

plotScores.pca = function(model, comp = c(1, 2), type = 'p', main = 'Scores', xlab = NULL, 
                          ylab = NULL, show.labels = F, show.legend = T,
                          show.axes = T, ...)
{
   # Makes a scores plot.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   comp: one or two numbers - which components to make the plot for
   #   type: type of the plot
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   show.labels: logical, show or not labels for the data objects
   #   show.legend: logical, show or not legend on the plot
   #   show.axes: logical, show or not coordinate axis lines crossing (0, 0)   
   
   ncomp = length(comp)
   legend = NULL
   
   if (ncomp == 1)
   {   
      # scores vs objects plot
      
      if (comp > model$ncomp || comp < 1)
         stop('Wrong number of components!')
      
      if (is.null(xlab))
         xlab = 'Objects'
      
      if (is.null(ylab))
         ylab = colnames(model$calres$scores)[comp]
      
      nobj.cal = nrow(model$calres$scores)      
      cdata = cbind(1:nobj.cal, model$calres$scores[, comp])      
      colnames(cdata) = c(xlab, ylab)
      rownames(cdata) = rownames(model$calres$scores)
      
      data = list(cdata = cdata)
      
      if (!is.null(model$testres))
      {
         nobj.test = nrow(model$testres$scores)
         tdata = cbind((nobj.cal + 1):(nobj.cal + nobj.test), model$testres$scores[, comp])      
         rownames(tdata) = rownames(model$testres$scores)
         data$tdata = tdata
         if (show.legend == T)
            legend = c('cal', 'test')
      }   
      
      mdaplotg(data, type = type, main = main, show.labels = show.labels, legend = legend, 
               xlab = xlab, ylab = ylab, ...)
   }
   else if (ncomp == 2)
   {
      # scores vs scores plot
      
      if (comp[1] > model$ncomp || comp[1] < 1 || comp[2] > model$ncomp || comp[2] < 1)
         stop('Wrong component numbers!')
      
      if (is.null(xlab))
         xlab = colnames(model$calres$scores)[comp[1]]
      
      if (is.null(ylab))
         ylab = colnames(model$calres$scores)[comp[2]]

      cdata = cbind(model$calres$scores[, comp[1]], model$calres$scores[, comp[2]])      
      rownames(cdata) = rownames(model$calres$scores)
      
      data = list(cdata = cdata)
      
      if (!is.null(model$testres))
      {
         tdata = cbind(model$testres$scores[, comp[1]], model$testres$scores[, comp[2]])      
         colnames(tdata) = c(xlab, ylab)
         rownames(tdata) = rownames(model$testres$scores)
         data$tdata = tdata
         if (show.legend == T)
            legend = c('cal', 'test')
      }   
      
      if (show.axes == T)
         show.lines = c(0, 0)      
      else
         show.lines = F

      mdaplotg(data, type = type, main = main, show.labels = show.labels, legend = legend, 
               show.lines = show.lines, xlab = xlab, ylab = ylab, ...)
      }
   else
   {
      stop('Wrong number of components!')
   }   
}  

plotResiduals.pca = function(model, ncomp = NULL, main = NULL, xlab = 'T2',
                             ylab = 'Q2', show.labels = F, show.legend = T, show.limits = T, ...)
{
   # Makes a residuals (T2 vs Q2) plot.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   ncomp: number of componnts to make the plot for
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   show.labels: logical, show or not labels for the data objects
   #   show.legend: logical, show or not legend on the plot
   #   show.limits: logical, show or not statistical limits on the plot   
   
   if (is.null(main))
   {
      if (is.null(ncomp))
         main = 'Residuals'
      else
         main = sprintf('Residuals (ncomp = %d)', ncomp)      
   }   
   
   if (is.null(ncomp))
      ncomp = model$ncomp.selected
   
   if (show.limits == T)
      show.lines = c(model$T2lim[1, ncomp], model$Q2lim[1, ncomp])
   else
      show.lines = F
   
   if (ncomp > model$ncomp || ncomp < 1)
      stop('Wrong number of components!')

   cdata = cbind(model$calres$T2[, ncomp], model$calres$Q2[, ncomp])
   rownames(cdata) = rownames(model$calres$scores)
   legend = 'cal'   
   data = list(cdata = cdata)

   if (!is.null(model$cvres))
   {
      cvdata = cbind(model$cvres$T2[, ncomp], model$cvres$Q2[, ncomp])      
      rownames(cvdata) = rownames(model$cvres$T2)      
      data$cvdata = cvdata
      legend = c(legend, 'cv')
   }      
   
   if (!is.null(model$testres))
   {
      tdata = cbind(model$testres$T2[, ncomp], model$testres$Q2[, ncomp])      
      rownames(tdata) = rownames(model$testres$scores)      
      data$tdata = tdata
      legend = c(legend, 'test')
   }      

   if (show.legend == F)
      legend = NULL
   
   mdaplotg(data, main = main, xlab = xlab, ylab = ylab,
            show.labels = show.labels, legend = legend, show.lines = show.lines, ...)
}  

plotLoadings.pca = function(model, comp = c(1, 2), type = NULL, main = 'Loadings', xlab = NULL, 
                            ylab = NULL, show.labels = T, show.legend = T,  show.axes = T, ...)
{
   # Makes a loadings plot.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   comp: one or two numbers - which components to make the plot for
   #   type: type of the plot ('p', 'b', 'l', 'h')
   #   main: main title for the plot
   #   xlab: label for x axis
   #   ylab: label for y axis
   #   show.labels: logical, show or not labels for the data objects
   #   show.legend: logical, show or not legend on the plot
   #   show.axes: logical, show or not coordinate axis lines crossing (0, 0)   
   
   ncomp = length(comp)
   
   if (ncomp == 2 && (type == 'p' || is.null(type)))
   {
      # scatter plot
      
      data = model$loadings[, c(comp[1], comp[2])]      
      mdaplot(data, show.labels = show.labels, main = main, xlab = xlab, ylab = ylab, 
              show.lines = c(0, 0), ...)      
   }  
   else if (ncomp < 1 | ncomp > 8 )
   {
      stop ('Number of components must be between 1 and 8!')
   }  
   else
   {
      # loadings vs objects
      
      if (is.null(type))
         type = 'b'
      
      if (is.null(xlab))
         xlab = 'Variables'
      
      if (is.null(ylab))
         ylab = 'Loadings'
      
      data = cbind(1:nrow(model$loadings), model$loadings[, comp, drop = F])            
      rownames(data) = rownames(model$loadings)
      
      if (show.legend == T)
         legend = colnames(data)[-1];

      mdaplotg(data, legend = legend, type = type, show.labels = show.labels, 
               main = main, ylab = ylab, xlab = xlab, ...)
   }   
}

plot.pca = function(model, comp = c(1, 2), show.labels = F, show.legend = T)
{   
   # Makes plots for PCA overview.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)
   #   comp: which components to show scores and loadings plots for
   #   show.labels: logical, show data objects labels on the plots or not
   #   show.legend: logical, show legend or not   
   
   par(mfrow = c(2, 2))
   plotScores(model, comp = comp, show.labels = show.labels, show.legend = show.legend)
   plotLoadings(model, comp = comp, show.labels = show.labels, show.legend = show.legend)
   plotResiduals(model, ncomp = model$ncomp.selected,  show.labels = show.labels, 
                 show.legend = show.legend, show.limits = T)
   plotCumVariance(model, show.legend = show.legend)
   par(mfrow = c(1, 1))
}

print.pca = function(model, ...)
{
   # Prints information about the PCA model.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)     
   
   cat('\nPCA model (class pca)\n')
   
   if (length(model$info) > 1)
   {
      cat('\nInfo:\n')
      cat(model$info)      
   }   
   
   cat('\n\nCall:\n')
   print(model$call)
   
   cat('\nMajor fields:\n')   
   cat('$loadings - matrix with loadings\n')
   cat('$eigenvals - eigenvalues for components\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - number of selected components\n')
   cat('$center - values for centering data\n')
   cat('$scale - values for scaling data\n')
   cat('$cv - number of segments for cross-validation\n')
   cat('$alpha - significance level for Q2 residuals\n')
   cat('$calres - results (scores, etc) for calibration set\n')
   
   if (!is.null(model$cvres))
   {
      cat('$cvres - results for cross-validation\n')      
   }   
   if (!is.null(model$testres))
   {
      cat('$testres - results for test set\n')      
   }    
}

summary.pca = function(model)
{
   # Shows summary for a PCA model.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   
   ncomp = model$ncomp
   
   cat('\nPCA model (class pca) summary\n')

   if (length(model$info) > 0)
      cat(sprintf('\nInfo:\n%s\n\n', model$info))
   
   data = cbind(round(model$eigenvals[1:model$ncomp], 3), 
                round(model$calres$expvar, 2),
                round(model$calres$cumexpvar, 2))
   
   colnames(data) = c('Eigvals', 'Expvar', 'Cumexpvar')
   show(data)
}
