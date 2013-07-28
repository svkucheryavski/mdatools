## class and methods for Principal Component Analysis based methods ##

pca = function(data, ncomp = 20, center = T, scale = F, cv = NULL, test.data = NULL, 
               alpha = 0.05, method = 'svd', info = '', ...)
{
   # Calibrate and validate a PCA model.
   #
   # Arguments:
   #   data: a matrix with data values
   #   ncomp: maximum number of components to calculate
   #   center: logical, mean center or not data values
   #   scale: logical, standardize or not data values
   #   cv: number of segments for random cross-validation (1 - for full CV)
   #   test.data: a matrix with data values for test set validation
   #   alpha: a significance level for Q2 residuals
   #   method: method to estimate principal component space (only SVD is supported so far)
   #   info: a short text with information about the model
   #
   # Returns:
   #   model: a PCA model (object of pca class)
      
   data = as.matrix(data)

   # check if data has missing values
   if (sum(is.na(data)) > 0)
   {
      warning('Data has missing values, will try to fix using pca.mvreplace.')
      data = pca.mvreplace(data, center = center, scale = scale)
   }   
   
   # correct maximum number of components
   ncomp = min(ncomp, ncol(data), nrow(data) - 1)

   # calibrate model  
   model = pca.cal(data, ncomp, center = center, scale = scale, method = method)
   model$ncomp = ncomp
   model$ncomp.selected = model$ncomp
   model$info = info
   model$alpha = alpha
   
   # apply model to calibration set
   model$calres = predict.pca(model, data)
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = pca.crossval(model, data, cv, center = center, scale = scale)
   
   # apply model to test set if provided
   if (!is.null(test.data))
      model$testres = predict.pca(model, test.data)
   
   # calculate and assign limit values for T2 and Q2 residuals
   lim = ldecomp.getResLimits(model$eigenvals, nrow(data), model$ncomp.selected, model$alpha)
   model$T2lim = lim$T2lim
   model$Q2lim = lim$Q2lim
   
   model$call = match.call()   
   class(model) = "pca"
   
   model
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
   lim = ldecomp.getResLimits(model$eigenvals, nrow(model$calres$scores), 
                              model$ncomp.selected, model$alpha)
   model$T2lim = lim$T2lim
   model$Q2lim = lim$Q2lim         
   
   model$calres$ncomp.selected = ncomp
   
   if (!is.null(model$testres))
      model$testres$ncomp.selected = ncomp

   if (!is.null(model$cvres))
      model$cvres$ncomp.selected = ncomp

   model
}

pca.mvreplace = function(data, center = T, scale = F, maxncomp = 7,
                         expvarlim = 0.95, covlim = 10^-6, maxiter = 100)
{
   # Makes missing values replacement with PCA based estimations.
   #
   # Arguments:
   #   data: a matrix with data values
   #   center: logical, mean center data or not
   #   scale: logical, standardize data or not
   #   maxncomp: maximum number of components to calculate for PCA estimation
   #   expvarlim: limit for explained variance to chose optimal components
   #   covlim: limit for covergence 
   #   maxiter: maximum number of iterations if covergence is not reached
   #
   # Returns:
   #   cdata: data matrix with replaced missing values 
   
   cdata = data
   mvidx = is.na(cdata)

   # calculate number of missing values for every variable
   # and make initial estimates with mean values
   for (i in 1:ncol(data))
   {
      mv = is.na(data[, i])
      
      if (sum(mv)/length(data[, i]) > 0.2)
         stop(sprintf('To many missing values in column #%d', i))
      
      cdata[mv, i] = mean(data[, i], na.rm = T)        
   }  
   
   # autoscale 
   cdata = scale(cdata, center = center, scale = scale)
   
   if (scale == T)
      gsd = attr(cdata, 'scaled:scale')
   
   if (center == T)
      gmean = attr(cdata, 'scaled:center');         
   
   data = cdata
   
   n = 1
   scoresp = 0
   scores = 1
   cond = 1
   while (cond > covlim && n < maxiter)
   {    
      n = n + 1
      
      # rescale data on every iteration
      cdata = scale(cdata, center = T, scale = F)
      lmean = attr(cdata, 'scaled:center')
      
      res = pca.svd(cdata, maxncomp)
      
      expvar = cumsum(res$eigenvals/sum(res$eigenvals))
      ncomp = min(which(expvar >= expvarlim), maxncomp)
            
      if (ncomp == 0)
         ncomp = 1
      if (ncomp == length(expvar))
         ncomp = ncomp - 1
      
      # get and trancate scores and loadings and reestimate the values
      scoresp = scores
      loadings = res$loadings[, 1:ncomp]      
      scores = cdata %*% loadings
      newdata = scores %*% t(loadings)   
      
      # remove centering
      newdata = sweep(newdata, 2L, lmean, '+', check.margin = F)

      cdata = data
      cdata[mvidx] = newdata[mvidx]
      
      if (n > 2)
      {
         # calculate difference between scores for convergence 
         ncompcond = min(ncol(scores), ncol(scoresp))
         cond = sum((scores[, 1:ncompcond] - scoresp[, 1:ncompcond])^2)
      }      
   }   
   
   # rescale the data back and return
   if (scale == T)
      cdata = sweep(cdata, 2L, gsd, '*', check.margin = F)
   
   if (center == T)
      cdata = sweep(cdata, 2L, gmean, '+', check.margin = F)

   cdata
}

pca.cal = function(data, ncomp, center = T, scale = F, method = 'svd')
{
   # Calibrates a PCA model.
   #
   # Arguments:
   #   data: a matrix with data values  
   #   ncomp: number of principal components to calculate
   #   center: logical, mean center the data values or not
   #   scale: logical, standardize the data values or not
   #   method: which method to use for computing principal component space
   #
   # Returns:
   #   model: a calibrated PCA model 
   
   
   data = prep.autoscale(data, center = center, scale = scale)
   model = pca.svd(data, ncomp)
   
   model$tnorm = sqrt(colSums(model$scores ^ 2)/(nrow(model$scores) - 1));   
   
   rownames(model$loadings) = colnames(data)
   colnames(model$loadings) = paste('Comp', 1:ncol(model$loadings))
   model$center = attr(data, 'prep:center')
   model$scale = attr(data, 'prep:scale')

   model
}  

pca.svd = function(data, ncomp = NULL)
{
   # Singulare Value Decomposition based PCA algorithm.
   #
   # Arguments:
   #   data: a matrix with data values (preprocessed)  
   #   ncomp: number of components to calculate
   #
   # Returns:
   #   res: a list with scores, loadings and eigenvalues of the components 
   
   if (is.null(ncomp)) 
      ncomp = min(ncol(data), nrow(data) - 1)
   else
      ncomp = min(ncomp, ncol(data), nrow(data) - 1)
   
   s = svd(data)
   loadings = s$v[, 1:ncomp]
      
   res = list(
      loadings = loadings,
      scores = data %*% loadings,
      eigenvals = (s$d^2)/(nrow(data) - 1)
   )
}

pca.nipals = function(data, ncomp)
{
   # NIPALS based PCA algorithm.
   #
   # Arguments:
   #   data: a matrix with data values (preprocessed)  
   #   ncomp: number of components to calculate
   #
   # Returns:
   #   res: a list with scores, loadings and eigenvalues of the components 
   
   
   nobj = nrow(data)
   nvar = ncol(data)   
   ncomp = min(ncomp, nobj - 1, nvar)
   
   scores = matrix(0, nrow = nobj, ncol = ncomp)
   loadings = matrix(0, nrow = nvar, ncol = ncomp)
   eigenvals = rep(0, ncomp);
   
   E = data
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

pca.pp = function(data, ncomp)
{
   # Projection Pursuite based PCA algorithm.
   #
   # Arguments:
   #   data: a matrix with data values (preprocessed)  
   #   ncomp: number of components to calculate
}

pca.crossval = function(model, data, cv, center = T, scale = F)
{
   # Cross-validates a PCA model
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   data: a matrix with data values
   #   cv: number of segments for cross-validation (1 for full CV)
   #   center: logical, mean center data values or not
   #   scale: logical, standardize data values or not
   #
   # Returns:
   #   res: results of cross-validation (object of class pcares) 
      
   ncomp = model$ncomp   
   nobj = nrow(data)
   nvar = ncol(data)
   
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
         datac = data[-ind, , drop = F]
         datat = data[ind, , drop = F]
         
         m = pca.cal(datac, ncomp, center, scale)               
         res = predict.pca(m, datat, cv = T)
         Q2[ind, ] = res$Q2
         T2[ind, ] = res$T2
      }
   }  
   
   rownames(Q2) = rownames(T2) = rownames(data)
   colnames(Q2) = colnames(T2) = colnames(model$scores)

   # in CV results there are no scores only residuals and variances
   res = pcares(NULL, NULL, NULL, model$calres$totvar, model$tnorm, model$ncomp.selected,
                T2, Q2)
}  

predict.pca = function(model, data, cv = F)
{
   # Applies PCA model to a data.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   data: a matrix with data values
   #   cv: logical, will prediction be used for cross-validation or not
   #
   # Returns:
   #   res: list with PCA results or (for CV) residual distances for each object 
      
   data = prep.autoscale(data, model$center, model$scale)
   scores = data %*% model$loadings
   residuals = data - scores %*% t(model$loadings)
   
   if (cv == F)
   {   
      totvar = sum(data^2)
      res = pcares(scores, model$loadings, residuals, totvar, model$tnorm, model$ncomp.selected)
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

plotScores.pca = function(model, comp = c(1, 2), main = 'Scores', xlab = NULL, ylab = NULL,
                          show.labels = F, show.legend = T,
                          show.axes = T, ...)
{
   # Makes a scores plot.
   #
   # Arguments:
   #   model: a PCA model (object of class pca)  
   #   comp: one or two numbers - which components to make the plot for
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
      
      mdaplotg(data, type = 'p', main = main, show.labels = show.labels, legend = legend, 
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

      mdaplotg(data, type = 'p', main = main, show.labels = show.labels, legend = legend, 
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
   
   
   if (show.limits == T && (is.null(ncomp) || ncomp == model$ncomp.selected))
      show.lines = c(model$T2lim, model$Q2lim)
   else
      show.lines = F

   if (is.null(ncomp))
      ncomp = model$ncomp.selected

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
   
   if (is.null(main))
      main = sprintf('Residuals (ncomp = %d)', ncomp)
   
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
      
      if (show.legend == T)
         legend = colnames(data)[-1];

      mdaplotg(data, legend = legend, type = type, 
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
   
   if (length(model$info) > 0)
   {
      cat('\nInfo:\n')
      cat(model$info)      
   }   
   
   cat('\n\nCall:\n')
   print(model$call)
   
   cat('\nMajor fields and methods:\n')   
   cat('$loadings - matrix with loadings\n')
   cat('$eigenvals - eigenvalues for components\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - selected number of components\n')
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
   cat('\nTry also: show(model$calres), summary(model) and plot(model)\n\n')
   
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
