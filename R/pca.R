# class and methods for Principal component analysis #
pca = function(data, ...) UseMethod("pca")

pca.default = function(data, ncomp = 20, center = T, scale = F, cv = NULL, 
                       test.data = NULL, alpha = 0.05, info = '', ...)
{
   data = as.matrix(data)

   # check if data has missing values
   if (sum(is.na(data)) > 0)
   {
      warning('Data has missing values, will try to fix using pca.mvreplace.')
      data = pca.mvreplace(data, center = center, scale = scale)
   }   
   
   # calibrate model and select number of components 
   ncomp = min(ncomp, ncol(data), nrow(data) - 1)
   res = pca.cal(data, ncomp, center = center, scale = scale)
   model = list(
      loadings = res$loadings,
      eigenvals = res$eigenvals,
      tnorm = res$tnorm,
      center = res$center,
      scale = res$scale, 
      ncomp = ncol(res$loadings)
   )
   
   model$ncomp.selected = model$ncomp
   model$info = info
   model$alpha = alpha
   
   # apply model to calibration set
   model$calres = predict.pca(model, data)
   
   # do cross-validation if needed
   if (!is.null(cv))
      model$cvres = pca.crossval(model, data, cv)
   
   # apply model to test set if any
   if (!is.null(test.data))
      model$testres = predict.pca(model, test.data)
   
   # calculate and assign limit values for T2 and Q2 residuals
   lim = ldecomp.getResLimits(model$eigenvals, nrow(data), model$ncomp.selected, model$alpha)
   model$T2lim = lim$T2lim
   model$Q2lim = lim$Q2lim
   
   model$call = match.call()   
   class(model) = "pca"
   
   return (model)
}

selectCompNum.pca = function(model, ncomp)
{
   if (is.null(model))
      stop('Object with model is not specified!')   

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
      
   return (model)
}

pca.mvreplace = function(data, center = T, scale = F, maxncomp = 7,
                         expvarlim = 0.95, covlim = 10^-6, maxiter = 100)
{
   # initial estimates with mean values   
   cdata = data
   mvidx = is.na(cdata)
   
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
   
   # iterations
   n = 1
   scoresp = 0
   scores = 1
   cond = 1
   while (cond > covlim && n < maxiter)
   {    
      n = n + 1
      
      cdata = scale(cdata, center = T, scale = F)
      lmean = attr(cdata, 'scaled:center')
      
      res = pca.svd(cdata, maxncomp)
      
      expvar = cumsum(res$eigenvals/sum(res$eigenvals))
      ncomp = min(which(expvar >= expvarlim), maxncomp)
            
      if (ncomp == 0)
         ncomp = 1
      if (ncomp == length(expvar))
         ncomp = ncomp - 1
      
      scoresp = scores
      loadings = res$loadings[, 1:ncomp]      
      scores = cdata %*% loadings
      newdata = scores %*% t(loadings)   
      
      newdata = sweep(newdata, 2L, lmean, '+', check.margin = F)

      cdata = data
      cdata[mvidx] = newdata[mvidx]
      
      if (n > 2)
      {
         # calculate difference between scores
         ncompcond = min(ncol(scores), ncol(scoresp))
         cond = sum((scores[, 1:ncompcond] - scoresp[, 1:ncompcond])^2)
      }      
   }   
   
   # rescale the data back and return
   if (scale == T)
      cdata = sweep(cdata, 2L, gsd, '*', check.margin = F)
   
   if (center == T)
      cdata = sweep(cdata, 2L, gmean, '+', check.margin = F)
   
   return (cdata)
}

pca.cal = function(data, ncomp, center = T, scale = F)
{
   data = prep.autoscale(data, center = center, scale = scale)
   res = pca.svd(data, ncomp)
   
   res$tnorm = sqrt(colSums(res$scores ^ 2)/(nrow(res$scores) - 1));   
   
   rownames(res$loadings) = colnames(data)
   colnames(res$loadings) = paste('Comp', 1:ncol(res$loadings))
   res$center = attr(data, 'prep:center')
   res$scale = attr(data, 'prep:scale')
   
   return (res)
}  

pca.svd = function(data, ncomp)
{
   nobj = nrow(data)
   nvar = ncol(data)   
   ncomp = min(ncomp, nobj - 1, nvar)

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

pca.crossval = function(model, data, cv)
{
   scale = model$scale
   center = model$center
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
         
         m = pca.cal(datac, ncomp, model$center, model$scale)               
         res = predict.pca(m, datat, stripped = T)
         Q2[ind, ] = res$Q2
         T2[ind, ] = res$T2
      }
   }  
   
   rownames(Q2) = rownames(T2) = rownames(data)
   colnames(Q2) = colnames(T2) = colnames(model$scores)
   
   res = pcacvres(T2, Q2, model$calres$fullvar)
}  

predict.pca = function(model, data, stripped = F)
{
   data = prep.autoscale(data, model$center, model$scale)
   scores = data %*% model$loadings
   residuals = data - scores %*% t(model$loadings)
   
   if (stripped == F)
   {   
      fullvar = sum(data^2)
      res = pcares(scores, model$loadings, residuals, fullvar, model$tnorm, model$ncomp.selected)
   }   
   else
   {
      res = ldecomp.getDistances(scores, model$loadings, residuals, model$tnorm)   
   }
   
   return (res)
}  


plotCumVariance.pca = function(obj, show.labels = F, show.legend = T)
{
   legend = NULL
   
   cdata = cbind(0:length(obj$calres$cumexpvar), c(0, obj$calres$cumexpvar))
   colnames(cdata) = c('Components', 'Explained variance, %')
   rownames(cdata) = c(0, round(obj$calres$cumexpvar, 1))
   data = list(cdata = cdata)
   
   if (show.legend == T)
      legend  = 'cal'
   
   if (!is.null(obj$cvres))
   {
      cvdata = cbind(0:length(obj$cvres$cumexpvar), c(0, obj$cvres$cumexpvar))
      colnames(cvdata) = c('Components', 'Explained variance, %')
      rownames(cvdata) = c(0, round(obj$cvres$cumexpvar, 1))
      
      data$cvdata = cvdata
      if (show.legend == T)
         legend = c(legend, 'cv')
   }      

   if (!is.null(obj$testres))
   {
      tdata = cbind(0:length(obj$testres$cumexpvar), c(0, obj$testres$cumexpvar))
      colnames(tdata) = c('Components', 'Explained variance, %')
      rownames(tdata) = c(0, round(obj$testres$cumexpvar, 1))
      
      data$tdata = tdata
      if (show.legend == T)
         legend = c(legend, 'test')
   }      
   
   mdaplots.lineg(data, main = 'Cumulative variance', ylab = 'Explained variance, %', 
                  show.labels = show.labels, legend = legend, pch = 16, type = 'b')   
}

plotVariance.pca = function(obj, show.labels = F, show.legend = T)
{
   legend = NULL
   
   cdata = cbind(1:length(obj$calres$expvar), obj$calres$expvar)
   colnames(cdata) = c('Components', 'Explained variance, %')
   rownames(cdata) = round(obj$calres$expvar, 1)
   data = list(cdata = cdata)

   if (show.legend == T)
      legend  = 'cal'
   
   if (!is.null(obj$cvres))
   {
      cvdata = cbind(1:length(obj$cvres$expvar), obj$cvres$expvar)
      colnames(cvdata) = c('Components', 'Explained variance, %')
      rownames(cvdata) = round(obj$cvres$expvar, 1)      
      data$cvdata = cvdata
      
      if (show.legend == T)
         legend = c(legend, 'cv')
   }      
   
   if (!is.null(obj$testres))
   {
      tdata = cbind(1:length(obj$testres$expvar), obj$testres$expvar)
      colnames(tdata) = c('Components', 'Explained variance, %')
      rownames(tdata) = round(obj$testres$expvar, 1)      
      data$tdata = tdata
      
      if (show.legend == T)
         legend = c(legend, 'test')
   }      
   
   mdaplots.lineg(data, main = 'Variance', ylab = 'Explained variance, %', 
                  show.labels = show.labels, legend = legend, pch = 16, type = 'b')   
}

plotScores.pca = function(obj, comp = c(1, 2), 
                          show.labels = F, show.legend = T,
                          show.axes = T)
{
   legend = NULL;
   if (length(comp) == 1)
   {   
      nobj.cal = nrow(obj$calres$scores)
      
      # scores vs objects
      cdata = cbind(1:nobj.cal, obj$calres$scores[, comp])      
      colnames(cdata) = c('Objects', colnames(obj$calres$scores)[comp])
      rownames(cdata) = rownames(obj$calres$scores)
      
      data = list(cdata = cdata)
      
      if (!is.null(obj$testres))
      {
         nobj.test = nrow(obj$testres$scores)
         tdata = cbind((nobj.cal + 1):(nobj.cal + nobj.test), obj$testres$scores[, comp])      
         colnames(tdata) = c('Objects', colnames(obj$testres$scores)[comp])
         rownames(tdata) = rownames(obj$testres$scores)
         data$tdata = tdata
         if (show.legend == T)
            legend = c('cal', 'test')
      }   
      
      mdaplots.scatterg(data, main = 'Scores', 
                        show.labels = show.labels, 
                        legend = legend
      )
   }
   else if (length(comp) == 2)
   {
      # scores vs scores
      cdata = cbind(obj$calres$scores[, comp[1]], obj$calres$scores[, comp[2]])      
      colnames(cdata) = colnames(obj$calres$scores)[comp]
      rownames(cdata) = rownames(obj$calres$scores)
      
      data = list(cdata = cdata)
      
      if (!is.null(obj$testres))
      {
         tdata = cbind(obj$testres$scores[, comp[1]], obj$testres$scores[, comp[2]])      
         colnames(tdata) = colnames(obj$testres$scores)[comp]
         rownames(tdata) = rownames(obj$testres$scores)
         data$tdata = tdata
         if (show.legend == T)
            legend = c('cal', 'test')
      }   
      
      if (show.axes == T)
         show.lines = c(0, 0)      
      else
         show.lines = F

      mdaplots.scatterg(data, main = 'Scores', 
                        show.labels = show.labels, 
                        legend = legend, 
                        show.lines = show.lines)
      
   }
   else
   {
      stop('Wrong number of components!')
   }   
}  

plotResiduals.pca = function(obj, ncomp = NULL, show.labels = F, show.legend = T, show.limits = T)
{
   if (show.limits == T && (is.null(ncomp) || ncomp == obj$ncomp.selected))
      show.lines = c(obj$T2lim, obj$Q2lim)
   else
      show.lines = F

   if (is.null(ncomp))
      ncomp = obj$ncomp.selected
   
   cdata = cbind(obj$calres$T2[, ncomp], obj$calres$Q2[, ncomp])
   colnames(cdata) = c('T2', 'Q2')
   rownames(cdata) = rownames(obj$calres$scores)
   
   data = list(cdata = cdata)
   legend = NULL

   if (show.legend == T)
      legend = 'cal'

   if (!is.null(obj$cvres))
   {
      cvdata = cbind(obj$cvres$T2[, ncomp], obj$cvres$Q2[, ncomp])      
      colnames(cvdata) = c('T2', 'Q2')
      rownames(cvdata) = rownames(obj$cvres$T2)
      
      data$cvdata = cvdata
      if (show.legend == T)
         legend = c(legend, 'cv')
   }      
   
   if (!is.null(obj$testres))
   {
      tdata = cbind(obj$testres$T2[, ncomp], obj$testres$Q2[, ncomp])      
      colnames(tdata) = c('T2', 'Q2')
      rownames(tdata) = rownames(obj$testres$scores)
      
      data$tdata = tdata
      if (show.legend == T)
         legend = c(legend, 'test')
   }      

   mdaplots.scatterg(data, main = sprintf('Residuals (ncomp = %d)', ncomp),
                     show.labels = show.labels,
                     legend = legend,
                     show.lines = show.lines)
}  

plotLoadings.pca = function(obj, comp = c(1, 2), show.labels = T, 
                            show.legend = T, type = 'p')
{
   ncomp = length(comp)
   
   if (ncomp == 2 && type == 'p')
   {
      # scatter plot
      data =obj$loadings[, c(comp[1], comp[2])]      
      mdaplots.scatter(data, show.labels = show.labels, main = 'Loadings', show.lines = c(0, 0))      
   }  
   else if (ncomp < 1 | ncomp > 8 )
   {
      stop ('Number of components must be between 1 and 8!')
   }  
   else
   {
      if (type == 'p')
         type = 'l'
      
      # line plot
      data = cbind(1:nrow(obj$loadings), obj$loadings[, comp, drop = F])            
      mdaplots.line(data, show.legend = show.legend, type = type, main = 'Loadings', 
                    ylab = 'Loadings', xlab = 'Variables')
   }   
}

plot.pca = function(obj, comp = c(1, 2), show.labels = F, show.legend = T)
{   
   par(mfrow = c(2, 2))
   plotScores(obj, comp = comp, show.labels = show.labels, 
              show.legend = show.legend)
   plotLoadings(obj, comp = comp, show.labels = show.labels, 
                show.legend = show.legend)
   plotResiduals(obj, ncomp = obj$ncomp.selected, 
                 show.labels = show.labels, show.legend = show.legend, show.limits = T)
   plotCumVariance(obj, show.legend = show.legend)
   par(mfrow = c(1, 1))
}

print.pca = function(model, ...)
{
   cat('\nPCA model (class pca)\n')
   
   if (length(model$info) > 0)
   {
      cat('\nInfo:\n')
      cat(model$info)      
   }   
   
   cat('\nCall:\n')
   print(model$call)
   
   cat('\nMajor fields and methods:\n')   
   cat('$loadings - matrix with loadings\n')
   cat('$eigenvals - eigenvalues for components\n')
   cat('$ncomp - number of calculated components\n')
   cat('$ncomp.selected - selected number of components\n')
   cat('$center - values for centering data\n')
   cat('$scale - values for scaling data\n')
   cat('$cv - number of segments for cross-validation\n')
   cat('$alpha - significance level for Q2 residuals\n\n')
   cat('$calres - results (scores, etc) for calibration set\n')
   
   if (!is.null(model$cvres))
   {
      cat('$cvres - results for cross-validation\n')      
   }   
   if (!is.null(model$testres))
   {
      cat('$testres - results for test set\n')      
   }   
   cat('\nTry also: show(model$calres), summary(model) and plot(model)\n')
   
}

summary.pca = function(model)
{
   ncomp = model$ncomp.selected
   cat('\nPCA model (class pca) summary\n')

   if (length(model$info) > 0)
      cat(sprintf('\nInfo:\n%s\n\n', model$info))
   
   data = cbind(round(model$eigenvals[1:model$ncomp], 3), 
                round(model$calres$expvar, 2),
                round(model$calres$cumexpvar, 2))
   
   colnames(data) = c('Eigvals', 'Exp. var', 'Cum. exp. var')
   show(data)
}

