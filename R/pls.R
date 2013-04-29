# class and methods for Partial Least Squares regression #
pls = function(X, y, ...) UseMethod("pls")

pls.default = function(X, y, ncomp = 12, cv = 0, autoscale = 1, Xt = NULL, yt = NULL, ...)
{
   X = as.matrix(X)
   y = as.matrix(y)
   
   ncomp = min(ncol(X), nrow(X) - 1, ncomp)
   
   # build a model and apply to calibration set
   model = pls.cal(X, y, ncomp, autoscale)
   model$calres = predict.pls(model, X, y)
   
   # do cross-validation if needed
   if (cv > 0)
      model$cvres = pls.crossval(model, X, y, cv)    
   
   # do test set validation if provided
   if (!is.null(Xt) && !is.null(yt))
   {
      Xt = as.matrix(Xt)
      yt = as.matrix(yt)
      model$testres = predict.pls(model, Xt, yt)
   }
   
   model$ncomp.selected = ncomp
   model$call = match.call()
   
   class(model) = "pls"
   
   model
}

pls.cal = function(X, y, ncomp, autoscale)
{   
   X = as.matrix(X)
   y = as.matrix(y)
   
   if (autoscale > 0)
   {
      # find mean values for X and y
      mX = apply(X, 2, mean)
      my = apply(y, 2, mean)
      
      if (autoscale == 2)
      {
         # calculate stadnard deviations for X variables
         sdX = apply(X, 2, sd)
         sdy = apply(y, 2, sd)
      }
      else
      {
         # use vector with ones if no standardization is needed 
         sdX = rep(1, ncol(X))         
         sdy = rep(1, ncol(y))
      }   
      
      # autoscale X and y
      X = scale(X, center = mX, scale = sdX)
      #y = scale(y, center = my, scale = sdy)
      y = y - my
   }
   
   
   # do SIMPLS
   model = pls.simpls(X, y, ncomp)
   
   model$autoscale = autoscale
   
   if (autoscale > 0)
   {   
      model$mX = mX
      model$sdX = sdX
      model$my = my
      model$sdy = sdy
   }
   
   model$ncomp = ncomp
   
   return (model)
}

## SIMPLS algorithm ###
pls.simpls = function(X, y, ncomp, stripped = FALSE)
{
   X = as.matrix(X)
   y = as.matrix(y)
   
   objnames = rownames(X);
   prednames = colnames(X);
   respnames = colnames(y);
   
   nobj = dim(X)[1]
   npred = dim(X)[2]
   nresp = 1
   
   V = R = matrix(0, nrow = npred, ncol = ncomp)
   tQ = matrix(0, nrow = ncomp, ncol = nresp)
   B = array(0, dim = c(npred, nresp, ncomp))
   
   if (!stripped) {
      P = R
      U = TT = matrix(0, nrow = nobj, ncol = ncomp)
   }
   
   S = crossprod(X, y)
   for (a in 1:ncomp) {
      q.a = 1
      r.a = S %*% q.a
      t.a = X %*% r.a
      t.a = t.a - mean(t.a)
      tnorm = sqrt(c(crossprod(t.a)))
      
      t.a = t.a/tnorm
      r.a = r.a/tnorm
      p.a = crossprod(X, t.a)
      q.a = crossprod(y, t.a)
      v.a = p.a
      if (a > 1) {
         v.a = v.a - V %*% crossprod(V, p.a)
      }
      v.a = v.a/sqrt(c(crossprod(v.a)))
      S = S - v.a %*% crossprod(v.a, S)
      R[, a] = r.a
      tQ[a, ] = q.a
      V[, a] = v.a
      B[, , a] = R[, 1:a, drop = FALSE] %*% tQ[1:a, , drop = FALSE]
      if (!stripped) {
         u.a = y %*% q.a
         if (a > 1) 
            u.a = u.a - TT %*% crossprod(TT, u.a)
         P[, a] = p.a
         TT[, a] = t.a
         U[, a] = u.a
      }
   }
   
   B = B[, 1, ]
   
   if (stripped) {
      list(coeffs = B)
   }
   else {
      
      lvnames = paste("Comp", 1:ncomp)
      ncompnames = paste(1:ncomp, " components")
      rownames(B) = prednames
      colnames(B) = ncompnames
      rownames(TT) = rownames(U) = objnames
      colnames(TT) = colnames(U) = lvnames
      list(
         coeffs = B, 
         xloadings = P, 
         yloadings = t(tQ), 
         weights = R
      )
   }
}   

pls.crossval = function(model, X, y, cv)
{
   X = as.matrix(X)
   y = as.matrix(y)
   
   autoscale = model$autoscale
   ncomp = model$ncomp
   nobj = nrow(X)
   nvar = ncol(X)
   
   # get matrix with indices for cv segments
   idx = crossval(nobj, cv)
   
   seglen = ncol(idx);
   
   yp = matrix(0, nrow = nobj, ncol = ncomp)
   
   # loop over segments
   for (i in 1:nrow(idx))
   {
      ind = na.exclude(idx[i,])
      
      if (length(ind) > 0)
      {   
         Xc = X[-ind, , drop = F]
         yc = y[-ind, , drop = F]
         Xt = X[ind, , drop = F]
         yt = y[ind, , drop = F]
         
         m = pls.cal(Xc, yc, ncomp, autoscale)               
         res = predict.pls(m, Xt, stripped = T)
         yp[ind, ] = res$yp         
      }
   }  
   
   res = plsresult(yp, y)
}

## select optimal ncomp for the model ##
pls.selectncomp = function(model, ncomp)
{
   if (ncomp <= model$ncomp && ncomp > 0)
   {   
      model$ncomp.selected = ncomp;
      model$calres$ncomp.selected = ncomp
      
      if (!is.null(model$cvres)) 
         model$cvres$ncomp.selected = ncomp
      
      if (!is.null(model$testres)) 
         model$testres$ncomp.selected = ncomp
   }
   
   return (model)
}   

predict.pls = function(model, X, y = NULL, stripped = FALSE)
{   
   
   if (model$autoscale > 0)
      X = scale(X, center = model$mX, scale = model$sdX) 
   
   yp = X %*% as.matrix(model$coeffs)
   
   if (model$autoscale > 0) 
      yp = yp + model$my   
   
   if (stripped == FALSE)
   {   
      xscores = X %*% (model$weights %*% solve(t(model$xloadings) %*% model$weights))  
      
      if (!is.null(y))
      {   
         yy = y - model$my
         yscores = as.matrix(yy) %*% model$yloadings   
      }
      
      rownames(xscores) = rownames(yscores) = rownames(X)
      colnames(xscores) = colnames(yscores) = paste("LV", 1:model$ncomp)
      
      res = plsresult(yp, y, X = X,
                      xscores = xscores, 
                      yscores = yscores,
                      xloadings = model$xloadings,
                      yloadings = model$yloadings
      )
   }      
   else
   {
      res = list(yp = yp)
   }   
}  


plotRMSE.pls = function(model, 
                        main = 'RMSE', xlab = 'ncomps', ylab = 'RMSE', 
                        type = 'b',
                        show.legend = T)
{
   ncomp = model$ncomp
   
   legend = c('cal')
   cdata = cbind(1:ncomp, model$calres$rmse)   
   data = list(cdata = cdata)
   
   if (!is.null(model$cvres)) 
   { 
      cvdata = cbind(1:ncomp, model$cvres$rmse)
      data$cvdata = cvdata
      legend = c(legend, 'cv')
   }   
   
   if (!is.null(model$testres)) 
   { 
      testdata = cbind(1:ncomp, model$testres$rmse)
      data$testdata = testdata
      legend = c(legend, 'test')
   }     
   
   if (show.legend == F)
      legend = NULL
   
   mdaplots.lineg(data, legend = legend, type = type, 
                  main = main, xlab = xlab, ylab = ylab)
}

plotXYScores.pls = function(model, ncomp = 1, main = 'XY Scores',
                            show.labels = F, show.legend = T)
{
   main = sprintf('XY scores (ncomp = %d)', ncomp)
   
   cdata = cbind(model$calres$xscores[, ncomp], model$calres$yscores[, ncomp])
   colnames(cdata) = c('X scores', 'Y scores')
   data = list(cdata = cdata)
   legend = c('cal')
      
   if (!is.null(model$testres)) 
   { 
      tdata = cbind(model$testres$xscores[, ncomp], model$testres$yscores[, ncomp])
      data$tdata = tdata
      legend = c(legend, 'test')      
   }   
   
   if (show.legend == F)
      legend = NULL
   
   mdaplots.scatterg(data, legend = legend, show.labels = show.labels, main = main)   
}  

## plot with measured vs predicted y values ##
plotPredictions.pls = function(model, ncomp = 0, main = 'Predictions', 
                                xlab = 'y, measured',
                                ylab = 'y, predicted', 
                                show.labels = F,
                                show.legend = T)
{
   if (ncomp == 0) 
      ncomp = model$ncomp.selected
   
   if (ncomp > model$ncomp) 
   { 
      warning(sprintf('\nChosen ncomp is larger than model has. Use %d instead.', model$ncomp))
      ncomp = model$ncomp
   }
   
   cdata = cbind(model$calres$y, model$calres$yp[, ncomp])
   colnames(cdata) = c('y, measured', 'y, predicted')
   rownames(cdata) = rownames(model$calres$yp)
   legend = c('cal')
   data = list(cdata = cdata)
   
   if (!is.null(model$cvres)) 
   { 
      cvdata = cbind(model$cvres$y, model$cvres$yp[, ncomp])
      colnames(cvdata) = c('y, measured', 'y, predicted')
      rownames(cvdata) = rownames(model$cvres$yp)
      legend = c(legend, 'cv')
      data$cvdata = cvdata
   }   
   
   if (!is.null(model$testres)) 
   { 
      testdata = cbind(model$testres$y, model$testres$yp[, ncomp])
      colnames(testdata) = c('y, measured', 'y, predicted')
      rownames(testdata) = rownames(model$testres$yp)
      legend = c(legend, 'test')
      data$testdata = testdata
   }   
   
   if (show.legend == F)
      legend = NULL
   
   mdaplots.scatterg(data, legend = legend, show.labels = show.labels, 
                     main = main)      
}

plotYResiduals.pls = function(model, ncomp = 0, main = 'Y residuals', 
                               xlab = 'y residuals',
                               ylab = 'y values', 
                               show.labels = F,
                               show.legend = T)
{
   if (ncomp == 0) 
      ncomp = model$ncomp.selected
   
   if (ncomp > model$ncomp) 
   { 
      warning(sprintf('\nChosen ncomp is larger than model has. Use %d instead.', model$ncomp))
      ncomp = model$ncomp
   }
   
   cdata = cbind(model$calres$y, model$calres$yp[, ncomp] - model$calres$y)
   colnames(cdata) = c('y values', 'y residuals')
   rownames(cdata) = rownames(model$calres$yp)
   legend = c('cal')
   data = list(cdata = cdata)
   
   if (!is.null(model$cvres)) 
   { 
      cvdata = cbind(model$cvres$y, model$cvres$yp[, ncomp] - model$cvres$y)
      colnames(cvdata) = c('y values', 'y residuals')
      rownames(cvdata) = rownames(model$cvres$yp)
      legend = c(legend, 'cv')
      data$cvdata = cvdata
   }   
   
   if (!is.null(model$testres)) 
   { 
      testdata = cbind(model$testres$y, model$testres$yp[, ncomp] - model$testres$y)
      colnames(testdata) = c('y values', 'y residuals')
      rownames(testdata) = rownames(model$testres$yp)
      legend = c(legend, 'test')
      data$testdata = testdata
   }   
   
   if (show.legend == F)
      legend = NULL
   
   mdaplots.scatterg(data, legend = legend, show.labels = show.labels, 
                     main = sprintf('Y residuals (ncomp = %d)', ncomp))      
}

plotRegcoeffs.pls = function(model, main = 'Regression coefficients', 
                             show.labels = F, ncomp = 0)
{
   if (ncomp == 0)
      ncomp = model$ncomp.selected
   
   if (ncomp > model$ncomp) 
   { 
      warning(sprintf('\nChosen ncomp is larger than model has. Use %d instead.', model$ncomp))
      ncomp = model$ncomp
   }
   
   data = cbind(1:nrow(model$coeffs), model$coeffs[, ncomp])
   colnames(data) = c('Variables', 'Coefficients')
   rownames(data) = rownames(model$xloadings)
   mdaplots.line(data, type = 'l', main = main, show.labels = show.labels)
}


plotXResiduals.pls = function(model, ncomp = NULL, show.labels = F, show.legend = T)
{
   if (is.null(ncomp))
      ncomp = model$ncomp.selected
   cdata = cbind(model$calres$T2[, ncomp], model$calres$Q2[, ncomp])

   colnames(cdata) = c('T2', 'Q2')
   rownames(cdata) = rownames(model$calres$scores)
   
   data = list(cdata = cdata)
   legend = NULL
   
   if (!is.null(model$testres))
   {
      tdata = cbind(model$testres$T2[, ncomp], model$testres$Q2[, ncomp])      
      colnames(tdata) = c('T2', 'Q2')
      rownames(tdata) = rownames(model$testres$scores)
      
      data$tdata = tdata
      if (show.legend == T)
         legend = c('cal', 'test')
   }      
   
   mdaplots.scatterg(data, main = sprintf('X residuals (ncomp = %d)', ncomp),
                     show.labels = show.labels,
                     legend = legend)
} 

## makes a plot with regression results ##
plot.pls = function(model, show.legend = T, show.labels = F)
{
   par(mfrow = c(2, 2))      
   plotXYScores(model, show.labels = show.labels, show.legend = show.legend)   
   plotRegcoeffs(model, show.labels = show.labels)   
   plotRMSE(model, show.legend = show.legend)   
   plotPredictions(model, show.labels = show.labels, show.legend = show.legend)   
   par(mfrow = c(1, 1))
}


## show summary for a model ##
summary.pls = function(model)
{
   ncomp = model$ncomp.selected
   cat('\nPLS model (class pls) summary\n')
   cat('\nPerformance and validation:\n')
   cat(sprintf('Selected LVs: %d\n\n', ncomp))
   cat('     ')
   cat(sprintf('%6s\t', colnames(as.matrix(model$calres))))
   cat('\n')
   cat('Cal:  ')
   cat(sprintf('%.4f\t', as.matrix(model$calres)[ncomp, , drop = F]))
   cat('\n')
   
   if (!is.null(model$cvres))
   {
      cat('CV:   ')
      cat(sprintf('%.4f\t', as.matrix(model$cvres)[ncomp, , drop = F]))      
      cat('\n')
   }
   
   if (!is.null(model$testres))
   {
      cat('Test: ')
      cat(sprintf('%.4f\t', as.matrix(model$testres)[ncomp, , drop = F]))      
      cat('\n')
   }   
   
}

## print information about a model ##
print.pls = function(model, ...)
{
   cat('\nPLS model (class pls)\n')
   cat('\nCall:\n')
   print(model$call)
   cat('\nMajor fields:\n')
   cat('$coeffs - vector with regression coefficients\n')
   cat('$xloadings - vector with X loadings\n')
   cat('$yloadings - vector with Y loadings\n')
   cat('$calres - results for calibration set\n')
   if (!is.null(model$cvres))
   {
      cat('$cvres - results for cross-validation\n')      
   }   
   if (!is.null(model$testres))
   {
      cat('$testres - results for test set\n')      
   }   
   cat('\nTry summary(model) and plot(model) to see the results of validation\n')
   
}
