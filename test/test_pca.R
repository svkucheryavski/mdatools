library(mdatools)

do_people = T

if (do_people == T) {
   data(People)
   data = people
   data[4, 4] = NA
   pcamodel = pca(data, ncomp = 8, scale = T, cv = 1, info = 'My first model')
   pcamodel = selectCompNum(pcamodel, 5)
}  else
{  
   ncobj = 200
   ntobj = 100
   nvar = 2000
   ncomp = 10
   values = rnorm((ncobj + ntobj) * nvar, 2, 2)
   data = matrix(values[1:(ncobj * nvar)], nrow = ncobj, ncol = nvar)
   tdata = matrix(values[(ncobj * nvar + 1):length(values)], nrow = ntobj, ncol = nvar)

   t1 = system.time({
      pcamodel = pca(data, ncomp = ncomp, scale = T, cv = 1, test.data = tdata)
   })
   t2 = system.time({
      pcamodel = selectCompNum(pcamodel, 8)
   })
   show(t1)
   show(t2)
}

show(pcamodel)
summary(pcamodel)

cat('\n1. Check variance plot for cv results\n')
par(mfrow = c(2, 2))
plotVariance(pcamodel$cvres)
plotVariance(pcamodel$cvres, show.labels = T)
plotCumVariance(pcamodel$cvres)
plotCumVariance(pcamodel$cvres, show.labels = T)
readline('Press Enter to continue...')

cat('\n2. Check residual plot for cv results\n')
par(mfrow = c(2, 2))
plotResiduals(pcamodel$cvres)
plotResiduals(pcamodel$cvres, show.labels = T)
plotResiduals(pcamodel$cvres, ncomp = 2)
plotResiduals(pcamodel$cvres, ncomp = 2, show.labels = T)
readline('Press Enter to continue...')

cat('\n3. print and summary for cv results\n')
print(pcamodel$cvres)
summary(pcamodel$cvres)
readline('Press Enter to continue...')

cat('\n4. Check scores plots for cal results\n')
par(mfrow = c(2, 2))
plotScores(pcamodel$calres)
plotScores(pcamodel$calres, comp = 1, cgroup = data[, 1])
plotScores(pcamodel$calres, c(1, 3), show.labels = T)
plotScores(pcamodel$calres, c(1, 2), cgroup = data[, 1], 
           show.labels = T, show.colorbar = F, show.axes = F)
readline('Press Enter to continue...')

cat('\n5. Check residuals plots for cal results\n')
par(mfrow = c(2, 2))
plotResiduals(pcamodel$calres)
plotResiduals(pcamodel$calres, ncomp = 3, show.labels = T)
plotResiduals(pcamodel$calres, ncomp = 3, cgroup = data[, 1], show.colorbar = T)
plotResiduals(pcamodel$calres, ncomp = 3, cgroup = data[, 2], show.labels = T)
readline('Press Enter to continue...')

cat('\n6. Check variance plot for cal results\n')
par(mfrow = c(2, 2))
plotVariance(pcamodel$calres)
plotVariance(pcamodel$calres, show.labels = T)
plotCumVariance(pcamodel$calres)
plotCumVariance(pcamodel$calres, show.labels = T)
readline('Press Enter to continue...')

cat('\n7. print and summary for cal results\n')
print(pcamodel$calres)
summary(pcamodel$calres)
readline('Press Enter to continue...')

cat('\n8. Check scores plots for pca model\n')
par(mfrow = c(2, 2))
plotScores(pcamodel, 1, show.labels = T)
plotScores(pcamodel)
plotScores(pcamodel, c(1, 3), show.labels = T, show.legend = F)
plotScores(pcamodel, c(1, 2), show.axes = F)
readline('Press Enter to continue...')

cat('\n9. Check loadings plots for pls model\n')
par(mfrow = c(2, 2))
plotLoadings(pcamodel)
plotLoadings(pcamodel, c(1, 3), show.labels = F)
plotLoadings(pcamodel, c(1, 3), type = 'b')
plotLoadings(pcamodel, 1:4, show.legend = T)
readline('Press Enter to continue...')

cat('\n10. Check residuals plots for pca model\n')
par(mfrow = c(2, 2))
plotResiduals(pcamodel)
plotResiduals(pcamodel, ncomp = 2)
plotResiduals(pcamodel, show.labels = T)
plotResiduals(pcamodel, ncomp = 1, show.legend = F)
readline('Press Enter to continue...')

cat('\n11. Check variance plots for pca model\n')
par(mfrow = c(2, 2))
plotVariance(pcamodel)
plotVariance(pcamodel, show.legend = F, show.labels = T)
plotCumVariance(pcamodel)
plotCumVariance(pcamodel, show.legend = F, show.labels = T)
readline('Press Enter to continue...')

cat('\n12. Check summary plot for pca model\n')
par(mfrow = c(1, 1))
plot(pcamodel, show.labels = F)
readline('Press Enter to continue...')




