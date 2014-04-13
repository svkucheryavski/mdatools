library(mdatools)

do_people = T

if (do_people == T) 
{
   data(people)
   data = people
   data[4, 4] = NA
   pcamodel = pca(data, ncomp = 8, scale = T, cv = 1, info = 'My first model')
   pcamodel = selectCompNum(pcamodel, 5)
   gpch = c(16, 17)
   glty = c(1, 2)
} else
{  
   ncobj = 60
   ntobj = 30
   nvar = 10
   ncomp = 6
   values = rnorm((ncobj + ntobj) * nvar, 2, 2)
   data = matrix(values[1:(ncobj * nvar)], nrow = ncobj, ncol = nvar)
   tdata = matrix(values[(ncobj * nvar + 1):length(values)], nrow = ntobj, ncol = nvar)
   gpch = c(16, 17, 19)
   glty = c(1, 2, 3)
   t1 = system.time({
      pcamodel = pca(data, ncomp = ncomp, scale = T, cv = 1, test.data = tdata)
   })
   t2 = system.time({
      pcamodel = selectCompNum(pcamodel, 4)
   })
   show(t1)
   show(t2)
}

show(pcamodel)
summary(pcamodel)

cat('\n1. Check variance plot for cv results\n')
par(mfrow = c(2, 2))
plotVariance(pcamodel$cvres)
plotVariance(pcamodel$cvres, show.labels = T, col = 'red', pch = 17, lty = 2)
plotCumVariance(pcamodel$cvres, colmap = 'gray')
plotCumVariance(pcamodel$cvres, main = 'My variance', xlab = 'PCs', show.labels = T, type = 'h')
readline('Press Enter to continue...')

cat('\n2. Check residual plot for cv results\n')
par(mfrow = c(2, 2))
plotResiduals(pcamodel$cvres)
plotResiduals(pcamodel$cvres, main = 'My residuals', xlab = 'Hott', ylab = 'Res', show.labels = T)
plotResiduals(pcamodel$cvres, ncomp = 2, col = 'red', pch = 17)
plotResiduals(pcamodel$cvres, ncomp = 2, show.labels = T, cgroup = data[, 1])
readline('Press Enter to continue...')

cat('\n3. print and summary for cv results\n')
print(pcamodel$cvres)
summary(pcamodel$cvres)
readline('Press Enter to continue...')

cat('\n4. Check scores plots for cal results\n')
par(mfrow = c(2, 2))
plotScores(pcamodel$calres)
plotScores(pcamodel$calres, c(1, 3), show.labels = T)
plotScores(pcamodel$calres, c(1, 2), cgroup = data[, 1], 
           show.labels = T, show.colorbar = F, show.axes = F)
plotScores(pcamodel$calres, comp = 1, cgroup = data[, 1], show.labels = T,
           colmap = c('red', 'green'), pch = 17)
readline('Press Enter to continue...')


cat('\n5. Check residuals plots for cal results\n')
par(mfrow = c(2, 2))
plotResiduals(pcamodel$calres)
plotResiduals(pcamodel$calres, ncomp = 3, show.labels = T)
plotResiduals(pcamodel$calres, ncomp = 3, cgroup = data[, 1], colmap = 'gray')
plotResiduals(pcamodel$calres, ncomp = 3, cgroup = data[, 2], show.labels = T)
readline('Press Enter to continue...')

cat('\n6. Check variance plot for cal results\n')
par(mfrow = c(2, 2))
plotVariance(pcamodel$calres)
plotVariance(pcamodel$calres, show.labels = T)
plotCumVariance(pcamodel$calres, col = 'red', pch = 17, lty = 2, show.labels = T)
plotCumVariance(pcamodel$calres, colmap = 'gray')
readline('Press Enter to continue...')

cat('\n7. print and summary for cal results\n')
print(pcamodel$calres)
summary(pcamodel$calres)
readline('Press Enter to continue...')

cat('\n8. Check scores plots for pca model\n')
par(mfrow = c(2, 2))
plotScores(pcamodel, 1, show.labels = T)
plotScores(pcamodel, colmap = c('red', 'green'), pch = 17)
plotScores(pcamodel, c(1, 3), show.labels = T, show.legend = F)
plotScores(pcamodel, c(1, 2), show.axes = F, main = 'My scores', xlab = 'PC1', ylab = 'PC2')
readline('Press Enter to continue...')

cat('\n9. Check loadings plots for pls model\n')
par(mfrow = c(2, 2))
plotLoadings(pcamodel)
plotLoadings(pcamodel, c(1), type = 'h')
plotLoadings(pcamodel, c(1, 3), type = 'h', lty = c(1, 2), pch = c(16, 17), show.labels = T)
plotLoadings(pcamodel, 1:4, show.legend = T, colmap = 'gray', lty = c(2, 2, 1, 1),
             legend.position = 'top')
readline('Press Enter to continue...')

cat('\n10. Check residuals plots for pca model\n')
par(mfrow = c(2, 2))
plotResiduals(pcamodel, main = 'My residuals', legend.position = 'top')
plotResiduals(pcamodel, colmap = 'gray', pch = gpch, ncomp = 2)
plotResiduals(pcamodel, show.labels = T, colmap = c('green', 'orange'), legend.position = 'topleft')
plotResiduals(pcamodel, ncomp = 1, show.legend = F)
readline('Press Enter to continue...')

cat('\n11. Check variance plots for pca model\n')
par(mfrow = c(2, 2))
plotVariance(pcamodel, pch = gpch, show.legend = T)
plotVariance(pcamodel, show.legend = F, show.labels = T)
plotCumVariance(pcamodel, colmap = c('red', 'blue'), pch = gpch, lty = glty)
plotCumVariance(pcamodel, show.legend = T, show.labels = T, type = 'h', legend.position = 'topleft')
readline('Press Enter to continue...')

cat('\n12. Check summary plot for pca model\n')
par(mfrow = c(1, 1))
plot(pcamodel, show.labels = T)
readline('Press Enter to continue...')

cat('\n13. Check summary plot printing\n')
png(filename = '../../PCA summary.png', width = 2400, height = 1800, pointsize = 9, res = 300)
plot.new()
plot(pcamodel, show.labels = T)
dev.off()

