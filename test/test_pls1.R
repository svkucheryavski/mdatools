library('mdatools')
data(simdata)
ny = 3

conc.c = simdata$conc.c

y = simdata$conc.c[, ny]
x = simdata$spectra.c

y.t = simdata$conc.t[, ny]
x.t = simdata$spectra.t

x = prep.savgol(x, 3, 1, 0)
x.t = prep.savgol(x.t, 3, 1, 0)

cat('\n\n### A. Checking valibration and prediction ###\n')
readline('Press Enter to continue...')

cat('\n1. Calibaration and CV, nobj << nvar\n')
plsmodel = pls(x[1:9, ], y[1:9], cv = 10)
summary(plsmodel)
plot(plsmodel)
readline('Press Enter to continue...')

cat('\n2. Calibaration and CV, nvar << nobj\n')
plsmodel = pls(x[, 1:9], y, cv = 10)
summary(plsmodel)
plot(plsmodel)
readline('Press Enter to continue...')

cat('\n3. Calibaration and CV, small nvar and nobj\n')
plsmodel = pls(x[1:9, 1:9], y[1:9], cv = 1)
summary(plsmodel)
plot(plsmodel)
readline('Press Enter to continue...')

cat('\n4. Calibaration and CV plus prediction for small set\n')
plsmodel = pls(x, y, cv = 1, x.test = x.t[1:4, ], y.test = y.t[1:4])
res = predict(plsmodel, x.t[1:10, ], y.t[1:10])
summary(res)
plot(res)
readline('Press Enter to continue...')

cat('\n5. Calibaration and prediction, small nvar and nobj\n')
plsmodel = pls(x[1:9, 1:4], y[1:9], cv = 1)
summary(plsmodel)
plot(plsmodel)
res = predict(plsmodel, x.t[1:9, 1:4], y[1:9])
summary(res)
plot(res)
readline('Press Enter to continue...')

cat('\n\n### B. Checking calibration results ###\n')
readline('Press Enter to continue...')

plsmodel = pls(x, y, cv = 2, x.test = x.t, y.test = y.t, ncomp = 7)
plsmodel = selectCompNum(plsmodel, 3)

plspred = predict(plsmodel, x.t)


cat('\n1. Print and summary\n')
print(plsmodel$calres)
summary(plsmodel$calres)
summary(plsmodel$calres, ncomp = 4)
readline('Press Enter to continue...')

cat('\n2. Prediction plots\n')
par(mfrow = c(2, 2))
plotPredictions(plsmodel$calres)
plotPredictions(plsmodel$calres, colmap = 'gray')
plotPredictions(plsmodel$calres, col = 'blue', 
                xlab = 'C2, m', 
                ylab = 'C2, p')
plot.regres(plsmodel$calres, cgroup = conc.c[, ny], 
            colmap = c('red', 'green'), 
            show.labels = T,
            show.line = F)
readline('Press Enter to continue...')

cat('\n3. Y residuals plots\n')
par(mfrow = c(2, 2))
plotYResiduals(plsmodel$calres, cgroup = y)
plotYResiduals(plsmodel$calres, colmap = 'gray')
plotYResiduals(plsmodel$calres, show.line = T, show.labels = T)
plotYResiduals(plsmodel$calres, ncomp = 1, col = 'red', show.line = F, show.labels = T)
readline('Press Enter to continue...')

cat('\n4. RMSE plots\n')
par(mfrow = c(2, 2))
plotRMSE(plsmodel$calres)
plotRMSE(plsmodel$calres, colmap = 'gray')
plotRMSE(plsmodel$calres, col = 'red', show.labels = T)
plotRMSE(plsmodel$calres, col = 'red', type = 'h', show.labels = T)
readline('Press Enter to continue...')

cat('\n5. X scores plots\n')
par(mfrow = c(2, 2))
plotXScores(plsmodel$calres)
plotXScores(plsmodel$calres, c(1, 3))
plotXScores(plsmodel$calres, colmap = c('red', 'blue'), cgroup = y)
plotXScores(plsmodel$calres, col = 'orange', show.labels = T)
readline('Press Enter to continue...')

cat('\n6. XY scores plots\n')
par(mfrow = c(2, 2))
plotXYScores(plsmodel$calres)
plotXYScores(plsmodel$calres, 2)
plotXYScores(plsmodel$calres, 2, colmap = c('red', 'blue'), cgroup = y)
plotXYScores(plsmodel$calres, 3, col = 'orange', show.labels = T)
readline('Press Enter to continue...')

cat('\n7. X residuals plots\n')
par(mfrow = c(2, 2))
plotXResiduals(plsmodel$calres)
plotXResiduals(plsmodel$calres, ncomp = 2, col = 'red', show.labels = T)
plotXResiduals(plsmodel$calres, cgroup = y)
plotXResiduals(plsmodel$calres, ncomp = 3)
readline('Press Enter to continue...')

cat('\n8. X variance plots\n')
par(mfrow = c(2, 2))
plotXVariance(plsmodel$calres)
plotXVariance(plsmodel$calres, type = 'h', colmap = 'gray', show.labels = T)
plotXCumVariance(plsmodel$calres, col = 'red', show.labels = T)
plotXCumVariance(plsmodel$calres, col = 'red', show.labels = T, type = 'h')
readline('Press Enter to continue...')

cat('\n9. Y variance plots\n')
par(mfrow = c(2, 2))
plotYVariance(plsmodel$calres)
plotYVariance(plsmodel$calres, type = 'h', colmap = 'gray', show.labels = T)
plotYCumVariance(plsmodel$calres, col = 'red', show.labels = T)
plotYCumVariance(plsmodel$calres, type = 'h')
readline('Press Enter to continue...')

cat('\n10. Summary plots\n')
plot(plsmodel$calres)
readline('Press Enter to continue...')

cat('\n\n### C. Checking cross-validation results ###')
readline('Press Enter to continue...')

cat('\n1. Print and summary\n')
print(plsmodel$cvres)
summary(plsmodel$cvres)
summary(plsmodel$cvres, ncomp = 4)
readline('Press Enter to continue...')

cat('\n2. Prediction plots\n')
par(mfrow = c(2, 2))
plotPredictions(plsmodel$cvres)
plotPredictions(plsmodel$cvres, colmap = 'gray')
plotPredictions(plsmodel$cvres, col = 'blue', 
                xlab = 'C2, measured', 
                ylab = 'C2, predicted')
plot.regres(plsmodel$cvres, cgroup = y, 
            colmap = c('red', 'green'), 
            show.labels = T,
            show.line = F)
readline('Press Enter to continue...')

cat('\n3. Y residuals plots\n')
par(mfrow = c(2, 2))
plotYResiduals(plsmodel$cvres, cgroup = y)
plotYResiduals(plsmodel$cvres, colmap = 'gray')
plotYResiduals(plsmodel$cvres, show.line = T, show.labels = T)
plotYResiduals(plsmodel$cvres, ncomp = 1, col = 'red', show.line = F, show.labels = T)
readline('Press Enter to continue...')

cat('\n4. RMSE plots\n')
par(mfrow = c(2, 2))
plotRMSE(plsmodel$cvres)
plotRMSE(plsmodel$cvres, colmap = 'gray')
plotRMSE(plsmodel$cvres, col = 'red', show.labels = T)
plotRMSE(plsmodel$cvres, col = 'red', type = 'h')
readline('Press Enter to continue...')

cat('\n5. X scores plots\n')
par(mfrow = c(1, 1))
plotXScores(plsmodel$cvres)
readline('Press Enter to continue...')

cat('\n6. XY scores plots\n')
par(mfrow = c(1, 1))
plotXYScores(plsmodel$cvres)
readline('Press Enter to continue...')

cat('\n7. X residuals plots\n')
par(mfrow = c(2, 2))
plotXResiduals(plsmodel$cvres)
plotXResiduals(plsmodel$cvres, ncomp = 2, col = 'red', show.labels = T)
plotXResiduals(plsmodel$cvres, cgroup = y)
plotXResiduals(plsmodel$cvres, ncomp = 3)
readline('Press Enter to continue...')

cat('\n8. X variance plots\n')
par(mfrow = c(2, 2))
plotXVariance(plsmodel$cvres)
plotXVariance(plsmodel$cvres, type = 'h', colmap = 'gray', show.labels = T)
plotXCumVariance(plsmodel$cvres, col = 'red', show.labels = T)
plotXCumVariance(plsmodel$cvres, col = 'red', show.labels = T, type = 'h')
readline('Press Enter to continue...')

cat('\n9. Y variance plots\n')
par(mfrow = c(2, 2))
plotYVariance(plsmodel$cvres)
plotYVariance(plsmodel$cvres, type = 'h', colmap = 'gray', show.labels = T)
plotYCumVariance(plsmodel$cvres, col = 'red', show.labels = T)
plotYCumVariance(plsmodel$cvres, type = 'h')
readline('Press Enter to continue...')

cat('\n10. Check summary plot for cal results\n')
plot(plsmodel$calres)
readline('Press Enter to continue...')

cat('\n\n### D. Checking prediction results ###\n')
readline('Press Enter to continue...')

cat('\n1. Print and summary\n')
print(plspred)
summary(plspred)
readline('Press Enter to continue...')

cat('\n2. Prediction plots\n')
par(mfrow = c(2, 2))
plotPredictions(plspred)
plotPredictions(plspred, ncomp = 2)
plotPredictions(plspred, ncomp = 2, cgroup = y.t)
plotPredictions(plspred, ncomp = 2, col = 'red', pch = 18, show.labels = T)
readline('Press Enter to continue...')

cat('\n3. Y residuals plots\n')
plotYResiduals(plspred)
readline('Press Enter to continue...')

cat('\n\n### E. Checking plsmodel methods ###\n')
readline('Press Enter to continue...')

cat('\n1. Print and summary\n')
print(plsmodel)
summary(plsmodel)
summary(plsmodel, ncomp = 6)
readline('Press Enter to continue...')

cat('\n2. Prediction plots\n')
par(mfrow = c(2, 2))
plotPredictions(plsmodel)
plotPredictions(plsmodel, ncomp = 2, colmap = 'gray', show.legend = F)
plotPredictions(plsmodel, colmap = c('blue', 'red', 'orange'), 
                xlab = 'C2, measured', 
                ylab = 'C2, predicted',
                main = 'Super predictions',
                ncomp = 1,
                pch = c(16, 17, 18),
                legend.position = 'bottomright')
plotPredictions(plsmodel, colmap = c('red', 'green', 'brown'), 
                show.labels = T,
                show.lines = F)
readline('Press Enter to continue...')

cat('\n3. Y residuals plots\n')
par(mfrow = c(2, 2))
plotYResiduals(plsmodel)
plotYResiduals(plsmodel, ncomp = 1, colmap = 'gray')
plotYResiduals(plsmodel, show.line = T, show.labels = T, ncomp = 2, 
               xlab = 'Values', ylab = 'Residuals', main = 'Plot residuals')
plotYResiduals(plsmodel, ncomp = 1, colmap = c('red', 'blue'), 
               show.line = F, show.labels = T)
readline('Press Enter to continue...')

cat('\n4. Regression coefficients plots\n')
par(mfrow = c(2, 2))
plotRegcoeffs(plsmodel)
plotRegcoeffs(plsmodel, ncomp = 1, colmap = 'gray')
plotRegcoeffs(plsmodel, col = 'red', show.labels = T, show.line = F)
plotRegcoeffs(plsmodel, col = 'red', type = 'h')
readline('Press Enter to continue...')

cat('\n4. RMSE plots\n')
par(mfrow = c(2, 2))
plotRMSE(plsmodel)
plotRMSE(plsmodel, colmap = 'gray', show.legend = F)
plotRMSE(plsmodel, colmap = c('red', 'blue'), type = 'l', lty = c(1, 2, 3))
plotRMSE(plsmodel, type = 'h', show.grid = F)
readline('Press Enter to continue...')

cat('\n5. X scores plots\n')
par(mfrow = c(2, 2))
plotXScores(plsmodel)
plotXScores(plsmodel, c(1, 3), show.axes = F, show.legend = F)
plotXScores(plsmodel, colmap = c('red', 'red'), pch = c(15, 18))
plotXScores(plsmodel, colmap = 'gray', show.labels = T)
readline('Press Enter to continue...')

cat('\n6. XY scores plots\n')
par(mfrow = c(2, 2))
plotXYScores(plsmodel)
plotXYScores(plsmodel, 2, colmap = 'gray', show.legend = F, show.axes = F)
plotXYScores(plsmodel, 3, colmap = c('red', 'blue'), pch = c(16, 15))
plotXYScores(plsmodel, 3,  show.labels = T)
readline('Press Enter to continue...')

cat('\n7. X residuals plots\n')
par(mfrow = c(2, 2))
plotXResiduals(plsmodel)
plotXResiduals(plsmodel, ncomp = 2, colmap = c('red', 'orange'), show.labels = T)
plotXResiduals(plsmodel, colmap = 'gray')
plotXResiduals(plsmodel, ncomp = 1, legend.position = 'top',
               xlab = 'T res', ylab = 'Q res', main = 'Residuals for X')
readline('Press Enter to continue...')

cat('\n8. X variance plots\n')
par(mfrow = c(2, 2))
plotXVariance(plsmodel)
plotXVariance(plsmodel, type = 'h', colmap = 'gray', show.labels = T)
plotXCumVariance(plsmodel, col = 'red', show.labels = T)
plotXCumVariance(plsmodel, col = 'red', show.labels = T, type = 'h')
readline('Press Enter to continue...')

cat('\n9. Y variance plots\n')
par(mfrow = c(2, 2))
plotYVariance(plsmodel)
plotYVariance(plsmodel, type = 'h', colmap = 'gray', show.labels = T)
plotYCumVariance(plsmodel, col = 'red', show.labels = T)
plotYCumVariance(plsmodel, type = 'h')
readline('Press Enter to continue...')

cat('\n10. Summary plots\n')
par(mfrow = c(1, 2))
plot(plsmodel)
readline('Press Enter to continue...')

cat('\n11. Selectivity ratio plots\n')
par(mfrow = c(2, 2))
plotSelectivityRatio(plsmodel)
plotSelectivityRatio(plsmodel, ncomp = 2)
plotSelectivityRatio(plsmodel, ncomp = 4, type = 'h')
plotSelectivityRatio(plsmodel, ncomp = 1, color = 'red', lty = 2)
readline('Press Enter to continue...')

cat('\n12. VIP scores plots\n')
par(mfrow = c(2, 2))
plotVIPScores(plsmodel)
plotVIPScores(plsmodel, ncomp = 2)
plotVIPScores(plsmodel, ncomp = 4, type = 'h')
plotVIPScores(plsmodel, ncomp = 1, col = 'red', lty = 2)
readline('Press Enter to continue...')

cat('\n13. Jack-knifing regcoeffs plots\n')
plsmodel = pls(x, y, coeffs.ci='jk')
plsmodel = selectCompNum(plsmodel, 3)
par(mfrow = c(2, 2))
plotRegcoeffs(plsmodel, show.ci = T)
plotRegcoeffs(plsmodel, show.ci = T, ncomp = 2)
plotRegcoeffs(plsmodel, show.ci = T, ncomp = 2, type = 'h')
plotRegcoeffs(plsmodel, show.ci = T, ncomp = 1, col = c('red', 'green', 'green'), lwd = 2)
readline('Press Enter to continue...')
