library('mdatools')
data(simdata)

ny = 3

y = simdata$conc.c[, ny]
x = simdata$spectra.c

y.t = simdata$conc.t[, ny]
x.t = simdata$spectra.t

model = pls(x, y, ncomp = 10)
model = selectCompNum(model, 3)

# 1. Y and X variance
par(mfrow = c(2, 2))
plotXVariance(model)
plotXCumVariance(model)
plotYVariance(model)
plotYCumVariance(model)
par(mfrow = c(1, 1))

# 1. RMSE, predictions and residuals
par(mfrow = c(2, 2))
plotRMSE(model)
plotPredictions(model, ncomp = 3)
plotXResiduals(model, ncomp = 3)
plotYResiduals(model, ncomp = 10)
par(mfrow = c(1, 1))

# N. Selectivity ratio
par(mfrow = c(2, 2))
plotSelectivityRatio(model, ncomp = 1)
plotSelectivityRatio(model, ncomp = 2)
plotSelectivityRatio(model)
plotSelectivityRatio(model, ncomp = 5)
par(mfrow = c(1, 1))

# N. VIP scores
par(mfrow = c(2, 2))
plotVIPScores.pls(model, ncomp = 1)
plotVIPScores.pls(model, ncomp = 2)
plotVIPScores.pls(model)
plotVIPScores.pls(model, ncomp = 5)
par(mfrow = c(1, 1))

