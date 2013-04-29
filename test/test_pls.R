library('mdatools')

data(Simdata)

nobj = nrow(spectra)
nvar = ncol(spectra)

varstep = 1
yvar = 1

X = spectra[1:70, seq(1, nvar, varstep)]
y = conc[1:70, yvar]
# y[5] = y[5] * 2
Xt = spectra[71:nobj, seq(1, nvar, varstep)]
yt = conc[71:nobj, yvar]

plsmodel = pls(X, y, Xt = Xt, yt = yt, autoscale = 1, cv = 0)
plsmodel = pls.selectncomp(plsmodel, 2)

summary(plsmodel)
plot(plsmodel)
#par(mfrow = c(1, 2))
#plotYResiduals(plsmodel, show.labels = T)
#plotXResiduals(plsmodel, show.labels = T)
#par(mfrow = c(2, 2))
#plotRMSE(plsmodel)
plotXYScores(plsmodel, 1, show.labels = T)
#plotRegcoeffs(plsmodel)
plotPredictions(plsmodel, show.labels = T)

#par(mfrow = c(1, 1))
##plotResiduals(plsmodel, ncomp = 1)
#plotXLoadings(plsmodel, 2)
#plotWeights(plsmodel)
#plotYVariance(plsmodel)
