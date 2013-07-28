library('mdatools')
data(Simdata)

pspectra.c = prep.savgol(prep.snv(spectra.c), 5, 3, 0)
pspectra.t = prep.savgol(prep.snv(spectra.t), 5, 3, 0)

odata = list(
      cal = cbind(210:359, t(spectra.c)),
      val = cbind(210:359, t(spectra.t))
)

pdata = list(
   cal = cbind(210:359, t(pspectra.c)),
   val = cbind(210:359, t(pspectra.t))
)

par(mfrow = c(2, 1))
mdaplotg(odata, single.x = T, type = 'l')
mdaplotg(pdata, single.x = T, type = 'l')
par(mfrow = c(1, 1))
