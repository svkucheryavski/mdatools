library(mdatools)
data(Simdata)

spectra = spectra.c
conc = conc.c
wavelength = 201:350

spectra1 = spectra[conc[, 1] < 0.3, ]
spectra2 = spectra[conc[, 1] >= 0.3 & conc[, 1] < 0.65, ]
spectra3 = spectra[conc[, 1] >= 0.65, ]

# make three series for scatter plots
datas1 = cbind(spectra[, 100], spectra[, 95])
colnames(datas1) = c('100 nm', '95 nm')

datas2 = cbind(spectra[, 100], spectra[, 110])
colnames(datas2) = c('100 nm', '110 nm')

datas3 = cbind(spectra[, 100], spectra[, 115])
colnames(datas3) = c('100 nm', '115 nm')

# make three series for line plots
datal = cbind(wavelength, t(spectra))
colnames(datal)[1] = 'Wavelength'
colnames(datal)[2] = 'Absorbance'

datal1 = cbind(wavelength, t(spectra1))
colnames(datal1)[1] = 'Wavelength'
colnames(datal1)[2] = 'Absorbance'

datal2 = cbind(wavelength, t(spectra2))
colnames(datal2)[1] = 'Wavelength'
colnames(datal2)[2] = 'Absorbance'

datal3 = cbind(wavelength, t(spectra3))
colnames(datal2)[1] = 'Wavelength'
colnames(datal2)[2] = 'Absorbance'

# make three groups for scatterg plots using list, 
# matrix with one x column and multiple x columns 
datags1 = list(
   first = datas1,
   second = datas2,
   third = datas3
)

datags2 = cbind(datas1, datas2, datas3)

datags3 = cbind(datas1, datas2[, 2], datas3[, 2])

legendgs = c('95nm', '110nm', '115nm')

# make three groups for lineg plots using list, 
# matrix with one x column and multiple x columns 
datagl1 = list(
   first = datal1,
   second = datal2,
   third = datal3
)

datagl2 = cbind(datal1[, 1:2], datal2[, 1:2], datal3[, 1:2])

datagl3 = cbind(datal1[, 1:2], datal2[, 2], datal3[, 3])

legendgl = c('< 0.3', '< 0.65', '< 1')

### Test conventional plots ###

cat('\n1. scatter plots - labels and lines\n')
#png(filename = '~/Desktop/Fig2.png', width = 2150, height = 1600, pointsize = 9, res = 300, bg = 'white')
par(mfrow = c(2, 2))
mdaplot(datas1)
mdaplot(datas1, show.labels = T)
mdaplot(datas2, show.labels = T, show.lines = c(0.1, 0.03))
mdaplot(datas2, show.labels = T, show.lines = c(NA, 0.03))
#dev.off()
readline('Press Enter to continue...')

cat('\n2. scatter plots - cgroup and colormaps\n')
#png(filename = '~/Desktop/Fig2.png', width = 2150, height = 1600, pointsize = 9, res = 300, bg = 'white')
par(mfrow = c(2, 2))
mdaplot(datas1, cgroup = conc[, 1])
mdaplot(datas2, cgroup = conc[, 1], colmap = 'gray')
mdaplot(datas3, cgroup = conc[, 1], colmap = c('red', 'green'))
mdaplot(datas1, cgroup = conc[, 1], colmap = rainbow(4))
#dev.off()
res = tryCatch( mdaplot(datas1, cgroup = conc[, 1], colmap = c('#100000', '#400000', 'fignya')), 
          error = function(e) FALSE)
cat(sprintf('If color value for colmap is not valid, will function works: %s', res))
readline('Press Enter to continue...')

cat('\n3. Test line plots - default settings\n')
#png(filename = '~/Desktop/Fig2.png', width = 2150, height = 1600, pointsize = 9, res = 300, bg = 'white')
par(mfrow = c(2, 2))
mdaplot(datal1, type = 'l', col = 'red')
mdaplot(datal, type = 'l', cgroup = conc[, 1])
mdaplot(datal1[, 1:5], type = 'l', lty = 3, lwd = 3)
mdaplot(datal2, type = 'l', cgroup = conc[, 1], colmap = 'gray')
#dev.off()
readline('Press Enter to continue...')

cat('\n4. Test linescatter plots - default settings\n')
#png(filename = '~/Desktop/Fig2.png', width = 2150, height = 1600, pointsize = 9, res = 300, bg = 'white')
par(mfrow = c(2, 2))
mdaplot(datal1[30:40, 1:2], type = 'b', col = 'red', show.labels = T)
mdaplot(datal[30:40, 1:5], type = 'b', cgroup = c(1, 2, 3))
mdaplot(datal1[30:40, 1:5], type = 'b', lty = 2, pch = 17)
mdaplot(datal2[30:40, 1:5], type = 'b', colmap = 'gray')
#dev.off()
readline('Press Enter to continue...')
### Test group plots ###

cat('\n3. Test scatter group plots - cgroup and colormaps\n')
par(mfrow = c(2, 2))
mdaplotg(datags1)
mdaplotg(datags2, single.x = F, legend = legendgs, colmap = c('red', 'blue'))
mdaplotg(datags3, legend = legendgs, xlab = 'wave', ylab = 'y wave', colmap = 'gray')
mdaplotg(datags3, legend = legendgs, 
                  xlab = 'x wave', 
                  ylab = 'y wave',
                  main = 'Spectra',
                  pch = c(15, 16, 17), 
                  col = rainbow(10),
                  show.grid = F,
                  show.lines = c(0.04, 0.04))
readline('Press Enter to continue...')

cat('\n4. Test scatter group plots - legend position\n')
par(mfrow = c(3, 2))
mdaplotg(datags1, legend = legendgs, legend.position = 'topright')
mdaplotg(datags1, legend = legendgs, legend.position = 'top')
mdaplotg(datags1, legend = legendgs, legend.position = 'topleft')
mdaplotg(datags1, legend = legendgs, legend.position = 'bottomleft')
mdaplotg(datags1, legend = legendgs, legend.position = 'bottom')
mdaplotg(datags1, legend = legendgs, legend.position = 'bottomright')
readline('Press Enter to continue...')

cat('\n6. Test line group plots - data organizing\n')
par(mfrow = c(2, 2))
mdaplotg(datagl1, type = 'l', legend = legendgl, legend.position = 'top')
mdaplotg(datagl2, type = 'l', legend = legendgl, single.x = F)
mdaplotg(datagl3, type = 'l', legend = legendgl)
mdaplotg(datagl1, type = 'l', legend = NULL)
readline('Press Enter to continue...')

cat('\n7. Test different markers and line types for groups\n')
par(mfrow = c(2, 2))
mdaplotg(datagl1, type = 'l', legend = legendgl, show.grid = F, 
         colmap = rainbow(4), xlab = 'X variable')
mdaplotg(datagl3, type = 'l', lty = c(1, 2, 3), legend = legendgl, show.grid = T)
mdaplotg(datagl3, type = c('l', 'p', 'b'), lty = c(1, 2, 3), pch = c(16, 17, 18),
         colmap = c('#000000', '#000000'), legend = legendgl)
mdaplotg(datagl3[30:40, ], type = 'b', lty = c(3, 2, 1), colmap = c('red', 'orange', 'pink'), 
               pch = c(18, 17, 16), legend = legendgl, show.labels = T)
readline('Press Enter to continue...')


