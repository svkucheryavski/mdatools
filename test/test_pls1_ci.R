# people data
data(people)

x = people[, -4]
y = people[, 4]

m = pls(x, y, scale = T, coeffs.ci = 'jk', coeffs.alpha = 0.01)
m = selectCompNum(m, 2)
par(mfrow = c(2, 2))
plotRegcoeffs(m, show.ci = T, show.labels = T)
plotRegcoeffs(m, show.ci = T, type = 'h', show.labels = T)
plotRegcoeffs(m, show.ci = T, type = 'b')
plotRegcoeffs(m, show.ci = T, type = 'l')
par(mfrow = c(1, 1))

show(m$coeffs$t.values[, 3, 1])
show(round(m$coeffs$p.values[, 3, 1], 3))


# spectral data
data(simdata)
ny = 3

y = simdata$conc.c[, ny]
x = simdata$spectra.c

m = pls(x, y, coeffs.ci = 'jk', coeffs.alpha = 0.01, cv = 15)
m = selectCompNum(m, 3)
par(mfrow = c(2, 1))
plotRegcoeffs(m, show.ci = T)
plotRegcoeffs(m, show.ci = T, type = 'h', col = c('green'))
par(mfrow = c(1, 1))
plotRegcoeffs(m, show.ci = F, colmap = c('red', 'green'))
