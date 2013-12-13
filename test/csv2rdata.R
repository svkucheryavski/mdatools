data = read.csv('data/Simdata.csv', row.names = 1, header = T)
data = as.matrix(data)
data = data[, 1:(ncol(data) - 1)]
spectra.c = data[1:100, -(1:3)]
spectra.t = data[101:150, -(1:3)]
conc.c = data[1:100, (1:3)]
conc.t = data[101:150, (1:3)]

save(spectra.c, spectra.t, conc.c, conc.t, file = 'data/Simdata.rdata')