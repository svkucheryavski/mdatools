library(mdatools)
data(people)

data = people
mvdata = data

# add missing values
mvdata[c(1, 11, 17), 1] = NA
mvdata[3, 4] = NA
mvdata[c(5, 7, 9, 23), 8] = NA
mvdata[c(25, 27, 31, 32), 2] = NA
mvdata[c(4, 2), 12] = NA

# make PCA models
mvdata2 = pca.mvreplace(mvdata, center = T, scale = T)
#res = prep.autoscale(data, scale = T)
#data = res$x

show(cbind(data[, 1], mvdata[, 1], mvdata2[, 1]))

model1 = pca(data, scale = T)
model2 = pca(mvdata2, scale = T)
par(mfrow = c(2, 2))
plotScores(model1, show.labels = T)
plotScores(model2, show.labels = T)
plotLoadings(model1, show.labels = T)
plotLoadings(model2, show.labels = T)
par(mfrow = c(1, 1))

a = matrix(c(1:10, 2*(1:10), 3 * (1:10)), ncol = 3)
a[5, 1] = NA
a[8, 3] = NA


ap = pca.mvreplace(a, center = T, scale = F)

show(cbind(a, round(ap, 2)))



