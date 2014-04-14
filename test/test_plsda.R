library(mdatools)

calset = iris[seq(1, nrow(iris), 2), 1:4] # calibration set, 3 classes
calset.c = iris[seq(1, nrow(iris), 2), 5]

testset = iris[seq(2, nrow(iris), 2), 1:4] # test set, 3 classes
testset.c = iris[seq(2, nrow(iris), 2), 5] # test set, 3 classes

cc = as.numeric(calset.c)
ct = as.numeric(testset.c)
model = plsda(calset, cc, ncomp = 3, cv = 1, info = 'IRIS data example')
plot(model)
readline('Press enter to continue...')

cc = as.numeric(calset.c)
cc[cc == 1] = 10
cc[cc == 3] = 30
ct = as.numeric(testset.c)
ct[ct == 1] = 10
ct[ct == 3] = 30
model = plsda(calset, cc, ncomp = 3, cv = 1, info = 'IRIS data example')
plot(model)
readline('Press enter to continue...')

cc = calset.c
ct = testset.c
model = plsda(calset, cc, ncomp = 3, cv = 1, info = 'IRIS data example')
plot(model)
readline('Press enter to continue...')

res = predict(model, testset, testset.c)
#res = predict(model, testset[1:50, ], testset.c[1:50, ])
#res = predict(model, tetset)

cat('1. Show print and summary\n')
print(model)
summary(model)
summary(model$calres)
readline('Press enter to continue...')
summary(model, ncomp = 2)
summary(model$calres, nc = 2)
readline('Press enter to continue...')

cat('2. Show print and summary\n')
plot(model)
readline('Press enter to continue...')

cat('3. Show performance for model\n')
par(mfrow = c(2, 2))
plotPerformance(model)
plotPerformance(model, param = 'sensitivity')
plotSensitivity(model, nc = 1)
plotMisclassified(model, nc = 2)
readline('Press enter to continue...')

cat('2. Show predictions model\n')
par(mfrow = c(2, 2))
plotPredictions(model)
plotPredictions(model, res = 'cvres', ncomp = 2, nc = 2)
plotPredictions.pls(model)
plotPredictions.pls(model, ncomp = 2, ny = 2)
readline('Press enter to continue...')
