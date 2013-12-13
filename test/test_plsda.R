library(mdatools)

calset = iris[seq(1, nrow(iris), 2), 1:4] # calibration set, 3 classes
calset.c = iris[seq(1, nrow(iris), 2), 5]

testset = iris[seq(2, nrow(iris), 2), 1:4] # test set, 3 classes
testset.c = iris[seq(2, nrow(iris), 2), 5] # test set, 3 classes

model = plsda(calset, calset.c, ncomp = 3, cv = 1, info = 'IRIS data example')

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
