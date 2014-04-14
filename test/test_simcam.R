library(mdatools)

calset = iris[seq(1, nrow(iris), 2), ] # calibration set, 3 classes
tstset = iris[seq(2, nrow(iris), 2), ] # test set, 3 classes
tstset2c = iris[seq(2, nrow(iris) - 50, 2), ] # first 2 classes

se = calset[calset[, 5] == 'setosa', 1:4]
ve = calset[calset[, 5] == 'versicolor', 1:4]
vi = calset[calset[, 5] == 'virginica', 1:4]

# calibration data with reference
cal.data = calset[, 1:4]
cal.c = calset[, 5]

# test data with reference
test.data = tstset[, 1:4]
test.c = tstset[, 5]

# test data with reference for two classes
test2c.data = tstset2c[, 1:4]
test2c.c = tstset2c[, 5]

# test data with reference for with unknown classes
test2cu.data = tstset2c[, 1:4]
test2cu.c = c(matrix('setosa', 1, 25), matrix('none', 1, 10), matrix('Ar', 1, 15))

# make individual models
semodel = simca(se, 'setosa', ncomp = 4, alpha = 0.01, cv = 1)
semodel = selectCompNum(semodel, 1)

# make individual models
vimodel = simca(vi, 'virginica')
vimodel = selectCompNum(vimodel, 3)

# make individual models
vemodel = simca(ve, 'versicolor', ncomp = 4, cv = 5, x.test = test.data, c.test = test.c)
vemodel = selectCompNum(vemodel, 3)

# make group models
model = simcam(list(semodel, vemodel, vimodel))
model2 = simcam(list(semodel, vimodel))

### Make predictions and show plots for 3 classes model ###

pred = predict(model, test.data, test.c)
cpred = predict(model, cal.data, cal.c)
vepred = predict(model, ve)
pred2c = predict(model, test2c.data, test2c.c)
pred2cu = predict(model, test2cu.data, test2cu.c)

cat('1. Show print and summary\n')
print(model)
summary(model)
plot(model)
readline('Press enter to continue...')

cat('2. Show predictions for results\n')
showPredictions(pred)
showPredictions(vepred)
showPredictions(pred2c)
showPredictions(pred2cu)

cat('3. Show predictions for results\n')
par(mfrow = c(3, 2))
plotPredictions(pred)
plotPredictions(cpred)
plotPredictions(vepred)
plotPredictions(pred2c)
plotPredictions(pred2cu)
readline('Press enter to continue...')

cat('4. Show sicrimination power plot')
par(mfrow = c(2, 2))
plotDiscriminationPower(model)
plotDiscriminationPower(model, nc = c(1, 3))
plotDiscriminationPower(model, nc = c(2, 3))
plotDiscriminationPower(model, nc = c(2, 3), show.labels = T, type = 'b')
readline('Press enter to continue...')

cat('5. Show model distance plots')
par(mfrow = c(2, 2))
plotModelDistance(model)
plotModelDistance(model, nc = 2, show.labels = T)
plotModelDistance(model, nc = 3, show.labels = T)
plotModelDistance(model, nc = 3, show.labels = T, type = 'b', col = 'red')
readline('Press enter to continue...')

cat('6. Show Coomans plots')
par(mfrow = c(2, 2))
plotCooman(model)
plotCooman(model, nc = c(1, 3))
plotCooman(model, nc = c(2, 3), show.labels = T)
plotCooman(model, nc = c(2, 3), show.labels = T, show.limits = F)
readline('Press enter to continue...')
