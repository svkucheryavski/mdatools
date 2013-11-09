library(mdatools)

calset = iris[seq(1, nrow(iris), 2), ]
tstset = iris[seq(2, nrow(iris), 2), ]

se = calset[calset[, 5] == 'setosa', 1:4]
ve = calset[calset[, 5] == 'versicolor', 1:4]
vi = calset[calset[, 5] == 'virginica', 1:4]

test.data = tstset[, 1:4]
test.c = c(matrix(1, 1, 25), matrix(2, 1, 25), matrix(3, 1, 25))


# Select which model to calculate
option = 3

if (option == 1)
{
   # Setosa model with cross-validation
   model = simca(se, 'Se', ncomp = 4, alpha = 0.01, cv = 1)
   model = selectCompNum(model, 1)
   pred = predict(model, test.data, test.c == 1)   
} else if (option == 2) {  
   # Virginica model with test set validation (no CV)
   model = simca(vi, 'Vi', test.data = test.data[test.c == 3, ])
   model = selectCompNum(model, 3)
   pred = predict(model, test.data, test.c == 3)   
} else {
   # Versicolor model with 5 segments CV and test set
   model = simca(ve, 'Ve', ncomp = 4, cv = 5, test.data = test.data, test.c = test.c == 2)
   model = selectCompNum(model, 3)
   pred = predict(model, test.data, test.c == 2)   
}

cat('1. Show print and summary')
print(model)
summary(model)
plot(model)
readline('Press enter to continue...')

cat('2. Show predictions plot for model')
par(mfrow = c(2, 1))
plotPredictions(model);
plotPredictions(model, ncomp = 1, type = 'p', show.labels = T)
readline('Press enter to continue...')

cat('3. Show predictions plot for results')
par(mfrow = c(2, 1))
plotPredictions(model$calres);
plotPredictions(pred, ncomp = 1, type = 'p', show.labels = T)
readline('Press enter to continue...')

cat('4. Show residuals plot for model and results')
par(mfrow = c(2, 2))
plotResiduals(model)
plotResiduals(model, ncomp = 1, show.labels = T)
plotResiduals(pred, ncomp = 3)
plotResiduals(model$calres, ncomp = 1, show.labels = T)
readline('Press enter to continue...')

cat('5. Show sensitivity plot for model and results')
par(mfrow = c(2, 2))
plotSensitivity(model, legend.position = 'bottomright')
plotSensitivity(model, type = 'b', show.labels = T, legend.position = 'bottomright')
plotSensitivity(model$calres, legend.position = 'bottomright')
plotPerformance(pred, type = 'b', show.labels = T, legend.position = 'bottomright')
readline('Press enter to continue...')

cat('6. Show modelling power and loadings plot')
par(mfrow = c(2, 2))
plotModellingPower(model, show.labels = T)
plotModellingPower(model, ncomp = 1, type = 'b', show.labels = T)
plotLoadings(model, comp = 1, type = 'h', show.labels = T)
plotLoadings(model, comp = c(1, 2), show.labels = T, type = 'p')
readline('Press enter to continue...')

cat('7. Show scores and variance plots')
par(mfrow = c(2, 2))
plotScores(model, show.labels = T)
plotScores(model, comp = 1, type = 'b', show.labels = T)
plotVariance(model, type = 'h', show.labels = T)
plotCumVariance(model, show.labels = T)
readline('Press enter to continue...')
