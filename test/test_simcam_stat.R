## Prepare data ##
se = iris[iris[, 5] == 'setosa', 1:4]
se.c = matrix('Se', ncol = 1, nrow = nrow(se))

ve = iris[iris[, 5] == 'versicolor', 1:4]
ve.c = matrix('Ve', ncol = 1, nrow = nrow(ve))

vi = iris[iris[, 5] == 'virginica', 1:4]
vi.c = matrix('Vi', ncol = 1, nrow = nrow(vi))


## Make individual modelsi ##
semodel = simca(se, 'Se', scale = T, cv = 1)
semodel = selectCompNum(semodel, 1)

vimodel = simca(vi, 'Vi', scale = T, cv = 1)
vimodel = selectCompNum(vimodel, 1)

vemodel = simca(ve, 'Ve', scale = T, cv = 1)
vemodel = selectCompNum(vemodel, 1)

## Select case ##
#  1 - 3 classes in model and 3 same classes in prediction data with reference values
#  2 - 3 classes in model and 3 same classes in prediction data but without reference values provided
#  3 - 2 classes in model and 3 classes in prediction data with reference values
#  4 - 3 classes in model and 2 classes in prediction data with reference values
case = 4 

if (case == 1)
{
   model = simcam(list(semodel, vemodel, vimodel))
   calset = rbind(se, ve, vi)
   calset.c = rbind(se.c, ve.c, vi.c)

   res = predict(model, calset, calset.c)
} else if (case == 2)
{
   model = simcam(list(semodel, vemodel, vimodel))
   calset = rbind(se, ve, vi)
   calset.c = rbind(se.c, ve.c, vi.c)

   res = predict(model, calset)
} else if (case == 3)
{
   model = simcam(list(semodel, vimodel))
   calset = rbind(se, ve, vi)
   calset.c = rbind(se.c, ve.c, vi.c)

   res = predict(model, calset, calset.c)
} else if (case == 4)
{
   model = simcam(list(semodel, vemodel, vimodel))
   calset = rbind(se, vi)
   calset.c = rbind(se.c, vi.c)

   res = predict(model, calset, calset.c)
}

## Show plots and results ##

cat('Show plots for SIMCAM model\n')
par(mfrow = c(2, 2))
plotDiscriminationPower(model)
plotDiscriminationPower(model, nc = c(1, 2), show.labels = T)
plotModelDistance(model)
plotModelDistance(model, nc = 2, show.labels = T)
readline('Press Enter key to continue...')


cat('Show summary and print values\n')
print(res)
summary(res)
plot(res)

cat('Show predictions plots\n')
par(mfrow = c(2, 2))
plotPredictions(res)
plotPredictions(res, nc = 2)
plotPredictions(res, nc = 2, show.labels = T)
plotPredictions(res, nc = c(1, 2) )
readline('Press Enter key to continue...')

cat('Show residuals and Coomans plots\n')
par(mfrow = c(2, 2))
plotResiduals(res)
plotResiduals(res, nc = 2, show.labels = T)
plotCooman(res)
plotCooman(res, nc = c(2, 1), show.labels = T)
readline('Press Enter key to continue...')


