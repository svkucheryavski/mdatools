# new tests on top

## prepare datasets 
data(iris)
cal.ind = c(1:25, 51:75, 101:125)
val.ind = c(26:50, 76:100, 126:150)

Xc = iris[cal.ind, 1:4]
Xv = iris[val.ind, 1:4]

cc.all = iris[cal.ind, 5]
cc.vir = cc.all == 'virginica'
cv.all = iris[val.ind, 5]
cv.vir = cv.all == 'virginica'

# Tests written for 0.9.1

context('PLS-DA: calibration and cross-validation')

## calibrate multiclass model with full CV
m.all = plsda(Xc, cc.all, 3, cv = 1)

## calibrate one class model with full CV
m.vir = plsda(Xc, cc.vir, 3, cv = 1, classname = 'virginica')

# apply models to test set
res11 = predict(m.all, Xv, cv.all)
#res12 = predict(m.all, Xv, cv.vir) # this should give an error
res21 = predict(m.vir, Xv, cv.all)
res22 = predict(m.vir, Xv, cv.vir)
