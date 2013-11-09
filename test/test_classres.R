## Define dta

objnames = c(
   'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 
   'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 
   'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8' 
)
classnames = c('A', 'B', 'C')
lvnames = c('LV1', 'LV2', 'LV3')

c.ref = matrix(
   c(
      1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1      
      ), ncol = 3)

c.pred = array(
   c(
      1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      
            
      0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
      
      
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0,      
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1
      
      ), dim = c(24, 3, 3))

p.pred = array(
   c(
      1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
            
      0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
            
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0,      
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1      
   ), dim = c(24, 3, 3))

c.ref = c.ref * 2 - 1;
c.pred = c.pred * 2 - 1
p.pred[p.pred == 0] = runif(length(p.pred[p.pred == 0]), -1.2, -0.0001)
p.pred[p.pred == 1] = runif(length(p.pred[p.pred == 1]), 0.0001, 1.4)

dimnames(c.ref)[[1]] = objnames;
dimnames(c.ref)[[2]] = classnames;

dimnames(c.pred)[[1]] = objnames;
dimnames(c.pred)[[2]] = lvnames;
dimnames(c.pred)[[3]] = classnames;

dimnames(p.pred)[[1]] = objnames;
dimnames(p.pred)[[2]] = lvnames;
dimnames(p.pred)[[3]] = classnames;

## 1. one class one component
lc.ref = c.ref[ , 1, drop = F]
lc.pred = c.pred[ , 1, 1, drop = F]
lp.pred = p.pred[ , 1, 1, drop = F]

# with predictions and reference classes
res = classres(lc.pred, c.ref = lc.ref)
print(res)
summary(res)
par(mfrow = c(1, 1))
plotPredictions(res, show.labels = T)
readline('Press enter to continue...')

# with predictions, probability and reference classes
res = classres(lc.pred, c.ref = lc.ref, p.pred = lp.pred)
print(res)
summary(res)
par(mfrow = c(1, 1))
plotPredictions(res, show.labels = T)
readline('Press enter to continue...')


# with predictions and probability no reference
res = classres(lc.pred, p.pred = lp.pred)
print(res)
summary(res)
par(mfrow = c(1, 1))
plotPredictions(res, show.labels = T)
readline('Press enter to continue...')

# with predictions, no probability and reference
res = classres(lc.pred)
print(res)
summary(res)
par(mfrow = c(1, 1))
plotPredictions(res, show.labels = T)
readline('Press enter to continue...')


## 2. many classes one component
lc.ref = c.ref
lc.pred = c.pred[ , 1, , drop = F]
lp.pred = p.pred[ , 1, , drop = F]

# with predictions and reference classes
res = classres(lc.pred, c.ref = lc.ref)
print(res)
summary(res)
par(mfrow = c(3, 1))
plotPredictions(res, show.labels = T)
plotPredictions(res, nc = 2, show.labels = T)
plotPredictions(res, nc = 3, show.labels = T)
readline('Press enter to continue...')

# with predictions, probability and reference classes
res = classres(lc.pred, c.ref = lc.ref, p.pred = lp.pred)
print(res)
summary(res)
par(mfrow = c(3, 1))
plotPredictions(res, show.labels = T)
plotPredictions(res, nc = 2, show.labels = T)
plotPredictions(res, nc = 3, show.labels = T)
readline('Press enter to continue...')

# with predictions and probability no reference
res = classres(lc.pred, p.pred = lp.pred)
print(res)
summary(res)
par(mfrow = c(3, 1))
plotPredictions(res, show.labels = T)
plotPredictions(res, nc = 2, show.labels = T)
plotPredictions(res, nc = 3, show.labels = T)
readline('Press enter to continue...')

# with predictions, no probability and reference
res = classres(lc.pred)
print(res)
summary(res)
par(mfrow = c(3, 1))
plotPredictions(res, show.labels = T)
plotPredictions(res, nc = 2, show.labels = T)
plotPredictions(res, nc = 3, show.labels = T)
readline('Press enter to continue...')

## 3. one class many components
lc.ref = c.ref[ , 1, drop = F]
lc.pred = c.pred[ , , 1, drop = F]
lp.pred = p.pred[ , , 1, drop = F]

# with predictions and reference classes
res = classres(lc.pred, c.ref = lc.ref)
print(res)
summary(res)
par(mfrow = c(3, 1))
plotPredictions(res, show.labels = T)
plotPredictions(res, ncomp = 2, show.labels = T)
plotPredictions(res, ncomp = 3, show.labels = T)
readline('Press enter to continue...')

# with predictions, probability and reference classes
res = classres(lc.pred, c.ref = lc.ref, p.pred = lp.pred)
print(res)
summary(res)
par(mfrow = c(3, 1))
plotPredictions(res, show.labels = T)
plotPredictions(res, ncomp = 2, show.labels = T)
plotPredictions(res, ncomp = 3, show.labels = T)
readline('Press enter to continue...')

# with predictions and probability no reference
res = classres(lc.pred, p.pred = lp.pred)
print(res)
summary(res)
par(mfrow = c(3, 1))
plotPredictions(res, show.labels = T)
plotPredictions(res, ncomp = 2, show.labels = T)
plotPredictions(res, ncomp = 3, show.labels = T)
readline('Press enter to continue...')

# with predictions, no probability and reference
res = classres(lc.pred)
print(res)
summary(res)
par(mfrow = c(3, 1))
plotPredictions(res, show.labels = T)
plotPredictions(res, ncomp = 2, show.labels = T)
plotPredictions(res, ncomp = 3, show.labels = T)
readline('Press enter to continue...')


## 4. many classes many components
lc.ref = c.ref
lc.pred = c.pred
lp.pred = p.pred

# with predictions and reference classes
res = classres(lc.pred, c.ref = lc.ref)
print(res)
summary(res)
par(mfcol = c(3, 2))
plotPredictions(res, show.labels = T)
plotPredictions(res, ncomp = 2, show.labels = T)
plotPredictions(res, ncomp = 3, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 1, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 2, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 3, show.labels = T)
readline('Press enter to continue...')

# with predictions, probability and reference classes
res = classres(lc.pred, c.ref = lc.ref, p.pred = lp.pred)
print(res)
summary(res)
par(mfcol = c(3, 2))
plotPredictions(res, show.labels = T)
plotPredictions(res, ncomp = 2, show.labels = T)
plotPredictions(res, ncomp = 3, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 1, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 2, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 3, show.labels = T)
readline('Press enter to continue...')

# with predictions and probability no reference
res = classres(lc.pred, p.pred = lp.pred)
print(res)
summary(res)
par(mfcol = c(3, 2))
plotPredictions(res, show.labels = T)
plotPredictions(res, ncomp = 2, show.labels = T)
plotPredictions(res, ncomp = 3, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 1, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 2, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 3, show.labels = T)
readline('Press enter to continue...')

# with predictions, no probability and reference
res = classres(lc.pred)
print(res)
summary(res)
par(mfcol = c(3, 2))
plotPredictions(res, show.labels = T)
plotPredictions(res, ncomp = 2, show.labels = T)
plotPredictions(res, ncomp = 3, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 1, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 2, show.labels = T)
plotPredictions(res, nc = 3, ncomp = 3, show.labels = T)
readline('Press enter to continue...')
