# 1. One response case
context('Regres: calculation of performance statistics for one response')

## number of objects, components and y-variables
nobj = 30
ncomp = 3
nresp = 1

## create x and y matrices
x = matrix(rnorm(nobj), nrow = nobj, ncol = nresp)
y = x * 2 + 1 + rnorm(nobj, 0, 0.1)

## make nr x nc x nv array for yp with best performance at nc = 2
yp = array(
   c(
      y + rnorm(nobj, 0, 0.50),
      y + rnorm(nobj, 0, 0.25),
      y + rnorm(nobj, 0, 0.75) + 0.1 # added small bias
   ),
   dim = c(nobj, ncomp, nresp)
)
dimnames(yp) = list(
   paste0('O', 1:nobj),
   paste0('Comp ', 1:ncomp),
   paste0('Y', 1:nresp)
)

## compute R2 manually - hardcore: no loops, no apply functions

ytot = sum((y - mean(y))^2)
r2 = 1 - c(
   sum((y - yp[, 1, ])^2), 
   sum((y - yp[, 2, ])^2), 
   sum((y - yp[, 3, ])^2) 
) / ytot
r2 = matrix(r2, nrow = nresp, ncol = ncomp)

## compute RMSE manually
rmse = sqrt(
   c(
      sum((y - yp[, 1, ])^2), 
      sum((y - yp[, 2, ])^2), 
      sum((y - yp[, 3, ])^2) 
   ) / nobj
)
rmse = matrix(rmse, nrow = nresp, ncol = ncomp)

## compute bias manually
bias = c(
   mean(y - yp[, 1, ]), 
   mean(y - yp[, 2, ]), 
   mean(y - yp[, 3, ])
)
bias = matrix(bias, nrow = nresp, ncol = ncomp)

## create regres object
res = regres(yp, y, 2)

## run tests
test_that("R2 is computed correctly", {
   expect_equal(dim(res$r2), c(nresp, ncomp))
   expect_lte(sum(abs(res$r2 - r2)), 2 * .Machine$double.eps)
})

test_that("RMSE is computed correctly", {
   expect_equal(dim(res$rmse), c(nresp, ncomp))
   expect_lte(sum(abs(res$rmse - rmse)), 2 * .Machine$double.eps)
})

test_that("Bias is computed correctly", {
   expect_equal(dim(res$bias), c(nresp, ncomp))
   expect_lte(sum(abs(res$bias - bias)), 2 * .Machine$double.eps)
})


# 2. Several responses test
context('Regres: calculation of performance statistics for several responses')

## number of objects, components and y-variables
nobj = 30
ncomp = 2
nresp = 3

## create x and y matrices
x = matrix(rnorm(nobj), nrow = nobj, ncol = 1)
y = matrix(
   c(
      x *  2 + 1 + rnorm(nobj, 0, 0.1),
      x * -2 - 1 + rnorm(nobj, 0, 0.2),
      x *  3 - 1 + rnorm(nobj, 0, 0.2)
   ),
   nrow = nobj, ncol = nresp
)

## make nr x nc x nv array for yp with best performance at nc = 2
yp = array(
   c(
      y[, 1] + rnorm(nobj, 0, 0.75) + 1,
      y[, 1] + rnorm(nobj, 0, 0.25),
      y[, 2] + rnorm(nobj, 0, 0.50) + 1,
      y[, 2] + rnorm(nobj, 0, 0.75),
      y[, 3] + rnorm(nobj, 0, 0.75) + 1,
      y[, 3] + rnorm(nobj, 0, 0.25)
   ),
   dim = c(nobj, ncomp, nresp)
)

dimnames(yp) = list(
   paste0('O', 1:nobj),
   paste0('Comp ', 1:ncomp),
   paste0('Y', 1:nresp)
)

## compute R2 manually
ytot = c( 
   sum((y[, 1] - mean(y[, 1]))^2),
   sum((y[, 2] - mean(y[, 2]))^2),
   sum((y[, 3] - mean(y[, 3]))^2)
)

r2 = 1 - c(
   sum((y[, 1] - yp[, 1, 1])^2) / ytot[1], 
   sum((y[, 2] - yp[, 1, 2])^2) / ytot[2], 
   sum((y[, 3] - yp[, 1, 3])^2) / ytot[3], 
   sum((y[, 1] - yp[, 2, 1])^2) / ytot[1], 
   sum((y[, 2] - yp[, 2, 2])^2) / ytot[2], 
   sum((y[, 3] - yp[, 2, 3])^2) / ytot[3]
) 

r2 = matrix(r2, nrow = nresp, ncol = ncomp)

## compute RMSE manually
rmse = sqrt(
   c(
      sum((y[, 1] - yp[, 1, 1])^2), 
      sum((y[, 2] - yp[, 1, 2])^2), 
      sum((y[, 3] - yp[, 1, 3])^2), 
      sum((y[, 1] - yp[, 2, 1])^2), 
      sum((y[, 2] - yp[, 2, 2])^2), 
      sum((y[, 3] - yp[, 2, 3])^2)
   ) / nobj
)
rmse = matrix(rmse, nrow = nresp, ncol = ncomp)

## compute bias manually
bias = c(
   mean(y[, 1] - yp[, 1, 1]), 
   mean(y[, 2] - yp[, 1, 2]), 
   mean(y[, 3] - yp[, 1, 3]), 
   mean(y[, 1] - yp[, 2, 1]), 
   mean(y[, 2] - yp[, 2, 2]), 
   mean(y[, 3] - yp[, 2, 3])
)
bias = matrix(bias, nrow = nresp, ncol = ncomp)

## create regres object
res = regres(yp, y, 2)

## run tests
test_that("R2 is computed correctly", {
   expect_equal(dim(res$r2), c(nresp, ncomp))
   expect_lte(sum(abs(res$r2 - r2)), 2 * .Machine$double.eps)
})

test_that("RMSE is computed correctly", {
   expect_equal(dim(res$rmse), c(nresp, ncomp))
   expect_lte(sum(abs(res$rmse - rmse)), 2 * .Machine$double.eps)
})

test_that("Bias is computed correctly", {
   expect_equal(dim(res$bias), c(nresp, ncomp))
   expect_lte(sum(abs(res$bias - bias)), 2 * .Machine$double.eps)
})
