# new tests on top

# Tests written for 0.9.1

## test that maximum number of components is estimated correctly when
## number of rows is small and CV is used

context('PCA: estimation of maximum number of components')

data(people)
x.cal = people[1:5, ]

m01 = pca(x.cal)
m02 = pca(x.cal, cv = 1)
m03 = pca(x.cal, cv = 2)

test_that("Estimation of maximum number of components (w/wo CV) is correct", {
   expect_equal(m01$ncomp, 4)
   expect_equal(m02$ncomp, 3)
   expect_equal(m03$ncomp, 1)
})

m01 = pca(x.cal, ncomp = 5)
m02 = pca(x.cal, ncomp = 5, cv = 1)
m03 = pca(x.cal, ncomp = 5, cv = 2)

test_that("Adjustment of number of components (w/wo CV) when proposed value too large is correct", {
   expect_equal(m01$ncomp, 4)
   expect_equal(m02$ncomp, 3)
   expect_equal(m03$ncomp, 1)
})

m01 = pca(x.cal, ncomp = 1)
m02 = pca(x.cal, ncomp = 1, cv = 1)
m03 = pca(x.cal, ncomp = 1, cv = 2)

test_that("No adjustment of number of components (w/wo CV) when it is small enough", {
   expect_equal(m01$ncomp, 1)
   expect_equal(m02$ncomp, 1)
   expect_equal(m03$ncomp, 1)
})

# Tests written for 0.9.0
data(iris)

## compute SVD decomposition for Setosa samples of Iris dataset
x.cal1 = scale(iris[1:50, 1:4], center = T, scale = F)
x.cal2 = scale(iris[1:50, 1:4], center = T, scale = T)
loads1 = svd(x.cal1)$v
loads2 = svd(x.cal2)$v

## check loadings obtained using conventional algorithms
context('PCA: loadings by conventional algorithms (max nPCs)')

scores1 = x.cal1 %*% loads1
scores2 = x.cal2 %*% loads2

## centered only
m01 = pca(x.cal1)
m11 = pca(x.cal1, method = 'svd')
m21 = pca(x.cal1, method = 'nipals')

## scaled
m02 = pca(x.cal2, scale = T)
m12 = pca(x.cal2, scale = T, method = 'svd')
m22 = pca(x.cal2, scale = T, method = 'nipals')

## test for loading size
test_that("Loadings have correct size (centered)", {
   expect_equal(nrow(m01$loadings), 4)
   expect_equal(nrow(m11$loadings), 4)
   expect_equal(nrow(m21$loadings), 4)
   expect_equal(ncol(m01$loadings), 4)
   expect_equal(ncol(m11$loadings), 4)
   expect_equal(ncol(m21$loadings), 4)
})

test_that("Loadings have correct size (autoscaled)", {
   expect_equal(nrow(m02$loadings), 4)
   expect_equal(nrow(m12$loadings), 4)
   expect_equal(nrow(m22$loadings), 4)
   expect_equal(ncol(m02$loadings), 4)
   expect_equal(ncol(m12$loadings), 4)
   expect_equal(ncol(m22$loadings), 4)
})

# test for loading values
test_that("Loadings are computed correctly (centered)", {
   expect_lt(sum(abs(m01$loadings) - abs(loads1)), 0.001)
   expect_lt(sum(abs(m11$loadings) - abs(loads1)), 0.001)
   expect_lt(sum(abs(m21$loadings) - abs(loads1)), 0.001)
})

test_that("Loadings are computed correctly (autoscaled)", {
   expect_lt(sum(abs(m02$loadings) - abs(loads2)), 0.001)
   expect_lt(sum(abs(m12$loadings) - abs(loads2)), 0.001)
   expect_lt(sum(abs(m22$loadings) - abs(loads2)), 0.001)
})

context('PCA: loadings by conventional algorithms (2 PCs)')

ncomp = 2
m01 = pca(x.cal1, ncomp)
m11 = pca(x.cal1, ncomp, method = 'svd')
m21 = pca(x.cal1, ncomp, method = 'nipals')

test_that("Loadings have correct size", {
   expect_equal(nrow(m01$loadings), 4)
   expect_equal(nrow(m11$loadings), 4)
   expect_equal(nrow(m21$loadings), 4)
   expect_equal(ncol(m01$loadings), 2)
   expect_equal(ncol(m11$loadings), 2)
   expect_equal(ncol(m21$loadings), 2)
})

test_that("Loadings are computed correctly", {
   expect_lt(sum(abs(m01$loadings) - abs(loads1[, 1:2])), 0.001)
   expect_lt(sum(abs(m11$loadings) - abs(loads1[, 1:2])), 0.001)
   expect_lt(sum(abs(m21$loadings) - abs(loads1[, 1:2])), 0.001)
})

