data(iris)

# compute SVD decomposition for Setosa samples of Iris dataset
x.cal = iris[1:50, 1:4]
x.cal = scale(x.cal, center = T, scale = F)
loads = svd(x.cal)$v

# check loadings obtained using conventional algorithms

context('Conventional algorithms (max nPCs)')

m0 = pca(x.cal)
m1 = pca(x.cal, method = 'svd')
m2 = pca(x.cal, method = 'nipals')

test_that("Loadings have correct size (auto)", {
   expect_equal(nrow(m0$loadings), 4)
   expect_equal(nrow(m1$loadings), 4)
   expect_equal(nrow(m2$loadings), 4)
   expect_equal(ncol(m0$loadings), 4)
   expect_equal(ncol(m1$loadings), 4)
   expect_equal(ncol(m2$loadings), 4)
})

test_that("Loadings are computed correctly", {
   expect_equal(sum(m0$loadings - loads), 0)
   expect_equal(sum(m1$loadings - loads), 0)
   expect_equal(sum(m2$loadings - loads), 0)
})

context('Conventional algorithms (2 PCs)')

ncomp = 2
m0 = pca(x.cal, ncomp)
m1 = pca(x.cal, ncomp, method = 'svd')
m2 = pca(x.cal, ncomp, method = 'nipals')

test_that("Loadings have correct size", {
   expect_equal(nrow(m0$loadings), 4)
   expect_equal(nrow(m1$loadings), 4)
   expect_equal(nrow(m2$loadings), 4)
   expect_equal(ncol(m0$loadings), 2)
   expect_equal(ncol(m1$loadings), 2)
   expect_equal(ncol(m2$loadings), 2)
})

test_that("Loadings are computed correctly", {
   expect_equal(sum(m0$loadings - loads[, 1:2]), 0)
   expect_equal(sum(m1$loadings - loads[, 1:2]), 0)
   expect_equal(sum(m2$loadings - loads[, 1:2]), 0)
})

