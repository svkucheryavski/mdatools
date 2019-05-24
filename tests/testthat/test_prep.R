# new tests on top

context('Preprocessing: autoscale')

data(simdata)

# normal spectra
X1 = simdata$spectra.c
p11 = prep.autoscale(X1)
p12 = prep.autoscale(X1, center = T, scale = F)
p13 = prep.autoscale(X1, center = T, scale = T)

# spectra with constant and small variation variables
X2 = simdata$spectra.c
X2[, 20:50] = matrix(apply(X2[, 20:50], 2, mean), nrow = nrow(X2), ncol = 31, byrow = T) # constant
X2[, 51:60] = X2[, 51:60] * 0.001 + 0.1 # small CV (around 0.17%)
p21 = prep.autoscale(X2)
p22 = prep.autoscale(X2, center = T, scale = F)
p23 = prep.autoscale(X2, center = T, scale = T)
p24 = prep.autoscale(X2, center = T, scale = T, max.cov = 0.2)
p25 = prep.autoscale(-X2, center = T, scale = T, max.cov = 0.2)

test_that("Default autoscaling is correct", {
   expect_equal(attr(p11, 'prep:center'), apply(X1, 2, mean))
   expect_equal(attr(p11, 'prep:scale'), FALSE)
   expect_equal(p11, p12)
   expect_equal(p11, scale(X1, center = T, scale = F), check.attributes = FALSE)
})

test_that("Full autoscaling is correct", {
   expect_equal(attr(p13, 'prep:center'), apply(X1, 2, mean))
   expect_equal(attr(p13, 'prep:scale'), apply(X1, 2, sd))
   expect_equal(p13, scale(X1, center = T, scale = T), check.attributes = FALSE)
})

test_that("Default autoscaling for data with constant variables is correct", {
   expect_equal(attr(p21, 'prep:center'), apply(X2, 2, mean))
   expect_equal(attr(p21, 'prep:scale'), FALSE)
   expect_equal(p21, p22)
   expect_equal(p21, scale(X2, center = T, scale = F), check.attributes = FALSE)
})

test_that("Full autoscaling for data with constant variables is correct", {
   expect_equal(attr(p23, 'prep:center'), apply(X2, 2, mean))
   expect_equal(attr(p23, 'prep:scale')[-(20:50)], apply(X2, 2, sd)[-(20:50)])
   expect_equal(attr(p23, 'prep:scale')[(20:50)], rep(1, 31), check.attributes = FALSE)
})

test_that("Full autoscaling with limits for max.cov is correct", {
   expect_equal(attr(p24, 'prep:scale')[-(20:60)], apply(X2, 2, sd)[-(20:60)])
   expect_equal(attr(p24, 'prep:scale')[(20:60)], rep(1, 41), check.attributes = FALSE)
})

test_that("Full autoscaling with limits for max.cov is correct for negative data", {
   expect_equal(attr(p25, 'prep:scale')[-(20:60)], apply(-X2, 2, sd)[-(20:60)])
   expect_equal(attr(p25, 'prep:scale')[(20:60)], rep(1, 41), check.attributes = FALSE)
})
