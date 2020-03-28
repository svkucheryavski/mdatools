# new tests on top

context("prep: autoscale")

data(simdata)

# normal spectra
X1 <- simdata$spectra.c
p11 <- prep.autoscale(X1)
p12 <- prep.autoscale(X1, center = T, scale = F)
p13 <- prep.autoscale(X1, center = T, scale = T)

# spectra with constant and small variation variables
X2 <- simdata$spectra.c
X2[, 20:50] <- matrix(apply(X2[, 20:50], 2, mean), nrow = nrow(X2), ncol = 31, byrow = T) # constant
X2[, 51:60] <- X2[, 51:60] * 0.001 + 0.1 # small CV (around 0.17%)
p21 <- prep.autoscale(X2)
p22 <- prep.autoscale(X2, center = T, scale = F)
p23 <- prep.autoscale(X2, center = T, scale = T)
p24 <- prep.autoscale(X2, center = T, scale = T, max.cov = 0.2)
p25 <- prep.autoscale(-X2, center = T, scale = T, max.cov = 0.2)

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


context("prep: snv, msc, norm")

test_that("SNV works correctly", {
   spectra <- simdata$spectra.c
   expect_silent(pspectra <- prep.snv(spectra))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))
   expect_equivalent(apply(pspectra, 1, mean), rep(0, nrow(spectra)))
   expect_equivalent(apply(pspectra, 1, sd), rep(1, nrow(spectra)))
})

test_that("MSC works correctly", {
   spectra <- simdata$spectra.c
   mn <- apply(spectra, 2, mean)
   md <- apply(spectra, 2, median)

   spectra <- simdata$spectra.c
   expect_silent(pspectra1 <- prep.msc(spectra))
   expect_silent(pspectra2 <- prep.msc(spectra, mspectrum = mn))
   expect_silent(pspectra3 <- prep.msc(spectra, mspectrum = md))
   expect_silent(pspectra4 <- prep.msc(spectra, mspectrum = matrix(md)))

   expect_equivalent(attr(pspectra1, "mspectrum"), mn)
   expect_equivalent(attr(pspectra2, "mspectrum"), mn)
   expect_equivalent(attr(pspectra3, "mspectrum"), md)
   expect_equivalent(attr(pspectra4, "mspectrum"), md)

   expect_equal(pspectra1, pspectra2)
   expect_equal(pspectra3, pspectra4)

   expect_equal(mda.getattr(spectra), mda.getattr(pspectra1))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra2))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra3))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra4))

   # wrong scenario
   expect_error(prep.msc(spectra, mspectrum = c(1, 2, 3)))
   expect_error(prep.msc(spectra, mspectrum = matrix(1:10, ncol = 5)))
})

test_that("Normalization to unit area works correctly", {
   spectra <- simdata$spectra.c
   expect_silent(pspectra <- prep.norm(spectra, "area"))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))
   expect_equivalent(apply(pspectra, 1, function(x) sum(abs(x))), rep(1, nrow(pspectra)))
})

test_that("Normalization to unit length works correctly", {
   spectra <- simdata$spectra.c
   expect_silent(pspectra <- prep.norm(spectra, "length"))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))
   expect_equivalent(apply(pspectra, 1, function(x) sqrt(sum(x^2))), rep(1, nrow(pspectra)))
})

context("prep: savgol")

test_that("SavGol smoothing works correctly", {
   spectra <- simdata$spectra.c
   expect_silent(pspectra <- prep.savgol(spectra, 3, 1, 0))
   expect_silent(pspectra <- prep.savgol(spectra, 3, 1, 1))
   expect_silent(pspectra <- prep.savgol(spectra, 3, 2, 2))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))

   # wrong values for filter width
   expect_error(pspectra <- prep.savgol(spectra,  1, 1, 0))
   expect_error(pspectra <- prep.savgol(spectra, -1, 1, 0))
   expect_error(pspectra <- prep.savgol(spectra,  2, 1, 0))

   # wrong values for pdorder
   expect_error(pspectra <- prep.savgol(spectra, 3, 0, 0))
   expect_error(pspectra <- prep.savgol(spectra, 3, -1, 0))
   expect_error(pspectra <- prep.savgol(spectra, 3, 10, 0))

   # wrong values for ddorder
   expect_error(pspectra <- prep.savgol(spectra, 3, 1, -1))
   expect_error(pspectra <- prep.savgol(spectra, 3, 1,  2))
   expect_error(pspectra <- prep.savgol(spectra, 3, 2,  4))

})
