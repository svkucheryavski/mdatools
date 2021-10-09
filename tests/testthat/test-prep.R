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


context("prep: snv, msc, norm, km")

test_that("SNV works correctly", {
   spectra <- simdata$spectra.c
   expect_silent(pspectra <- prep.snv(spectra))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))
   expect_equivalent(apply(pspectra, 1, mean), rep(0, nrow(spectra)))
   expect_equivalent(apply(pspectra, 1, sd), rep(1, nrow(spectra)))

   spectra <- simdata$spectra.c
   expect_silent(pspectra <- prep.norm(spectra, type = "snv"))
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

test_that("Normalization to unit sum works correctly", {
   spectra <- simdata$spectra.c
   expect_silent(pspectra <- prep.norm(spectra, "sum"))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))
   expect_equivalent(apply(pspectra, 1, function(x) sum(x)), rep(1, nrow(pspectra)))
})

test_that("Normalization to IS peak works (one value)", {
   spectra <- simdata$spectra.c

   expect_error(prep.norm(spectra, "is"))
   expect_error(prep.norm(spectra, "is", col.ind = -1))
   expect_error(prep.norm(spectra, "is", col.ind = 1000))
   expect_error(prep.norm(spectra, "is", col.ind = c(1, 1000)))
   expect_error(prep.norm(spectra, "is", col.ind = c(-1, 100)))


   col.ind <- 100
   expect_silent(pspectra <- prep.norm(spectra, type = "is", col.ind))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))
   expect_equivalent(apply(pspectra[, col.ind, drop = FALSE], 1, function(x) sum(x)), rep(1, nrow(pspectra)))

   col.ind <- 100:110
   expect_silent(pspectra <- prep.norm(spectra, type = "is", col.ind))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))
   expect_equivalent(apply(pspectra[, col.ind, drop = FALSE], 1, function(x) sum(x)), rep(1, nrow(pspectra)))

   col.ind <- 1:ncol(spectra)
   pspectra1 <- prep.norm(spectra, type = "is", col.ind)
   pspectra2 <- prep.norm(spectra, type = "sum")
   expect_equivalent(pspectra1, pspectra2)
})

test_that("PQN normalization works correctly", {
   data(simdata)
   spectra <- simdata$spectra.c
   ref.spectrum1 <- apply(simdata$spectra.c[c(1, 20, 30, 80, 80, 100), ], 2, mean)
   ref.spectrum2 <- apply(spectra, 2, mean)

   expect_error(prep.norm(spectra, "pqn", ref.spectrum = 1))
   expect_error(prep.norm(spectra, "pqn", ref.spectrum = ref.spectrum1[, 1:10, drop = FALSE]))

   # manual preprocessing of spectra
   ref.spectrum1 <- ref.spectrum1 / sum(abs(ref.spectrum1))
   ref.spectrum2 <- ref.spectrum2 / sum(abs(ref.spectrum2))
   pspectra2 <- pspectra1 <- matrix(0, nrow(spectra), ncol(spectra))

   for (i in seq_len(nrow(spectra))) {
      s <- spectra[i, ] / sum(abs(spectra[i, ]))

      q1 <- s / ref.spectrum1
      pspectra1[i, ] <- s / median(q1)

      q2 <- s / ref.spectrum2
      pspectra2[i, ] <- s / median(q2)
   }

   par(mfrow = c(2, 2))
   mdaplot(prep.norm(spectra, type = "sum"), type = "l")
   mdaplot(pspectra1, type = "l")
   mdaplot(pspectra2, type = "l")
   mdaplot(prep.norm(spectra, type = "pqn"), type = "l")

   expect_equivalent(prep.norm(spectra, type = "pqn", ref.spectrum = ref.spectrum1), pspectra1)
   expect_equivalent(prep.norm(spectra, type = "pqn", ref.spectrum = ref.spectrum2), pspectra2)
   expect_equivalent(prep.norm(spectra, type = "pqn"), pspectra2)
})

test_that("Kubelka-Munk works correctly", {
   spectra <- simdata$spectra.c
   spectra <- spectra - min(spectra) + 0.01 * max(spectra)
   expect_silent(pspectra <- prep.ref2km(spectra))
   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))
   expect_equal(dim(spectra), dim(pspectra))
   expect_true(all(pspectra > 0))

   spectra[10, 10] <- 0
   expect_error(prep.km(spectra))

   spectra[10, 10] <- -1
   expect_error(prep.km(spectra))
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

   # check that the derivatives are correct based on sin(x) and cos(x)
   x <- pi * seq(-2.4, 2.4, by = 0.1)
   y <- matrix(sin(x), nrow = 1)
   zeros <- which((abs(y) < 10^-12))
   ones <- which((abs(y) > 0.999999))

   # first derivative should have zeros where original function has maximums
   dy2 <- prep.savgol(y, 3, 1, 1)
   expect_equal(dim(dy2), dim(y))
   expect_equivalent(which((abs(dy2) < 10^-12)), ones)

   # second derivative should have zeros where the original function has zeros
   dy3 <- prep.savgol(y, 3, 2, 2)
   expect_equal(dim(dy3), dim(y))
   expect_equivalent(which((abs(dy3) < 10^-12)), zeros)

})

context("prep: alsbasecorr")

test_that("ALS baseline correction works correctly", {
   data(carbs)

   spectra <- mda.t(carbs$S)
   nvar <- ncol(spectra)

   spoiled_spectra <- spectra + rbind(
      dnorm(seq_len(nvar), nvar/2, nvar/5) * 100 * max(spectra[1, ]),
      dnorm(seq_len(nvar), nvar/3, nvar/5) * 500 * max(spectra[2, ]),
      dnorm(seq_len(nvar), nvar/1.5, nvar/5) * 100 * max(spectra[3, ])
   )

   corrected_spectra1 <- prep.alsbasecorr(spoiled_spectra)
   corrected_spectra2 <- prep.alsbasecorr(spoiled_spectra, plambda = 3, p = 0.05)
   corrected_spectra3 <- prep.alsbasecorr(spoiled_spectra, plambda = 4, p = 0.01)

   err_before <- sum((spectra - spoiled_spectra)^2) / sum(spectra^2)
   err_after1 <- sum((spectra - corrected_spectra1)^2) / sum(spectra^2)
   err_after2 <- sum((spectra - corrected_spectra2)^2) / sum(spectra^2)
   err_after3 <- sum((spectra - corrected_spectra3)^2) / sum(spectra^2)

   # preprocessing improves the spectra
   expect_true(err_after1 < err_before)
   expect_true(err_after2 < err_before)
   expect_true(err_after3 < err_before)

   # third case is the best
   expect_true(err_after3 < err_after1)
   expect_true(err_after3 < err_after2)
})

context("prep: trasnform")

test_that("Transformation works correctly", {
   x <- -5:5
   y <- cbind(exp(x), x^2)
   colnames(y) <- c("Y1", "Y2")
   attr(y, "yaxis.values") <- x
   attr(y, "yaxis.name") <- "Time to sleep"

   # errors

   # normal behaviour no extra parameters
   yp <- prep.transform(y, log)
   expect_equal(dim(yp), dim(y))
   expect_equal(colnames(yp), colnames(y))
   expect_equal(rownames(yp), rownames(y))
   expect_equal(mda.getattr(yp), mda.getattr(y))
   expect_equivalent(yp, log(y))

   # normal behaviour extra parameters
   yp <- prep.transform(y, log, 2)
   expect_equal(dim(yp), dim(y))
   expect_equal(colnames(yp), colnames(y))
   expect_equal(rownames(yp), rownames(y))
   expect_equal(mda.getattr(yp), mda.getattr(y))
   expect_equivalent(yp, log2(y))

   # normal behaviour extra parameters user defined function
   yp <- prep.transform(y, function(x, p) x^p, 1.25)
   expect_equal(dim(yp), dim(y))
   expect_equal(colnames(yp), colnames(y))
   expect_equal(rownames(yp), rownames(y))
   expect_equal(mda.getattr(yp), mda.getattr(y))
   expect_equivalent(yp, y^1.25)

});

context("prep: varsel")

test_that("Variable selection works correctly", {

   checkprep <- function(x, xp, ind) {
      expect_equal(ncol(xp), length(ind))
      expect_equal(nrow(xp), nrow(x))
      expect_equal(attr(xp, "yaxis.name"), attr(x, "yaxis.name"))
      expect_equal(attr(xp, "yaxis.values"), attr(x, "yaxis.values"))
      expect_equal(attr(xp, "xaxis.name"), attr(x, "xaxis.name"))
      expect_equal(attr(xp, "xaxis.values"), attr(x, "xaxis.values")[ind])
      expect_equal(attr(xp, "exclrows"), attr(x, "exclrows"))
   }

   x <- simdata$spectra.c
   attr(x, "xaxis.values") <- simdata$wavelength
   attr(x, "xaxis.name") <- "Wavelength, nm"
   attr(x, "yaxis.values") <- seq_len(nrow(x)) * 10
   attr(x, "yaxis.name") <- "Time, s"
   x <- mda.exclrows(x, c(1, 20, 30, 70))

   ind <- 21:140
   expect_silent(xp <- prep.varsel(x, ind))
   checkprep(x, xp, ind)

   ind <- c(21:30, 130:140)
   expect_silent(xp <- prep.varsel(x, ind))
   checkprep(x, xp, ind)

   ind <- rep(FALSE, ncol(x))
   ind[21:140] <- TRUE
   expect_silent(xp <- prep.varsel(x, ind))
   checkprep(x, xp, which(ind))

   ind <- rep(FALSE, ncol(x))
   ind[c(11:20, 130:140)] <- TRUE
   expect_silent(xp <- prep.varsel(x, ind))
   checkprep(x, xp, which(ind))

   # errors
   x <- mda.exclcols(x, c(1, 150))
   expect_error(xp <- prep.varsel(x, 20:140))

});


context("prep: combine methods together")

test_that("List of available methods is shown corectly", {
   expect_output(prep.list())
})

test_that("Errors are raised when necessary", {
   expect_error(prep("varsel", var.ind = 50:150))
   expect_error(prep("snv120"))
})

test_that("Method works with one preprocessing method in the list", {

   x <- simdata$spectra.c
   attr(x, "xaxis.values") <- simdata$wavelength
   attr(x, "xaxis.name") <- "Wavelength, nm"
   attr(x, "yaxis.values") <- seq_len(nrow(x)) * 10
   attr(x, "yaxis.name") <- "Time, s"
   x <- mda.exclrows(x, c(1, 20, 30, 70))

   p <- list(
      prep("savgol", list(width = 11, porder = 2, dorder = 1))
   )

   px1 <- employ.prep(p, x)
   px2 <- prep.savgol(x, width = 11, porder = 2, dorder = 1)

   expect_equal(px1, px2)
   expect_equal(mda.getattr(px1), mda.getattr(x))
   expect_equal(mda.getattr(px2), mda.getattr(x))
})

test_that("Method works with several preprocessing methods in the list", {

   x <- simdata$spectra.c
   attr(x, "xaxis.values") <- simdata$wavelength
   attr(x, "xaxis.name") <- "Wavelength, nm"
   attr(x, "yaxis.values") <- seq_len(nrow(x)) * 10
   attr(x, "yaxis.name") <- "Time, s"
   x <- mda.exclrows(x, c(1, 20, 30, 70))

   p <- list(
      prep("savgol", list(width = 11, porder = 2, dorder = 1)),
      prep("snv"),
      prep("autoscale", list(center = TRUE, scale = TRUE))
   )

   px1 <- employ.prep(p, x)
   px2 <- x
   px2 <- prep.savgol(px2, width = 11, porder = 2, dorder = 1)
   px2 <- prep.snv(px2)
   px2 <- prep.autoscale(px2, center = TRUE, scale = TRUE)

   expect_equal(px1, px2)
   expect_equal(mda.getattr(px1), mda.getattr(x))
   expect_equal(mda.getattr(px2), mda.getattr(x))
})

test_that("Method works with preprocessing methods and varable selection", {

   x <- simdata$spectra.c
   attr(x, "xaxis.values") <- simdata$wavelength
   attr(x, "xaxis.name") <- "Wavelength, nm"
   attr(x, "yaxis.values") <- seq_len(nrow(x)) * 10
   attr(x, "yaxis.name") <- "Time, s"
   x <- mda.exclrows(x, c(1, 20, 30, 70))

   p <- list(
      prep("savgol", list(width = 11, porder = 2, dorder = 1)),
      prep("snv"),
      prep("autoscale", list(center = TRUE, scale = TRUE)),
      prep("varsel", list(var.ind = 50:130))
   )

   px1 <- employ.prep(p, x)
   px2 <- x
   px2 <- prep.savgol(px2, width = 11, porder = 2, dorder = 1)
   px2 <- prep.snv(px2)
   px2 <- prep.autoscale(px2, center = TRUE, scale = TRUE)
   px2 <- mda.subset(px2, select = 50:130)

   expect_equal(px1, px2)
   expect_equal(mda.getattr(px1), mda.getattr(px2))
})


test_that("Method works with user defined preprocessing method", {

   x <- simdata$spectra.c
   attr(x, "xaxis.values") <- simdata$wavelength
   attr(x, "xaxis.name") <- "Wavelength, nm"
   attr(x, "yaxis.values") <- seq_len(nrow(x)) * 10
   attr(x, "yaxis.name") <- "Time, s"
   x <- mda.exclrows(x, c(1, 20, 30, 70))

   p <- list(
      prep("r2a", method = function(data) log(1/abs(data))),
      prep("savgol", list(width = 11, porder = 2, dorder = 1)),
      prep("snv")
   )

   px1 <- employ.prep(p, x)
   px2 <- log(1/abs(x))
   px2 <- prep.savgol(px2, width = 11, porder = 2, dorder = 1)
   px2 <- prep.snv(px2)

   expect_equal(px1, px2)
   expect_equal(mda.getattr(px1), mda.getattr(px2))
})
