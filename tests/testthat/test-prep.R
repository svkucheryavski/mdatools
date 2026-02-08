# new tests on top

context("prep: checks")

data(simdata)

# spectra as data frame or vector — should show error
Xe <- as.data.frame(simdata$spectra.c)
xe <- 1:10

test_that("Preprocessing methods raise error if data is not a matrix", {

   expect_error(prep.snv((Xe)))
   expect_error(prep.snv((xe)))

   expect_error(prep.msc((Xe)))
   expect_error(prep.msc((xe)))

   expect_error(prep.norm((Xe)))
   expect_error(prep.norm((xe)))

   expect_error(prep.km((Xe)))
   expect_error(prep.km((xe)))

   expect_error(prep.savgol((Xe)))
   expect_error(prep.savgol((xe)))

   expect_error(prep.alsbasecorr((Xe)))
   expect_error(prep.alsbasecorr((xe)))

   expect_error(prep.varsel((Xe)))
   expect_error(prep.varsel((xe)))
})

context("prep: autoscale")

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

   # par(mfrow = c(2, 2))
   # mdaplot(prep.norm(spectra, type = "sum"), type = "l")
   # mdaplot(pspectra1, type = "l")
   # mdaplot(pspectra2, type = "l")
   # mdaplot(prep.norm(spectra, type = "pqn"), type = "l")

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

   # check for small numeric sets
   x = matrix(c(1, 1, 1, 3, 4, 7, 4, 3, 1, 1, 1.0), nrow = 1)
   y = matrix(c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5.0), nrow = 1)
   z = matrix(c(4, 6, 4, 6, 4, 6, 4, 6, 4, 6.0), nrow = 1)

   # no derivartive
   expect_equivalent(prep.savgol(x, width = 5, porder = 1, dorder = 0), c(0.4, 1.2, 2.0, 3.2, 3.8, 4.2, 3.8, 3.2, 2.0, 1.2, 0.4))
   expect_equivalent(prep.savgol(y, width = 3, porder = 1, dorder = 0), c(5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0))
   expect_equivalent(prep.savgol(z, width = 5, porder = 1, dorder = 0), c(4.8, 4.8, 4.8, 5.2, 4.8, 5.2, 4.8, 5.2, 5.2, 5.2))

   # first derivartive
   expect_equivalent(prep.savgol(x, width = 3, porder = 1, dorder = 1), c(0.0, 0.0, 1.0, 1.5, 2.0, 0.0, -2.0, -1.5, -1.0, 0.0, 0.0))
   expect_equivalent(prep.savgol(y, width = 3, porder = 1, dorder = 1), c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
   expect_equivalent(prep.savgol(z, width = 3, porder = 1, dorder = 1), c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

   # second derivartive
   expect_equivalent(prep.savgol(x, width = 3, porder = 2, dorder = 2), c(0.0, 0.0, 2.0, -1.0, 2.0, -6.0, 2.0, -1.0, 2.0, 0.0, 0.0))
   expect_equivalent(prep.savgol(y, width = 3, porder = 2, dorder = 2), c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
   expect_equivalent(prep.savgol(z, width = 3, porder = 2, dorder = 2), c(-4.0, -4.0, 4.0, -4.0, 4.0, -4.0, 4.0, -4.0, 4.0, 4.0))


   y1 <- prep.savgol(
      rbind(
         matrix(c(1, 1, 1, 3, 4, 7, 4, 3, 1, 1, 1.0), nrow = 1),
         matrix(c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5.0), nrow = 1)
      ), width = 5, porder = 1, dorder = 0
   )

   y2 <- rbind(
      prep.savgol(matrix(c(1, 1, 1, 3, 4, 7, 4, 3, 1, 1, 1.0), nrow = 1), width = 5, porder = 1, dorder = 0),
      prep.savgol(matrix(c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5.0), nrow = 1), width = 5, porder = 1, dorder = 0)
   )

   expect_equivalent(y1, y2)
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
      prep("savgol", width = 11, porder = 2, dorder = 1)
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
      prep("savgol", width = 11, porder = 2, dorder = 1),
      prep("norm", type="snv"),
      prep("center", type = "median"),
      prep("scale", type = "iqr")
   )

   p <- prep.fit(p, x)

   px1 <- prep.apply(p, x)

   px2 <- x
   px2 <- prep.savgol(px2, width = 11, porder = 2, dorder = 1)
   px2 <- prep.norm(px2, type = "snv")
   px2 <- prep.center(px2, type="median")
   px2 <- prep.scale(px2, type="iqr")

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
      prep("savgol", width = 11, porder = 2, dorder = 1),
      prep("norm", type = "snv"),
      prep("center"),
      prep("scale"),
      prep("varsel", var.ind = 50:130)
   )

   p <- prep.fit(p, x)

   px1 <- prep.apply(p, x)
   px2 <- x
   px2 <- prep.savgol(px2, width = 11, porder = 2, dorder = 1)
   px2 <- prep.norm(px2, type = "snv")
   px2 <- prep.center(px2)
   px2 <- prep.scale(px2)
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
      prep("savgol", width = 11, porder = 2, dorder = 1),
      prep("norm", type = "snv")
   )

   px1 <- prep.apply(p, x)
   px2 <- x
   px2 <- prep.savgol(px2, width = 11, porder = 2, dorder = 1)
   px2 <- prep.norm(px2, "snv")

   expect_equal(px1, px2)
   expect_equal(mda.getattr(px1), mda.getattr(px2))
})


context("prep: spikes")

test_that("Removing spikes works correctly", {
   data(carbs)
   spectra <- carbs$D
   sspectra <- spectra
   sspectra[1,  110] <- sspectra[1,  110] * 10
   sspectra[10, 220] <- sspectra[10, 220] * 10
   sspectra[15, 330] <- sspectra[15, 330] * 10
   sspectra[20, 450] <- sspectra[20, 450] * 10

   expect_silent(pspectra <- prep.spikes(sspectra, width = 5, threshold = 6))
   ind2 <- which(abs(spectra - pspectra) / abs(spectra) > 0.5)
   ind1 <- which(abs(spectra - sspectra) / abs(spectra) > 0.5)

   expect_equal(mda.getattr(spectra), mda.getattr(pspectra))
   expect_equal(length(ind2), 0)
   expect_equal(length(ind1), 4)

})

context("prep: center")

do_test_center <- function(type, f) {
   data(people)

   # full scale
   xp <- prep.center(people, type = type)
   v <- apply(people, 2, f)
   expect_equal(attr(xp, "prep:center"), v)
   expect_equivalent(xp, scale(people, center = v, scale = FALSE))

   # parameters only
   p <- prep.center.params(people, type = type)
   expect_true(is.list(p))
   expect_equal(length(p), 2)
   expect_equal(names(p), c("type", "center"))
   expect_equal(p[[1]], type)
   expect_equal(p[[2]], v)

   # with outliers
   o <- c(5, 20)
   v <- apply(people[-o, ], 2, f)
   people <- mda.exclrows(people, o)
   xp <- prep.center(people, type = type)
   expect_equal(attr(xp, "prep:center"), v)
   expect_equivalent(xp, scale(people, center = v, scale = FALSE))

   p <- prep.center.params(people, type = type)
   expect_equal(p[[2]], v)
}

test_that("Centering works correctly", {
   do_test_center("mean", mean)
   do_test_center("median", median)
})

test_that("Manual centering works correctly", {
   data(people)
   center <- rep(10, ncol(people))
   xp <- prep.center(people, center = center)
   expect_equal(attr(xp, "prep:center"), center)
   expect_equivalent(xp, scale(people, center = center, scale = FALSE))
})



context("prep: scale")

do_test_scale <- function(type, f) {
   data(people)

   # full scale
   xp <- prep.scale(people, type = type)
   v <- apply(people, 2, f)
   expect_equal(attr(xp, "prep:scale"), v)
   expect_equivalent(xp, scale(people, center = FALSE, scale = v))

   # parameters only
   p <- prep.scale.params(people, type = type)
   expect_true(is.list(p))
   expect_equal(length(p), 2)
   expect_equal(names(p), c("type", "scale"))
   expect_equal(p[[1]], type)
   expect_equal(p[[2]], v)

   # with outliers
   o <- c(5, 20)
   v <- apply(people[-o, ], 2, f)
   people <- mda.exclrows(people, o)
   xp <- prep.scale(people, type = type)
   expect_equal(attr(xp, "prep:scale"), v)
   expect_equivalent(xp, scale(people, center = FALSE, scale = v))

   p <- prep.scale.params(people, type = type)
   expect_equal(p[[2]], v)
}

test_that("scaling works correctly", {
   pareto <- function(x) sqrt(sd(x))
   rng <- function(x) max(x)- min(x)

   do_test_scale("sd", sd)
   do_test_scale("iqr", IQR)
   do_test_scale("range", rng)
   do_test_scale("pareto", pareto)
})


context("prep: emsc")

test_that("emsc works correctly (p = 0)", {
   x <- c(-3, -2, -1, 0, 2, 3, 4)^2
   data <- rbind(x, x + 10, x - 10, x + 5, x - 5, x * 2 - 3, x * (-2) - 3)
   mspectrum <- apply(data, 2, mean)

   pdata1 <- prep.emsc(data)
   pdata2 <- prep.emsc(data, mspectrum = mspectrum)
   pdata3 <- prep.emsc(data, degree = 0)
   pdata4 <- prep.emsc(data, degree = 0, mspectrum = mspectrum)

   # conventional MSC
   pdata5 <- matrix(0, nrow(data), ncol(data))
   for (i in seq_len(nrow(data))) {
      coef <- coef(lm(data[i, ] ~ mspectrum))
      pdata5[i, ] <- (data[i, ] - coef[1]) / coef[2]
   }

   expect_equivalent(pdata1, pdata5)
   expect_equivalent(pdata2, pdata5)
   expect_equivalent(pdata3, pdata5)
   expect_equivalent(pdata4, pdata5)

})

test_that("emsc works correctly (p > 0)", {
   x <- c(-3, -2, -1, 0, 2, 3, 4)^2
   data <- rbind(x^1.1 + 5, x^1.4 - 5, 2*x^1.1 + 5, 2*x^1.4 - 5,x ^ 1.2 - 3, x^1.3 + 3, x ^ 1.5 - 13, x^1.23 + 13)
   mspectrum <- apply(data, 2, mean)

   pdata1a <- prep.emsc(data, 1)
   pdata1b <- prep.emsc(data, 4)

   pdata2a <- prep.emsc(data, 1, mspectrum = x)
   pdata2b <- prep.emsc(data, 4, mspectrum = x)

   ref1 <- matrix(mspectrum, nrow = nrow(data), ncol = ncol(data), byrow = TRUE)
   ref2 <- matrix(x, nrow = nrow(data), ncol = ncol(data), byrow = TRUE)

   err1a <- sum((pdata1a - ref1)^2) / sum(ref1^2)
   err1b <- sum((pdata1b - ref1)^2) / sum(ref1^2)
   err2a <- sum((pdata2a - ref2)^2) / sum(ref2^2)
   err2b <- sum((pdata2b - ref2)^2) / sum(ref2^2)

   expect_true(err1a < 0.01)
   expect_true(err1b < 0.001)
   expect_true(err2a < 0.01)
   expect_true(err2b < 0.001)

})



context("prep: prep.fit and writeJSON")

test_that("fitting prep model works correctly (case Tecator)", {

   dc <- read.csv2('./tecator_train.csv', row.names = 1, check.names = FALSE)
   Xc <- as.matrix(dc[, -1])

   dt <- read.csv2('./tecator_test.csv', row.names = 1, check.names = FALSE)
   Xt <- as.matrix(dt[, -1])

   dcp.ref <- read.csv2('./tecator_train-preprocessed.csv', row.names = 1, check.names = FALSE)
   Xcp.ref <- as.matrix(dcp.ref[, -1])

   dtp.ref <- read.csv2('./tecator_test-preprocessed.csv', row.names = 1, check.names = FALSE)
   Xtp.ref <- as.matrix(dtp.ref[, -1])

   p <- list(
      prep("spikes", width = 5, threshold = 5),
      prep("savgol", width = 7, porder = 2, dorder = 2),
      prep("varsel", var.ind = 11:90),
      prep("norm", type = "snv"),
      prep("emsc", degree = 2),
      prep("center", type = "mean"),
      prep("scale", type = "pareto")
   )


   pm  <- prep.fit(p, Xc)
   Xcp <- prep.apply(pm, Xc)
   Xtp <- prep.apply(pm, Xt)

   expect_equivalent(Xcp, Xcp.ref)
   expect_equivalent(Xtp, Xtp.ref)
})

test_that("fitting prep model works correctly (case HS-)", {

   dc <- read.csv('./hs_train.csv', row.names = 1, check.names = FALSE)
   Xc <- as.matrix(dc[, -1])

   dcp.ref <- read.csv('./hs_train-preprocessed.csv', row.names = 1, check.names = FALSE)
   Xcp.ref <- as.matrix(dcp.ref[, -1])

   p <- list(
      prep("varsel", var.ind = 121:1755),
      prep("spikes", width = 5, threshold = 10),
      prep("savgol", width = 5, porder = 1, dorder = 0),
      prep("alsbasecorr", plambda = 2.5, p = 0.008),
      prep("varsel", var.ind = 81:1580),
      prep("norm", type = "is", col.ind = 976),
      prep("center", type = "mean")
   )

   pm  <- prep.fit(p, Xc)
   Xcp <- prep.apply(pm, Xc)
   expect_equivalent(Xcp, Xcp.ref)
})



testCase <- function(p, Xc, jsonFilename) {


   # fit pm model
   pm  <- prep.fit(p, Xc)
   #print(pm)
   # test JSON strings
   json1 <- prep.asjson(pm)
   fileConn <- file(jsonFilename)
   json2 <- readLines(fileConn, warn = FALSE)
   close(fileConn)

   expect_equal(extractStringArray(json1, "info"), extractStringArray(json2, "info"))
   expect_equal(extractValue(json1, "npred"), extractValue(json2, "npred"))


   # save it to JSON file and then load it back
   prep.writeJSON(pm, './preprocessing-tmp.json')
   pm1 <- prep.readJSON('./preprocessing-tmp.json')

   # load model developed in web-app
   pm2 <- prep.readJSON(jsonFilename)

   # all models should be the same
   for (n in seq_along(pm)) {
      if (!is.list(pm[[n]])) next
      testList(pm[[n]]$params, pm1[[n]]$params)
      testList(pm[[n]]$params, pm2[[n]]$params)
   }
}

test_that("JSON() and readJSON() methods work correctly.", {

   dc <- read.csv2('./tecator_train.csv', row.names = 1, check.names = FALSE)
   Xc <- as.matrix(dc[, -1])

   # one method
   p <- list(
      prep("emsc", degree = 2)
   )
   testCase(p, Xc, "./preprocessing-model-1.json")

   # two methods
   p <- list(
      prep("norm", type = "snv"),
      prep("emsc", degree = 2)
   )
   testCase(p, Xc, "./preprocessing-model-2.json")

   # three methods
   p <- list(
      prep("spikes", width = 5, threshold = 5),
      prep("norm", type = "snv"),
      prep("emsc", degree = 2)
   )
   testCase(p, Xc, "./preprocessing-model-3.json")

   # four methods
   p <- list(
      prep("spikes", width = 5, threshold = 5),
      prep("savgol", width = 7, porder = 2, dorder = 2),
      prep("norm", type = "snv"),
      prep("emsc", degree = 2)
   )
   testCase(p, Xc, "./preprocessing-model-4.json")

   # five methods
   p <- list(
      prep("spikes", width = 5, threshold = 5),
      prep("savgol", width = 7, porder = 2, dorder = 2),
      prep("varsel", var.ind = 11:90),
      prep("norm", type = "snv"),
      prep("emsc", degree = 2),
      prep("center", type = "mean"),
      prep("scale", type = "pareto")
   )
   testCase(p, Xc, "./preprocessing-model-full-tecator.json")


   # five methods only scale
   p <- list(
      prep("spikes", width = 5, threshold = 5),
      prep("savgol", width = 7, porder = 2, dorder = 2),
      prep("varsel", var.ind = 11:90),
      prep("norm", type = "snv"),
      prep("emsc", degree = 2),
      prep("scale", type = "range")
   )
   testCase(p, Xc, "./preprocessing-model-full-tecator-scale.json")

   # five methods only center
   p <- list(
      prep("spikes", width = 5, threshold = 5),
      prep("savgol", width = 7, porder = 2, dorder = 2),
      prep("varsel", var.ind = 11:90),
      prep("norm", type = "snv"),
      prep("emsc", degree = 2),
      prep("center", type = "median")
   )
   testCase(p, Xc, "./preprocessing-model-full-tecator-center.json")

   # full hs- case

   dc <- read.csv('./hs_train.csv', row.names = 1, check.names = FALSE)
   Xc <- as.matrix(dc[, -1])

   p <- list(
      prep("varsel", var.ind = 121:1755),
      prep("spikes", width = 5, threshold = 10),
      prep("savgol", width = 5, porder = 1, dorder = 0),
      prep("alsbasecorr", plambda = 2.5, p = 0.008),
      prep("varsel", var.ind = 81:1580),
      prep("norm", type = "is", col.ind = 976),
      prep("center", type = "mean")
   )

   testCase(p, Xc, "./preprocessing-model-full-hs.json")

})

