######################################
# Tests for pls class methods        #
######################################

setup({
   pdf(file = tempfile("mdatools-test-pls-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-pls-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

#################################################
# Block 1. Tests pls.run() using random values  #
#################################################

context("pls: test basic static methods")
x <- matrix(rnorm(1000), ncol = 100)
y1 <- matrix(rnorm(10), ncol = 1)
y2 <- matrix(rnorm(20), ncol = 2)

test_that("pls.run() works correctly", {

   # wrong values for method and ncomp raise error
   expect_error(pls.run(x, y1, 5, method = "wrong"))
   expect_error(pls.run(x, y1, 0))
   expect_error(pls.run(x, y1, 1000))

   # correct values give no errors or warnings
   ncomp <- 5
   expect_silent(pls.run(x, y1, ncomp, method = "simpls"))
   expect_silent(pls.run(x, y1, ncomp))
   expect_silent(pls.run(x, y2, ncomp, method = "simpls"))
   expect_silent(pls.run(x, y2, ncomp))

   # outcome has correct dimensions for different number of components
   for (ncomp in c(1, 5, 9)) {
      m1 <- pls.run(x, y1, ncomp)
      expect_equal(dim(m1$coeffs), c(ncol(x), ncomp, ncol(y1)))
      expect_equal(dim(m1$weights), c(ncol(x), ncomp))
      expect_equal(dim(m1$xloadings), c(ncol(x), ncomp))
      expect_equal(dim(m1$yloadings), c(ncol(y1), ncomp))
      expect_equal(m1$ncomp, ncomp)

      m2 <- pls.run(x, y2, ncomp)
      expect_equal(dim(m2$coeffs), c(ncol(x), ncomp, ncol(y2)))
      expect_equal(dim(m2$weights), c(ncol(x), ncomp))
      expect_equal(dim(m2$xloadings), c(ncol(x), ncomp))
      expect_equal(dim(m2$yloadings), c(ncol(y2), ncomp))
      expect_equal(m2$ncomp, ncomp)

   }
})

############################################################
# Block 2. Tests pls.cal() using people and spectral data  #
############################################################

context("pls: test pls.cal() method for different settings")

# mock some data
datasets <- list()

## small data no names no extra arguments
data(people)
xc <- people[, -4, drop = F]
yc <- people[, 4, drop = F]
datasets[["people"]] <- list(xc = xc, yc = yc, center = T, scale = T, ncomp = 8)

## spectral data with extra attributes, three responses, hidden values and test set
data(simdata)
xc <- simdata$spectra.c
yc <- simdata$conc.c
xc <- mda.exclrows(xc, c(1, 10, 20, 30, 40))
xc <- mda.exclcols(xc, c(1, 100:110))
yc <- mda.exclrows(yc, c(1, 10, 20, 30, 40))
attr(xc, "name") <- "Spectra, cal"
attr(xc, "xaxis.name") <- "Wavelength, nm"
attr(xc, "xaxis.values") <- simdata$wavelength

xt <- simdata$spectra.t
yt <- simdata$conc.t
xt <- mda.exclrows(xt, c(15, 35))
yt <- mda.exclrows(yt, c(15, 35))
attr(xt, "name") <- "Spectra, val"
attr(xt, "xaxis.name") <- "Wavelength, nm"
attr(xt, "xaxis.values") <- simdata$wavelength

datasets[["spectra"]] <- list(xc = xc, yc = yc, xt = xt, yt = yt, center = T, scale = F, ncomp = 6)

for (i in seq_along(datasets)) {

   d <- datasets[[i]]
   name <- names(datasets)[i]
   xaxis.name.x <- if (is.null(attr(d$xc, "xaxis.name"))) "Predictors" else attr(d$xc, "xaxis.name")
   xaxis.name.y <- if (is.null(attr(d$yc, "xaxis.name"))) "Responses" else attr(d$yc, "xaxis.name")

   for (ncomp in 1:d$ncomp) {
      obj <- pls.cal(d$xc, d$yc, ncomp = ncomp, center = d$center, scale = d$scale)

      test_that("dimension and attributes for x-loadings are correct", {
         expect_equal(dim(obj$xloadings), c(ncol(d$xc), ncomp))
         expect_equal(attr(obj$xloadings, "exclrows"), attr(d$xc, "exclcols"))
         expect_equal(attr(obj$xloadings, "xaxis.name"), "Components")
         expect_equal(attr(obj$xloadings, "yaxis.values"), attr(d$xc, "xaxis.values"))
         expect_equal(attr(obj$xloadings, "yaxis.name"), xaxis.name.x)
      })

      test_that("dimension and attributes for y-loadings are correct", {
         expect_equal(dim(obj$yloadings), c(ncol(d$yc), ncomp))
         expect_equal(attr(obj$yloadings, "exclrows"), attr(d$yc, "exclcols"))
         expect_equal(attr(obj$yloadings, "xaxis.name"), "Components")
         expect_equal(attr(obj$yloadings, "yaxis.values"), attr(d$yc, "xaxis.values"))
         expect_equal(attr(obj$yloadings, "yaxis.name"), xaxis.name.y)
      })

      test_that("dimension and attributes for weights are correct", {
         expect_equal(dim(obj$weights), c(ncol(d$xc), ncomp))
         expect_equal(attr(obj$weights, "exclrows"), attr(d$xc, "exclcols"))
         expect_equal(attr(obj$weights, "xaxis.name"), "Components")
         expect_equal(attr(obj$weights, "yaxis.values"), attr(d$xc, "xaxis.values"))
         expect_equal(attr(obj$weights, "yaxis.name"), xaxis.name.x)
      })

      test_that("dimension and attributes for coeffs are correct", {
         expect_equal(class(obj$coeffs), "regcoeffs")
         expect_equal(dim(obj$coeffs$values), c(ncol(d$xc), ncomp, ncol(d$yc)))
         expect_equal(attr(obj$coeffs$values, "exclrows"), attr(d$xc, "exclcols"))
         expect_equal(attr(obj$coeffs$values, "yaxis.values"), attr(d$xc, "xaxis.values"))
         expect_equal(attr(obj$coeffs$values, "yaxis.name"), xaxis.name.x)
      })

      test_that("other model parameters are correct", {
         expect_equal(obj$ncomp, ncomp)
         expect_equal(obj$exclrows, attr(d$xc, "exclrows"))
         expect_equal(obj$exclcols, attr(d$xc, "exclcols"))
         expect_equal(class(obj), c("pls", "regmodel"))
      })

      if (length(attr(d$xc, "exclcols")) > 0) {
         exclcols <- attr(d$xc, "exclcols")
         expect_equivalent(obj$xloadings[exclcols, , drop = F], matrix(0, length(exclcols), ncomp))
         expect_equivalent(obj$weights[exclcols, , drop = F], matrix(0, length(exclcols), ncomp))
         expect_equivalent(
            obj$coeffs$values[exclcols, , , drop = F],
            array(0, dim = c(length(exclcols), ncomp, ncol(d$yc)))
         )
      }
   }

   test_that("pls: test pls.cal() method estimates ncomp correctly", {
      exclrows <- attr(d$xc, "exclrows")
      exclcols <- attr(d$xc, "exclcols")
      nrows <- nrow(d$xc) - length(exclrows)

      obj <- pls.cal(d$xc, d$yc, ncomp = 50, center = d$center, scale = d$scale)
      exp_ncomp <- min(nrows - 1, ncol(d$xc) - length(exclcols))
      expect_lte(obj$ncomp, exp_ncomp)

      obj <- pls.cal(d$xc, d$yc, ncomp = 50, center = d$center, scale = d$scale, cv = 1)
      exp_ncomp <- min(nrows - 2, ncol(d$xc) - length(exclcols))
      expect_lte(obj$ncomp, exp_ncomp)

      obj <- pls.cal(d$xc, d$yc, ncomp = 50, center = d$center, scale = d$scale, cv = 2)
      exp_ncomp <- min(nrows - 1 - nrows/2, ncol(d$xc) - length(exclcols))
      expect_lte(obj$ncomp, exp_ncomp)

   })
}

################################################################
# Block 3. Test pls.predict() using people and spectral data   #
################################################################

context("pls: test pls.predict() method for different settings")

for (i in seq_along(datasets)) {

   d <- datasets[[i]]
   name <- names(datasets)[i]
   xaxis.name.x <- if (is.null(attr(d$xc, "xaxis.name"))) "Predictors" else attr(d$xc, "xaxis.name")
   xaxis.name.y <- if (is.null(attr(d$yc, "xaxis.name"))) "Responses" else attr(d$yc, "xaxis.name")
   yaxis.name.x <- if (is.null(attr(d$xc, "yaxis.name"))) "Objects" else attr(d$xc, "yaxis.name")
   yaxis.name.y <- if (is.null(attr(d$yc, "yaxis.name"))) "Objects" else attr(d$yc, "yaxis.name")

   for (ncomp in 1:d$ncomp) {
      obj <- pls.cal(d$xc, d$yc, ncomp = ncomp, center = d$center, scale = d$scale)
      res <- predict(obj, d$xc, d$yc)

      test_that("dimension and attributes for xdecomp are correct", {
         # scores
         expect_equal( dim(res$xdecomp$scores), c(nrow(d$xc), ncomp))
         expect_equal(attr(res$xdecomp$scores, "exclrows"), attr(d$xc, "exclrows"))
         expect_equal(attr(res$xdecomp$scores, "yaxis.values"), attr(d$xc, "yaxis.values"))
         expect_equal(attr(res$xdecomp$scores, "yaxis.name"), yaxis.name.x)

         # othogonal distance
         expect_equal( dim(res$xdecomp$Q), c(nrow(d$xc), ncomp))
         expect_equal(attr(res$xdecomp$Q, "exclrows"), attr(d$xc, "exclrows"))
         expect_equal(attr(res$xdecomp$Q, "yaxis.values"), attr(d$xc, "yaxis.values"))
         expect_equal(attr(res$xdecomp$Q, "yaxis.name"), yaxis.name.x)

         # score distance
         expect_equal( dim(res$xdecomp$T2), c(nrow(d$xc), ncomp))
         expect_equal(attr(res$xdecomp$T2, "exclrows"), attr(d$xc, "exclrows"))
         expect_equal(attr(res$xdecomp$T2, "yaxis.values"), attr(d$xc, "yaxis.values"))
         expect_equal(attr(res$xdecomp$T2, "yaxis.name"), yaxis.name.x)

         # explained variance
         expect_lte(sum(res$xdecomp$expvar), 100)
      })

      test_that("dimension and attributes for ydecomp are correct", {
         # scores
         expect_equal( dim(res$ydecomp$scores), c(nrow(d$yc), ncomp))
         expect_equal(attr(res$ydecomp$scores, "exclrows"), attr(d$xc, "exclrows"))
         expect_equal(attr(res$ydecomp$scores, "yaxis.values"), attr(d$xc, "yaxis.values"))
         expect_equal(attr(res$ydecomp$scores, "yaxis.name"), yaxis.name.x)

         # othogonal distance
         expect_equal( dim(res$ydecomp$Q), c(nrow(d$yc), ncomp))
         expect_equal(attr(res$ydecomp$Q, "exclrows"), attr(d$xc, "exclrows"))
         expect_equal(attr(res$ydecomp$Q, "yaxis.values"), attr(d$xc, "yaxis.values"))
         expect_equal(attr(res$ydecomp$Q, "yaxis.name"), yaxis.name.x)

         # othogonal distance
         expect_equal( dim(res$ydecomp$T2), c(nrow(d$yc), ncomp))
         expect_equal(attr(res$ydecomp$T2, "exclrows"), attr(d$xc, "exclrows"))
         expect_equal(attr(res$ydecomp$T2, "yaxis.values"), attr(d$xc, "yaxis.values"))
         expect_equal(attr(res$ydecomp$T2, "yaxis.name"), yaxis.name.x)

         # explained variance
         expect_lte(sum(res$ydecomp$expvar), 100)
      })

      test_that("dimension and attributes for y.pred are correct", {
         expect_equal( dim(res$y.pred), c(nrow(d$yc), ncomp, ncol(d$yc)))
         expect_equal(attr(res$y.pred, "exclrows"), attr(d$xc, "exclrows"))
         expect_equal(attr(res$y.pred, "yaxis.values"), attr(d$xc, "yaxis.values"))
         expect_equal(attr(res$y.pred, "yaxis.name"), yaxis.name.x)
      })
   }
}

########################################################################
# Block 4. Tests predict.pls() by comparing with plsregres() outcomes  #
########################################################################

context("pls: compare pedict.pls() outcome with known values")

test_that("predictions for people data are correct", {
   # model
   f <- system.file("testdata", "pls-xloadings.csv", package = "mdatools")
   xloadings <- read.csv(f, header = FALSE)

   f <- system.file("testdata", "pls-yloadings.csv", package = "mdatools")
   yloadings <- read.csv(f, header = FALSE)

   f <- system.file("testdata", "pls-weights.csv", package = "mdatools")
   weights <- read.csv(f, header = FALSE)

   f <- system.file("testdata", "pls-coeffs.csv", package = "mdatools")
   coeffs <- read.csv(f, header = FALSE)

   xc <- people[, -4]
   yc <- people[, 4, drop = F]
   obj <- pls.cal(xc, yc, ncomp = 4, center = TRUE, scale = TRUE)
   expect_equivalent(obj$xloadings, as.matrix(xloadings), tolerance = 10^-5)
   expect_equivalent(obj$yloadings, as.matrix(yloadings), tolerance = 10^-5)
   expect_equivalent(obj$weights, as.matrix(weights), tolerance = 10^-5)
   expect_equivalent(obj$coeffs$values[, 4, ], as.matrix(coeffs), tolerance = 10^-5)

   # predictions
   f <- system.file("testdata", "pls-xscores.csv", package = "mdatools")
   xscores <- read.csv(f, header = FALSE)
   print(f)

   f <- system.file("testdata", "pls-yscores.csv", package = "mdatools")
   yscores <- read.csv(f, header = FALSE)
   print(f)

   f <- system.file("testdata", "pls-xres.csv", package = "mdatools")
   xresid <- read.csv(f, header = FALSE)
   print(f)

   f <- system.file("testdata", "pls-yres.csv", package = "mdatools")
   yresid <- read.csv(f, header = FALSE)
   print(f)

   f <- system.file("testdata", "pls-expvar.csv", package = "mdatools")
   expvar <- read.csv(f, header = FALSE)
   print(f)

   res <- predict(obj, xc, yc)
   expect_equivalent(res$xdecomp$scores, as.matrix(xscores), tolerance = 10^-4)
   expect_equivalent(res$ydecomp$scores, as.matrix(yscores), tolerance = 10^-4)
   expect_equivalent(res$xdecomp$residuals, as.matrix(xresid), tolerance = 10^-3)
   expect_equivalent(res$ydecomp$residuals, as.matrix(yresid), tolerance = 10^-3)
   expect_equivalent(res$xdecomp$expvar, as.matrix(expvar[1, ] * 100), tolerance = 10^-3)
   expect_equivalent(res$ydecomp$expvar, as.matrix(expvar[2, ] * 100), tolerance = 10^-3)
})

#####################################
# Block 4. Tests pls() constructor  #
#####################################

context("pls: test constructor")

for (i in seq_along(datasets)) {

   d <- datasets[[i]]
   name <- names(datasets)[i]

   test_that("test constructor with minimum settings", {
      expect_warning(m <- pls(d$xc, d$yc))
      expect_equal(m$calres, m$res[["cal"]])
      expect_equal(class(m$calres), c("plsres", "regres"))
      expect_equal(m$ncomp, m$ncomp.selected)
      expect_null(m$res[["cv"]])
      expect_null(m$res[["test"]])
   })


   test_that("test constructor with ncomp specified by user", {
      for (ncomp in seq_len(d$ncomp)) {
         if (ncomp > 1) {
            expect_warning(m <- pls(d$xc, d$yc, ncomp = ncomp))
         } else {
            expect_silent(m <- pls(d$xc, d$yc, ncomp = ncomp))
         }
         expect_equal(m$calres, m$res[["cal"]])
         expect_equal(class(m$calres), c("plsres", "regres"))
         expect_equal(m$ncomp, ncomp)
         expect_equal(m$ncomp, m$ncomp.selected)
         expect_null(m$res[["cv"]])
         expect_null(m$res[["test"]])
      }
   })

   expect_warning(m <- pls(d$xc, d$yc, ncomp = ncomp))
   fprintf("\n\nCalibration (%s): -------\n", name)
   #summary(m)
   #summary(m$calres)

   test_that("test constructor with cross-validation", {
      for (cv in list(1, list("rand", 4), list("rand", 4, 4), list("ven", 8))) {

         m <- pls(d$xc, d$yc, ncomp = d$ncomp, cv = cv)

         # check calibration results
         expect_equal(m$calres, m$res[["cal"]])
         expect_equal(class(m$calres), c("plsres", "regres"))
         expect_equal(m$calres$ncomp, m$ncomp)
         expect_equal(m$calres$ncomp.selected, m$ncomp.selected)
         expect_equal(class(m$calres$xdecomp), "ldecomp")
         expect_equal(class(m$calres$ydecomp), "ldecomp")
         expect_gte(max(m$calres$r2), 0.9)
         expect_lte(max(m$calres$r2), 1.0)

         # check cross-validation results
         expect_equal(m$cvres, m$res[["cv"]])
         expect_equal(class(m$cvres), c("plsres", "regres"))
         expect_equal(m$cvres$ncomp, m$ncomp)
         expect_equal(m$cvres$ncomp.selected, m$ncomp.selected)
         expect_equal(length(m$cvres$rmse), m$ncomp * ncol(d$yc))
         expect_null(m$cvres$xdecomp)
         expect_null(m$cvres$ydecomp)
         expect_gte(max(m$cvres$r2), 0.9)
         expect_lte(max(m$cvres$r2), 1.0)

         # check test set results (must be null)
         expect_null(m$res[["test"]])

         # summary
         fprintf("\n\nCross-validation (%s): -------\n", name)
         #summary(m)
         #summary(m$calres)
         #summary(m$cvres)
      }
   })
}

test_that("test constructor with test set", {

   d <- datasets[["spectra"]]
   m <- pls(d$xc, d$yc, ncomp = d$ncomp, x.test = d$xt, y.test = d$yt)

   # check calibration results
   expect_equal(m$calres, m$res[["cal"]])
   expect_equal(class(m$calres), c("plsres", "regres"))
   expect_equal(m$calres$ncomp, m$ncomp)
   expect_equal(m$calres$ncomp.selected, m$ncomp.selected)
   expect_equal(class(m$calres$xdecomp), "ldecomp")
   expect_equal(class(m$calres$ydecomp), "ldecomp")
   for (i in seq_len(ncol(d$yc))) {
      expect_gte(max(m$calres$r2[i, ]), 0.9)
      expect_lte(max(m$calres$r2[i, ]), 1.0)
   }

   # check cross-validation results
   expect_equal(m$testres, m$res[["test"]])
   expect_equal(class(m$testres), c("plsres", "regres"))
   expect_equal(m$testres$ncomp, m$ncomp)
   expect_equal(m$testres$ncomp.selected, m$ncomp.selected)
   expect_equal(class(m$testres$xdecomp), "ldecomp")
   expect_equal(class(m$testres$ydecomp), "ldecomp")
   expect_equal(length(m$testres$rmse), m$ncomp * ncol(d$yc))
   for (i in seq_len(ncol(d$yc))) {
      expect_gte(max(m$testres$r2[i, ]), 0.9)
      expect_lte(max(m$testres$r2[i, ]), 1.0)
   }

   # check test set results (must be null)
   expect_null(m$res[["cv"]])

   cat("\n\nTest set: -------\n")
   summary(m)
   summary(m$calres)
   summary(m$testres)

   m <- pls(d$xc, d$yc, ncomp = d$ncomp, x.test = d$xt, y.test = d$yt, cv = 1)
   cat("\n\nTest set: -------\n")
   summary(m)

})


######################################
# Block 5. Tests vipscores() method  #
######################################

context("pls: test vipscores()")

for (i in seq_along(datasets)) {

   d <- datasets[[i]]
   name <- names(datasets)[i]
   m <- pls(d$xc, d$yc, d$ncomp, center = d$center, scale = d$scale, cv = 10)

   xaxis.name <- if (is.null(attr(d$xc, "xaxis.name"))) "Predictors" else attr(d$xc, "xaxis.name")

   #returns correct outcome
   expect_silent(v <- vipscores(m))
   expect_equal(dim(v), c(ncol(d$xc), ncol(d$yc)))
   expect_equal(rownames(v), colnames(d$xc))
   expect_equal(colnames(v), colnames(d$yc))
   expect_equal(attr(v, "yaxis.values"), attr(d$xc, "xaxis.values"))
   expect_equal(attr(v, "yaxis.name"), xaxis.name)
   expect_equal(attr(v, "exclrows"), attr(d$xc, "exclcols"))

   # can work with different number of ny and components
   expect_silent(v <- vipscores(m, ncomp = 1))
   expect_silent(v <- vipscores(m, ncomp = m$ncomp))

   # raise error if ny or ncomp are not correct
   expect_error(vipscores(m, ncomp = 0))
   expect_error(vipscores(m, ncomp = 100))
   expect_error(vipscores(m, ncomp = 1:3))

   # check deprecated method
   expect_warning(v1 <- getVIPScores(m))
   expect_warning(v2 <- getVIPScores(m, ncomp = m$ncomp))
   expect_equal(v1, vipscores(m))
   expect_equal(v2, vipscores(m, ncomp = m$ncomp))

   # send output to file for visual checking
   print(v)
}

test_that("vipscores for people data (A = 4) identical to once computed in MATLAB", {
   d <- datasets[[1]]
   m <- pls(d$xc, d$yc, d$ncomp, center = d$center, scale = d$scale, cv = 10)

   f <- system.file("testdata", "pls-vipscores.csv", package = "mdatools")
   vip <- read.csv(f, header = FALSE)
   expect_equivalent(vipscores(m, ncomp = 4), as.matrix(vip), tolerance = 10^-4)
})


######################################
# Block 6. Tests selratio()  method  #
######################################

context("pls: test selratio()")

for (i in seq_along(datasets)) {

   d <- datasets[[i]]
   name <- names(datasets)[i]
   m <- pls(d$xc, d$yc, d$ncomp, center = d$center, scale = d$scale, cv = 10)

   xaxis.name <- if (is.null(attr(d$xc, "xaxis.name"))) "Predictors" else attr(d$xc, "xaxis.name")

   #returns correct outcome
   selratio(m)
   expect_silent(v <- selratio(m))
   expect_equal(dim(v), c(ncol(d$xc), ncol(d$yc)))
   expect_equal(rownames(v), colnames(d$xc))
   expect_equal(colnames(v), colnames(d$yc))
   expect_equal(attr(v, "yaxis.values"), attr(d$xc, "xaxis.values"))
   expect_equal(attr(v, "yaxis.name"), xaxis.name)
   expect_equal(attr(v, "exclrows"), attr(d$xc, "exclcols"))

   # can work with different number of ny and components
   expect_silent(v <- selratio(m, ncomp = 1))
   expect_silent(v <- selratio(m, ncomp = m$ncomp))

   # raise error if ny or ncomp are not correct
   expect_error(selratio(m, ncomp = 0))
   expect_error(selratio(m, ncomp = 100))
   expect_error(selratio(m, ncomp = 1:3))

   # check deprecated method
   expect_warning(v1 <- getSelectivityRatio(m))
   expect_warning(v2 <- getSelectivityRatio(m, ncomp = m$ncomp))
   expect_equal(v1, selratio(m))
   expect_equal(v2, selratio(m, ncomp = m$ncomp))

   # send output to file for visual checking
   print(v)
}

#test_that("vipscores for people data (A = 4) identical to once computed in MATLAB", {
#   d <- datasets[[1]]
#   m <- pls(d$xc, d$yc, d$ncomp, center = d$center, scale = d$scale, cv = 10)
#
#   vip <- as.matrix(read.csv("../matlab/pls-vipscores.csv", header = FALSE))
#   expect_equivalent(vipscores(m, ncomp = 4), vip, tolerance = 10^-4)
#})

#########################################
# Block 7. Tests for outlier detection  #
#########################################

context("pls: xy residuals and categorization")

test_that("XY-residual limits are computed correctly", {

   # test for people data
   d <- datasets[[1]]
   m1 <- pls(d$xc, d$yc, d$ncomp, center = d$center, scale = d$scale, cv = 1)
   m2 <- pls(d$xc, d$yc, d$ncomp, center = d$center, scale = d$scale, cv = 1, lim.type = "ddrobust")
   m3 <- pls(d$xc, d$yc, d$ncomp, center = d$center, scale = d$scale, cv = 1, lim.type = "chisq")

   expect_null(m3$Zlim)
   expect_equal(dim(m1$Zlim), c(4, d$ncomp))
   expect_equal(dim(m2$Zlim), c(4, d$ncomp))

   c1 <- categorize(m1, m1$res$cal)
   c2 <- categorize(m2, m2$res$cal)

   par(mfrow = c(1, 2))
   plotXYResiduals(m1, show.labels = T)
   plotXYResiduals(m2, show.labels = T)

   expect_error(c3 <- categorize(m3, m1$res$cal))
   expect_equal(sum(c1 == "extreme"), 1)
   expect_equal(sum(c1 == "outlier"), 0)

   expect_equal(sum(c2 == "extreme"), 4)
   expect_equal(sum(c2 == "outlier"), 1)

   # make two outliers and test again
   d <- datasets[[1]]
   d$yc[9] <- 25
   d$xc[1, 1] <- 115

   m1 <- pls(d$xc, d$yc, 4, center = d$center, scale = d$scale, cv = 1)
   m2 <- pls(d$xc, d$yc, 4, center = d$center, scale = d$scale, cv = 1, lim.type = "ddrobust")

   c1 <- categorize(m1, m1$res$cal, ncomp = 2)
   c2 <- categorize(m2, m1$res$cal, ncomp = 2)

   par(mfrow = c(1, 2))
   plotXYResiduals(m1, show.labels = T, ncomp = 2)
   plotXYResiduals(m2, show.labels = T, ncomp = 2)

   expect_equal(sum(c1 == "extreme"), 3)
   expect_equal(sum(c1 == "outlier"), 0)

   expect_equal(sum(c2 == "extreme"), 2)
   expect_equal(sum(c2 == "outlier"), 2)

})