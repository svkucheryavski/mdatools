######################################
# Tests for pls class methods        #
######################################

setup({
   pdf(file = "dump/mdatools-test-pls.pdf")
   sink("dump/mdatools-test-pls.txt", append = FALSE, split = FALSE)
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
   xloadings <- read.csv("auxfiles/pls-xloadings.csv", header = FALSE)
   yloadings <- read.csv("auxfiles/pls-yloadings.csv", header = FALSE)
   weights <- read.csv("auxfiles/pls-weights.csv", header = FALSE)
   coeffs <- read.csv("auxfiles/pls-coeffs.csv", header = FALSE)

   xc <- people[, -4]
   yc <- people[, 4, drop = F]
   obj <- pls.cal(xc, yc, ncomp = 4, center = TRUE, scale = TRUE)
   expect_equivalent(obj$xloadings, as.matrix(xloadings), tolerance = 10^-5)
   expect_equivalent(obj$yloadings, as.matrix(yloadings), tolerance = 10^-5)
   expect_equivalent(obj$weights, as.matrix(weights), tolerance = 10^-5)
   expect_equivalent(obj$coeffs$values[, 4, ], as.matrix(coeffs), tolerance = 10^-5)

   # predictions
   xscores <- read.csv("auxfiles/pls-xscores.csv", header = FALSE)
   yscores <- read.csv("auxfiles/pls-yscores.csv", header = FALSE)
   xresid  <- read.csv("auxfiles/pls-xres.csv", header = FALSE)
   yresid  <- read.csv("auxfiles/pls-yres.csv", header = FALSE)
   expvar  <- read.csv("auxfiles/pls-expvar.csv", header = FALSE)

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
   summary(m)
   summary(m$calres)

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

   vip <- read.csv("auxfiles/pls-vipscores.csv", header = FALSE)
   expect_equivalent(vipscores(m, ncomp = 4), as.matrix(vip), tolerance = 10^-4)
})

test_that("vipscores for people data (A = 4) identical to once computed in MATLAB for PLS2", {
   data(people)
   X <- people[, -c(4, 6)]
   Y <- people[,  c(4, 6)]
   m <- pls(X, Y, 8, center = TRUE, scale = TRUE, cv = 1)

   vip <- read.csv("auxfiles/pls2-vipscores.csv", header = FALSE)[[1]]
   #dim(vip) <- c(ncol(X), 2)

   expect_equivalent(vipscores(m, ncomp = 4), matrix(vip, ncol = 2), tolerance = 10^-4)
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

test_that("PLS gives results comparable to other software", {


   # read model parameters and calibration scores made in PLS_Toolbox
   xscores <-   as.matrix(read.delim("auxfiles/plstlbx-people-xscores.csv", sep = " ", header = FALSE))
   yscores <-   as.matrix(read.delim("auxfiles/plstlbx-people-yscores.csv", sep = " ", header = FALSE))
   weights <-   as.matrix(read.delim("auxfiles/plstlbx-people-weights.csv", sep = " ", header = FALSE))
   xloadings <- as.matrix(read.delim("auxfiles/plstlbx-people-xloadings.csv", sep = " ", header = FALSE))
   yloadings <- c(5.3643, 1.0338, 0.4675, 0.3567)
   coeffs <- c(0.2078, 0.2647, 0.0073, 0.0722, -0.0016, 0.1829, 0.1420, -0.1984, 0.2153, 0.0151, -0.0405)

   # make a model with manual cv splits
   data(people)
   X <- people[, -4]
   y <- people[, 4, drop = FALSE]
   m <- pls(X, y, 4, scale = TRUE, cv = rep(1:10, length.out = nrow(X)))

   # compare main model parameters
   # here we re-normalize results from PLS_Toolbox
   xnorm <- sqrt(colSums(xscores^2))
   expect_equivalent(m$xloadings, xloadings %*% diag(xnorm), tolerance = 10^-3)
   expect_equivalent(m$res$cal$xdecomp$scores, xscores %*% diag(1/xnorm), tolerance = 10^-3)
   expect_equivalent(m$weights, weights %*% diag(1/xnorm), tolerance = 10^-3)

   expect_equivalent(m$yloadings, yloadings, tolerance = 10^-4)
   expect_equivalent(m$coeffs$values[, 4, 1], coeffs, tolerance = 10^-4)

   # check selectivity ratio
   # here we change the result from PLS toolbox a bit for large values as they add
   # given portion of x-variance when compute SR
   sr <- c(24.8, 21.7, 2.2359, 0.1179, 0.1552, 0.9740, 0.0076, 5.9018, 10.0, 0.0256, 0.0138)
   expect_equivalent(selratio(m, 4), sr, tolerance = 10^-1)

   # compare calibration results
   ypred <-  as.matrix(read.delim("auxfiles/plstlbx-people-ypred.csv", sep = " ", header = FALSE))
   xqdist <- as.matrix(read.delim("auxfiles/plstlbx-people-xqdist.csv", sep = " ", header = FALSE))
   xhdist <- as.matrix(read.delim("auxfiles/plstlbx-people-xhdist.csv", sep = " ", header = FALSE))
   yqdist <- as.matrix(read.delim("auxfiles/plstlbx-people-yqdist.csv", sep = " ", header = FALSE))
   rmsec <- c(1.0273, 0.7404, 0.6668, 0.6198)
   r2c <- c(0.9283, 0.9627, 0.9698, 0.9739)

   expect_equivalent(m$res$cal$xdecomp$T2[, 4], xhdist, tolerance = 10^-3)
   expect_equivalent(m$res$cal$xdecomp$Q[, 4], xqdist, tolerance = 10^-3)

   expect_equivalent(m$res$cal$ydecomp$Q[, 4], yqdist, tolerance = 10^-3)
   expect_equivalent(m$res$cal$ydecomp$scores, yscores, tolerance = 10^-4)
   expect_equivalent(m$res$cal$y.pred[, 4, 1], ypred, tolerance = 10^-4)

   expect_equivalent(m$res$cal$rmse, rmsec, tolerance = 10^-3)
   expect_equivalent(m$res$cal$r2, r2c, tolerance = 10^-3)


   # compare cross-validation results
   rmsecv <- c(1.1044, 0.8673, 0.8627, 0.8186)
   r2cv <- c(0.9173, 0.9489, 0.9498, 0.9545)
   expect_equivalent(m$res$cv$rmse, rmsecv, tolerance = 10^-3)
   expect_equivalent(m$res$cv$r2, r2cv, tolerance = 10^-3)
})


###########################################
# Block 8: testing JSON methods in PLS    #
###########################################

context("pls: testing JSON methods")

# expected vector for case 1: People (Shoesize is response), A = 5, center, scale
v.exp1 <- c(1,11,5,39.90625,173.125,64.46875,0,34.4375,27437.5,249.5,131.59375,0,81.5,0,115.125,3.896726293499663,10.057095072152924,15.191220671291395,1.016001016001524,9.517174628177514,8929.608234732721,90.59694504569622,49.49673131925301,1.016001016001524,7.317676345451078,1.016001016001524,12.164862143392414,-5.353200818268451,-5.414932764105852,4.8588992441366,-1.9837488956598004,-2.559398676321643,-4.425206288227119,1.5130612621713437,5.066417074564808,-5.193996368200587,1.758374392095087,0.5721361976664939,-0.9383004890009083,-0.6297524284595666,0.7454867418184986,0.6938069030899766,2.583122841623663,2.7503970924682357,-4.818715118060619,1.5421145822609588,-0.7896256745597711,-4.976471531850274,-0.22215404578300857,0.5168729619723487,0.03279725544254994,1.8997770822859827,-0.22954686265778762,0.0523099593965144,0.7809987245698705,0.4347331447717977,0.8765740013734368,0.7806227880341323,-1.0144621116160533,-4.2161600925816245,-0.18280256298618333,-0.5566358689282596,-1.4148376310352047,-0.34783331634116615,-0.375534485287149,-0.6727319101878396,-0.7805979540763109,-0.1642152073772969,-0.368498450755976,0.19667639054109168,-3.4835064352648946,-0.27980897031530133,0.49351359405100237,-0.17778635658665923,-1.7784249496414082,-2.045398518427989,0.6895546763711853,-0.9767201384735575,-0.9467024401267798,-0.6000513760117836,-0.20518092550430247,-0.2706606509668198,0.03225806451612902,0.03225806451612903,0.03225806451612902,0.032258064516129004,0.03225806451612903,5.343432980475502,3.231320948658815,2.45273783899078,1.9525754085543192,1.6244179535093726,11,4,3,2,1,5.555402075957791,3.1326839251404173,2.2512148971963244,1.7040524426157482,1.4785368769701401,22,6,5,3,9,0.9687500000000006,1.9375000000000002,2.90625,3.875000000000001,4.843750000000001,3,8,11,9,7,1.1706531603048451,1.915232854779818,2.760792914543984,3.6362575078957726,4.081956330384955,4,19,17,27,25,0.06950473825612581,0.036106824973159535,0.029277447672190437,0.0253005494182464,0.022930272365119995,1,1,1,1,1,0.05624687662837766,0.03939893980553692,0.02397824904440559,0.025208672030409933,0.021738646785584577,1,1,1,1,1,-0.03312651550807127,-0.03340998085162388,0.02837038012165645,-0.011134467572355478,-0.012218672244370633,-0.024027389360909347,0.003062913972119423,0.03174111835393245,-0.03215550881658082,0.005337002201103687,0.004143648487157818,-0.02254805107823774,-0.02145238426527492,-0.0028845749015263067,0.007345718469255104,0.039012034308196386,0.027537919066542755,-0.08030069370019638,0.026695426538091253,-0.022072515936632023,-0.06886616128460553,0.010707954546391237,0.023015118380384985,0.053170119441992535,0.1685401966819274,0.007797598705645983,-0.0267098795722126,0.07731599267768807,0.08439341073117022,0.016981871250269037,0.039794700777265046,-0.033422451824542995,-0.11679583650856114,0.011116973357864239,-0.10789587693351882,-0.217903433097581,-0.045938511388928384,0.04008716045220307,-0.12977545847951552,-0.10073803319498592,0.023824735677888658,-0.003996854915263689,0.033085893335042735,-0.13283828643194703,-0.20093393775314733,0.20539079381270028,0.1040611309782249,0.10223912470636444,-0.24663220147993942,0.19929154894222675,-0.016789009546850434,-0.2416219423533065,-0.24075429510131735,0.005487028769566186,-0.04571804154673186,-5.364312479321461,-1.0337955431587753,0.4674826987505736,-0.35673623887460504,0.2754067277693107,28.775848375803932,0.076678217365749,0.00814532974267895,0.003846059124526113,0.0019809217806203293)

# expected vector for case 2: People, A = 5, center, scale, excluded outliers: 1, 18
v.exp2 <- c(1,11,5,39.43333333333333,171.83333333333334,62.9,0.06666666666666667,34.233333333333334,27216.666666666668,242.16666666666666,130.6,0.06666666666666667,80.5,0,115.5,3.5300027677368577,8.855480801734313,14.312642108916544,1.0148325268098497,9.409032099927774,8391.994680208258,87.3574043338859,50.339947803186945,1.0148325268098497,6.366669675632676,1.0170952554312156,12.232658379701558,5.167036534594553,5.228030829190363,-4.728882404241247,1.7482431875408604,2.3821674959524066,4.110017964842157,-1.7244044385613206,-4.959047772220314,5.019700028929209,-1.7394937100283532,-0.10381546729864008,-0.9239374541025428,-0.5439467373861087,0.69206823131495,-0.15632661160136416,1.849442924832116,3.0201957950983007,-4.668643933803133,1.6346263466879356,-0.5153580191925352,-4.8417823479112325,0.16720868427653526,0.38630008334584987,0.143213505757451,1.4283003305027493,-0.9387133635029506,-0.6531523154642436,0.9010337218756459,0.010806545879369911,0.45061586555408656,0.6525403414273143,-1.1649005460660424,-4.379234007350946,0.16439174881513277,-0.37061015448166934,-1.6388466833964428,-2.165361785879482,-1.9587896034611723,-0.26143133574529526,-1.2561330361143428,0.12337205008863622,0.24208224815951424,0.08243570075402926,-2.2964200523635943,-0.5917651400335773,-0.7263178802406306,0.15274331856615986,4.35815005333933,3.8902076009714484,-1.1150528339135881,0.5657950174794601,-0.030805033853458424,-0.7080611306579396,0.3177169191631606,-1.9617632033928485,0.03448275862068965,0.034482758620689655,0.03448275862068966,0.03448275862068966,0.034482758620689655,5.372801484645785,3.292814647914144,2.4430447712145447,1.8305039349899273,0.46238838573820135,9,5,3,3,4,5.718482199998322,3.314958528234614,2.3849666134787815,1.5971405374968821,0.4772422677319428,16,6,3,5,2,0.9666666666666666,1.9333333333333327,2.8999999999999995,3.8666666666666667,4.833333333333331,4,12,27,9,11,1.2562488181064069,2.019711961263013,2.7517947059553127,3.3876081168678143,4.520022926668442,4,21,40,15,14,0.07660898519947587,0.04314451722070581,0.0357266085477844,0.029951238741401775,0.028751555307275647,1,1,1,1,1,0.06276186997334957,0.04707339179297959,0.021906878794163878,0.022818557929763526,0.03145819674262024,1,1,2,1,1,0.034341739394014484,0.03481460747842039,-0.029810550693857057,0.011382057735719016,0.012749569811072062,0.02356075561749924,-0.0048558057834932455,-0.03397662144299859,0.03289874375207857,-0.0058888035482204315,-0.002088127113734265,-0.02212756948502009,-0.023321528810950707,-0.0021285057377206642,-0.004206509478170949,0.03241403926367183,0.0343122800202062,-0.08391457241444747,0.03529765867018251,-0.015086414239753309,-0.07095758604654538,0.019770314391770877,0.008662611217330654,0.06877000555071348,0.15595770168571824,-0.0031847345171661405,-0.042661447747177964,0.08165954425509157,0.08338751412755738,-0.024092930165792788,0.006747247165468441,-0.04641919258671079,-0.13954434807513186,0.03032917688239172,-0.13163503707586777,-0.23211872003165498,-0.07574226629677394,0.017625689798944216,-0.1317622044481019,-0.13643297980575894,0.07344972964741295,0.05439341511088944,0.050302303575883944,-0.08889344722930773,0.03692976230648283,-0.0983653799716417,-0.07013129306582033,0.05597300942802837,0.1376947299402081,-0.09616317477107651,-0.02587458183105663,0.07300600951630035,0.06296486894835229,0.02543073070709012,-0.06508751070107231,5.167371715293544,-1.001965088894372,0.4717385506693983,-0.416246434449021,-0.1897116312296237,26.70173044401576,0.07956244961670857,0.009932340676981424,0.0064034806435255345,0.001115131188092191)

# expected vector for case 3: People, A = 5, no center, no scale, excluded outliers: 1, 18, preprocessing: median centering, pareto scaling
v.exp3 <- c(1,11,5,39.43333333333333,0.280035236167282,0.23789351514354326,-0.9264875492561486,0.07606840502256676,-30.383137623926928,-0.0891597605201777,0.5778664440307266,-0.9264875492561486,1.4802973661668754e-17,0,0.42887510298886755,1,1,1,1,1,1,1,1,1,1,1,1,-4.369072431556267,-7.01665462442142,1.4946687246357,-14.398358374198796,-493.1584702035098,-20.001383840408693,10.772491032136443,1.9250920072320379,-3.447010733436216,2.240761191209107,0.3261108795607687,-12.750119250486966,-16.14658799493031,3.595646265686174,4.062234805217164,12.727329907745364,-44.162197811833984,20.000959959585447,3.3006823925062996,-11.087305829754664,2.8681067342080206,0.6479397901855253,-7.46309696037699,-8.482546147761553,2.009138815544004,-4.510451418811089,-0.7428804971795553,11.23933153140222,-29.519478299918504,3.079060381211827,-5.3591338200516025,-3.530728927094598,-0.3626279622073626,3.528979790842163,4.179652238026205,-1.1569612462579992,1.8150191535194387,-0.06493049051489663,-6.935946670022737,-8.023915357505674,-1.5019539623975513,2.679030644964855,1.2015852046168392,-6.96174501432534,-0.937131269093622,-3.0099582783114998,0.9856753077658834,-1.1189975210743357,0.038223604169188796,2.3728264313356653,2.70852472527337,0.9494297970863727,-0.15113029322930815,-1.0720531287955617,-17.21237193930726,0.03448275862068965,0.03448275862068964,0.03448275862068964,0.034482758620689655,0.03448275862068966,163.7860042700048,60.202670038449874,20.165845574994552,13.286058693638937,2.503760339808396,4,5,5,3,3,154.16872835859564,54.87418685650802,23.58187359375835,13.28199199735387,2.4376521814535357,6,3,5,2,2,0.9666666666666667,1.9333333333333333,2.9,3.866666666666668,4.833333333333334,2,5,11,9,13,0.9313036909798958,2.124833003703574,2.7988684002385025,3.8374724349512515,4.869125548994928,2,5,9,12,17,10.38044120576084,3.91285138980401,0.8114504041268893,0.44229622444365974,0.41807628547112935,2,2,1,1,1,11.02974994574613,6.223934113817969,0.8654873301013054,0.32040611678597053,0.3918379582467818,3,1,1,1,1,-0.00017527567908591244,-0.00022589920690953231,0.00005150631569498668,-0.00005988065698751396,-0.00200318670199435,-0.0003776872920908968,0.000059089591248803866,0.00005870440329211213,-0.00014237337189983303,0.000010185941506544799,0.000012525949417072818,-0.007663893149108617,-0.009600850632009246,0.0022100377692002614,-0.00004265194085886927,0.0008663261365817327,-0.014401468128627472,0.0007277481523796677,0.002474680651367761,-0.006245584961659324,0.000048857536112559575,0.0005449248579029361,-0.012484501991185776,-0.015387698816691522,0.003693647429258069,-0.004773039431238275,0.000002701575658878823,0.0015957688553433222,-0.022102805916113567,0.005018761966723686,-0.009323029848274805,-0.0033913858099069313,0.0012042497503043244,0.019384753029772033,0.032834813502611006,-0.009460749033859937,0.018007884912950526,-0.000594609363138521,-0.03682558229269142,-0.03605414186488826,-0.01409643324951626,0.01125803765490238,0.008832833223418692,-0.02068771997702553,0.005911441078292689,-0.02661405946358808,0.009285260279085038,-0.026689970251149445,0.0009690219350314108,0.009408065345737019,0.014719350165296342,0.018104917081882974,0.016319776114622337,-0.006690318172931741,-0.04661127164517689,-7.067774083390152,-13.929382415552562,-9.645829646552624,3.327855975022992,-0.8524072789314314,622.465678142225,2083.5445601643346,376.61341522177963,9.29638887793663,0.33245340715231547)

# expected vector for case 4: People, A = 5, no center, no scale, excluded outliers: 1, 18, preprocessing: median centering, pareto scaling, excluded variables: 2, 7, 8
v.exp4 <- c(1,8,5,39.43333333333333,0.280035236167282,0.23789351514354326,-0.9264875492561486,-30.383137623926928,-0.0891597605201777,0.5778664440307266,-0.9264875492561486,1.4802973661668754e-17,1,1,1,1,1,1,1,1,1,-4.36955024274344,-7.018142796581476,1.4949997409444327,-493.1556709001425,-20.013465713827593,10.782512593145526,1.924904121530639,-3.448084094666105,-12.747693666082228,-16.158698146609392,3.600278783876413,12.829256858584078,-44.154007557121254,19.986200996070963,3.30301182077693,-11.084070877113357,-7.393702516719928,-8.451795755514754,2.0006107048299233,-0.8602159449864306,11.173425164041335,-29.631844029966498,3.0504174440649465,-5.341686322204416,3.747870774400014,5.151131195727044,-1.464515383728583,-0.006101844099083653,-7.630798573020698,-8.124585905208384,-1.7022244101896449,2.9380031519080387,-0.9841027677223171,1.8665919800459605,-0.887211675928654,-0.011818061864496077,0.41569159626551166,0.2989545469954821,-1.0899023607314926,-3.3328843779148167,0.03448275862068966,0.034482758620689655,0.03448275862068965,0.03448275862068965,0.03448275862068965,149.0314739746212,46.23313442897461,7.180721998189546,1.2309572003473739,0.6376876383889428,4,3,2,3,1,154.76434221997815,46.33296243531861,6.153705594937575,1.184183795894525,0.5124310517232593,5,2,1,2,2,0.9666666666666668,1.933333333333334,2.9000000000000004,3.8666666666666667,4.833333333333333,2,5,11,10,13,0.9315661954745584,2.129968374580221,2.882325802183152,3.851172444095112,4.531213324896707,2,5,10,13,11,10.380479989501993,3.919343015164548,0.8799781195760324,0.5216337404425826,0.46489623966586896,2,2,1,2,2,11.030058995220564,6.107678400192306,0.9930031461686298,0.6220967739368077,0.5172830724182558,3,1,1,1,1,-0.0001754316317881644,-0.00022610020223265894,0.00005155214378224652,-0.0020049690506973796,-0.00037802334187324543,0.000059142166605968534,0.00005875663592572815,-0.00014250004955524022,-0.007667172751819764,-0.009604784604697101,0.0022109454011157996,0.0008651556655924671,-0.0144054216871509,0.0007260885933826242,0.002475808990892457,-0.006248120839942676,-0.012693190889394756,-0.015620717152203662,0.003750080049319831,-0.00015090798563771005,0.0020546695899885453,-0.022859549946688266,0.005115846578605996,-0.00946539382776875,0.023591347549080616,0.03753139591646412,-0.010719923584646355,-0.00000898729812141286,-0.04309491991858265,-0.03768426112404348,-0.016366951510305865,0.013504441349440703,-0.04816331073821331,0.09570677050085444,-0.03522314104763689,-0.0007321735827161008,0.009952502323749323,0.00267658439155183,-0.14418522099284975,-0.174207263212145,-7.06769177183095,-13.922431871987138,-9.548871497075211,3.2787697958231585,1.3046551357686353,622.4511797098855,2081.4735435995567,369.69249022948435,9.78626522825556,0.9185025957914107)

compareResults <- function(rr, rjs) {

   expect_equal(abs(rjs$xdecomp$scores), abs(rr$xdecomp$scores))
   expect_equal(rjs$xdecomp$Q, rr$xdecomp$Q)
   expect_equal(rjs$xdecomp$T2, rr$xdecomp$T2)
   expect_equal(rjs$xdecomp$expvar, rr$xdecomp$expvar)
   expect_equal(rjs$xdecomp$cumexpvar, rr$xdecomp$cumexpvar)

   expect_equal(abs(rjs$ydecomp$scores), abs(rr$ydecomp$scores))
   expect_equal(rjs$ydecomp$Q, rr$ydecomp$Q)
   expect_equal(rjs$ydecomp$T2, rr$ydecomp$T2)
   expect_equal(rjs$ydecomp$expvar, rr$ydecomp$expvar)
   expect_equal(rjs$ydecomp$cumexpvar, rr$ydecomp$cumexpvar)

   rr.xres <- mda.purgeCols(rr$xdecomp$residuals)
   rjs.xres <- mda.purgeCols(rjs$xdecomp$residuals)
   rr.yres <- mda.purgeCols(rr$ydecomp$residuals)
   rjs.yres <- mda.purgeCols(rjs$ydecomp$residuals)

   attr(rr.xres, "exclcols") <- NULL
   attr(rr.xres, "prep:center") <- NULL
   attr(rr.xres, "prep:scale") <- NULL

   attr(rr.yres, "exclcols") <- NULL
   attr(rr.yres, "prep:center") <- NULL
   attr(rr.yres, "prep:scale") <- NULL

   attr(rjs.xres, "exclcols") <- NULL
   attr(rjs.xres, "prep:center") <- NULL
   attr(rjs.xres, "prep:scale") <- NULL

   attr(rjs.yres, "exclcols") <- NULL
   attr(rjs.yres, "prep:center") <- NULL
   attr(rjs.yres, "prep:scale") <- NULL

   expect_equal(rr.xres, rjs.xres)
   expect_equal(rr.yres, rjs.yres)
   expect_equal(rr$y.pred, rjs$y.pred)
}

compareModels <- function(mr, mjs) {

   if (is.null(attr(mr$xloadings, "yaxis.values"))) {
      attr(mjs$xloadings, "yaxis.values") <- NULL
   }


   expect_equal(mjs$ncomp, mr$ncomp)
   expect_equal(mjs$ncomp.selected, mr$ncomp.selected)
   expect_equal(mjs$method, mr$method)
   expect_equal(mjs$center, mr$center)
   expect_equal(mjs$scale, mr$scale)

   attr(mjs$weights, "yaxis.values") <- NULL
   attr(mjs$coeffs$values, "yaxis.values") <- NULL

   expect_equal(abs(mjs$yloadings), abs(mda.purge(mr$yloadings)))
   expect_equal(abs(mjs$xloadings), abs(mda.purge(mr$xloadings)))
   expect_equal(abs(mjs$weights), abs(mda.purge(mr$weights)))

   if (length(mr$exclcols) > 0 && dim(mjs$coeffs$values)[1] < dim(mr$coeffs$values)[1]) {
      attr(mr$coeffs$values, "exclcols") <- NULL
      expect_equal(mjs$coeffs$values, mr$coeffs$values[-mr$exclcols, , , drop = FALSE], check.attributes = FALSE)
   } else {
      expect_equal(mjs$coeffs$values, mr$coeffs$values)
   }

   expect_equal(mr$xeigenvals, mjs$xeigenvals)
   expect_equal(mr$yeigenvals, mjs$yeigenvals)


   # we need to round DoF in R as they are rounded in web version
   mr$limParams$Q$moments$Nu <- round(mr$limParams$Q$moments$Nu)
   mr$limParams$Q$robust$Nu <- round(mr$limParams$Q$robust$Nu)
   mr$limParams$T2$moments$Nu <- round(mr$limParams$T2$moments$Nu)
   mr$limParams$T2$robust$Nu <- round(mr$limParams$T2$robust$Nu)
   mr$limParams$Z$moments$Nu <- round(mr$limParams$Z$moments$Nu)
   mr$limParams$Z$robust$Nu <- round(mr$limParams$Z$robust$Nu)

   testList(mr$limParams, mjs$limParams)

   # same for Qlim/T2lim
   expect_equal(mr$Qlim, mjs$Qlim)
   expect_equal(mr$T2lim, mjs$T2lim)
   expect_equal(mr$Zlim, mjs$Zlim)

   # prep
   ncr <- ncol(mr$xloadings)
   ncjs <- ncol(mjs$xloadings)
   if (!is.null(mr$prep) && ncr == ncjs) {
      for (n in seq_along(mr$prep)) {
         if (!is.list(mr$prep[[n]])) next
         testList(mr$prep[[n]]$params, mjs$prep[[n]]$params)
      }
   }
}

test_that("JSON methods work correctly", {
   data(people)
   X <- people[, -4]
   y <- people[, 4, drop = FALSE]

   # case 1
   m <- pls(X, y, 5, scale = TRUE, cv = 1)
   m <- selectCompNum(m, 4)

   v <- asvector(m)
   expect_equivalent(abs(v), abs(v.exp1), tolerance = 0.00001)

   writeJSON(m, "dump/pls-model-1-r.json")
   m1 <- readJSON("dump/pls-model-1-r.json")
   m2 <- readJSON("jsonfiles/pls-model-1.json")
   compareModels(m, m2)
   compareModels(m, m1)

   r  <- predict(m, X, y)
   r1 <- predict(m1, X, y)
   r2 <- predict(m2, X, y)
   compareResults(r, r1)
   compareResults(r, r2)


   # case 2
   m <- pls(X, y, 5, scale = TRUE, cv = 1, exclrows = c(1, 18))
   m <- selectCompNum(m, 3)
   v <- asvector(m)

   expect_equivalent(abs(v), abs(v.exp2), tolerance = 0.0001)

   writeJSON(m, "dump/pls-model-2-r.json")
   m1 <- readJSON("dump/pls-model-2-r.json")
   m2 <- readJSON("jsonfiles/pls-model-2.json")
   compareModels(m, m1)
   compareModels(m, m2)

   r  <- predict(m, X, y)
   r1 <- predict(m1, X, y)
   r2 <- predict(m2, X, y)
   compareResults(r, r1)
   compareResults(r, r2)

   # case 3
   p <- list(
      prep("center", type = "median"),
      prep("scale", type = "pareto")
   )

   m <- pls(X, y, 5, center = TRUE, scale = FALSE, exclrows = c(1, 18), prep = p, cv = 1)
   m <- selectCompNum(m, 3)
   v <- asvector(m)
   expect_equivalent(abs(v), abs(v.exp3), tolerance = 0.0001)

   writeJSON(m, "dump/pls-model-3-r.json")
   m1 <- readJSON("dump/pls-model-3-r.json")
   m2 <- readJSON("jsonfiles/pls-model-3.json")
   compareModels(m, m2)
   compareModels(m, m1)

   r  <- predict(m, X, y)
   r1 <- predict(m1, X, y)
   r2 <- predict(m2, X, y)
   compareResults(r, r1)
   compareResults(r, r2)


   # case 4
   p <- list(
      prep("center", type = "median"),
      prep("scale", type = "pareto")
   )

   m <- pls(X, y, 5, center = TRUE, scale = FALSE, exclrows = c(1, 18), exclcols = c(4, 10, 11), prep = p, cv = 1)
   m <- selectCompNum(m, 4)
   v <- asvector(m)
   expect_equal(abs(v), abs(v.exp4), tolerance = 0.0001)

   writeJSON(m, "dump/pls-model-4-r.json")
   m1 <- readJSON("dump/pls-model-4-r.json")
   m2 <- readJSON("jsonfiles/pls-model-4.json")
   compareModels(m, m2)
   compareModels(m, m1)

   r  <- predict(m, X, y)
   r1 <- predict(m1, X, y)
   r2 <- predict(m2, X, y)
   compareResults(r1, r2)

   plot(r1)
   plot(r2)
   plot(r)
   summary(m1)
   summary(m2)
   summary(m)
   summary(m1$prep)
   summary(m2$prep)
   summary(m$prep)

})
