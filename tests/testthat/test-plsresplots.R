######################################
# Tests for plsres plotting methods  #
######################################

setup({
   pdf(file = tempfile("mdatools-test-plsresplots-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-plsresplots-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

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
attr(xc, "xaxis.values") <- simdata$Wavelength

xt <- simdata$spectra.t
yt <- simdata$conc.t
xt <- mda.exclrows(xt, c(15, 35))
yt <- mda.exclrows(yt, c(15, 35))
attr(xt, "name") <- "Spectra, val"
attr(xt, "xaxis.name") <- "Wavelength, nm"
attr(xt, "xaxis.values") <- simdata$Wavelength

datasets[["spectra"]] <- list(xc = xc, yc = yc, xt = xt, yt = yt, center = T, scale = F, ncomp = 6)

for (i in seq_along(datasets)) {

   d <- datasets[[i]]
   name <- names(datasets)[i]
   context(sprintf("plsres: test plots for %s", name))

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("PLS results - plots - ", name), pos = 4)

   m <- pls(d$xc, d$yc, ncomp = d$ncomp, center = d$center, scale = d$scale, cv = 10)
   r <- m$calres

   test_that("X-variance plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXVariance(r))
      expect_silent(plotXVariance(r, col = "red", type = "h", show.labels = T, main = "Myvar"))
      expect_silent(plotXCumVariance(r))
      expect_silent(plotXCumVariance(r, col = "red", type = "h", show.labels = T, main = "Mycumvar"))
   })

   test_that("X-variance plot can return plot data", {
      expect_equivalent(plotXVariance(r, show.plot = F), r$xdecomp$expvar)
      expect_equivalent(plotXCumVariance(r, show.plot = F), r$xdecomp$cumexpvar)
      expect_null(plotXVariance(m$cvres, show.plot = F))
      expect_null(plotXCumVariance(m$cvres, show.plot = F))
   })

   test_that("Y-variance plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotYVariance(r))
      expect_silent(plotYVariance(r, col = "red", type = "h", show.labels = T, main = "Myvar"))
      expect_silent(plotYCumVariance(r))
      expect_silent(plotYCumVariance(r, col = "red", type = "h", show.labels = T, main = "Mycumvar"))
   })

   test_that("Y-variance ploy can return plot data", {
      expect_equivalent(plotYVariance(r, show.plot = F), r$ydecomp$expvar)
      expect_equivalent(plotYCumVariance(r, show.plot = F), r$ydecomp$cumexpvar)
      expect_null(plotYVariance(m$cvres, show.plot = F))
      expect_null(plotYCumVariance(m$cvres, show.plot = F))
   })

   test_that("X-scores plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXScores(r))
      expect_silent(plotXScores(r, c(1, 3), show.labels = T, show.excluded = T, col = "red"))
      expect_silent(plotXScores(r, c(1, 3), type = "l", show.labels = T, show.excluded = T))
      expect_silent(plotXScores(r, c(1, 3), type = "h", show.labels = T, show.excluded = T))
   })

   test_that("X-scores plot can return plot data", {
      expect_equivalent(plotXScores(r, show.plot = F), r$xdecomp$scores[, 1:2])
      expect_equivalent(plotXScores(r, c(1, 3), show.plot = F), r$xdecomp$scores[, c(1, 3)])
      expect_equivalent(plotXScores(r, 1, show.plot = F), r$xdecomp$scores[, 1])
      expect_null(plotXScores(m$cvres, show.plot = F))
   })

   test_that("X-residuals plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXResiduals(r))
      expect_silent(plotXResiduals(r, 1, show.labels = T, show.excluded = T, col = "red"))
   })

   test_that("X-residuals plot can return plot data", {
      pd_exp <- cbind(r$xdecomp$T2[, r$ncomp.selected], r$xdecomp$Q[, r$ncomp.selected])
      expect_equivalent(plotXResiduals(r, norm = F, show.plot = F), pd_exp)
      expect_null(plotXResiduals(m$cvres, show.plot = F))

      pd_exp <- cbind(r$xdecomp$T2[, 1], r$xdecomp$Q[, 1])
      expect_equivalent(plotXResiduals(r, ncomp = 1, norm = F, show.plot = F), pd_exp)
      expect_null(plotXResiduals(m$cvres, ncomp = 1, show.plot = F))
   })

   test_that("XY-scores plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXYScores(r))
      expect_silent(plotXYScores(r, 2, show.labels = T, show.excluded = T, col = "red"))
   })

   test_that("X-residuals plot can return plot data", {
      pd_exp <- cbind(r$xdecomp$scores[, 1], r$ydecomp$scores[, 1])
      expect_equivalent(plotXYScores(r, show.plot = F), pd_exp)
      expect_null(plotXYScores(m$cvres, show.plot = F))

      pd_exp <- cbind(r$xdecomp$scores[, 2], r$ydecomp$scores[, 2])
      expect_equivalent(plotXYScores(r, 2, show.plot = F), pd_exp)
      expect_null(plotXYScores(m$cvres, 2, show.plot = F))
   })

   test_that("RMSE plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotRMSE(r))
      expect_silent(plotRMSE(r, ny = ncol(d$yc)))
      expect_silent(plotRMSE(r, col = "red", type = "h", show.labels = T, main = "My RMSE"))
   })

   test_that("RMSE plot can return plot data", {
      expect_equivalent(plotRMSE(r, show.plot = F), r$rmse[1, ])
      expect_equivalent(plotRMSE(r, ny = ncol(d$yc), show.plot = F), r$rmse[ncol(d$yc), ])
   })

   test_that("Y-residuals plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotYResiduals(r))
      expect_silent(plotYResiduals(r, ncomp = 1, show.labels = T, show.excluded = T, col = "red"))
      expect_silent(plotYResiduals(r, ny = ncol(d$yc)))
      expect_silent(plotYResiduals(r, ny = ncol(d$yc), ncomp = 1))
   })

   test_that("Y-residuals plot can return plot data", {
      pd_exp <- cbind(r$y.ref[, 1], r$y.ref[, 1] - r$y.pred[, r$ncomp.selected, 1])
      expect_equivalent(plotYResiduals(r, show.plot = F), pd_exp)

      pd_exp <- cbind(r$y.ref[, ncol(d$yc)], r$y.ref[, ncol(d$yc)] - r$y.pred[, 1, ncol(d$yc)])
      expect_equivalent(plotYResiduals(r, ncomp = 1, ny = ncol(d$yc), show.plot = F), pd_exp)
   })

   test_that("predictions plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotPredictions(r))
      expect_silent(plotPredictions(r, ncomp = 1, ny = ncol(d$yc)))
      expect_silent(plotPredictions(r, show.labels = T, show.excluded = T))
      expect_silent(plotPredictions(r, ncomp = 1, ny = ncol(d$yc), col = "red", show.line = F))
   })

   cat("\n\n output test -----\n")
   print(r)
   summary(r)
}
