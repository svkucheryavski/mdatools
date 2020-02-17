######################################
# Tests for pls plotting methods     #
######################################

setup({
   pdf(file = tempfile("mdatools-test-plsplots-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-plsplots-", fileext = ".txt"), append = FALSE, split = FALSE)
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
   context(sprintf("pls: test plots for %s", name))

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("PLS - plots - ", name), pos = 4)

   m1 <- pls(d$xc, d$yc, ncomp = d$ncomp, center = d$center, scale = d$scale, cv = list("ven", 8))
   m2 <- pls(d$xc, d$yc, ncomp = d$ncomp, center = d$center, scale = d$scale, x.test = d$xt, y.test = d$yt, cv = 1)

   test_that("Variance (X) plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXVariance(m1))
      expect_silent(plotXVariance(m2))
      expect_silent(plotXVariance(m1, type = "h", show.labels = T, legend.position = "bottom"))
      expect_silent(plotXVariance(m2, type = "h", show.labels = T, legend.position = "bottom"))
   })

   test_that("Cumulative variance (X) plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXCumVariance(m1))
      expect_silent(plotXCumVariance(m2))
      expect_silent(plotXCumVariance(m1, type = "h", show.labels = T, legend.position = "bottom"))
      expect_silent(plotXCumVariance(m2, type = "h", show.labels = T, legend.position = "bottom"))
   })

   test_that("Variance (Y) plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotYVariance(m1))
      expect_silent(plotYVariance(m2))
      expect_silent(plotYVariance(m1, type = "h", show.labels = T, legend.position = "bottom"))
      expect_silent(plotYVariance(m2, type = "h", show.labels = T, legend.position = "bottom"))
   })

   test_that("Cumulative variance (Y) plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotYCumVariance(m1))
      expect_silent(plotYCumVariance(m2))
      expect_silent(plotYCumVariance(m1, type = "h", show.labels = T, legend.position = "bottom"))
      expect_silent(plotYCumVariance(m2, type = "h", show.labels = T, legend.position = "bottom"))
   })

   test_that("X-scores plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXScores(m1))
      expect_silent(plotXScores(m2))
      expect_silent(plotXScores(m1, c(1, 3), show.labels = T, show.excluded = T, legend.position = "top"))
      expect_silent(plotXScores(m2, c(1, 3), show.labels = T, show.excluded = T, legend.position = "top"))
   })

   # X-distance
   test_name <- "X-residual distance plot works fine"
   test_that(test_name, {
      par(mfrow = c(2, 2))
      expect_silent(plotXResiduals(m1))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotXResiduals(m2))
      expect_silent(plotXResiduals(m1, 3, show.labels = T, show.excluded = T))
      expect_silent(plotXResiduals(m2, 3, show.labels = T, show.excluded = T))
   })

   test_name <- "X-residual distance plot understands two values for show.limits"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotXResiduals(m1, show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotXResiduals(m1, show.limits = c(FALSE, TRUE)))
      expect_silent(plotXResiduals(m2, show.limits = c(TRUE, FALSE)))
      expect_silent(plotXResiduals(m2, show.limits = c(TRUE, TRUE)))
   })

   test_name <- "X-residual distance plot: show.limits works fine with log transform"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotXResiduals(m1, log = TRUE, show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotXResiduals(m1, log = TRUE, show.limits = c(FALSE, TRUE)))
      expect_silent(plotXResiduals(m2, log = TRUE, show.limits = c(TRUE, FALSE)))
      expect_silent(plotXResiduals(m2, log = TRUE, show.limits = c(TRUE, TRUE)))
   })

   test_name <- "X-residual distance plot: show.limits works fine with none normalized values"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotXResiduals(m1, norm = FALSE, show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotXResiduals(m1, norm = FALSE, show.limits = c(FALSE, TRUE)))
      expect_silent(plotXResiduals(m2, norm = FALSE, show.limits = c(TRUE, FALSE)))
      expect_silent(plotXResiduals(m2, norm = FALSE, show.limits = c(TRUE, TRUE)))
   })

   # XY-distance
   test_name <- "XY-residual distance plot works fine"
   test_that(test_name, {
      par(mfrow = c(2, 2))
      expect_silent(plotXYResiduals(m1))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotXYResiduals(m2))
      expect_silent(plotXYResiduals(m1, 3, show.labels = T, show.excluded = T))
      expect_silent(plotXYResiduals(m2, 3, show.labels = T, show.excluded = T))
   })

   test_name <- "XY-residual distance plot understands two values for show.limits"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotXYResiduals(m1, show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotXYResiduals(m1, show.limits = c(FALSE, TRUE)))
      expect_silent(plotXYResiduals(m2, show.limits = c(TRUE, FALSE)))
      expect_silent(plotXYResiduals(m2, show.limits = c(TRUE, TRUE)))
   })

   test_name <- "XY-residual distance plot and two values for show.limits works fine with category"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotXYResiduals(m1, cgroup = "categories", show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotXYResiduals(m1, cgroup = "categories", show.limits = c(FALSE, TRUE)))
      expect_silent(plotXYResiduals(m2, cgroup = "categories", show.limits = c(TRUE, FALSE)))
      expect_silent(plotXYResiduals(m2, cgroup = "categories", show.limits = c(TRUE, TRUE)))
   })

   test_name <- "XY-residual distance plot: show.limits works fine with log transform"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotXYResiduals(m1, log = TRUE, show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotXYResiduals(m1, log = TRUE, show.limits = c(FALSE, TRUE)))
      expect_silent(plotXYResiduals(m2, log = TRUE, show.limits = c(TRUE, FALSE)))
      expect_silent(plotXYResiduals(m2, log = TRUE, show.limits = c(TRUE, TRUE)))
   })

   test_name <- "XY-residual distance plot: show.limits works fine with none normalized values"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotXYResiduals(m1, norm = FALSE, show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotXYResiduals(m1, norm = FALSE, show.limits = c(FALSE, TRUE)))
      expect_silent(plotXYResiduals(m2, norm = FALSE, show.limits = c(TRUE, FALSE)))
      expect_silent(plotXYResiduals(m2, norm = FALSE, show.limits = c(TRUE, TRUE)))
   })


   # XY-scores
   test_that("XY-scores plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXYScores(m1))
      expect_silent(plotXYScores(m2))
      expect_silent(plotXYScores(m1, 3, show.labels = T, show.excluded = T, legend.position = "top"))
      expect_silent(plotXYScores(m2, 3, show.labels = T, show.excluded = T, legend.position = "top"))
   })

   test_that("X-loadings plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXLoadings(m1))
      expect_silent(plotXLoadings(m2))
      expect_silent(plotXLoadings(m1, c(1, 3), type = "h", show.labels = T, legend.position = "top"))
      expect_silent(plotXLoadings(m2, c(1, 3), type = "h", show.labels = T, legend.position = "top"))
   })

   test_that("XY-loadings plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotXYLoadings(m1))
      expect_silent(plotXYLoadings(m2))
      expect_silent(plotXYLoadings(m1, c(1, 3), show.labels = T, legend.position = "top"))
      expect_silent(plotXYLoadings(m2, c(1, 3), show.labels = T, legend.position = "top"))
   })

   test_that("Weights plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotWeights(m1))
      expect_silent(plotWeights(m2))
      expect_silent(plotWeights(m1, c(1, 3), type = "p", show.labels = T))
      expect_silent(plotWeights(m2, c(1, 3), type = "p", show.labels = T))
   })

   test_that("VIP scores plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotVIPScores(m1))
      expect_silent(plotVIPScores(m2))
      expect_silent(plotVIPScores(m1, ny = 1, type = "l", show.labels = T))
      expect_silent(plotVIPScores(m2, ny = ncol(d$yc), ncomp = m2$ncomp, type = "h", show.labels = T))
   })

   test_that("Selectivity ratio plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotSelectivityRatio(m1))
      expect_silent(plotSelectivityRatio(m2))
      expect_silent(plotSelectivityRatio(m1, ny = 1, type = "l", show.labels = T))
      expect_silent(plotSelectivityRatio(m2, ny = ncol(d$yc), ncomp = m2$ncomp, type = "h", show.labels = T))
   })

   test_that("Regression coefficients  plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotRegcoeffs(m1))
      expect_silent(plotRegcoeffs(m2))
      expect_silent(plotRegcoeffs(m1, ny = ncol(d$yc), ncomp = 2, type = "h", show.labels = T))
      expect_silent(plotRegcoeffs(m2, ny = ncol(d$yc), ncomp = 2, type = "h", show.labels = T))
   })

   test_that("RMSE plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotRMSE(m1))
      expect_silent(plotRMSE(m2))
      expect_silent(plotRMSE(m1, ny = ncol(d$yc), type = "h", show.labels = T, legend.position = "top"))
      expect_silent(plotRMSE(m2, ny = ncol(d$yc), type = "h", show.labels = T, legend.position = "top"))
   })

   test_that("Y-residuals plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotYResiduals(m1))
      expect_silent(plotYResiduals(m2))
      expect_silent(plotYResiduals(m1, ny = ncol(d$yc), ncomp = 2, show.labels = T, show.excluded = T))
      expect_silent(plotYResiduals(m2, ny = ncol(d$yc), ncomp = 2, show.labels = T, show.excluded = T))
   })

   test_that("predictions plot works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotPredictions(m1))
      expect_silent(plotPredictions(m2))
      expect_silent(plotPredictions(m1, ny = ncol(d$yc), ncomp = 2, show.labels = T, show.excluded = T))
      expect_silent(plotPredictions(m2, ny = ncol(d$yc), ncomp = 2, show.labels = T, show.excluded = T))
   })


   test_that("overall model plot works fine", {
      expect_silent(plot(m1))
      expect_silent(plot(m2, ny = ncol(d$yc), ncomp = 2))
   })
}
