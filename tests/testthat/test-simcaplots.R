###############################################
# Tests for simca and simcares class methods  #
###############################################

setup({
   pdf(file = tempfile("mdatools-test-simcaplots-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-simcaplots-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

## prepare datasets
data(iris)
ind.test <- seq(2, nrow(iris), by = 2)

x.cal <- iris[-ind.test, 1:4]
x.cal1 <- x.cal[1:25, ]
x.cal1 <- mda.exclrows(x.cal1, c(1, 10))
x.cal2 <- x.cal[26:50, ]
x.cal2 <- mda.exclrows(x.cal2, c(11, 21))
x.cal3 <- x.cal[51:75, ]
x.cal3 <- mda.exclrows(x.cal3, c(6, 16))

classnames <- levels(iris[, 5])
x.test <- iris[ind.test, 1:4]
c.test <- iris[ind.test, 5]

## create models
m1 <- simca(x.cal1, classnames[1], 4, cv = 1, scale = F)
m1 <- selectCompNum(m1, 2)
m2 <- simca(x.cal2, classnames[2], 4, cv = 5, scale = T, lim.type = "chisq", alpha = 0.01)
m2 <- selectCompNum(m2, 3)
m3 <- simca(x.cal3, classnames[3], 4, cv = list("rand", 5, 10), x.test = x.test,
   c.test = c.test, scale = T, lim.type = "ddrobust", alpha = 0.10)
m3 <- selectCompNum(m3, 3)
m3 <- setDistanceLimits(m3, lim.type = "ddrobust", alpha = 0.05)


models <- list("se" = m1, "ve" = m2, "vi" = m3)
for (i in seq_along(models)) {

   m <- models[[i]]
   name <- names(models)[i]
   calres <- list("cal" = m$res[["cal"]])

   context(sprintf("simca: test pcamodel related part (model %s)", name))

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("SIMCA - pca plots - ", name), pos = 4)

   test_that("loadings plot works well", {
      par(mfrow = c(2, 2))
      expect_silent(plotLoadings(m))
      expect_silent(plotLoadings(m, comp = c(1, 2), type = "h", show.labels = T))
      expect_silent(plotLoadings(m, type = "l", show.labels = T, labels = "values"))
      expect_silent(plotLoadings(m, type = "b"))
   })

   test_that("variance plot works well", {
      par(mfrow = c(2, 2))
      expect_silent(plotVariance(m))
      expect_silent(plotVariance(m, type = "h", show.labels = T))
      expect_silent(plotCumVariance(m))
      expect_silent(plotCumVariance(m, type = "h", show.labels = T))
   })

   if (m$lim.type != "jm") {
      test_that("DoF plot works well", {
         par(mfrow = c(2, 2))
         expect_silent(plotQDoF(m))
         expect_silent(plotT2DoF(m, type = "l", show.labels = T))
         expect_silent(plotDistDoF(m))
         expect_silent(plotDistDoF(m, type = "l", show.labels = T))
      })
   }

   test_that("scores plot works well", {
      par(mfrow = c(2, 2))
      expect_silent(plotScores(m))
      expect_silent(plotScores(m, comp = c(1, 3), show.labels = T))
      expect_silent(plotScores(m, type = "h", show.labels = T, labels = "values", res = calres))
      expect_silent(plotScores(m, type = "b", res = calres))
   })

   test_that("extreme plot works well", {
      par(mfrow = c(2, 2))
      expect_silent(plotExtreme(m))
      expect_silent(plotExtreme(m, comp = 2, show.labels = T))
      expect_silent(plotExtreme(m, comp = 1:3))
      expect_silent(plotExtreme(m, comp = 1:2, col = c("red", "green")))
   })

   test_that("residuals plot works well", {
      res <- list("cal" = m$res[["cal"]])
      par(mfrow = c(2, 2))
      expect_silent(plotResiduals(m))
      expect_silent(plotResiduals(m, ncomp = 3))
      if (m$lim.type != "chisq") {
         expect_silent(plotResiduals(m, ncomp = 3, cgroup = "categories", res = calres, log = T))
      } else {
         expect_silent(plotResiduals(m, ncomp = 3, cgroup = "categories", res = calres))
      }
      expect_silent(plotResiduals(m, ncomp = 2, show.labels = T))
   })

   context(sprintf("simca: test classmodel related part (model %s)", name))

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("SIMCA - classification plots - ", name), pos = 4)

   par(mfrow = c(2, 2))
   test_that("performance plot works correctly", {
      expect_silent(plotPerformance(m))
      expect_silent(plotPerformance(m, type = "h"))
      expect_error(plotSpecificity(m))
      expect_silent(plotSensitivity(m))
   })

   # classification results
   par(mfrow = c(2, 2))
   test_that("predictions and probabilities plot works correctly", {
      expect_silent(plotPredictions(m))
      expect_silent(plotPredictions(m, ncomp = 1, pch = 1))
      expect_silent(plotProbabilities(m$calres))
      expect_silent(plotProbabilities(m$calres, ncomp = 1, pch = 1))
   })

   # just output to check in txt file
   fprintf("\nSummary and print methods for model: %s\n", name)
   cat("-------------------------------\n")
   print(m)
   summary(m)

   context(sprintf("simca: new prdictions (model %s)", name))

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("SIMCA - results for predictions - ", name), pos = 4)

   res <- predict(m, x.test, c.test)

    test_that("variance plot works well", {
      par(mfrow = c(2, 2))
      expect_silent(plotVariance(res))
      expect_silent(plotVariance(res, type = "h", show.labels = T))
      expect_silent(plotCumVariance(res))
      expect_silent(plotCumVariance(res, type = "h", show.labels = T))
   })

   test_that("scores plot works well", {
      par(mfrow = c(2, 2))
      expect_silent(plotScores(res))
      expect_silent(plotScores(res, comp = c(1, 3), show.labels = T))
      expect_silent(plotScores(res, type = "h", show.labels = T, labels = "values"))
      expect_silent(plotScores(res, type = "b", res = calres))
   })

   test_that("residuals plot works well", {
      par(mfrow = c(2, 2))
      expect_silent(plotResiduals(res))
      expect_silent(plotResiduals(res, ncomp = 3))
      expect_silent(plotResiduals(m, cgroup = "categories", res = list("new" = res)))
      expect_silent(plotResiduals(m, cgroup = c.test, res = list("new" = res)))
   })

   par(mfrow = c(2, 2))
   test_that("performance plot works correctly", {
      expect_silent(plotPerformance(res))
      expect_silent(plotPerformance(res, type = "h"))
      expect_silent(plotSpecificity(res))
      expect_silent(plotSensitivity(res))
   })

   # classification results
   par(mfrow = c(2, 2))
   test_that("prediction plot works correctly", {
      expect_silent(plotPredictions(res))
      expect_silent(plotPredictions(res, ncomp = 1, pch = 1))
      expect_silent(plotProbabilities(res))
      expect_silent(plotProbabilities(res, ncomp = 1, pch = 1))
   })

   # just output to check in txt file
   fprintf("\nSummary and print methods for result: %s\n", name)
   cat("-------------------------------\n")
   print(res)
   summary(res)
   showPredictions(res)
}

# code for manual comparison  with DDSIMCA GUI - not uses as real tests

if (FALSE) {

   ## prepare datasets
   rm(list = ls())
   data(iris)
   ind.test <- seq(2, nrow(iris), by = 2)

   x.cal <- iris[-ind.test, 1:4]
   x.cal1 <- x.cal[1:25, ]
   x.cal2 <- x.cal[26:50, ]
   x.cal3 <- x.cal[51:75, ]

   classnames <- levels(iris[, 5])
   x.test <- iris[ind.test, 1:4]
   c.test <- iris[ind.test, 5]

   ## test for ddmoments
   m <- simca(x.cal1, classnames[1], 3, scale = F, lim.type = "ddmoments")
   plotResiduals(m, show.labels = T, labels = "indices", cgroup = "categories", log = T)
   plotExtreme(m)
   summary.pca(m)
   summary(m)

   r <- predict(m, x.test, c.test)
   plotResiduals(m, res = list(new = r), show.labels = T, labels = "indices", cgroup = c.test, log = T)
   plotExtreme(m, res = r)
   summary(r)

   ## test for ddrobust
   m <- simca(x.cal1, classnames[1], 3, scale = F, lim.type = "ddrobust")
   plotResiduals(m, show.labels = T, labels = "indices", cgroup = "categories", log = T)
   plotExtreme(m)
   summary.pca(m)
   summary(m)

   r <- predict(m, x.test, c.test)
   plotResiduals(m, res = list(new = r), show.labels = T, labels = "indices", cgroup = c.test, log = T)
   plotExtreme(m, res = r)
   summary(r)

   ## test for chisq
   m <- simca(x.cal1, classnames[1], 3, scale = F, lim.type = "chisq")
   plotResiduals(m, show.labels = T, labels = "indices", cgroup = "categories", log = T)
   plotExtreme(m)
   summary.pca(m)
   summary(m)

   r <- predict(m, x.test, c.test)
   plotResiduals(m, res = list(new = r), show.labels = T, labels = "indices", cgroup = c.test, log = T)
   plotExtreme(m, res = r)
   summary(r)

}
