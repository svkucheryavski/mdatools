#####################################################
# Tests for plotting methods of classres() class    #
#####################################################

setup({
   pdf(file = tempfile("mdatools-test-classresplots-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-classresplots-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# mock data

## reference values
c.ref <- as.factor(c(rep("C2", 10), rep("C1", 10), rep("C3", 10)))

## predicted values
nrows <- 30
ncomp <- 3
nclasses <- 2
n <- nrows * ncomp * nclasses
classnames <- paste0("C", seq_len(nclasses))

set.seed(42)
c.pred <- array((rnorm(n) > 0) * 2 - 1, dim = c(nrows, ncomp, nclasses))
dimnames(c.pred)[[1]] <- paste0("O", 1:nrows)
dimnames(c.pred)[[2]] <- paste0("Comp ", seq_len(ncomp))
dimnames(c.pred)[[3]] <- classnames

## excluded values
excluded_rows <- c(1, 7, 15, 25)
# 4. same as #2 but with excluded rows

## result with excluded values
c.predexcl <- c.pred
attr(c.predexcl, "exclrows") <- excluded_rows
res1 <- classres(c.predexcl, c.ref, ncomp.selected = 2)

## result with excluded values and one component
c.predexcl <- c.pred[, 1, , drop = FALSE]
attr(c.predexcl, "exclrows") <- excluded_rows
res2 <- classres(c.predexcl, c.ref, ncomp.selected = 1)

## result with excluded values and one class
c.predexcl <- c.pred[, , 1, drop = FALSE]
attr(c.predexcl, "exclrows") <- excluded_rows
res3 <- classres(c.predexcl, c.ref, ncomp.selected = 2)

## result with excluded values, one component and one class
c.predexcl <- c.pred[, 1, 1, drop = FALSE]
attr(c.predexcl, "exclrows") <- excluded_rows
res4 <- classres(c.predexcl, c.ref, ncomp.selected = 1)

res <- list(
   "exclrows" = res1,
   "exclrows + one comp" = res2,
   "exclrows + one class" = res3,
   "exclrows + one comp + one class" = res3
)

for (i in seq_along(res)) {

   nc <- res[[i]]$nclasses
   name <- names(res)[[i]]

   context(sprintf("classres: performance plots for %s", name))


   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("Classres - performance plots - ", name), pos = 4)

   # performance
   par(mfrow = c(2, 2))
   test_that("performance plot works correctly", {
      expect_silent(plotPerformance(res[[i]]))
      expect_silent(plotPerformance(res[[i]], nc = nc))
      expect_silent(plotPerformance(res[[i]], nc = nc, type = "b"))
      expect_silent(plotPerformance(res[[i]], nc = nc, type = "l", col = "red", lty = 1:3))
   })

   test_that("performance plot returns plot data if show.plot = FALSE", {
      pd <- plotPerformance(res[[i]], nc = nc, show.plot = F)
      expect_equal(dim(pd), c(3, res[[i]]$ncomp))
      expect_equivalent(pd,
         rbind(
            res[[i]]$sensitivity[nc,],
            res[[i]]$specificity[nc,],
            res[[i]]$misclassified[nc,]
         )
      )
   })

   # sensitivity
   par(mfrow = c(2, 2))
   test_that("sensitivity plot works correctly", {
      expect_silent(plotSensitivity(res[[i]]))
      expect_silent(plotSensitivity(res[[i]], nc = nc))
      expect_silent(plotSensitivity(res[[i]], nc = nc, type = "b"))
      expect_silent(plotSensitivity(res[[i]], nc = nc, type = "l", col = "red", lty = 2))
   })

   test_that("sensitivity plot returns plot data if show.plot = FALSE", {
      pd <- plotSensitivity(res[[i]], nc = nc, show.plot = F)
      expect_equal(dim(pd), c(1, res[[i]]$ncomp))
      expect_equivalent(pd, rbind(res[[i]]$sensitivity[nc,]))
   })

   # specificity
   par(mfrow = c(2, 2))
   test_that("specificity plot works correctly", {
      expect_silent(plotSpecificity(res[[i]]))
      expect_silent(plotSpecificity(res[[i]], nc = nc))
      expect_silent(plotSpecificity(res[[i]], nc = nc, type = "b"))
      expect_silent(plotSpecificity(res[[i]], nc = nc, type = "l", col = "red", lty = 2))
   })

   test_that("specificity plot returns plot data if show.plot = FALSE", {
      pd <- plotSpecificity(res[[i]], nc = nc, show.plot = F)
      expect_equal(dim(pd), c(1, res[[i]]$ncomp))
      expect_equivalent(pd, rbind(res[[i]]$specificity[nc,]))
   })

   # misclassified
   par(mfrow = c(2, 2))
   test_that("misclassified plot works correctly", {
      expect_silent(plotMisclassified(res[[i]]))
      expect_silent(plotMisclassified(res[[i]], nc = nc))
      expect_silent(plotMisclassified(res[[i]], nc = nc, type = "b"))
      expect_silent(plotMisclassified(res[[i]], nc = nc, type = "l", col = "red", lty = 2))
   })

   test_that("misclassified plot returns plot data if show.plot = FALSE", {
      pd <- plotMisclassified(res[[i]], nc = nc, show.plot = F)
      expect_equal(dim(pd), c(1, res[[i]]$ncomp))
      expect_equivalent(pd, rbind(res[[i]]$misclassified[nc,]))
   })



   context(sprintf("classres: prediction plots for %s", name))

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("Classres - prediction plots - ", name), pos = 4)

   # classification results
   par(mfrow = c(2, 2))
   test_that("prediction plot works correctly", {
      expect_silent(plotPredictions(res[[i]]))
      expect_silent(plotPredictions(res[[i]], nc = 1))
      expect_silent(plotPredictions(res[[i]], nc = 1, ncomp = 1))
      expect_silent(plotPredictions(res[[i]], col = c("red", "orange", "brown"), pch = 1:3))
   })

   # classification results - labels
   par(mfrow = c(2, 2))
   test_that("prediction plot works correctly", {
      expect_silent(plotPredictions(res[[i]]))
      expect_silent(plotPredictions(res[[i]], nc = 1, show.excluded = T, show.labels = T, labels = "values"))
      expect_silent(plotPredictions(res[[i]], show.labels = T, labels = "indices"))
      expect_silent(plotPredictions(res[[i]], show.labels = T, labels = "names"))
   })

   test_that("prediction plot returns plot data if show.plot = FALSE", {
      pd <- plotPredictions(res[[i]], show.plot = F)
      expect_equal(ncol(pd), 2)
      expect_equal(min(pd[, 2]), 1)
      expect_equal(max(pd[, 2]), nc + 1)
   })

   # just output to check in txt file
   fprintf("\nResult object: %s\n", names(res)[[i]])
   cat("-------------------------------\n")
   print(res[[i]])
   summary(res[[i]])
   showPredictions(res[[i]])
   print(getConfusionMatrix(res[[i]], 1))

}
