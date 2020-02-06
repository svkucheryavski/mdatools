#######################################1
# Tests for classmodel class methods  #
#######################################

setup({
   pdf(file = tempfile("mdatools-test-classmodel-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-test-classmodel-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# mock data

## reference values
c.ref.cal <- as.factor(c(rep("C2", 10), rep("C1", 10), rep("C3", 10)))
c.ref.test <- as.factor(c(rep("C2", 5), rep("C1", 5), rep("C3", 5)))

## excluded values
excluded_rows <- c(1, 7, 15, 25)

## predicted values
nrows.cal <- 30
nrows.test <- 15
ncomp <- 3
nclasses <- 2
n.cal <- nrows.cal * ncomp * nclasses
n.test <- nrows.test * ncomp * nclasses
classnames <- paste0("C", seq_len(nclasses))

set.seed(42)
c.pred.cal <- array((rnorm(n.cal) > 0) * 2 - 1, dim = c(nrows.cal, ncomp, nclasses))
dimnames(c.pred.cal)[[1]] <- paste0("O", seq_len(nrows.cal))
dimnames(c.pred.cal)[[2]] <- paste0("Comp ", seq_len(ncomp))
dimnames(c.pred.cal)[[3]] <- classnames
attr(c.pred.cal, "exclrows") <- excluded_rows
attr(c.ref.cal, "exclrows") <- excluded_rows

c.pred.cv <- array((rnorm(n.cal) > 0) * 2 - 1, dim = c(nrows.cal, ncomp, nclasses))
dimnames(c.pred.cv)[[1]] <- paste0("O", seq_len(nrows.cal))
dimnames(c.pred.cv)[[2]] <- paste0("Comp ", seq_len(ncomp))
dimnames(c.pred.cv)[[3]] <- classnames
attr(c.pred.cv, "exclrows") <- excluded_rows

c.pred.test <- array((rnorm(n.test) > 0) * 2 - 1, dim = c(nrows.test, ncomp, nclasses))
dimnames(c.pred.test)[[1]] <- paste0("T", seq_len(nrows.test))
dimnames(c.pred.test)[[2]] <- paste0("Comp ", seq_len(ncomp))
dimnames(c.pred.test)[[3]] <- classnames

m3 <- m2 <- m1 <- list("nclasses" = 2, "ncomp.selected" = 1, "classnames" = classnames)
class(m3) <- class(m2) <- class(m1) <- "classmodel"

m1$name <- "cal + cv + test"
m1$res <- list(
   "cal" = classres(c.pred.cal, c.ref.cal),
   "cv" = classres(c.pred.cv, c.ref.cal),
   "test" = classres(c.pred.test, c.ref.test)
)

m2$name <- "cal + cv"
m2$res <- list(
   "cal" = classres(c.pred.cal, c.ref.cal),
   "cv" = classres(c.pred.cv, c.ref.cal)
)

m3$name <- "cal"
m3$res <- list(
   "cal" = classres(c.pred.cal, c.ref.cal)
)

for (m in list(m1, m2, m3)) {

   context(sprintf("classmodel: plots for model with %s", m$name))

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("Classmodel - plots - model with ", m$name), pos = 4)

   test_that("performance plots works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotPerformance(m))
      expect_silent(plotPerformance(m, type = "b"))
      expect_silent(plotSensitivity(m))
      expect_silent(plotSpecificity(m))
   })

   test_that("prediction plots works fine", {
      par(mfrow = c(2, 2))
      expect_silent(plotPredictions(m))
      expect_silent(plotPredictions(m, col = c("red", "orange", "green"), show.labels = T))
      expect_silent(plotPredictions(m, res.name = "cal"))
      expect_silent(plotPredictions(m, res.name = "cal", show.excluded = T, show.labels = T))
   })
}

## test output sinked to external file
# TODO: implement print and summary for classmodel
print(m1)
print(m2)
print(m3)
print(summary(m1))
print(summary(m2))
print(summary(m3))

