####################################
# Tests for all bugs found in 2022 #
####################################

setup({
   pdf(file = tempfile("mdatools-test-classmodel-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-test-classmodel-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

context("test for bug #107:")
data(people)

test_that("bug is fixed in PCA", {
   m <- pca(people, 1)
   expect_silent(plotScores(m))
   expect_silent(plotLoadings(m))
   expect_silent(plotScores(m, 1))
   expect_silent(plotLoadings(m, 1))

   expect_error(plotScores(m, 0))
   expect_error(plotScores(m, c(1, 3)))
   expect_error(plotScores(m$res$cal, 0))
   expect_error(plotScores(m$res$cal, c(1, 3)))

   expect_error(plotLoadings(m, 0))
   expect_error(plotLoadings(m, c(1, 3)))
   expect_error(plotLoadings(m$res$cal, 0))
   expect_error(plotLoadings(m$res$cal, c(1, 3)))

})

test_that("bug is fixed in PLS", {
   m <- pls(people[, -4], people[, 4], 1)
   expect_silent(plotXScores(m))
   expect_silent(plotXLoadings(m))

   expect_error(plotXScores(m, 0))
   expect_error(plotXScores(m, c(1, 3)))
   expect_error(plotXScores(m$res$cal, 0))
   expect_error(plotXScores(m$res$cal, c(1, 3)))

   expect_error(plotXLoadings(m, 0))
   expect_error(plotXLoadings(m, c(1, 3)))
   expect_error(plotXLoadings(m$res$cal, 0))
   expect_error(plotXLoadings(m$res$cal, c(1, 3)))

})

context("test for bug in getRegcoeffs():")

test_that("the bug is fixed in PLS1", {

   X <- people[, -4]
   y <- people[, 4, drop = FALSE]

   # without centering and scaling
   sx <- rep(1, ncol(X))
   mx <- rep(0, ncol(X))
   sy <- 1
   my <- 0
   s <- sy / sx

   m <- pls(X, y, 10, cv = 1, scale = FALSE, center = FALSE)
   b1 <- getRegcoeffs(m)
   b2 <- m$coeffs$values[, m$ncomp.selected, 1]
   expect_equivalent(b1,  c(my - sum(s * b2 * mx), s * b2))

   # without centering
   sx <- apply(X, 2, sd)
   mx <- rep(0, ncol(X))
   sy <- as.numeric(apply(y, 2, sd))
   my <- 0
   s <- sy / sx

   m <- pls(X, y, 10, cv = 1, scale = TRUE, center = FALSE)
   b1 <- getRegcoeffs(m)
   b2 <- m$coeffs$values[, m$ncomp.selected, 1]
   expect_equivalent(b1,  c(my - sum(s * b2 * mx), s * b2))


   # without scaling
   sx <- rep(1, ncol(X))
   mx <- apply(X, 2, mean)
   sy <- 1
   my <- as.numeric(apply(y, 2, mean))
   s <- sy / sx

   m <- pls(X, y, 10, cv = 1, scale = FALSE)
   b1 <- getRegcoeffs(m)
   b2 <- m$coeffs$values[, m$ncomp.selected, 1]
   expect_equivalent(b1,  c(my - sum(s * b2 * mx), s * b2))

   # with centering and scaling
   sx <- apply(X, 2, sd)
   mx <- apply(X, 2, mean)
   sy <- as.numeric(apply(y, 2, sd))
   my <- as.numeric(apply(y, 2, mean))
   s <- sy / sx

   m <- pls(X, y, 10, cv = 1, scale = TRUE)
   b1 <- getRegcoeffs(m)
   b2 <- m$coeffs$values[, m$ncomp.selected, 1]
   expect_equivalent(b1,  c(my - sum(s * b2 * mx), s * b2))
})

test_that("the bug is fixed in PLS2", {

   testit <- function(m, sx, mx, sy, my) {

      # for ny = 1
      s <- sy[1] / sx
      b1 <- getRegcoeffs(m, ny = 1)
      b2 <- m$coeffs$values[, m$ncomp.selected, 1]
      expect_equivalent(b1,  c(my[1] - sum(s * b2 * mx), s * b2))

      # for ny = 2
      s <- sy[2] / sx
      b1 <- getRegcoeffs(m, ny = 2)
      b2 <- m$coeffs$values[, m$ncomp.selected, 2]
      expect_equivalent(b1,  c(my[2] - sum(s * b2 * mx), s * b2))
   }

   X <- people[, -c(4, 6)]
   y <- people[, c(4, 6), drop = FALSE]

   # without centering and scaling
   sx <- rep(1, ncol(X))
   mx <- rep(0, ncol(X))
   sy <- rep(1, ncol(y))
   my <- rep(0, ncol(y))
   m <- pls(X, y, 10, cv = 1, scale = FALSE, center = FALSE)
   testit(m, sx, mx, sy, my)

   # without centering
   sx <- apply(X, 2, sd)
   mx <- rep(0, ncol(X))
   sy <- as.numeric(apply(y, 2, sd))
   my <- rep(0, ncol(y))
   m <- pls(X, y, 10, cv = 1, scale = TRUE, center = FALSE)
   testit(m, sx, mx, sy, my)

   # without scaling
   sx <- rep(1, ncol(X))
   mx <- apply(X, 2, mean)
   sy <- rep(1, ncol(y))
   my <- as.numeric(apply(y, 2, mean))
   m <- pls(X, y, 10, cv = 1, scale = FALSE, center = TRUE)
   testit(m, sx, mx, sy, my)

   # with centering and scaling
   sx <- apply(X, 2, sd)
   mx <- apply(X, 2, mean)
   sy <- as.numeric(apply(y, 2, sd))
   my <- as.numeric(apply(y, 2, mean))
   m <- pls(X, y, 10, cv = 1, scale = TRUE, center = TRUE)
   testit(m, sx, mx, sy, my)

})
