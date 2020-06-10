######################################
# Tests for bugs discovered for PLS  #
######################################

setup({
   pdf(file = tempfile("mdatools-test-pls-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-pls-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

context("pls: fixes of bugs")

## 9.06.2020
## bug #88 related to cross-validation in PLS


data(people)
people <- as.data.frame(people)
X <- people[, -4]
y <- people[,  4, drop = FALSE]
test_that("bug 88 is fixed", {
   expect_error(expect_warning({ set.seed(4); m <- pls(X, y, cv = 8)}))
   expect_silent({
      for (i in 1:10) {
         set.seed(i)
         for (j in 1:10) m <- pls(X, y, ncomp = 10, cv = 8)
      }
   })
})

## 26.03.2020
## bug #85 related to using y as data frame

test_that("bug #85 is fixed", {
   expect_silent({
      data(people)
      people <- as.data.frame(people)
      X <- people[, -4]
      y <- people[,  4, drop = FALSE]
      m <- pls(X, y, ncomp = 5, cv = 1)
      m <- pls(X, y, cv = 1)
      m <- pls(X, y, ncomp = 10, cv = 8)
   })
})


## bug #86 related to small eigenvalues problem in SIMPLS

test_that("bug #86 is fixed", {
   data(simdata)
   X <- simdata$spectra.c
   y <- simdata$conc.c

   # no cross-validation - just warning
   expect_warning(pls(X, y, ncomp = 80))

   # too many components - should show both warning and error (because of cross-validation)
   expect_error(expect_watning(pls(X, y, ncomp = 80, cv = 1)))
   expect_error(expect_watning(pls(X, y, ncomp = 80, cv = 4)))
   expect_error(expect_watning(pls(X, y, ncomp = 80, cv = list("rand", 4, 4))))
   expect_error(expect_watning(pls(X, y, ncomp = 80, cv = list("ven", 4))))

   # number of components limied - should be ok
   expect_silent(pls(X, y, ncomp = 30, cv = 1))
   expect_silent(pls(X, y, ncomp = 30, cv = 4))
   expect_silent(pls(X, y, ncomp = 30, cv = list("rand", 4, 4)))
   expect_silent(pls(X, y, ncomp = 30, cv = list("ven", 4)))
})
