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

context("pls: fixes of bugs for 0.10.3")

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
      m <- pls(X, y, cv = 8)
   })
})

#data(simdata)
#X <- simdata$spectra.c
#y <- simdata$conc.c
#pls(X, y, ncomp = 80, cv = 1)
#pls(X, y, ncomp = 80, cv = 2)

test_that("bug #86 is fixed", {
   data(simdata)
   X <- simdata$spectra.c
   y <- simdata$conc.c

   # too many components - should show warning
   expect_warning(pls(X, y, ncomp = 80))
   expect_warning(pls(X, y, ncomp = 80, cv = 1))
   expect_warning(pls(X, y, ncomp = 80, cv = 4))
   expect_warning(pls(X, y, ncomp = 80, cv = list("rand", 4, 4)))
   expect_warning(pls(X, y, ncomp = 80, cv = list("ven", 4)))

   # number of components limied - should be ok
   expect_silent(pls(X, y, ncomp = 50, cv = 1))
   expect_silent(pls(X, y, ncomp = 50, cv = 4))
   expect_silent(pls(X, y, ncomp = 50, cv = list("rand", 4, 4)))
   expect_silent(pls(X, y, ncomp = 50, cv = list("ven", 4)))
})
