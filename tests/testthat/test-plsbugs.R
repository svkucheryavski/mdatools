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

context("fixes of PLS bugs for 0.10.3")

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

