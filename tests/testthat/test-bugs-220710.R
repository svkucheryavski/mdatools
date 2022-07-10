
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