#####################################################
# Tests for Procrustes Cross-Validation methods     #
#####################################################

setup({
   pdf(file = tempfile("mdatools-test-pcv-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-pcv-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# Simple tests for the PCV R implementation
I <- 100
J <- 50
A <- 10
K <- 4
set.seed(42)
x <- matrix(rnorm(I * J), I, J)

params <- list()
params[[1]] <- list(x = x)
params[[2]] <- list(x = x, ncomp = 1)
params[[3]] <- list(x = x, ncomp = A)
params[[4]] <- list(x = x, ncomp = A, nseg = 4)
params[[5]] <- list(x = x, ncomp = A, nseg = 10)
params[[6]] <- list(x = x, ncomp = A, nseg = nrow(X))
params[[7]] <- list(x = x, ncomp = A, nseg = 4, scale = TRUE)
params[[8]] <- list(x = x, ncomp = A, nseg = 10, scale = TRUE)
params[[9]] <- list(x = x, ncomp = A, nseg = nrow(X), scale = TRUE)

context("pcv: for PCA")
for (i in seq_along(params)) {
   test_that(paste0("pcv() for PCA works well with parameters set #", i), {
      x.pv <- do.call(pcv, params[[i]])
      expect_equal(nrow(x.pv), nrow(x))
      expect_equal(ncol(x.pv), ncol(x))
      expect_gt(ks.test(x, x.pv)$p.value, 0.01)
   })
}
