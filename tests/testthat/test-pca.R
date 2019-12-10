data(people)

# create several datasets with different layout
x1 <- people
x2 <- x1[1:2, , drop = F]
x3 <- x1[, 1:2, drop = F]

# same but scaled
x4 <- scale(people, center = TRUE, scale = TRUE)
x5 <- x4[1:2, , drop = F]
x6 <- x5[, 1:2, drop = F]

# ! large random data 1000 x 10 (make bigger later)
x7 <- cbind(
   rnorm(1000, 0, 50),
   rnorm(1000, 0, 10),
   rnorm(1000, 0, 7),
   rnorm(1000, 0, 5),
   rnorm(1000, 0, 3),
   rnorm(1000, 0, 2),
   rnorm(1000, 0, 1),
   rnorm(1000, 0, 0.1),
   rnorm(1000, 0, 0.01),
   rnorm(1000, 0, 0.001)
)

# combine all to a list
x <- list(x1, x2, x3, x4, x5, x6, x7)

context('PCA: testing different algorithms')

# testing function for results of one model
tf <- function(x, f) {

   ncomp_correct <- c(1, min(ncol(x), nrow(x) - 1), round(min(ncol(x), nrow(x) - 1)/2) + 1)
   ncomp_wrong <- c(0, 9999)

   for (ncomp in ncomp_correct) {
      res <- f(x, ncomp)
      residuals <- x - res$scores %*% t(res$loadings)

      expect_equal(nrow(res$loadings), ncol(x))
      expect_equal(ncol(res$loadings), ncomp)
      expect_equal(nrow(res$scores), nrow(x))
      expect_equal(ncol(res$scores), ncomp)
      expect_equal(length(res$eigenvals), ncomp)

      expect_equivalent(res$scores, x %*% res$loadings)
      expect_equivalent(res$eigenvals, colSums(res$scores^2)/(nrow(x) - 1))
   }

   for (ncomp in ncomp_wrong) {
      expect_error(f(x, ncomp))
   }
}

test_that("pca.svd works correctly", {
   for (x.c in x) {
      tf(x.c, pca.svd)
   }
})

test_that("pca.nipals works correctly", {
   for (x.c in x) {
      tf(x.c, pca.nipals)
   }
})

# testing function for results of one model
tf <- function(x, f1, f2) {

   ncomp_correct <- c(1, min(ncol(x), nrow(x) - 1), round(min(ncol(x), nrow(x) - 1)/2) + 1)

   for (ncomp in ncomp_correct) {
      res1 <- f1(x, ncomp)
      res2 <- f2(x, ncomp)

      expect_equivalent(abs(res1$scores), abs(res2$scores), tolerance = 10^-3)
      expect_equivalent(abs(res1$loadings), abs(res2$loadings), tolerance = 10^-3)
      expect_equivalent(res1$eigenvals, res2$eigenvals, tolerance = 10^-3)
   }
}

test_that("SVD and NIPALS gives similar results works correctly", {
   for (x.c in x) {
      tf(x.c, pca.nipals, pca.svd)
   }
})

context('PCA: testing pca.run()')

# testing function for results of one model
tf <- function(x, f1, method, rand, tol = 10^-6) {

   ncomp_correct <- c(1, min(ncol(x), nrow(x) - 1), round(min(ncol(x), nrow(x) - 1)/2) + 1)

   for (ncomp in ncomp_correct) {
      res1 <- f1(x, ncomp)
      res2 <- pca.run(x, ncomp, method, rand)

      expect_equivalent(abs(res1$scores), abs(res2$scores), tolerance = tol)
      expect_equivalent(abs(res1$loadings), abs(res2$loadings), tolerance = tol)
      expect_equivalent(res1$eigenvals, res2$eigenvals, tolerance = tol)
   }
}

test_that("pca.run works correctly without rand", {
   for (x.c in x) {
      tf(x.c, pca.svd, method = "svd", rand = NULL)
      tf(x.c, pca.nipals, method = "nipals", rand = NULL)
   }
})

test_that("pca.run works correctly with rand (svd)", {
   for (x.c in x[1:7]) {
      tf(x.c, pca.svd, method = "svd", rand = c(5, 0), tol = 10^-2)
      tf(x.c, pca.svd, method = "svd", rand = c(5, 1), tol = 10^-3)
      tf(x.c, pca.svd, method = "svd", rand = c(5, 2), tol = 10^-4)
   }
})

test_that("pca.run works correctly with rand (nipals)", {
   for (x.c in x) {
      tf(x.c, pca.nipals, method = "svd", rand = c(5, 0), tol = 10^-2)
      tf(x.c, pca.nipals, method = "svd", rand = c(5, 1), tol = 10^-3)
      tf(x.c, pca.nipals, method = "svd", rand = c(5, 2), tol = 10^-4)
   }
})

