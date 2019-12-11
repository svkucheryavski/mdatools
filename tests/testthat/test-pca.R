################################
# Tests for PCA class methods  #
################################

## create several datasets with different layout based on People data
data(people)

x1 <- people
x2 <- x1[1:2, , drop = F]
x3 <- x1[, 1:2, drop = F]

## same but scaled
x4 <- scale(people, center = TRUE, scale = TRUE)
x5 <- x4[1:2, , drop = F]
x6 <- x5[, 1:2, drop = F]

## large random data 1000 x 10
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

## combine all to a list
x <- list(x1, x2, x3, x4, x5, x6, x7)


########################################################
# Block 1: testing methods implementing PCA algorithms #
########################################################

context('PCA: testing different algorithms')

## testing function for testing results of one model
tf <- function(x, f) {

   # set of correct and incorrect number of components
   ncomp_correct <- c(1, min(ncol(x), nrow(x) - 1), round(min(ncol(x), nrow(x) - 1)/2) + 1)
   ncomp_wrong <- c(0, 9999)

   for (ncomp in ncomp_correct) {
      res <- f(x, ncomp)

      # check size of scores, loadings and eigenvalues
      expect_equal(dim(res$loadings), c(ncol(x), ncomp))
      expect_equal(dim(res$scores), c(nrow(x), ncomp))
      expect_equal(length(res$eigenvals), ncomp)

      # check computed values
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

# testing function for direct comparing two methods
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

######################################
# Block 2: testing pca.run() method  #
######################################

context('PCA: testing pca.run()')

## testing function for comparing pca.run() results with algorithms
tf <- function(x, f1, method, rand, tol = 10^-6) {

   ncomp_correct <- c(1, min(ncol(x), nrow(x) - 1), round(min(ncol(x), nrow(x) - 1)/2) + 1)

   for (ncomp in ncomp_correct) {
      res1 <- f1(x, ncomp)
      res2 <- pca.run(x, ncomp, method, rand)

      expect_equivalent(abs(res1$scores), abs(res2$scores), tolerance = tol)
      expect_equivalent(abs(res1$loadings), abs(res2$loadings), tolerance = tol)
      expect_equivalent(res1$eigenvals, res2$eigenvals, tolerance = tol)

      expect_equal(dim(res2$residuals), dim(x))
      expect_equivalent(res2$residuals, x - tcrossprod(res2$scores, res2$loadings))
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
      ## the bigger q or p in rand = c(p, q) the more precise results are
      tf(x.c, pca.svd, method = "svd", rand = c(5, 0), tol = 10^-2)
      tf(x.c, pca.svd, method = "svd", rand = c(5, 1), tol = 10^-3)
      tf(x.c, pca.svd, method = "svd", rand = c(5, 2), tol = 10^-4)
   }
})

test_that("pca.run works correctly with rand (nipals)", {
   for (x.c in x) {
      ## the bigger q or p in rand = c(p, q) the more precise results are
      tf(x.c, pca.nipals, method = "svd", rand = c(5, 0), tol = 10^-2)
      tf(x.c, pca.nipals, method = "svd", rand = c(5, 1), tol = 10^-3)
      tf(x.c, pca.nipals, method = "svd", rand = c(5, 2), tol = 10^-4)
   }
})

######################################
# Block 3: testing pca.cal() method  #
######################################

## setup two datasets - one without attributes and with
x1 <- people
x2 <- people
attr(x2, "yaxis.name") <- "Persons"
attr(x2, "yaxis.values") <- 1:32
attr(x2, "xaxis.name") <- "Measurements"
attr(x2, "xaxis.values") <- (1:12)/10

## number of components to test
ncomp_correct <- c(1, 5, 12)

tf <- function(x, ncomp, center, scale) {
   m <- pca.cal(x, ncomp, center = center, scale = scale, method = "svd")
   expect_equal(m$ncomp, ncomp)

   # excluded rows to compute scale and center correctly
   attrs <- attributes(x)
   if (length(attrs$exclrows) > 0) x <- x[-attrs$exclrows, , drop = FALSE]
   expect_equivalent(m$scale, (if (is.logical(scale) && scale == TRUE) apply(x, 2, sd) else scale))
   expect_equivalent(m$center, (if (is.logical(center) && center == TRUE) apply(x, 2, mean) else center))

   # excluded columns to compute other thinngs correctly
   xt <- prep.autoscale(x, center = center, scale = scale)
   if (length(attrs$exclcols) > 0) xt <- xt[, -attrs$exclcols, drop = F]
   mt <- pca.svd(xt, ncomp)

   expect_equivalent(m$eigenvals, mt$eigenvals)
   expect_equivalent(length(m$eigenvals), ncomp)
   expect_equivalent(names(m$eigenvals), colnames(m$loadings))

   expect_equivalent(dim(m$loadings), c(ncol(x), ncomp))
   expect_equivalent(rownames(m$loadings), colnames(x))
   expect_equivalent(attr(m$loadings, "yaxis.name"), attr(x, "xaxis.name"))
   expect_equivalent(attr(m$loadings, "yaxis.values"), attr(x, "xaxis.values"))

   if (length(attrs$exclcols) > 0) {
      expect_equivalent(m$loadings[-attrs$exclcols, , drop = FALSE], mt$loadings)
      expect_equivalent(m$loadings[attrs$exclcols, , drop = FALSE], matrix(0, length(attrs$exclcols), ncomp))
      expect_equivalent(attr(m$loadings, "exclrows"), attrs$exclcols)
   } else {
      expect_equivalent(m$loadings, mt$loadings)
   }

}

test_that("pca.cal works fine with full data", {

   for (x in list(x1, x2)) {
      for (ncomp in ncomp_correct) {
         tf(x, ncomp, center = FALSE, scale = FALSE)
         tf(x, ncomp, center = TRUE, scale = FALSE)
         tf(x, ncomp, center = TRUE, scale = TRUE)
         tf(x, ncomp, center = apply(x, 2, median), scale = TRUE)
         tf(x, ncomp, center = apply(x, 2, median), scale = apply(x, 2, IQR))
      }
   }
})

## exclude several columns and rows and redefine components
x1 <- mda.exclrows(x1, c(1, 10, 20))
x2 <- mda.exclcols(x1, c(1, 6, 12))
ncomp_correct <- c(1, 5, 9)

test_that("pca.cal works fine with excluded data", {

   for (x in list(x1, x2)) {
      for (ncomp in ncomp_correct) {
         tf(x, ncomp, center = FALSE, scale = FALSE)
         tf(x, ncomp, center = TRUE, scale = FALSE)
         tf(x, ncomp, center = TRUE, scale = TRUE)
         tf(x, ncomp, center = apply(x, 2, median), scale = TRUE)
         tf(x, ncomp, center = apply(x, 2, median), scale = apply(x, 2, IQR))
      }
   }
})
