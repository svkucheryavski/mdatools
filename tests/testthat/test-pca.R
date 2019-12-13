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

context('PCA: testing pca.cal()')


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
   expect_equivalent(m$scale,
      (if (is.logical(scale) && scale == TRUE) apply(x, 2, sd) else scale))
   expect_equivalent(m$center,
      (if (is.logical(center) && center == TRUE) apply(x, 2, mean) else center))

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
      expect_equivalent(m$loadings[attrs$exclcols, , drop = FALSE],
      matrix(0, length(attrs$exclcols), ncomp))
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

test_that("pca.cal calculates limit parameters correctly (full data)", {
   for (x in list(x1, x2)) {
      m <- pca.cal(x, 11, center = TRUE, scale = TRUE, method = "svd")

      # parameters computed using DDSimca Toolbox
      # using people data, autoscaled

      momentsParamsQ <- list(
         u0 = c(5.39623603, 3.22376502, 1.65661907, 0.68981822, 0.38111636, 0.22105099, 0.12476396,
            0.07148856, 0.04489746, 0.02151256, 0.00678829),
         Nu = c(10.00000000, 4.00000000, 4.00000000, 3.00000000, 3.00000000, 3.00000000, 5.0000000,
            4.00000000, 2.00000000, 2.00000000, 1.00000000),
         nobj = nrow(x)
      )

      momentsParamsT2 <- list(
         u0 = c(0.96875000, 1.93750000, 2.90625000, 3.87500000, 4.84375000, 5.81250000, 6.78125000,
            7.75000000, 8.71875000, 9.68750000, 10.65625000),
         Nu = c(2.00000000, 10.00000000, 13.00000000, 14.00000000, 13.00000000, 12.00000000, 8.000,
            9.00000000, 10.00000000, 12.00000000, 15.00000000),
         nobj = nrow(x)
      )

      robustParamsQ <- list(
         u0 = c(5.59398247, 3.14412554, 1.68373781, 0.63464055, 0.33473045, 0.21529782, 0.13614951,
            0.06862920, 0.04666029, 0.02054372, 0.00374660),
         Nu = c(18.00000000, 10.00000000, 6.00000000, 3.00000000, 4.00000000, 3.00000000, 5.0000000,
            3.00000000, 4.00000000, 2.00000000, 1.00000000),
         nobj = nrow(x)
      )

      robustParamsT2 <- list(
         u0 = c(1.23877706, 1.84445603, 2.67812507, 3.91673370, 4.51759976, 5.66332668, 5.96116787,
            7.19886896, 8.06726598, 8.95166844, 10.01296469),
         Nu = c(3.00000000, 22.00000000, 17.00000000, 11.00000000, 15.00000000, 18.000000, 16.00000,
            18.00000000, 18.00000000, 10.00000000, 14.00000000),
         nobj = nrow(x)
      )

      # check parameters for Q residuals (moments)
      expect_equal(m$limParams$Q$moments$u0, momentsParamsQ$u0, tolerance = 10^-5)
      expect_equal(round(m$limParams$Q$moments$Nu), momentsParamsQ$Nu)
      expect_equal(m$limParams$Q$moments$nobj, momentsParamsQ$nobj)

      # check parameters for Q residuals (robust)
      expect_equal(m$limParams$Q$robust$u0, robustParamsQ$u0, tolerance = 10^-5)
      expect_equal(round(m$limParams$Q$robust$Nu), robustParamsQ$Nu)
      expect_equal(m$limParams$Q$robust$nobj, robustParamsQ$nobj)

      # check parameters for T2 residuals (moments)
      expect_equal(m$limParams$T2$moments$u0, momentsParamsT2$u0, tolerance = 10^-5)
      expect_equal(round(m$limParams$T2$moments$Nu), momentsParamsT2$Nu)
      expect_equal(m$limParams$T2$moments$nobj, momentsParamsT2$nobj)

      # check parameters for T2 residuals (robust)
      expect_equal(m$limParams$T2$robust$u0, robustParamsT2$u0, tolerance = 10^-5)
      expect_equal(round(m$limParams$T2$robust$Nu), robustParamsT2$Nu)
      expect_equal(m$limParams$T2$robust$nobj, robustParamsT2$nobj)

   }
})

## exclude several columns and rows and redefine components
x1 <- mda.exclrows(x1, c(1, 10, 20))
x2 <- mda.exclcols(x1, c(3, 12))
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

test_that("pca.cal calculates limit parameters correctly (excluded data)", {
   for (x in list(x2)) {
      m <- pca.cal(x, 9, center = TRUE, scale = TRUE, method = "svd")

      # parameters computed using DDSimca Toolbox
      # using people data, autoscaled

      momentsParamsQ <- list(
         u0 = c(4.43617757, 2.19049660, 0.44027710, 0.25151771, 0.14625433, 0.08237345, 0.05117557,
            0.02675517, 0.00723308),
         Nu = c(7.00000000, 3.00000000, 5.00000000, 3.00000000, 4.00000000, 3.00000000, 3.00000000,
            2.00000000, 1.00000000),
         nobj = nrow(x) - length(attr(x, "exclrows"))
      )

      momentsParamsT2 <- list(
         u0 = c(0.96551724, 1.93103448, 2.89655172, 3.86206897, 4.82758621, 5.79310345, 6.75862069,
            7.72413793, 8.68965517),
         Nu = c(3.00000000, 13.00000000, 15.00000000, 15.00000000, 14.00000000, 14.000000, 20.00000,
            14.00000000, 17.00000000),
         nobj = nrow(x) - length(attr(x, "exclrows"))
      )

      robustParamsQ <- list(
         u0 = c(4.62537076, 2.40195599, 0.43916799, 0.23975306, 0.13645627, 0.07560735, 0.04351790,
            0.02600730, 0.00787551),
         Nu = c(16.00000000, 3.00000000, 4.00000000, 3.00000000, 5.00000000, 4.00000000, 2.0000000,
            2.00000000, 1.00000000),
         nobj = nrow(x) - length(attr(x, "exclrows"))
      )

      robustParamsT2 <- list(
         u0 = c(1.39125723, 2.06314350, 2.61209799, 3.71494609, 5.03000803, 5.67830732, 6.11739113,
            7.31069105, 8.12330118),
         Nu = c(3.00000000, 19.00000000, 25.00000000, 12.00000000, 9.00000000, 12.0000, 37.0000000,
            21.00000000, 15.00000000),
         nobj = nrow(x) - length(attr(x, "exclrows"))
      )

      # check parameters for Q residuals (moments)
      expect_equal(m$limParams$Q$moments$u0, momentsParamsQ$u0, tolerance = 10^-5)
      expect_equal(round(m$limParams$Q$moments$Nu), momentsParamsQ$Nu)
      expect_equal(m$limParams$Q$moments$nobj, momentsParamsQ$nobj)

      # check parameters for Q residuals (robust)
      expect_equal(m$limParams$Q$robust$u0, robustParamsQ$u0, tolerance = 10^-5)
      expect_equal(round(m$limParams$Q$robust$Nu), robustParamsQ$Nu)
      expect_equal(m$limParams$Q$robust$nobj, robustParamsQ$nobj)

      # check parameters for T2 residuals (moments)
      expect_equal(m$limParams$T2$moments$u0, momentsParamsT2$u0, tolerance = 10^-5)
      expect_equal(round(m$limParams$T2$moments$Nu), momentsParamsT2$Nu)
      expect_equal(m$limParams$T2$moments$nobj, momentsParamsT2$nobj)

      # check parameters for T2 residuals (robust)
      expect_equal(m$limParams$T2$robust$u0, robustParamsT2$u0, tolerance = 10^-5)
      expect_equal(round(m$limParams$T2$robust$Nu), robustParamsT2$Nu)
      expect_equal(m$limParams$T2$robust$nobj, robustParamsT2$nobj)

   }
})

