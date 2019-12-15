#####################################################
# Tests for basic functionality of ldecomp() class  #
#####################################################

# Prepare 3 datasets

data(people)

## plain data, no names no attributes
X1 <- list()
X1$data <- people
rownames(X1$data) <- NULL
colnames(X1$data) <- NULL
X1$exp_exclrows <- NULL
X1$exp_exclcols <- NULL

## same data but with names and attributes
X2 <- list()
X2$data <- people
X2$exp_exclrows <- NULL
X2$exp_exclcols <- NULL
attr(X2$data, "name") <- "People"
attr(X2$data, "xaxis.name") <- "Properties"
attr(X2$data, "xaxis.values") <- 12:1
attr(X2$data, "yaxis.name") <- "Values"
attr(X2$data, "yaxis.values") <- 32:1

## for full PCA model of People data we know some results - use them for double check
Q_exp <- NULL
T2_exp <- NULL
expvar_exp <- NULL

## data with names and attributes and excluded rows
X3 <- X2
X3$exp_exclrows <- c(1, 10, 20)
X3$data <- mda.exclrows(X3$data, X3$exp_exclrows)

## data with names and attributes, excluded rows and columns
X4 <- X3
X4$data <- mda.exclcols(X4$data, c("Hairleng", "IQ"))

# function to get scores, loadings and residuals from data
getPCARes <- function(X, ncomp) {
   rows_excluded <- attr(X, "exclrows")
   cols_excluded <- attr(X, "exclcols")

   if (is.null(rownames(X))) {
      rownames(X) <- paste0('O', 1:nrow(X))
   }

   if (is.null(colnames(X))) {
      colnames(X) <- paste0('X', 1:ncol(X))
   }

   rownames <- rownames(X)
   colnames <- colnames(X)

   X_cal <- X

   # remove excluded rows from calibration data
   if (length(rows_excluded) > 0) {
      X_cal <- X_cal[-rows_excluded, , drop = FALSE]
   }

   # get mean and center and do autoscaling of the calibration data
   m <- apply(X_cal, 2, mean)
   s <- apply(X_cal, 2, sd)
   X_cal <- scale(X_cal, center = m, scale = s)

   # remove excluded columns
   if (length(cols_excluded) > 0) {
      X_cal <- X_cal[, -cols_excluded, drop = FALSE]
   }

   # find loadings
   loadings_visible <- svd(X_cal)$v[, 1:ncomp, drop = F]
   loadings <- matrix(0, nrow = ncol(X), ncol = ncomp)
   if (length(cols_excluded) > 0) {
      loadings[-cols_excluded, ] <- loadings_visible
   } else {
      loadings <- loadings_visible
   }

   # eigenvalues using only visible rows
   eigenvals <- colSums((X_cal %*% loadings_visible)^2) / (nrow(X_cal) - 1)

   X <- scale(X, center = m, scale = s)
   scores <- X %*% loadings
   residuals <- X - tcrossprod(scores, loadings)

   scores <- mda.setattr(scores, mda.getattr(X), type = 'row')
   residuals <- mda.setattr(residuals, mda.getattr(X))

   attr(loadings, "exclrows") <- attr(X, "exclcols")
   attr(loadings, "yaxis.name") <- attr(X, "xaxis.name")
   attr(loadings, "yaxis.values") <- attr(X, "xaxis.values")

   rownames(scores) <- rownames(residuals) <- rownames
   rownames(loadings) <- colnames(residuals) <- colnames
   colnames(scores) <- colnames(loadings) <- paste("Comp", 1:ncomp)

   return(list(scores = scores, loadings = loadings, residuals = residuals,
      eigenvals = eigenvals, totvar = sum(X_cal^2)))
}

##################################################
# Block 1. Computing distances and tnorm values  #
##################################################

context("ldecomp: computing distances and tnorm")

## shortcut for testing function
tf <- function(X) {

   for (A in c(1, 5, ncol(X$data) - length(attr(X$data, "exclcols")))) {
      m <- getPCARes(X$data, A)
      res <- ldecomp.getDistances(m$scores, m$loadings, m$residuals, m$eigenvals)

      # check dimensions
      expect_equal(dim(res$Q), c(32, A))
      expect_equal(dim(res$T2), c(32, A))

      # check values
      for (a in 1:A) {
         E <- m$residuals
         if (a < A) {
            E <- E + m$scores[, (a + 1):A, drop = F] %*% t(m$loadings[, (a + 1):A, drop = F])
         }
         U <- scale(m$scores[, 1:a], center = FALSE, scale = sqrt(m$eigenvals[1:a]))

         expect_equal(res$Q[, a], rowSums(E^2))
         expect_equal(res$T2[, a], rowSums(U^2))
      }


      # check names
      expect_equal(rownames(res$Q), rownames(m$scores))
      expect_equal(colnames(res$Q), colnames(m$loadings))
      expect_equal(rownames(res$T2), rownames(m$scores))
      expect_equal(colnames(res$T2), colnames(m$loadings))

      # check attributes
      expect_equal(attr(res$Q, "exclrows"), X$exp_exclrows)
      expect_equal(attr(res$Q, "yaxis.name"), attr(X$data, "yaxis.name"))
      expect_equal(attr(res$Q, "yaxis.values"), attr(X$data, "yaxis.values"))
      expect_equal(attr(res$T2, "exclrows"), X$exp_exclrows)
      expect_equal(attr(res$T2, "yaxis.name"), attr(X$data, "yaxis.name"))
      expect_equal(attr(res$T2, "yaxis.values"), attr(X$data, "yaxis.values"))
   }
}


test_that("Q and T2 are correct for plain data", { tf(X1) })
test_that("Q and T2 are correct for data attributes", { tf(X2) })
test_that("Q and T2 are correct for data with excluded rows", { tf(X3) })
test_that("Q and T2 are correct for data with excluded columns", { tf(X3) })

#######################################
# Block 2. Computing variance values  #
#######################################

context("ldecomp: computing variance")

## shortcut for testing function
tf <- function(X) {

   for (A in c(1, 5, ncol(X$data) - length(attr(X$data, "exclcols")))) {
      m <- getPCARes(X$data, A)
      dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals, m$eigenvals)
      res <- ldecomp.getVariances(m$scores, m$loadings, m$residuals, dist$Q)

      # check dimensions
      expect_equal(length(res$expvar), A)
      expect_equal(length(res$cumexpvar), A)

      # check values
      for (a in 1:A) {
         E <- m$residuals
         if (a < A) {
            E <- E + m$scores[, (a + 1):A, drop = F] %*% t(m$loadings[, (a + 1):A, drop = F])
         }

         # excluded rows should not be taken into account for calculation of variance
         rows_excluded <- attr(m$scores, "exclrows")
         if (length(rows_excluded) > 0) {
            E <- E[-rows_excluded, ]
         }

         expect_equivalent(res$cumexpvar[a], 100 * (1 - sum(E^2) / m$totvar))
      }

      # check names
      expect_equal(names(res$expvar), colnames(m$loadings))
      expect_equal(names(res$cumexpvar), colnames(m$loadings))
   }
}

test_that("explained variance is correct for plain data", { tf(X1) })
test_that("explained variance is correct for data attributes", { tf(X2) })
test_that("explained variance is correct for data with excluded rows", { tf(X3) })
test_that("explained variance is correct for data with excluded columns", { tf(X3) })


