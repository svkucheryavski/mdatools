#####################################################
# Tests for plotting methods of ldecomp() class     #
#####################################################

setup({
   pdf(file = tempfile("mdatools-test-ldecompplots-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-ldecompplots-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# function to get scores, loadings and residuals from data
getPCARes <- function(X, ncomp) {
   rows_excluded <- attr(X, "exclrows")
   cols_excluded <- attr(X, "exclcols")

   if (is.null(rownames(X))) {
      rownames(X) <- paste0("O", 1:nrow(X))
   }

   if (is.null(colnames(X))) {
      colnames(X) <- paste0("X", 1:ncol(X))
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

   scores <- mda.setattr(scores, mda.getattr(X), type = "row")
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

data(people)
x <- people
x <- prep.autoscale(x, center = TRUE, scale = TRUE)
m <- getPCARes(x, 5)
obj <- ldecomp(m$scores, m$loadings, m$residuals, m$eigenvals)

context("ldecomp: plots")

test_that("scores plot can return plot data.", {
   pd <- plotScores(obj, c(1, 3), show.plot = FALSE)
   expect_true("matrix" %in% class(pd))
   expect_equivalent(pd[, 1], obj$scores[, 1])
   expect_equivalent(pd[, 2], obj$scores[, 3])
})

test_that("scores plot works fine with different attributes.", {
   expect_silent({
      par(mfrow = c(2, 2))
      plotScores(obj)
      plotScores(obj, c(1, 3), show.labels = T)
      plotScores(obj, c(1, 2), cgroup = x[, 1], show.labels = T, show.axes = F)
      plotScores(obj, 1, show.labels = T, cgroup = x[, "Sex"], colmap = c("red", "green"), pch = 17)
   })

   expect_silent({
      par(mfrow = c(2, 2))
      plotScores(obj, type = "l")
      plotScores(obj, c(1, 3), type = "h", show.labels = T)
      plotScores(obj, c(1, 2), type = "b", show.labels = T, show.axes = F)
      plotScores(obj, 1, type = "l", show.labels = T, colmap = c("red", "green"), pch = 17)
   })

})

test_that("scores plot works fine with cgroup and convex hulls.", {

   cgroup <- interaction(
      factor(x[, "Sex"], labels = c("M", "F")),
      factor(x[, "Region"], labels = c("S", "M"))
   )

   expect_silent({
      par(mfrow = c(2, 2))
      p <- plotScores(obj, cgroup = cgroup)
      plotConvexHull(p)
      p <- plotScores(obj, cgroup = cgroup, colmap = c("red", "green"))
      plotConvexHull(p, opacity = 0.25)
      p <- plotScores(obj, cgroup = cgroup)
      plotConfidenceEllipse(p)
      p <- plotScores(obj, cgroup = cgroup, colmap = c("red", "green"))
      plotConfidenceEllipse(p, conf.level = 0.99, opacity = 0.25)

   })
})

test_that("scores plot works fine with Hotelling ellipse.", {

   expect_silent({
      par(mfrow = c(2, 2))
      p <- plotScores(obj, xlim = c(-8, 8), ylim = c(-8, 8))
      plotHotellingEllipse(p)
      p <- plotScores(obj, c(2, 3), xlim = c(-8, 8), ylim = c(-8, 8))
      plotHotellingEllipse(p)
      p <- plotScores(obj, xlim = c(-8, 8), ylim = c(-8, 8))
      plotHotellingEllipse(p, conf.lim = 0.90, col = "red", lty = 1, lwd = 0.5)
      p <- plotScores(obj, c(2, 3), xlim = c(-8, 8), ylim = c(-8, 8))
      plotHotellingEllipse(p, conf.lim = 0.90, col = "red", lty = 1, lwd = 0.5)

   })
})

test_that("residuals plot can return plot data.", {
   pd <- plotResiduals(obj, ncomp = 4, show.plot = FALSE)
   expect_true("matrix" %in% class(pd))
   expect_equivalent(pd[, 1], obj$T2[, 4])
   expect_equivalent(pd[, 2], obj$Q[, 4])
})

test_that("residuals plot works fine.", {
   expect_silent({
      par(mfrow = c(2, 2))
      plotResiduals(obj)
      plotResiduals(obj, 1, show.labels = T)
      plotResiduals(obj, 2, cgroup = x[, 1], show.labels = T)
      plotResiduals(obj, 3, cgroup = x[, 1], show.labels = T, colmap = c("red", "green"), pch = 17)
   })
})

test_that("varance plot can return plot data.", {
   pd <- plotVariance(obj, show.plot = FALSE)
   expect_equal(class(pd), "numeric")
   expect_equivalent(pd, obj$expvar)
})

test_that("variance plot works fine.", {
   expect_silent({
      par(mfrow = c(2, 2))
      plotVariance(obj)
      plotVariance(obj, type = "h")
      plotVariance(obj, type = "h", col = "red")
      plotVariance(obj, type = "h", col = "red", show.labels = TRUE)
   })
})

test_that("cumulative varance plot can return plot data.", {
   pd <- plotCumVariance(obj, show.plot = FALSE)
   expect_equal(class(pd), "numeric")
   expect_equivalent(pd, obj$cumexpvar)
})

test_that("variance plot works fine.", {
   expect_silent({
      par(mfrow = c(2, 2))
      plotCumVariance(obj)
      plotCumVariance(obj, type = "h")
      plotCumVariance(obj, type = "h", col = "red", show.labels = TRUE)
      plotCumVariance(obj, type = "h", col = "red", show.labels = TRUE, labels = "names")
   })
})
