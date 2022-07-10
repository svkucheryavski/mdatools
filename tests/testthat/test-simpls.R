######################################
# Tests for SIMPLS algorithms     #
######################################

setup({
   pdf(file = tempfile("mdatools-test-simpls-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-simpls-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

#################################################
# Block 1. Tests the old algorithm              #
#################################################

context("simpls: PLS2 example from the paper")

# convert results produced by SIMPLS algorithm to form
# suitable for testing
simpls2res <- function (m, X, Y, A) {
   R <- m$weights
   P <- m$xloadings
   Q <- m$yloadings

   T <- X %*% R
   Xexp <- rep(0, A)
   Yexp <- rep(0, A)

   for (a in 1:A) {
      B <- tcrossprod(R[, a, drop = FALSE], Q[, a, drop = FALSE])
      Xhat <- tcrossprod(T[, a, drop = FALSE], P[, a, drop = FALSE])
      Yhat <- X %*% B
      Xexp[a] <- sum(Xhat^2)/sum(X^2) * 100
      Yexp[a] <- sum(Yhat^2)/sum(Y^2) * 100
   }

   rnorm = sqrt(colSums(R^2))
   R = R %*% diag(1 / rnorm)

   return(list(R = R, T = T, P = P, Q = Q, Xexp = Xexp, Yexp = Yexp))
}

# deterministic small data
A <- 3
X <- matrix(c(-4, -4, 4, 4, 2, -2, 2, -2, 1, -1, -1, 1), nrow = 4)
Y <- matrix(c(430, -436, -361, 367, -94, 12, -22, 104), nrow = 4)

# expected values
expected = list(
   R = matrix(c(0.0164, 0.1798, 0.9836, 0.4562, -0.7719, 0.4428, 0.3392, 0.7141, -0.6124), ncol = A),
   T = matrix(c(0.6088, -0.6712, -0.2662, 0.3286, -0.6018, -0.1489, -0.0333, 0.7840, -0.1311, -0.5266, 0.8234, -0.1657), ncol = A),
   Xexp = c(6.72, 50.77, 42.51),
   Yexp = c(90.15, 4.54, 5.31)
)

test_that("new algorithm works correctly", {
   res <- simpls2res(pls.simpls(X, Y, A), X, Y, A)

   expect_equivalent(abs(res$R), abs(expected$R), tolerance = 10^-4)
   expect_equivalent(abs(res$T), abs(expected$T), tolerance = 10^-4)
   expect_equivalent(res$Xexp, expected$Xexp, tolerance = 10^-4)
   expect_equivalent(res$Yexp, expected$Yexp, tolerance = 10^-4)
})

test_that("old algorithm works correctly", {
   res <- simpls2res(pls.simplsold(X, Y, A), X, Y, A)

   expect_equivalent(abs(res$R), abs(expected$R), tolerance = 10^-4)
   expect_equivalent(abs(res$T), abs(expected$T), tolerance = 10^-4)
   expect_equivalent(res$Xexp, expected$Xexp, tolerance = 10^-4)
   expect_equivalent(res$Yexp, expected$Yexp, tolerance = 10^-4)
})


test_that("new algorithm is more numerically stable", {
   Xr <- matrix(rnorm(30000 * 1000), 30000, 1000)
   Yr <- matrix(rnorm(30000 * 2), 30000, 2)

   expect_warning(pls.simplsold(Xr, Yr / 100000000, 50))
   expect_silent(pls.simpls(Xr, Yr / 100000000, 50))
})


test_that("new algorithm gives results comparable to other software", {

   # read model parameters made in PLS_Toolbox
   dataFolder = file.path(system.file("/inst/testdata/", package="mdatools"))
   weights <- as.matrix(read.delim(paste0(dataFolder, "/plstlbx-people-weights.csv"), sep = " ", header = FALSE))
   xloadings <- as.matrix(read.delim(paste0(dataFolder, "/plstlbx-people-xloadings.csv"), sep = " ", header = FALSE))
   xscores <- as.matrix(read.delim(paste0(dataFolder, "/plstlbx-people-xscores.csv"), sep = " ", header = FALSE))
   yscores <- as.matrix(read.delim(paste0(dataFolder, "/plstlbx-people-yscores.csv"), sep = " ", header = FALSE))
   yloadings <- c(5.3643, 1.0338, 0.4675, 0.3567)
   coeffs <- c(0.2078, 0.2647, 0.0073, 0.0722, -0.0016, 0.1829, 0.1420, -0.1984, 0.2153, 0.0151, -0.0405)

   # make a model
   data(people)
   X <- scale(people[, -4], center = TRUE, scale = TRUE)
   y <- scale(people[, 4, drop = FALSE], center = TRUE, scale = TRUE)
   m <- pls.simpls(X, y, 4)

   # here we re-normalize results from PLS_Toolbox
   xnorm <- sqrt(colSums(xscores^2))
   expect_equivalent(m$xloadings, xloadings %*% diag(xnorm), tolerance = 10^-3)
   expect_equivalent(m$xscores, xscores %*% diag(1/xnorm), tolerance = 10^-3)
   expect_equivalent(m$weights, weights %*% diag(1/xnorm), tolerance = 10^-3)

   expect_equivalent(m$yscores, yscores, tolerance = 10^-4)
   expect_equivalent(m$yloadings, yloadings, tolerance = 10^-4)
   expect_equivalent(m$coeffs[, 4, 1], coeffs, tolerance = 10^-4)
})
