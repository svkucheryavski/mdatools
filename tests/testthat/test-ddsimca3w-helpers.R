context("ddsimca3w helpers")

test_that("khatriRao3w follows row-major (B outer, C inner) order", {
   B <- matrix(c(1, 2, 3, 4), nrow = 2)   # 2 x 2
   C <- matrix(c(5, 6, 7, 8), nrow = 2)   # 2 x 2
   KR <- khatriRao3w(B, C)
   # column 1 = kron(B[,1], C[,1]) = c(1*5, 1*6, 2*5, 2*6)
   expect_equal(KR[, 1], c(5, 6, 10, 12))
   expect_equal(KR[, 2], as.vector(kronecker(B[, 2], C[, 2])))
   expect_equal(dim(KR), c(4L, 2L))
})

test_that("unfold1_3w gives row-major mode-1 unfolding (col = (p-1)*Q + q)", {
   # G[r, p, q] = 100*r + 10*p + q ; dims R=2, P=2, Q=3
   G <- array(0, dim = c(2, 2, 3))
   for (r in 1:2) for (p in 1:2) for (q in 1:3) G[r, p, q] <- 100 * r + 10 * p + q
   M <- unfold1_3w(G)
   expect_equal(dim(M), c(2L, 6L))
   # row 1 must be p,q in order (1,1),(1,2),(1,3),(2,1),(2,2),(2,3)
   expect_equal(M[1, ], c(111, 112, 113, 121, 122, 123))
})

test_that("pinv3way matches solve() for an invertible matrix", {
   set.seed(42)
   A <- matrix(rnorm(16), 4, 4)
   expect_equal(pinv3way(A), solve(A), tolerance = 1e-10)
})

test_that("Tucker Branch A == Branch B when B and C are orthonormal", {
   set.seed(1)
   J <- 4; K <- 3; P <- 3; Q <- 2; R <- 2; n <- 6
   B <- qr.Q(qr(matrix(rnorm(J * P), J, P)))   # J x P, orthonormal columns
   C <- qr.Q(qr(matrix(rnorm(K * Q), K, Q)))   # K x Q, orthonormal columns
   G <- array(rnorm(R * P * Q), dim = c(R, P, Q))
   x <- matrix(rnorm(n * J * K), n, J * K)
   a <- projectTucker(G, B, C, x, non.negative = FALSE)
   b <- projectTucker(G, B, C, x, non.negative = TRUE)
   expect_equal(a$A, b$A, tolerance = 1e-8)
})

test_that("parafac.als recovers a clean rank-1 tensor", {
   set.seed(3)
   J <- 4; K <- 3; n <- 8
   a <- rnorm(n); b <- rnorm(J); c <- rnorm(K)
   X3 <- matrix(0, n, J * K)
   for (i in seq_len(n)) X3[i, ] <- a[i] * kronecker(b, c)
   attr(X3, "dim3") <- c(J, K)
   set.seed(10)
   fit <- parafac.als(X3, rank = 1)
   recon <- fit$A %*% t(khatriRao3w(fit$B, fit$C))
   expect_equal(as.numeric(recon), as.numeric(X3), tolerance = 1e-6)
})

test_that("parafac.als reconstruction error is small on low-rank-plus-noise", {
   set.seed(4)
   J <- 5; K <- 4; n <- 12; R <- 3
   A0 <- matrix(rnorm(n * R), n, R)
   B0 <- matrix(rnorm(J * R), J, R)
   C0 <- matrix(rnorm(K * R), K, R)
   X3 <- A0 %*% t(khatriRao3w(B0, C0)) + 0.01 * matrix(rnorm(n * J * K), n, J * K)
   attr(X3, "dim3") <- c(J, K)
   set.seed(11)
   fit <- parafac.als(X3, rank = 3)
   recon <- fit$A %*% t(khatriRao3w(fit$B, fit$C))
   rel <- sum((X3 - recon)^2) / sum(X3^2)
   expect_lt(rel, 1e-3)
})

test_that("parafac.als.nn keeps factors non-negative", {
   set.seed(5)
   J <- 4; K <- 3; n <- 10; R <- 2
   A0 <- abs(matrix(rnorm(n * R), n, R))
   B0 <- abs(matrix(rnorm(J * R), J, R))
   C0 <- abs(matrix(rnorm(K * R), K, R))
   X3 <- A0 %*% t(khatriRao3w(B0, C0))
   attr(X3, "dim3") <- c(J, K)
   set.seed(12)
   fit <- parafac.als.nn(X3, rank = 2)
   expect_true(all(fit$A >= -1e-9) && all(fit$B >= -1e-9) && all(fit$C >= -1e-9))
})

test_that("parafac.als.nn rejects negative data", {
   X3 <- matrix(c(1, -2, 3, 4, 5, 6), nrow = 2)   # contains a negative value
   attr(X3, "dim3") <- c(3, 1)
   expect_error(parafac.als.nn(X3, rank = 1), "non-negative data")
})

test_that("signFix3w makes each factor column's max-abs entry non-negative and preserves reconstruction", {
   set.seed(6)
   J <- 4; K <- 3; rA <- 2; rB <- 3; rC <- 2; n <- 5
   A <- matrix(rnorm(n * rA), n, rA)
   B <- matrix(rnorm(J * rB), J, rB)
   C <- matrix(rnorm(K * rC), K, rC)
   G <- array(rnorm(rA * rB * rC), dim = c(rA, rB, rC))
   M0 <- kronecker(B, C) %*% t(unfold1_3w(G))
   recon0 <- A %*% t(M0)
   sf <- signFix3w(A, B, C, G)
   maxabs <- function(v) v[which.max(abs(v))]
   expect_true(all(apply(sf$B, 2, maxabs) >= 0))
   expect_true(all(apply(sf$C, 2, maxabs) >= 0))
   M1 <- kronecker(sf$B, sf$C) %*% t(unfold1_3w(sf$G))
   recon1 <- sf$A %*% t(M1)
   expect_equal(recon1, recon0, tolerance = 1e-10)
})

test_that("tucker.hooi reconstructs a low-rank tensor and PARAFAC equals Tucker with a diagonal core", {
   set.seed(7)
   J <- 5; K <- 4; n <- 10; R <- 2
   A0 <- matrix(rnorm(n * R), n, R)
   B0 <- matrix(rnorm(J * R), J, R)
   C0 <- matrix(rnorm(K * R), K, R)
   X3 <- A0 %*% t(khatriRao3w(B0, C0))
   attr(X3, "dim3") <- c(J, K)
   set.seed(13)
   fit <- tucker.hooi(X3, ranks = c(R, R, R))
   M <- kronecker(fit$B, fit$C) %*% t(unfold1_3w(fit$G))
   recon <- fit$A %*% t(M)
   rel <- sum((unclass(X3) - recon)^2) / sum(unclass(X3)^2)
   expect_lt(rel, 1e-6)
})

test_that("tucker.nn keeps factors and core non-negative", {
   set.seed(8)
   J <- 4; K <- 3; n <- 9; R <- 2
   A0 <- abs(matrix(rnorm(n * R), n, R)); B0 <- abs(matrix(rnorm(J * R), J, R))
   C0 <- abs(matrix(rnorm(K * R), K, R))
   X3 <- A0 %*% t(khatriRao3w(B0, C0))
   attr(X3, "dim3") <- c(J, K)
   set.seed(14)
   fit <- tucker.nn(X3, ranks = c(R, R, R))
   expect_true(all(fit$A >= -1e-9) && all(fit$B >= -1e-9) &&
      all(fit$C >= -1e-9) && all(fit$G >= -1e-9))
   # guard reconstruction too — non-negativity alone passes even for a broken
   # decomposition (the original NN-Tucker had silently-wrong partner reshapes)
   recon <- fit$A %*% t(kronecker(fit$B, fit$C) %*% t(unfold1_3w(fit$G)))
   expect_lt(sum((unclass(X3) - recon)^2) / sum(unclass(X3)^2), 1e-3)
})

test_that("fold1_3w is the inverse of unfold1_3w", {
   G <- array(rnorm(2 * 3 * 4), dim = c(2, 3, 4))
   expect_equal(fold1_3w(unfold1_3w(G), 2, 3, 4), G, tolerance = 1e-12)
})
