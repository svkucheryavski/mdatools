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
