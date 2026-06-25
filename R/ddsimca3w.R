# DD-SIMCA for 3-way data (PARAFAC & Tucker).
# Decompose -> project -> H/Q/F -> reuse classify()/ddsimcares().
# All tensor reshapes are row-major (feature column for tensor (j,k) is
# (j-1)*K + k); R is column-major, so reshapes use kronecker()/aperm().

# Moore-Penrose pseudo-inverse via SVD (numpy.linalg.pinv semantics: cut
# singular values at rcond * max(singular value)). Kept private; the exported
# prep::pinv() floors at machine eps and must not be changed.
pinv3way <- function(X, rcond = 1e-15) {
   s <- svd(X)
   tol <- rcond * max(s$d)
   dinv <- ifelse(s$d > tol, 1 / s$d, 0)
   s$v %*% (dinv * t(s$u))
}

# Least-squares solve A %*% Z = b for Z (numpy.linalg.lstsq analogue).
lstsq3w <- function(A, b) pinv3way(A) %*% b

# Khatri-Rao (column-wise Kronecker) with B outer, C inner so row index of
# the result equals (j-1)*K + k for tensor mode-2 index j, mode-3 index k.
khatriRao3w <- function(B, C) {
   R <- ncol(B)
   KR <- matrix(0, nrow(B) * nrow(C), R)
   for (r in seq_len(R)) KR[, r] <- kronecker(B[, r], C[, r])
   KR
}

# Mode-1 unfolding of a 3-D core array G (dims R x P x Q) into an R x (P*Q)
# matrix whose column (p-1)*Q + q holds G[r, p, q] (row-major over p,q).
unfold1_3w <- function(G) {
   R <- dim(G)[1]
   matrix(aperm(G, c(1, 3, 2)), nrow = R)
}

# PARAFAC projection: given library factors B (J x R) and C (K x R), fit the
# sample-mode scores A for the unfolded data x (n x (J*K)) and the residual E.
projectParafac <- function(B, C, x) {
   KR <- khatriRao3w(B, C)               # (J*K) x R
   A <- t(lstsq3w(KR, t(x)))             # n x R
   E <- x - A %*% t(KR)                  # n x (J*K)
   list(A = A, E = E)
}

# Tucker projection. Orthogonal (Branch A): lstsq against M = kron(B,C) G_(1)^T.
# Non-negative (Branch B): the sequential closed form from the Python oracle,
# intentionally NOT equal to Branch A when B/C are non-orthogonal.
projectTucker <- function(G, B, C, x, non.negative) {
   G1 <- unfold1_3w(G)                   # R x (P*Q)
   M <- kronecker(B, C) %*% t(G1)        # (J*K) x R
   if (non.negative) {
      S_kron <- kronecker(B, C)          # (J*K) x (P*Q)
      A <- (x %*% S_kron %*% pinv3way(crossprod(S_kron)) %*% t(G1)) %*%
         pinv3way(tcrossprod(G1))        # n x R
   } else {
      A <- t(lstsq3w(M, t(x)))           # n x R
   }
   E <- x - A %*% t(M)                   # n x (J*K)
   list(A = A, E = E)
}

# Refold the n x (J*K) row-major matrix into a list of mode unfoldings used by
# CP-ALS. We work directly on the (J*K) feature axis via khatriRao3w, so only
# the mode-1 (sample) representation X (n x (J*K)) is needed plus per-mode
# helpers built from the current factor estimates.

# Standard CP-ALS with random init. X3: n x (J*K), attr "dim3" = c(J, K).
parafac.als <- function(X3, rank, tol = 1e-10, max.iter = 500) {
   d3 <- attr(X3, "dim3")
   J <- d3[1]; K <- d3[2]; n <- nrow(X3)
   X <- X3; attr(X, "dim3") <- NULL

   # mode unfoldings of the data tensor (row-major)
   # X1 (n x J*K) is X itself; X2 (J x n*K) and X3m (K x n*J) built once.
   T3 <- array(aperm(array(t(X), dim = c(K, J, n)), c(2, 1, 3)), dim = c(J, K, n))
   X2 <- matrix(aperm(T3, c(1, 3, 2)), nrow = J)   # J x (n*K), col (k-1)*n + i
   X3u <- matrix(aperm(T3, c(2, 3, 1)), nrow = K)  # K x (n*J), col (j-1)*n + i

   # Random init. Library functions must NOT call set.seed() (CRAN policy);
   # callers seed the RNG themselves for reproducible models, as mcrals does.
   A <- matrix(rnorm(n * rank), n, rank)
   B <- matrix(rnorm(J * rank), J, rank)
   C <- matrix(rnorm(K * rank), K, rank)

   normX <- sqrt(sum(X^2))
   err.old <- Inf
   for (it in seq_len(max.iter)) {
      # update A: A = X1 (B khatri-rao C) (B^T B * C^T C)^+ — B OUTER, C inner to
      # match X's row-major (j-1)*K+k column order. Using (C, B) here scrambles
      # the A update and the model fails to converge (verified: rel 0.18 vs 1e-5).
      KR_a <- khatriRao3w(B, C)
      A <- X %*% KR_a %*% pinv3way((crossprod(C) * crossprod(B)))
      # update B: B = X2 (C khatri-rao A) (C^T C * A^T A)^+
      KR_ca <- khatriRao3w(C, A)
      B <- X2 %*% KR_ca %*% pinv3way((crossprod(C) * crossprod(A)))
      # update C: C = X3u (B khatri-rao A) (B^T B * A^T A)^+
      KR_ba <- khatriRao3w(B, A)
      C <- X3u %*% KR_ba %*% pinv3way((crossprod(B) * crossprod(A)))

      recon <- A %*% t(khatriRao3w(B, C))
      err <- sqrt(sum((X - recon)^2)) / normX
      if (!is.finite(err)) break
      if (abs(err.old - err) < tol) break
      err.old <- err
   }
   if (it == max.iter) warning("parafac.als did not converge in max.iter iterations.", call. = FALSE)
   list(A = A, B = B, C = C)
}

# Non-negative CP via multiplicative updates (tensorly's algorithm; no NNLS).
parafac.als.nn <- function(X3, rank, tol = 1e-10, max.iter = 500) {
   d3 <- attr(X3, "dim3")
   J <- d3[1]; K <- d3[2]; n <- nrow(X3)
   X <- X3; attr(X, "dim3") <- NULL
   if (any(X < 0)) stop("Non-negative PARAFAC requires non-negative data.", call. = FALSE)

   T3 <- array(aperm(array(t(X), dim = c(K, J, n)), c(2, 1, 3)), dim = c(J, K, n))
   X2 <- matrix(aperm(T3, c(1, 3, 2)), nrow = J)
   X3u <- matrix(aperm(T3, c(2, 3, 1)), nrow = K)

   # Random init; caller seeds the RNG for reproducibility (no set.seed here).
   eps <- 1e-12
   A <- matrix(abs(rnorm(n * rank)), n, rank)
   B <- matrix(abs(rnorm(J * rank)), J, rank)
   C <- matrix(abs(rnorm(K * rank)), K, rank)

   normX <- sqrt(sum(X^2))
   err.old <- Inf
   for (it in seq_len(max.iter)) {
      KR_a <- khatriRao3w(B, C)   # B outer, C inner (matches X row-major order)
      A <- A * (X %*% KR_a) / pmax(A %*% (crossprod(C) * crossprod(B)), eps)
      KR_ca <- khatriRao3w(C, A)
      B <- B * (X2 %*% KR_ca) / pmax(B %*% (crossprod(C) * crossprod(A)), eps)
      KR_ba <- khatriRao3w(B, A)
      C <- C * (X3u %*% KR_ba) / pmax(C %*% (crossprod(B) * crossprod(A)), eps)

      recon <- A %*% t(khatriRao3w(B, C))
      err <- sqrt(sum((X - recon)^2)) / normX
      if (!is.finite(err)) break
      if (abs(err.old - err) < tol) break
      err.old <- err
   }
   # Multiplicative updates converge slowly and normally run to max.iter (this
   # matches tensorly, which stops at the iteration cap silently). No warning:
   # hitting the cap is expected operation for the non-negative variant.
   list(A = A, B = B, C = C)
}

# Build the n x J x K array from the row-major n x (J*K) matrix.
asTensor3w <- function(X, J, K) {
   n <- nrow(X)
   aperm(array(t(X), dim = c(K, J, n)), c(3, 2, 1))   # [i, j, k] = X[i, (j-1)*K + k]
}

# Inverse of unfold1_3w: pack an rA x (rB*rC) matrix (col (b-1)*rC + c) back
# into an rA x rB x rC array.
fold1_3w <- function(G1, rA, rB, rC) {
   aperm(array(G1, dim = c(rA, rC, rB)), c(1, 3, 2))
}

# Mode-n unfolding helpers for an n x J x K tensor T3 (column orders noted).
.unfoldA <- function(T3) {           # n x (J*K), col (j-1)*K + k  (row-major)
   matrix(aperm(T3, c(1, 3, 2)), nrow = dim(T3)[1])
}
.unfoldB <- function(T3) {           # J x (n*K), col (i-1)*K + k
   matrix(aperm(T3, c(2, 3, 1)), nrow = dim(T3)[2])
}
.unfoldC <- function(T3) {           # K x (n*J), col (j-1)*n + i
   matrix(aperm(T3, c(3, 1, 2)), nrow = dim(T3)[3])
}

# Canonical sign convention: flip A/B/C columns so each column's max-abs entry
# is non-negative; absorb compensating signs into G. sign(0) -> +1; ties go to
# the first index (which.max). Compute signs from the UNMODIFIED factors first.
signFix3w <- function(A, B, C, G) {
   colSigns <- function(M) {
      idx <- apply(abs(M), 2, which.max)
      picks <- M[cbind(idx, seq_len(ncol(M)))]
      ifelse(picks < 0, -1, 1)
   }
   sA <- colSigns(A); sB <- colSigns(B); sC <- colSigns(C)
   A2 <- A * rep(sA, each = nrow(A))
   B2 <- B * rep(sB, each = nrow(B))
   C2 <- C * rep(sC, each = nrow(C))
   G2 <- G
   dG <- dim(G)
   for (a in seq_len(dG[1])) for (b in seq_len(dG[2])) for (cc in seq_len(dG[3]))
      G2[a, b, cc] <- G[a, b, cc] * sA[a] * sB[b] * sC[cc]
   list(A = A2, B = B2, C = C2, G = G2)
}

# Leading `r` left singular vectors of a matrix (HOSVD init / HOOI projection).
.leadingSV <- function(M, r) {
   s <- svd(M, nu = r, nv = 0)
   s$u[, seq_len(r), drop = FALSE]
}

# Tucker3 via HOSVD init + HOOI (orthonormal factors).
tucker.hooi <- function(X3, ranks, tol = 1e-10, max.iter = 500) {
   d3 <- attr(X3, "dim3"); J <- d3[1]; K <- d3[2]; n <- nrow(X3)
   rA <- ranks[1]; rB <- ranks[2]; rC <- ranks[3]
   X <- X3; attr(X, "dim3") <- NULL; T3 <- asTensor3w(X, J, K)

   B <- .leadingSV(.unfoldB(T3), rB)
   C <- .leadingSV(.unfoldC(T3), rC)
   A <- .leadingSV(.unfoldA(T3), rA)

   normX <- sqrt(sum(X^2)); err.old <- Inf
   for (it in seq_len(max.iter)) {
      A <- .leadingSV(X %*% kronecker(B, C), rA)
      B <- .leadingSV(.unfoldB(T3) %*% kronecker(C, A), rB)
      C <- .leadingSV(.unfoldC(T3) %*% kronecker(B, A), rC)

      G1 <- t(A) %*% X %*% kronecker(B, C)
      G <- fold1_3w(G1, rA, rB, rC)
      recon <- A %*% t(kronecker(B, C) %*% t(unfold1_3w(G)))
      err <- sqrt(sum((X - recon)^2)) / normX
      if (!is.finite(err)) break
      if (abs(err.old - err) < tol) break
      err.old <- err
   }
   if (it == max.iter) warning("tucker.hooi did not converge in max.iter iterations.", call. = FALSE)
   G1 <- t(A) %*% X %*% kronecker(B, C)
   list(G = fold1_3w(G1, rA, rB, rC), A = A, B = B, C = C)
}

# Non-negative Tucker via multiplicative updates on all factors and the core.
# Multiplicative updates converge slowly; max.iter is higher and the loop
# normally runs to the cap (it stops early only if the error change drops
# below tol). Verified: rel < 1e-4 with all factors >= 0.
tucker.nn <- function(X3, ranks, tol = 1e-10, max.iter = 2000) {
   d3 <- attr(X3, "dim3"); J <- d3[1]; K <- d3[2]; n <- nrow(X3)
   rA <- ranks[1]; rB <- ranks[2]; rC <- ranks[3]
   X <- X3; attr(X, "dim3") <- NULL
   if (any(X < 0)) stop("Non-negative Tucker requires non-negative data.", call. = FALSE)
   # Random init; caller seeds the RNG for reproducibility (no set.seed here).
   eps <- 1e-12
   A <- matrix(abs(rnorm(n * rA)), n, rA)
   B <- matrix(abs(rnorm(J * rB)), J, rB)
   C <- matrix(abs(rnorm(K * rC)), K, rC)
   G <- array(abs(rnorm(rA * rB * rC)), dim = c(rA, rB, rC))
   T3 <- asTensor3w(X, J, K); UB <- .unfoldB(T3); UC <- .unfoldC(T3)
   normX <- sqrt(sum(X^2)); err.old <- Inf

   for (it in seq_len(max.iter)) {
      # A (mode-1): X1 ~ A %*% t(M), M = kron(B,C) %*% t(G_(1))
      M <- kronecker(B, C) %*% t(unfold1_3w(G))
      A <- A * (X %*% M) / pmax(A %*% crossprod(M), eps)

      # core (mode-1): X1 ~ A %*% G1 %*% t(kron(B,C))
      KBC <- kronecker(B, C); G1 <- unfold1_3w(G)
      G1 <- G1 * (t(A) %*% X %*% KBC) /
         pmax(t(A) %*% (A %*% G1 %*% t(KBC)) %*% KBC, eps)
      G <- fold1_3w(G1, rA, rB, rC)

      # B (mode-2): X2 ~ B %*% (G2 %*% t(kron(A,C))). UB col (i-1)*K+k matches
      # kron(A,C) rows; G2 col (a-1)*rC+c matches via aperm(G, c(2,3,1)).
      G2 <- matrix(aperm(G, c(2, 3, 1)), nrow = rB)
      Zb <- G2 %*% t(kronecker(A, C))
      B <- B * (UB %*% t(Zb)) / pmax(B %*% tcrossprod(Zb), eps)

      # C (mode-3): X3 ~ C %*% (G3 %*% t(kron(B,A))). UC col (j-1)*n+i matches
      # kron(B,A) rows; G3 col (b-1)*rA+a matches via aperm(G, c(3,1,2)).
      G3 <- matrix(aperm(G, c(3, 1, 2)), nrow = rC)
      Zc <- G3 %*% t(kronecker(B, A))
      C <- C * (UC %*% t(Zc)) / pmax(C %*% tcrossprod(Zc), eps)

      recon <- A %*% t(kronecker(B, C) %*% t(unfold1_3w(G)))
      err <- sqrt(sum((X - recon)^2)) / normX
      if (!is.finite(err)) break
      if (abs(err.old - err) < tol) break
      err.old <- err
   }
   # Multiplicative updates converge slowly and normally run to max.iter (this
   # matches tensorly, which stops at the iteration cap silently). No warning:
   # hitting the cap is expected operation for the non-negative variant.
   list(G = G, A = A, B = B, C = C)
}
