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
   R <- dim(G)[1]
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
