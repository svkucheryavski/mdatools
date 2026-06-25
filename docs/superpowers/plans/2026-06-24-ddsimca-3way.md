# DD-SIMCA for 3-way data (PARAFAC & Tucker) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `ddsimca.parafac()` and `ddsimca.tucker()` — DD-SIMCA one-class classifiers for 3-way data — that reuse the existing `classify()`/`ddsimcares()` machinery so the result class is identical across PCA, PARAFAC and Tucker variants.

**Architecture:** New 3-way models decompose the unfolded data (PARAFAC or Tucker), project to get scores `A` and residual `E`, compute the three distances `H` (Mahalanobis score distance), `Q` (orthogonal/residual distance) and the scaled full distance `F`, then hand `H`/`Q` plus a per-sub-model full-distance parameter set to the existing `classify()` → `ddsimcares()` pipeline. All new model/projection/decomposition code lives in one new file `R/ddsimca3w.R`; the only edits to shipped files are a backward-compatible `fp=NULL` argument on `classify()`, a `f0.model` field threaded through three result plot methods, an `ncol(T2)` keying fix in `as.matrix.ddsimcares`, and converting the three DoF plot functions to S3 generics.

**Tech Stack:** R (S3 OOP), base `stats`/`graphics` only — **zero new package dependencies**. testthat for tests. Python (`../../ddsimca-py`, with its `.venv`: numpy 2.4.0, scipy 1.16.3, tensorly 0.9.0) is a numerical oracle used ONLY to generate golden fixtures; it is not a runtime dependency.

## Global Constraints

- R >= 3.5.0; package uses S3 OOP throughout. Follow existing `R/ddsimca.R` style (4-space indent, `<-` assignment, `function(...)` spacing).
- **Zero new package dependencies.** Roll all decompositions in base R. Do not add to `DESCRIPTION` `Imports`.
- Max line length 120 chars; max cyclomatic complexity 20 (`.lintr`).
- **3-way scope exclusions (per package author):** NO excluded rows/columns, NO preprocessing (`prep`), NO `center`/`scale`. Constructors take a clean unfolded matrix and decompose as-is. PCV (`$res$pv`) and JSON export are deferred (not in this plan).
- Tensor unfolding is **row-major**: feature column index for tensor element `(j, k)` is `(j-1)*K + k` (1-indexed), i.e. mode-2 index `j` outer, mode-3 index `k` inner. R is column-major, so every reshape that must be row-major is built with `kronecker()` / `aperm()` (never bare `matrix()`).
- Model class: `c("ddsimca.parafac", "ddsimca3w", "ddsimca")` / `c("ddsimca.tucker", "ddsimca3w", "ddsimca")` — **no `"pca"`** (see spec "Why no `pca` class"). Result class unchanged: `c("ddsimcares", "pcares", "ldecomp")`.
- `classify(fp=NULL)` MUST be bit-identical to current behavior for PCA / scalar `ddsimca` (regression-tested).
- Commits: end message with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` is the harness default, **BUT this project's convention is NO co-author line** — omit it. Commit only what each task touches. Do not stage/commit unrelated working-tree changes already present (`DESCRIPTION`, `NEWS.md`, etc.).
- Run tests with `Rscript -e 'devtools::test()'`. Regenerate docs with `Rscript -e 'roxygen2::roxygenize()'`.
- Reference spec: `docs/superpowers/specs/2026-06-22-ddsimca-3way-design.md`. Python math reference: `../../ddsimca-py/ddsimca/ddsimca_parafac.py`, `ddsimca_tucker.py`, `_base.py`.

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `R/ddsimca3w.R` | **Create** | Everything 3-way: tensor/projection utilities, PARAFAC & Tucker decompositions, sign-fix, `project3w()`, `finalize3w()`, `ddsimca.parafac()`/`ddsimca.tucker()` constructors, `predict.ddsimca3w`, `selectCompNum.ddsimca3w`, `summary`/`print.ddsimca3w`, `plotFactors` + adapters, DoF `.ddsimca3w` methods. |
| `R/ddsimca.R` | Modify | `classify()` gains trailing `fp = NULL` storing `f0.model`; `plotDistances.ddsimca` f-branch uses `f0.model`. |
| `R/ddsimcares.R` | Modify | `f0.model` threaded into `plotExtremes`/`plotDistances`/`plotSelectivityArea`; `ncol(x$T2)` keying in `as.matrix`/`as.data.frame`/`plotAcceptance`/`plotDistances` bounds. |
| `R/pca.R` | Modify | `plotT2DoF`/`plotQDoF`/`plotDistDoF` converted from plain functions to S3 generics (`UseMethod`) with current bodies as `.pca` methods. |
| `R/defaults.R` | Modify | Add `plotFactors` generic; `plotT2DoF`/`plotQDoF`/`plotDistDoF` generic stubs live in `pca.R` (see Task 6 — keep them next to the `.pca` methods). |
| `NAMESPACE`, `man/*.Rd` | Regenerate | via `roxygen2::roxygenize()`. |
| `tests/testthat/test-ddsimca3w.R` | **Create** | Unit + integration + golden tests for all of the above. |
| `tests/testthat/test-ddsimca3w-helpers.R` | **Create** | Pure-function tests for utilities/decompositions (Tasks 1–3). |
| `../../ddsimca-py/tests/dump_3way_golden.py` | **Create** | Python script dumping per-sub-model factor matrices + result matrices as CSV (Task 7). |
| `tests/testthat/ddsimca3w/` | **Create** | Directory holding the golden CSV fixtures committed from the dump script. |

---

## Task 1: Tensor & projection utilities

**Files:**
- Create: `R/ddsimca3w.R` (start the file with these private helpers)
- Test: `tests/testthat/test-ddsimca3w-helpers.R`

**Interfaces:**
- Produces:
  - `pinv3way(X, rcond = 1e-15)` → matrix (Moore-Penrose pseudo-inverse, numpy semantics).
  - `lstsq3w(A, b)` → matrix, equals `pinv3way(A) %*% b`.
  - `khatriRao3w(B, C)` → `(nrow(B)*nrow(C)) x ncol(B)` matrix; column `r` is `kronecker(B[, r], C[, r])`.
  - `unfold1_3w(G)` → mode-1 unfolding of 3-D array `G` with dims `(R, P, Q)`: `R x (P*Q)` matrix, column index `(p-1)*Q + q` holds `G[r, p, q]`.
  - `projectParafac(B, C, x)` → list(`A` = `n x R` scores, `E` = `n x (J*K)` residual).
  - `projectTucker(G, B, C, x, non.negative)` → list(`A` = `n x R` scores, `E` = `n x (J*K)` residual). `R = dim(G)[1]`.

- [ ] **Step 1: Write the failing tests for the reshape/Khatri-Rao contract**

Add to `tests/testthat/test-ddsimca3w-helpers.R`:

```r
context("ddsimca3w helpers")

test_that("khatriRao3w follows row-major (B outer, C inner) order", {
   B <- matrix(c(1, 2, 3, 4), nrow = 2)   # 2 x 2
   C <- matrix(c(5, 6, 7, 8), nrow = 2)   # 2 x 2
   KR <- khatriRao3w(B, C)
   # column 1 = kron(B[,1], C[,1]) = c(1*5, 1*6, 2*5, 2*6)
   expect_equal(KR[, 1], c(5, 6, 10, 12))
   expect_equal(KR[, 2], kronecker(B[, 2], C[, 2]))
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w-helpers")'`
Expected: FAIL — could not find function "khatriRao3w" (etc.).

- [ ] **Step 3: Implement the utilities**

Create `R/ddsimca3w.R` beginning with a roxygen file banner comment and these helpers:

```r
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w-helpers")'`
Expected: PASS (3 tests).

- [ ] **Step 5: Add the Branch-A ≡ Branch-B invariant test**

Append to `tests/testthat/test-ddsimca3w-helpers.R`:

```r
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
```

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w-helpers")'`
Expected: PASS (4 tests).

- [ ] **Step 6: Commit**

```bash
git add R/ddsimca3w.R tests/testthat/test-ddsimca3w-helpers.R
git commit -m "Add ddsimca3w tensor and projection utilities"
```

---

## Task 2: PARAFAC decomposition

**Files:**
- Modify: `R/ddsimca3w.R` (append)
- Test: `tests/testthat/test-ddsimca3w-helpers.R`

**Interfaces:**
- Consumes: `khatriRao3w` (Task 1).
- Produces:
  - `parafac.als(X3, rank, tol = 1e-10, max.iter = 500)` → list(`A` = `n x rank`, `B` = `J x rank`, `C` = `K x rank`). `X3` is the `n x (J*K)` unfolded matrix; `dim` carried via `attr(X3, "dim3") = c(J, K)`. Standard CP-ALS, random init.
  - `parafac.als.nn(X3, rank, tol = 1e-10, max.iter = 500)` → same shape, non-negative (multiplicative updates). Requires `X3 >= 0`.

**Implementation note:** R's RNG differs from numpy/tensorly, so factors will NOT match the Python oracle bit-for-bit. Correctness is verified by math invariants (reconstruction error, rank-1 recovery), never by comparing factor values to Python. The DD-SIMCA numbers are validated separately against the oracle by feeding it fixed Python factors (Task 7/8).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-ddsimca3w-helpers.R`:

```r
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w-helpers")'`
Expected: FAIL — could not find function "parafac.als".

- [ ] **Step 3: Implement `parafac.als` and `parafac.als.nn`**

Append to `R/ddsimca3w.R`:

```r
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w-helpers")'`
Expected: PASS (7 tests). If the rank-1 recovery test is flaky, the ALS update or refold is wrong — debug the `T3`/`X2`/`X3u` construction before proceeding (this is the most error-prone reshape in the plan).

- [ ] **Step 5: Commit**

```bash
git add R/ddsimca3w.R tests/testthat/test-ddsimca3w-helpers.R
git commit -m "Add PARAFAC decomposition (ALS and non-negative)"
```

---

## Task 3: Tucker decomposition and sign-fix

**Files:**
- Modify: `R/ddsimca3w.R` (append)
- Test: `tests/testthat/test-ddsimca3w-helpers.R`

**Interfaces:**
- Consumes: `unfold1_3w`, `projectTucker`, `khatriRao3w` (Task 1).
- Produces:
  - `tucker.hooi(X3, ranks, tol = 1e-10, max.iter = 500)` → list(`G` = array `c(rA, rB, rC)`, `A` = `n x rA`, `B` = `J x rB`, `C` = `K x rC`). `ranks = c(rA, rB, rC)`. HOSVD init + HOOI.
  - `tucker.nn(X3, ranks, tol = 1e-10, max.iter = 2000)` → same, non-negative (multiplicative updates; higher cap, runs to it silently).
  - `signFix3w(A, B, C, G)` → list(A, B, C, G) with the canonical max-abs-positive sign convention applied to A/B/C columns and compensating signs absorbed into G (reconstruction-preserving).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-ddsimca3w-helpers.R`:

```r
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w-helpers")'`
Expected: FAIL — could not find function "signFix3w".

- [ ] **Step 3: Implement helpers + `signFix3w` + `tucker.hooi` + `tucker.nn`**

Append to `R/ddsimca3w.R`. All reshapes below were verified numerically (HOOI reconstructs a genuine distinct-rank Tucker tensor to ~1e-30; non-negative Tucker to rel < 1e-4 with all factors >= 0). The mode-2/mode-3 partner orders (`kron(A,C)` for the B update, `kron(B,A)` for the C update) and the core unfoldings `G2`/`G3` are the load-bearing, easy-to-break parts — do not "simplify" them.

```r
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
```

- [ ] **Step 3b: Add a fold/unfold round-trip assertion test**

Append to `tests/testthat/test-ddsimca3w-helpers.R`:

```r
test_that("fold1_3w is the inverse of unfold1_3w", {
   G <- array(rnorm(2 * 3 * 4), dim = c(2, 3, 4))
   expect_equal(fold1_3w(unfold1_3w(G), 2, 3, 4), G, tolerance = 1e-12)
})
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w-helpers")'`
Expected: PASS (all helper tests). Fix reshapes before moving on if `tucker.hooi` reconstruction or `signFix3w` round-trip fails.

- [ ] **Step 5: Commit**

```bash
git add R/ddsimca3w.R tests/testthat/test-ddsimca3w-helpers.R
git commit -m "Add Tucker decomposition (HOOI and non-negative) and sign-fix"
```

---

## Task 4: `classify(fp=NULL)` and `f0.model`

**Files:**
- Modify: `R/ddsimca.R:482-591` (the `classify` function)
- Test: `tests/testthat/test-ddsimca3w.R` (create) + existing `tests/testthat/test-ddsimca.R` (regression)

**Interfaces:**
- Consumes: nothing new.
- Produces: `classify(pcares, indices, numbers, qp, hp, nobj.cal, alpha, gamma, classname, c.ref = NULL, fp = NULL)`. When `fp = NULL`, behavior is identical to today (`Nf = Nh + Nq`, critical limits via plain `qchisq`), and each per-component outcome additionally carries `f0.model = Nf`. When `fp` is a list with `$u0`, `$Nu` (per component), `Nf = fp$Nu`, `f0.model = fp$u0`, and critical limits scale by `f0.model / Nf`.

- [ ] **Step 1: Write the failing regression + behavior tests**

Create `tests/testthat/test-ddsimca3w.R`:

```r
context("ddsimca3w: classify fp")

test_that("classify(fp=NULL) is bit-identical to the pre-fp behavior for a small PCA case", {
   data(people)
   x <- people[, -4]
   m <- ddsimca(x, "people", ncomp = 3, center = TRUE, scale = TRUE, pcv = NULL)
   v <- m$res$cal$simca$outcomes$moments$values[[3]]
   # f0.model must equal Nf for the PCA (sum-rule) path
   expect_equal(v$f0.model, v$Nf)
   # fce reproduces the plain chi-square critical limit
   expect_equal(v$fce, qchisq(1 - m$alpha, v$Nf))
})

test_that("classify with fp scales critical limits by f0.model/Nf", {
   data(people)
   x <- people[, -4]
   m <- ddsimca(x, "people", ncomp = 3, center = TRUE, scale = TRUE, pcv = NULL)
   r <- m$res$cal
   lp <- m$limParams
   ind <- r$simca$indices; num <- r$simca$numbers
   # craft a fake fp with Nf=10, f0=20 for every component
   ncomp <- ncol(r$T2)
   fp <- list(u0 = rep(20, ncomp), Nu = rep(10, ncomp))
   out <- classify(r, ind, num, lp$Q$moments, lp$T2$moments, m$nrows,
      m$alpha, m$gamma, "people", r$simca$c.ref, fp = fp)
   v <- out$values[[3]]
   expect_equal(v$Nf, 10)
   expect_equal(v$f0.model, 20)
   expect_equal(v$fce, 20 / 10 * qchisq(1 - m$alpha, 10))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: FAIL — `v$f0.model` is NULL (field does not exist yet).

- [ ] **Step 3: Edit `classify` in `R/ddsimca.R`**

Change the signature (line 482) to add the trailing argument. Wrap it across two lines — adding `, fp = NULL` to the current one-line signature pushes it past the 120-char lint limit:

```r
classify <- function(pcares, indices, numbers, qp, hp, nobj.cal, alpha = 0.05, gamma = 0.01,
   classname, c.ref = NULL, fp = NULL) {
```

Replace the distribution-parameter / critical-limit block (current lines 489-498) with:

```r
   # get distribution parameters
   q0 <- qp$u0
   Nq <- qp$Nu
   h0 <- hp$u0
   Nh <- hp$Nu

   # full-distance DoF and scale. fp = NULL keeps the basic DD-SIMCA sum rule
   # (Nf = Nh + Nq, scale f0.model = Nf) so PCA / scalar ddsimca is unchanged.
   # fp supplied (3-way data-driven F) uses the fitted (f0, Nf) directly.
   if (is.null(fp)) {
      Nf <- Nh + Nq
      f0.model <- Nf
   } else {
      Nf <- fp$Nu
      f0.model <- fp$u0
   }

   # compute critical limits (scale collapses to plain qchisq when f0.model == Nf)
   fCritE <- f0.model / Nf * qchisq(1 - alpha, Nf)
   fCritO <- f0.model / Nf * qchisq((1 - gamma)^(1 / nobj.cal), Nf)
```

Inside the per-component loop, add `f0.model` to the local scalars and into the `res` list. After `Nqa <- Nq[a]` (current line 531) add:

```r
      f0ma <- f0.model[a]
```

In the `res <- list(...)` add the field (next to `Nf = Nfa`):

```r
         f0.model = f0ma,
```

- [ ] **Step 4: Run the new tests AND the full ddsimca regression suite**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: PASS (3 tests).

Run: `Rscript -e 'devtools::test(filter = "ddsimca")'`
Expected: PASS — `test-ddsimca.R` unchanged (the new `f0.model` field is additive; all existing assertions hold).

- [ ] **Step 5: Commit**

```bash
git add R/ddsimca.R tests/testthat/test-ddsimca3w.R
git commit -m "Add fp argument and f0.model to classify for data-driven full distance"
```

---

## Task 5: Thread `f0.model` through result methods + `ncol(T2)` keying

**Files:**
- Modify: `R/ddsimcares.R` — `plotExtremes.ddsimcares:711-712`, `plotDistances.ddsimcares:1063,1087-1097`, `plotSelectivityArea.ddsimcares:467-468`, `as.matrix.ddsimcares:188,222`, `as.data.frame.ddsimcares:133`
- Modify: `R/ddsimca.R` — `plotDistances.ddsimca:1054-1055`
- Test: `tests/testthat/test-ddsimca3w.R` (append) + existing `test-ddsimcaplots.R`/`test-ddsimca.R` regression

**Interfaces:**
- Consumes: `f0.model` field on each outcome (Task 4).
- Produces: result plot methods that normalize/scale by `f0.model` (falling back to `Nf` when absent, so old objects and the PCA path are unchanged); `as.matrix`/bounds keyed off `ncol(x$T2)`.

- [ ] **Step 1: Write the failing test (T2-vs-scores width)**

Append to `tests/testthat/test-ddsimca3w.R`:

```r
test_that("as.matrix.ddsimcares keys component count off T2, not scores", {
   data(people)
   x <- people[, -4]
   m <- ddsimca(x, "people", ncomp = 4, center = TRUE, scale = TRUE, pcv = NULL)
   # simulate the 3-way shape: scores narrower than T2 (selected < total)
   r <- m$res$cal
   r$scores <- r$scores[, 1:2, drop = FALSE]   # pretend only 2 selected
   mm <- as.matrix(r)
   # rows must span all 4 sub-models (ncol of T2), not the 2 score columns
   expect_equal(nrow(mm), 4L)
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: FAIL — `nrow(mm)` is 2 (keyed off `ncol(scores)`).

- [ ] **Step 3: Edit `as.matrix.ddsimcares` (R/ddsimcares.R)**

Line 188: `ncomp <- ncol(x$scores)` → `ncomp <- ncol(x$T2)`
Line 222: `rownames(out) <- colnames(x$scores)` → `rownames(out) <- colnames(x$T2)`

Line 133 in `as.data.frame.ddsimcares`: `if (ncomp < 1 || ncomp > ncol(x$scores)) {` → `if (ncomp < 1 || ncomp > ncol(x$T2)) {`

Line 1063 in `plotDistances.ddsimcares`: `if (ncomp < 1 || ncomp > ncol(obj$scores)) {` → `if (ncomp < 1 || ncomp > ncol(obj$T2)) {`

In `plotAcceptance.ddsimcares` (R/ddsimcares.R:873) find the analogous `ncol(obj$scores)` bound check and change it to `ncol(obj$T2)` (if present; if the method uses `obj$ncomp.selected` only, no change needed).

- [ ] **Step 4: Edit the `f0.model` scaling in the three plot methods**

In `plotExtremes.ddsimcares` (R/ddsimcares.R) replace line 712:

```r
   f0m <- if (is.null(v$f0.model)) v$Nf else v$f0.model
   p <- 1 - pchisq(f * v$Nf / f0m, v$Nf)
```

In `plotDistances.ddsimcares` replace lines 1089-1090 and 1096:

```r
   f0m <- if (is.null(v$f0.model)) v$Nf else v$f0.model
   yle <- v$fce / f0m
   ylo <- v$fco / f0m
```
and (line ~1096, inside the `else` for `distance == "f"`):
```r
      y <- y / f0m
```

In `plotSelectivityArea.ddsimcares` replace lines 467-468:

```r
   f0m <- if (is.null(v$f0.model)) v$Nf else v$f0.model
   fcrit <- f0m / v$Nf * qchisq(1 - alpha, v$Nf)
   z <- fcrit / v$f0
```

In `plotDistances.ddsimca` (R/ddsimca.R) replace lines 1054-1055:

```r
   f0m <- if (is.null(v$f0.model)) v$Nf else v$f0.model
   yle <- v$fce / f0m
   ylo <- v$fco / f0m
```

> All five edits reduce to the current code when `f0.model == Nf` (the PCA path), so existing plots are unchanged. The scaling only differs for the 3-way data-driven F where `f0.model != Nf`; these diagnostic curves are not part of the golden set, so the choice to keep them self-consistent with the scaled decision boundary is intentional (it is more correct than the Python oracle, which leaves these particular curves unscaled).

- [ ] **Step 5: Run the new test + regression suites**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: PASS.

Run: `Rscript -e 'devtools::test(filter = "ddsimca")'`
Expected: PASS — existing `test-ddsimca.R` and `test-ddsimcaplots.R` unchanged.

- [ ] **Step 6: Commit**

```bash
git add R/ddsimcares.R R/ddsimca.R tests/testthat/test-ddsimca3w.R
git commit -m "Thread f0.model through DD-SIMCA result plots and key as.matrix off T2"
```

---

## Task 6: Convert DoF plots to S3 generics

**Files:**
- Modify: `R/pca.R:1897-1982` (the three DoF functions)
- Modify: `R/defaults.R` (no new generic needed here — see note)
- Test: existing `tests/testthat/test-pcaplots.R`, `test-simcaplots.R`, `test-ddsimca.R` (DoF calls); add a dispatch test to `test-ddsimca3w.R`

**Interfaces:**
- Produces: `plotT2DoF`, `plotQDoF`, `plotDistDoF` as S3 generics `function(obj, ...)` with the current bodies preserved as `*.pca` methods. `.ddsimca3w` methods are added in Task 10.

- [ ] **Step 1: Write the dispatch test**

Append to `tests/testthat/test-ddsimca3w.R`:

```r
test_that("DoF plot generics still dispatch for pca/ddsimca objects", {
   data(people)
   x <- people[, -4]
   m <- pca(x, ncomp = 3, center = TRUE, scale = TRUE, lim.type = "ddmoments")
   pdf(tempfile()); on.exit(dev.off())
   # plotT2DoF.pca draws (no show.plot arg); assert it dispatches without error
   expect_error(plotT2DoF(m), NA)
   expect_true(is.function(plotT2DoF))
})
```

- [ ] **Step 2: Run to verify current behavior (plain function, no dispatch)**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: PASS currently (it is a plain function) — this test mainly guards Step 3 from breaking dispatch. Proceed.

- [ ] **Step 3: Convert the three functions to generics + `.pca` methods in `R/pca.R`**

For `plotT2DoF` (line 1897), replace the single function with a generic stub and a `.pca` method (keep the existing roxygen block above; add `@export` to both, and change `@param obj`/`@rdname` as needed):

```r
#' @export
plotT2DoF <- function(obj, ...) {
   UseMethod("plotT2DoF")
}

#' @export
plotT2DoF.pca <- function(obj, type = "b", labels = "values", xticks = seq_len(obj$ncomp), ylab = "Nh", ...) {

   if (!(obj$lim.type %in% c("ddrobust", "ddmoments", "chisq"))) {
      stop("This plot can not be made for selected 'lim.type' method.", call. = FALSE)
   }

   plot_data <- mda.subset(obj$T2lim, subset = 4)
   attr(plot_data, "name") <- "Degrees of freedom"
   attr(plot_data, "xaxis.name") <- attr(obj$loadings, "xaxis.name")
   mdaplot(plot_data, xticks = xticks, labels = labels, type = type, ylab = ylab, ...)
}
```

Do the same for `plotQDoF` (line 1933 → generic + `plotQDoF.pca`, body unchanged, `ylab = "Nq"`) and `plotDistDoF` (line 1967 → generic + `plotDistDoF.pca`, body unchanged). The generic stubs are `function(obj, ...) UseMethod("...")`; all named formals (`type`, `labels`, `xticks`, `ylab`) move into the `.pca` methods only.

- [ ] **Step 4: Regenerate NAMESPACE and run regression**

Run: `Rscript -e 'roxygen2::roxygenize()'`

Confirm `NAMESPACE` now has `export(plotT2DoF)` + `S3method(plotT2DoF,pca)` (and the same for QDoF/DistDoF). The bare `export(plotT2DoF)` line is preserved (now exporting the generic).

Run: `Rscript -e 'devtools::test(filter = "pcaplots")'` and `Rscript -e 'devtools::test(filter = "simcaplots")'` and `Rscript -e 'devtools::test(filter = "ddsimca")'`
Expected: PASS — all existing DoF-plot calls dispatch to `.pca` unchanged.

- [ ] **Step 5: Commit**

```bash
git add R/pca.R NAMESPACE man/plotT2DoF.Rd man/plotQDoF.Rd man/plotDistDoF.Rd tests/testthat/test-ddsimca3w.R
git commit -m "Convert plotT2DoF/plotQDoF/plotDistDoF to S3 generics with pca methods"
```

---

## Task 7: Golden fixtures from the Python oracle

**Files:**
- Create: `../../ddsimca-py/tests/dump_3way_golden.py`
- Create: `tests/testthat/ddsimca3w/` (committed CSV outputs)

**Interfaces:**
- Produces, per variant (`parafac`, `tucker`) and limit type (`classic`, `robust`), CSV files the R golden test (Tasks 8/9) reads: the unfolded training matrix `X_train.csv`, per-sub-model factor matrices `B_<nc>.csv`, `C_<nc>.csv`, and (Tucker) `G_<nc>.csv` (flattened mode-1, row-major), plus the oracle result matrices `H.csv`, `Q.csv`, `F.csv` and the limit params `hParams_<lim>.csv`, `qParams_<lim>.csv`, `fParams_<lim>.csv` (each two rows: `u0`, `Nu`).

> Why factors-as-fixtures: R's ALS cannot reproduce tensorly's RNG, so we validate the DD-SIMCA delta (project → H/Q/F → limits) by feeding the SAME factor matrices into R and comparing to the oracle at 1e-9. R's own ALS is validated separately by invariants (Tasks 2/3).

- [ ] **Step 1: Write the dump script**

Create `../../ddsimca-py/tests/dump_3way_golden.py`:

```python
"""Dump per-sub-model factors + DD-SIMCA result matrices as CSV for the R
golden tests. Run inside the ddsimca-py venv:
    cd ../../ddsimca-py && .venv/bin/python tests/dump_3way_golden.py <out_dir>
"""
import sys, os
# make the sibling golden_parafac_setup.py importable when run from the repo root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from golden_parafac_setup import build_train_df, DIM, NCOMP
from ddsimca.ddsimca_parafac import ddsimca_parafac
from ddsimca.ddsimca_tucker import ddsimca_tucker


def _save(path, arr):
    np.savetxt(path, np.atleast_2d(arr), delimiter=",")


def dump_parafac(out):
    df = build_train_df("T")
    X = df.iloc[:, 1:].values.astype(float)
    m = ddsimca_parafac(df, dim=DIM, ncomp=NCOMP, nf_as_sum=False)
    _save(f"{out}/parafac_X_train.csv", X)
    for nc in range(NCOMP):
        _save(f"{out}/parafac_B_{nc + 1}.csv", m.models[nc]["B"])
        _save(f"{out}/parafac_C_{nc + 1}.csv", m.models[nc]["C"])
    H, Q, T, U, E = m._project(X)
    _save(f"{out}/parafac_H.csv", H)
    _save(f"{out}/parafac_Q.csv", Q)
    for lim in ("classic", "robust"):
        h0, Nh = m.hParams[lim]; q0, Nq = m.qParams[lim]; f0, Nf = m.fParams[lim]
        _save(f"{out}/parafac_hParams_{lim}.csv", np.vstack([h0, Nh]))
        _save(f"{out}/parafac_qParams_{lim}.csv", np.vstack([q0, Nq]))
        _save(f"{out}/parafac_fParams_{lim}.csv", np.vstack([f0, Nf]))


def dump_tucker(out):
    df = build_train_df("T")
    X = df.iloc[:, 1:].values.astype(float)
    rB, rC = DIM[0], DIM[1]
    m = ddsimca_tucker(df, dim=DIM, ncomp=(NCOMP, rB, rC), nf_as_sum=False)
    _save(f"{out}/tucker_X_train.csv", X)
    for nc in range(NCOMP):
        _save(f"{out}/tucker_B_{nc + 1}.csv", m.models[nc]["B"])
        _save(f"{out}/tucker_C_{nc + 1}.csv", m.models[nc]["C"])
        G = m.models[nc]["core"]
        _save(f"{out}/tucker_G_{nc + 1}.csv", G.reshape(G.shape[0], -1))  # row-major mode-1
    H, Q, T, U, E = m._project(X)
    _save(f"{out}/tucker_H.csv", H)
    _save(f"{out}/tucker_Q.csv", Q)
    for lim in ("classic", "robust"):
        h0, Nh = m.hParams[lim]; q0, Nq = m.qParams[lim]; f0, Nf = m.fParams[lim]
        _save(f"{out}/tucker_hParams_{lim}.csv", np.vstack([h0, Nh]))
        _save(f"{out}/tucker_qParams_{lim}.csv", np.vstack([q0, Nq]))
        _save(f"{out}/tucker_fParams_{lim}.csv", np.vstack([f0, Nf]))


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    dump_parafac(out)
    dump_tucker(out)
    print(f"wrote golden CSVs to {out}")
```

- [ ] **Step 2: Run the dump script into the fixtures directory**

Run:
```bash
cd ../../ddsimca-py && .venv/bin/python tests/dump_3way_golden.py \
  ../mdatools-r/mdatools/tests/testthat/ddsimca3w
```
Expected: `wrote golden CSVs to ...`. Verify the directory now contains `parafac_X_train.csv`, `parafac_B_1.csv` … `parafac_fParams_robust.csv`, and the `tucker_*` equivalents.

> If `ddsimca_parafac` errors because `m.models[nc]["B"]` shape differs from expectation, read the actual key names in `../../ddsimca-py/ddsimca/ddsimca_parafac.py:178-184` and adjust. Do NOT guess.

- [ ] **Step 3: Commit the fixtures and the script**

```bash
git add tests/testthat/ddsimca3w/
cd ../../ddsimca-py && git add tests/dump_3way_golden.py && git commit -m "Add 3-way DD-SIMCA golden dump script"
cd ../mdatools-r/mdatools && git commit -m "Add 3-way DD-SIMCA golden CSV fixtures"
```

> The `../../ddsimca-py` commit is in a separate repo; only commit it there if that repo is under version control and the user wants it tracked. Otherwise leave the script uncommitted and just keep the fixtures in mdatools.

---

## Task 8: PARAFAC model — `project3w`, `finalize3w`, constructor, predict

**Files:**
- Modify: `R/ddsimca3w.R` (append the model layer)
- Test: `tests/testthat/test-ddsimca3w.R` (append golden + end-to-end)

**Interfaces:**
- Consumes: decompositions (Task 2), projection helpers (Task 1), `classify` (Task 4), `ddsimcares` (`R/ddsimcares.R:93`), `ldecomp.getLimParams`/`ddmoments.param`/`ddrobust.param`/`clamp.dof` (`R/ldecomp.R`), golden CSVs (Task 7).
- Produces:
  - `project3w(model, x)` → list(`H`, `Q` = `n x ncomp`; `T`, `U` = `n x ncomp.selected`; `E` = `n x (J*K)`). Dispatches PARAFAC vs Tucker on `model$type`. Reads stored `model$models`, `model$non.negative`, `model$dim`, `model$ncomp`, `model$ncomp.selected`.
  - `finalize3w(H, Q, nf.as.sum)` → list(`T2`, `Q`, `F`) each `list(moments=list(u0,Nu), robust=list(u0,Nu))`, with F's `Nu` clamped via `clamp.dof`.
  - `ddsimca.parafac(x, classname, dim, ncomp, non.negative, nf.as.sum, alpha, gamma)` → model of class `c("ddsimca.parafac","ddsimca3w","ddsimca")`.
  - `predict.ddsimca3w(object, x, c.ref, alpha, gamma, ...)` → `ddsimcares`.

- [ ] **Step 1: Write the golden test (1e-9) for the PARAFAC projection + limits**

Append to `tests/testthat/test-ddsimca3w.R`:

```r
context("ddsimca3w: PARAFAC golden")

gpath <- function(f) test_path("ddsimca3w", f)
rdcsv <- function(f) as.matrix(read.csv(gpath(f), header = FALSE))

test_that("PARAFAC project3w + finalize3w match the Python oracle (1e-9) on fixed factors", {
   skip_if_not(file.exists(gpath("parafac_X_train.csv")))
   X <- rdcsv("parafac_X_train.csv")
   ncomp <- 3; J <- 5; K <- 4
   models <- lapply(seq_len(ncomp), function(nc) {
      B <- rdcsv(sprintf("parafac_B_%d.csv", nc))
      C <- rdcsv(sprintf("parafac_C_%d.csv", nc))
      pr <- projectParafac(B, C, X)
      A <- pr$A
      A_mean <- colMeans(A)
      A_c <- sweep(A, 2, A_mean)
      S <- crossprod(A_c) / max(nrow(A_c) - 1, 1)
      list(B = B, C = C, A_mean = A_mean, S_pinv = pinv3way(S),
           A_V = svd(A_c)$v, A_sigma = svd(A_c)$d)
   })
   model <- list(type = "parafac", models = models, non.negative = FALSE,
      dim = c(J, K), ncomp = ncomp, ncomp.selected = ncomp)
   pj <- project3w(model, X)
   expect_equal(unname(pj$H), unname(rdcsv("parafac_H.csv")), tolerance = 1e-9)
   expect_equal(unname(pj$Q), unname(rdcsv("parafac_Q.csv")), tolerance = 1e-9)

   lp <- finalize3w(pj$H, pj$Q, nf.as.sum = FALSE)
   pf_c <- rdcsv("parafac_fParams_classic.csv")
   expect_equal(unname(lp$F$moments$u0), unname(pf_c[1, ]), tolerance = 1e-9)
   expect_equal(unname(lp$F$moments$Nu), unname(pf_c[2, ]), tolerance = 1e-9)
   ph_c <- rdcsv("parafac_hParams_classic.csv")
   expect_equal(unname(lp$T2$moments$u0), unname(ph_c[1, ]), tolerance = 1e-9)
   expect_equal(unname(lp$T2$moments$Nu), unname(ph_c[2, ]), tolerance = 1e-9)
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: FAIL — `project3w`/`finalize3w` not found.

- [ ] **Step 3: Implement `project3w` and `finalize3w`**

Append to `R/ddsimca3w.R`:

```r
# Project x (n x (J*K)) through every stored sub-model. H/Q are n x ncomp
# (column nc = distances from the rank-nc sub-model; sub-models are NOT nested).
# T/U are n x ncomp.selected from the selected sub-model only; E is the residual
# of the final sub-model.
project3w <- function(model, x) {
   n <- nrow(x)
   ncomp <- model$ncomp
   H <- matrix(0, n, ncomp)
   Q <- matrix(0, n, ncomp)
   E <- matrix(0, n, ncol(x))
   Tsel <- NULL; Usel <- NULL

   for (nc in seq_len(ncomp)) {
      mdl <- model$models[[nc]]
      pr <- if (model$type == "parafac") {
         projectParafac(mdl$B, mdl$C, x)
      } else {
         projectTucker(mdl$G, mdl$B, mdl$C, x, model$non.negative)
      }
      A <- pr$A
      A_c <- sweep(A, 2, mdl$A_mean)
      H[, nc] <- rowSums((A_c %*% mdl$S_pinv) * A_c)
      Q[, nc] <- rowSums(pr$E^2)
      if (nc == ncomp) E <- pr$E

      if (nc == model$ncomp.selected) {
         Tt <- A_c %*% mdl$A_V
         sigma <- mdl$A_sigma
         Usel <- sweep(Tt, 2, ifelse(sigma > 1e-12, sigma, 1), "/")
         Usel[, sigma <= 1e-12] <- 0
         Tsel <- Tt
      }
   }
   list(H = H, Q = Q, T = Tsel, U = Usel, E = E)
}

# Build the limParams structure classify() consumes. T2 and Q come from
# ldecomp.getLimParams (clamps DoF). F is the scaled full distance: sum rule
# (nf.as.sum = TRUE) or data-driven (paper eq.7) with explicit clamp.dof on Nu.
finalize3w <- function(H, Q, nf.as.sum) {
   T2 <- ldecomp.getLimParams(H)
   Qp <- ldecomp.getLimParams(Q)

   if (nf.as.sum) {
      F <- list(
         moments = list(u0 = T2$moments$Nu + Qp$moments$Nu, Nu = T2$moments$Nu + Qp$moments$Nu),
         robust  = list(u0 = T2$robust$Nu  + Qp$robust$Nu,  Nu = T2$robust$Nu  + Qp$robust$Nu)
      )
   } else {
      Fc <- sweep(H, 2, T2$moments$u0, "/") %*% diag(T2$moments$Nu, ncol(H)) +
            sweep(Q, 2, Qp$moments$u0, "/") %*% diag(Qp$moments$Nu, ncol(Q))
      Fr <- sweep(H, 2, T2$robust$u0,  "/") %*% diag(T2$robust$Nu,  ncol(H)) +
            sweep(Q, 2, Qp$robust$u0,  "/") %*% diag(Qp$robust$Nu,  ncol(Q))
      mc <- ddmoments.param(Fc); rr <- ddrobust.param(Fr)
      # ddmoments.param/ddrobust.param do NOT clamp; only ldecomp.getLimParams
      # does. Clamp F's DoF to the package-wide integer [1,250] range here.
      mc$Nu <- clamp.dof(mc$Nu); rr$Nu <- clamp.dof(rr$Nu)
      F <- list(moments = mc, robust = rr)
   }
   list(T2 = T2, Q = Qp, F = F)
}
```

> The `sweep(..., "/") %*% diag(Nu)` form computes `(H/h0)*Nh + (Q/q0)*Nq` column-wise. For a single component `diag(Nu, 1)` must be forced to 1x1 — `diag(x, ncol)` with explicit `ncol` handles the scalar case (base `diag(scalar)` would build an identity of that size, a classic R trap). Verify with the golden test before trusting it.

- [ ] **Step 4: Run the golden test to verify projection + limits pass**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: PASS for the golden test. If `H`/`Q` mismatch, the projection or Mahalanobis is wrong; if only `F$...$Nu` mismatches, the clamp/sweep is wrong.

- [ ] **Step 5: Implement the constructor + predict + the result builder**

Append to `R/ddsimca3w.R`:

```r
#' DD-SIMCA with PARAFAC decomposition for 3-way data
#'
#' @param x numeric matrix (n x (J*K)), tensor unfolded row-major (column for
#'   element (j,k) is (j-1)*K + k).
#' @param classname name of the target class (<= 20 chars).
#' @param dim integer vector c(J, K), both >= 2.
#' @param ncomp number of PARAFAC components (<= min(J, K, n-1)).
#' @param non.negative logical, use non-negative PARAFAC.
#' @param nf.as.sum logical, full-distance DoF rule (FALSE = data-driven eq.7).
#' @param alpha significance level for extremes.
#' @param gamma significance level for outliers.
#'
#' @return model of class \code{c("ddsimca.parafac","ddsimca3w","ddsimca")}.
#' @export
ddsimca.parafac <- function(x, classname, dim,
   ncomp = min(dim[1], dim[2], nrow(x) - 1),
   non.negative = FALSE, nf.as.sum = FALSE, alpha = 0.05, gamma = 0.01) {

   x <- as.matrix(x)
   J <- dim[1]; K <- dim[2]
   if (J < 2 || K < 2 || ncol(x) != J * K) {
      stop("Argument 'dim' does not match the number of columns in 'x'.", call. = FALSE)
   }
   if (ncomp < 1 || ncomp > J || ncomp > K || ncomp > nrow(x) - 1) {
      stop("Wrong value for 'ncomp'.", call. = FALSE)
   }
   if (!is.character(classname) || nchar(classname) > 20) {
      stop("Argument 'classname' must be text of up to 20 symbols.", call. = FALSE)
   }

   X3 <- x; attr(X3, "dim3") <- c(J, K)
   models <- vector("list", ncomp)
   for (nc in seq_len(ncomp)) {
      fit <- if (non.negative) parafac.als.nn(X3, nc) else parafac.als(X3, nc)
      models[[nc]] <- buildSubModel3w(fit$B, fit$C, NULL, x, non.negative, "parafac")
   }

   model <- makeModel3w(x, models, "parafac", classname, c(J, K), ncomp,
      non.negative, nf.as.sum, alpha, gamma)
   model$call <- match.call()
   class(model) <- c("ddsimca.parafac", "ddsimca3w", "ddsimca")
   model$res[["cal"]] <- predict(model, x, rep(classname, nrow(x)))
   model$calres <- model$res[["cal"]]
   model
}

# Build one sub-model: project training data to get scores, then the H-stats
# (A_mean, S_pinv) and the SVD basis (A_V, A_sigma) for scores plots.
buildSubModel3w <- function(B, C, G, x, non.negative, type) {
   pr <- if (type == "parafac") projectParafac(B, C, x) else projectTucker(G, B, C, x, non.negative)
   A <- pr$A
   A_mean <- colMeans(A)
   A_c <- sweep(A, 2, A_mean)
   S <- crossprod(A_c) / max(nrow(A_c) - 1, 1)
   sv <- svd(A_c)
   list(B = B, C = C, G = G, A_mean = A_mean, S_pinv = pinv3way(S),
        A_V = sv$v, A_sigma = sv$d)
}

# Assemble the shared model state (limParams + metadata). Type-agnostic.
makeModel3w <- function(x, models, type, classname, dim3, ncomp, non.negative,
   nf.as.sum, alpha, gamma) {

   base <- list(type = type, models = models, non.negative = non.negative,
      dim = dim3, ncomp = ncomp, ncomp.selected = ncomp)
   pj <- project3w(base, x)
   lp <- finalize3w(pj$H, pj$Q, nf.as.sum)

   model <- base
   model$limParams <- lp
   model$tssX <- sum(x^2)        # total SS for cumexpvar (buildPcaresShape3w)
   model$Xcal <- x               # training data, reused by selectCompNum/setParams
   model$nrows <- nrow(x)
   model$nclasses <- 1
   model$classname <- classname
   model$alpha <- alpha
   model$gamma <- gamma
   model$nf.as.sum <- nf.as.sum
   model$limType <- "moments"
   model
}

#' Predictions for a 3-way DD-SIMCA model
#' @param object a 3-way DD-SIMCA model (class \code{ddsimca3w}).
#' @param x numeric matrix (n x (J*K)), unfolded row-major.
#' @param c.ref optional vector of reference class names.
#' @param alpha significance level for predictions.
#' @param gamma significance level for outliers.
#' @param ... ignored.
#' @return object of class \code{ddsimcares}.
#' @export
predict.ddsimca3w <- function(object, x, c.ref = NULL,
   alpha = object$alpha, gamma = object$gamma, ...) {

   x <- as.matrix(x)
   pj <- project3w(object, x)
   n <- nrow(x)

   # build the pcares-shaped object the result methods read
   res <- buildPcaresShape3w(object, pj)

   # indices/numbers (no excluded rows for 3-way)
   if (is.null(c.ref)) {
      indMembers <- NULL; indStrangers <- NULL; indUnknowns <- rep(TRUE, n)
      nMembers <- 0; nStrangers <- 0; nUnknowns <- n
   } else {
      indMembers <- c.ref == object$classname
      indStrangers <- !indMembers
      indUnknowns <- NULL
      nMembers <- sum(indMembers); nStrangers <- sum(indStrangers); nUnknowns <- 0
   }
   indices <- list(members = indMembers, strangers = indStrangers,
      unknown = indUnknowns, excluded = rep(FALSE, n))
   numbers <- list(members = nMembers, strangers = nStrangers,
      unknown = nUnknowns, excluded = 0)

   lp <- object$limParams
   outcomes <- list(
      moments = classify(res, indices, numbers, lp$Q$moments, lp$T2$moments,
         object$nrows, alpha, gamma, object$classname, c.ref, fp = lp$F$moments),
      robust = classify(res, indices, numbers, lp$Q$robust, lp$T2$robust,
         object$nrows, alpha, gamma, object$classname, c.ref, fp = lp$F$robust)
   )
   ddsimcares(res, outcomes, indices, numbers, alpha = alpha,
      classname = object$classname, c.ref = c.ref)
}

# Build the pcares-shaped result object (scores/T2/Q/residuals/expvar + dimnames)
# that classify() and the ddsimcares result methods consume.
buildPcaresShape3w <- function(object, pj) {
   ncomp <- object$ncomp
   compNames <- paste("Comp", seq_len(ncomp))
   H <- pj$H; Q <- pj$Q
   colnames(H) <- compNames; colnames(Q) <- compNames
   sel <- object$ncomp.selected
   scores <- pj$T
   if (!is.null(scores)) colnames(scores) <- paste("Comp", seq_len(ncol(scores)))

   # cumulative explained variance from per-sub-model residual SS (object$tssX
   # = sum(x^2) is stored on the model by makeModel3w). Non-nested sub-models can
   # make cumexpvar non-monotonic; that only affects the cosmetic Expvar columns.
   cumexpvar <- 100 * (1 - colSums(Q) / object$tssX)
   expvar <- c(cumexpvar[1], diff(cumexpvar))

   res <- list(
      scores = scores, T2 = H, Q = Q, residuals = pj$E,
      ncomp = ncomp, ncomp.selected = sel,
      expvar = expvar, cumexpvar = cumexpvar,
      U = pj$U
   )
   class(res) <- c("pcares", "ldecomp")
   res
}
```

- [ ] **Step 6: Write the end-to-end PARAFAC test**

Append to `tests/testthat/test-ddsimca3w.R`:

```r
test_that("ddsimca.parafac trains and predicts end-to-end", {
   set.seed(20)
   J <- 5; K <- 4; n <- 30; R <- 3
   A0 <- matrix(rnorm(n * R), n, R)
   B0 <- matrix(rnorm(J * R), J, R)
   C0 <- matrix(rnorm(K * R), K, R)
   x <- A0 %*% t(khatriRao3w(B0, C0)) + 0.01 * matrix(rnorm(n * J * K), n, J * K)
   m <- ddsimca.parafac(x, "T", dim = c(J, K), ncomp = 3)
   expect_s3_class(m, "ddsimca.parafac")
   expect_s3_class(m, "ddsimca3w")
   expect_s3_class(m$res$cal, "ddsimcares")
   expect_equal(ncol(m$res$cal$T2), 3L)
   # full-distance critical limit uses scaled F (f0.model != Nf in general)
   v <- m$res$cal$simca$outcomes$moments$values[[3]]
   expect_true(!is.null(v$f0.model))
   # predict on new data returns a ddsimcares
   r <- predict(m, x, rep("T", n))
   expect_s3_class(r, "ddsimcares")
})
```

- [ ] **Step 7: Run the full 3-way suite**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: PASS (golden + end-to-end).

- [ ] **Step 8: Commit**

```bash
git add R/ddsimca3w.R tests/testthat/test-ddsimca3w.R
git commit -m "Add ddsimca.parafac constructor, project3w/finalize3w, and predict"
```

---

## Task 9: Tucker model

**Files:**
- Modify: `R/ddsimca3w.R` (append the Tucker constructor)
- Test: `tests/testthat/test-ddsimca3w.R` (append golden + end-to-end)

**Interfaces:**
- Consumes: `tucker.hooi`/`tucker.nn`/`signFix3w` (Task 3), `project3w`/`finalize3w`/`makeModel3w`/`buildSubModel3w`/`predict.ddsimca3w` (Task 8).
- Produces: `ddsimca.tucker(x, classname, dim, ncomp, non.negative, nf.as.sum, alpha, gamma)` → model of class `c("ddsimca.tucker","ddsimca3w","ddsimca")`. `ncomp = c(R_A_max, R_B, R_C)`.

- [ ] **Step 1: Write the Tucker golden test**

Append to `tests/testthat/test-ddsimca3w.R`:

```r
context("ddsimca3w: Tucker golden")

test_that("Tucker project3w + finalize3w match the Python oracle (1e-9) on fixed factors", {
   skip_if_not(file.exists(gpath("tucker_X_train.csv")))
   X <- rdcsv("tucker_X_train.csv")
   ncomp <- 3; J <- 5; K <- 4; rB <- 5; rC <- 4
   models <- lapply(seq_len(ncomp), function(nc) {
      B <- rdcsv(sprintf("tucker_B_%d.csv", nc))
      C <- rdcsv(sprintf("tucker_C_%d.csv", nc))
      Gflat <- rdcsv(sprintf("tucker_G_%d.csv", nc))   # nc x (rB*rC), row-major
      G <- fold1_3w(Gflat, nc, rB, rC)
      buildSubModel3w(B, C, G, X, FALSE, "tucker")
   })
   model <- list(type = "tucker", models = models, non.negative = FALSE,
      dim = c(J, K), ncomp = ncomp, ncomp.selected = ncomp)
   pj <- project3w(model, X)
   expect_equal(unname(pj$H), unname(rdcsv("tucker_H.csv")), tolerance = 1e-9)
   expect_equal(unname(pj$Q), unname(rdcsv("tucker_Q.csv")), tolerance = 1e-9)

   lp <- finalize3w(pj$H, pj$Q, nf.as.sum = FALSE)
   tf_c <- rdcsv("tucker_fParams_classic.csv")
   expect_equal(unname(lp$F$moments$Nu), unname(tf_c[2, ]), tolerance = 1e-9)
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: FAIL only on the Tucker golden if the constructor isn't there yet — but `project3w`/`buildSubModel3w` already handle `type="tucker"` (Task 8), so this test may PASS already for the projection. If it passes, good — it validates the Tucker projection path. Proceed to add the constructor.

- [ ] **Step 3: Implement `ddsimca.tucker`**

Append to `R/ddsimca3w.R`:

```r
#' DD-SIMCA with Tucker3 decomposition for 3-way data
#'
#' @param x numeric matrix (n x (J*K)), unfolded row-major.
#' @param classname target class name (<= 20 chars).
#' @param dim integer vector c(J, K), both >= 2.
#' @param ncomp integer vector c(R_A_max, R_B, R_C): sample-mode max rank and the
#'   fixed mode-2 / mode-3 ranks. Requires R_B <= J, R_C <= K, R_A_max <= n-1,
#'   R_A_max <= R_B*R_C.
#' @param non.negative logical, use non-negative Tucker.
#' @param nf.as.sum logical, full-distance DoF rule (FALSE = data-driven eq.7).
#' @param alpha significance level for extremes.
#' @param gamma significance level for outliers.
#'
#' @return model of class \code{c("ddsimca.tucker","ddsimca3w","ddsimca")}.
#' @export
ddsimca.tucker <- function(x, classname, dim, ncomp,
   non.negative = FALSE, nf.as.sum = FALSE, alpha = 0.05, gamma = 0.01) {

   x <- as.matrix(x)
   J <- dim[1]; K <- dim[2]
   rA <- ncomp[1]; rB <- ncomp[2]; rC <- ncomp[3]
   if (J < 2 || K < 2 || ncol(x) != J * K) {
      stop("Argument 'dim' does not match the number of columns in 'x'.", call. = FALSE)
   }
   if (rA > nrow(x) - 1) stop("R_A_max cannot exceed n - 1.", call. = FALSE)
   if (rB > J) stop("R_B cannot exceed J.", call. = FALSE)
   if (rC > K) stop("R_C cannot exceed K.", call. = FALSE)
   if (rA > rB * rC) stop("R_A_max cannot exceed R_B * R_C.", call. = FALSE)
   if (!is.character(classname) || nchar(classname) > 20) {
      stop("Argument 'classname' must be text of up to 20 symbols.", call. = FALSE)
   }

   X3 <- x; attr(X3, "dim3") <- c(J, K)
   models <- vector("list", rA)
   for (nc in seq_len(rA)) {
      fit <- if (non.negative) tucker.nn(X3, c(nc, rB, rC)) else tucker.hooi(X3, c(nc, rB, rC))
      sf <- signFix3w(fit$A, fit$B, fit$C, fit$G)
      models[[nc]] <- buildSubModel3w(sf$B, sf$C, sf$G, x, non.negative, "tucker")
   }

   model <- makeModel3w(x, models, "tucker", classname, c(J, K), rA,
      non.negative, nf.as.sum, alpha, gamma)
   model$B.ncomp <- rB
   model$C.ncomp <- rC
   model$call <- match.call()
   class(model) <- c("ddsimca.tucker", "ddsimca3w", "ddsimca")
   model$res[["cal"]] <- predict(model, x, rep(classname, nrow(x)))
   model$calres <- model$res[["cal"]]
   model
}
```

- [ ] **Step 4: Write the end-to-end Tucker test**

Append to `tests/testthat/test-ddsimca3w.R`:

```r
test_that("ddsimca.tucker trains and predicts end-to-end", {
   set.seed(21)
   J <- 5; K <- 4; n <- 30; R <- 3
   A0 <- matrix(rnorm(n * R), n, R)
   B0 <- matrix(rnorm(J * R), J, R)
   C0 <- matrix(rnorm(K * R), K, R)
   x <- A0 %*% t(khatriRao3w(B0, C0)) + 0.01 * matrix(rnorm(n * J * K), n, J * K)
   m <- ddsimca.tucker(x, "T", dim = c(J, K), ncomp = c(3, 5, 4))
   expect_s3_class(m, "ddsimca.tucker")
   expect_s3_class(m$res$cal, "ddsimcares")
   expect_equal(ncol(m$res$cal$T2), 3L)
   expect_equal(m$B.ncomp, 5)
})
```

- [ ] **Step 5: Run the suite**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: PASS (Tucker golden + end-to-end).

- [ ] **Step 6: Commit**

```bash
git add R/ddsimca3w.R tests/testthat/test-ddsimca3w.R
git commit -m "Add ddsimca.tucker constructor"
```

---

## Task 10: S3 glue — selectCompNum, summary/print, plotFactors, adapters, DoF methods

**Files:**
- Modify: `R/ddsimca3w.R` (append) and `R/defaults.R` (add `plotFactors` generic)
- Test: `tests/testthat/test-ddsimca3w.R` (append)

**Interfaces:**
- Consumes: model + result objects (Tasks 8/9), `plotScores.ldecomp` (`R/ldecomp.R:145`), the DoF generics (Task 6).
- Produces:
  - `selectCompNum.ddsimca3w(obj, ncomp)` → model with updated `ncomp.selected` (model + each `$res`), `$scores` re-pointed to the selected sub-model.
  - `summary.ddsimca3w`, `print.ddsimca3w`.
  - `plotFactors` generic (in `defaults.R`) + `plotFactors.ddsimca3w(obj, comp, mode, ncomp, ...)`.
  - `plotScores.ddsimca3w`, `plotLoadings.ddsimca3w` adapters.
  - `plotT2DoF.ddsimca3w`, `plotQDoF.ddsimca3w`, `plotDistDoF.ddsimca3w` reading `limParams`.

- [ ] **Step 1: Write the tests**

Append to `tests/testthat/test-ddsimca3w.R`:

```r
context("ddsimca3w: S3 glue")

make_pf <- function() {
   set.seed(22)
   J <- 5; K <- 4; n <- 30; R <- 3
   A0 <- matrix(rnorm(n * R), n, R); B0 <- matrix(rnorm(J * R), J, R)
   C0 <- matrix(rnorm(K * R), K, R)
   x <- A0 %*% t(khatriRao3w(B0, C0)) + 0.01 * matrix(rnorm(n * J * K), n, J * K)
   ddsimca.parafac(x, "T", dim = c(J, K), ncomp = 3)
}

test_that("selectCompNum.ddsimca3w updates model and results", {
   m <- make_pf()
   m <- selectCompNum(m, 2)
   expect_equal(m$ncomp.selected, 2)
   expect_equal(m$res$cal$ncomp.selected, 2)
})

test_that("summary and print do not error", {
   m <- make_pf()
   expect_output(summary(m))
   expect_output(print(m))
})

test_that("plotFactors returns factor values for a PARAFAC model", {
   m <- make_pf()
   pdf(tempfile()); on.exit(dev.off())
   expect_error(plotFactors(m, comp = 1, mode = "B"), NA)
})

test_that("DoF methods read Nh/Nq/Nf from limParams", {
   m <- make_pf()
   pdf(tempfile()); on.exit(dev.off())
   expect_error(plotT2DoF(m), NA)
   expect_error(plotDistDoF(m), NA)
})

test_that("a PCA-only method has no applicable method on the 3-way model", {
   m <- make_pf()
   expect_error(plotEigenvalues(m), "no applicable method|not.*found")
})
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: FAIL — `selectCompNum`/`plotFactors`/DoF methods not defined for the class.

- [ ] **Step 3: Add the `plotFactors` generic to `R/defaults.R`**

Add next to the other generic definitions (follow the file's `UseMethod` style):

```r
#' Plot factors for a 3-way decomposition model
#'
#' @param obj a model object.
#' @param ... arguments passed to methods.
#'
#' @export
plotFactors <- function(obj, ...) {
   UseMethod("plotFactors")
}
```

- [ ] **Step 4: Implement the S3 methods in `R/ddsimca3w.R`**

Append:

```r
#' @export
selectCompNum.ddsimca3w <- function(obj, ncomp, ...) {
   if (ncomp < 1 || ncomp > obj$ncomp) stop("Wrong value for 'ncomp'.", call. = FALSE)
   obj$ncomp.selected <- ncomp
   # sub-models are not nested, so the selected sub-model's scores must be
   # re-projected; simplest correct path is to rebuild the cal result. The
   # training matrix is on the model (obj$Xcal, set by makeModel3w).
   obj$res$cal <- predict(obj, obj$Xcal, obj$res$cal$simca$c.ref)
   obj$res$cal$ncomp.selected <- ncomp
   obj$calres <- obj$res$cal
   obj
}

#' Update significance levels for a 3-way DD-SIMCA model
#'
#' @param obj a 3-way DD-SIMCA model (class \code{ddsimca3w}).
#' @param alpha significance level for extremes.
#' @param gamma significance level for outliers.
#' @param ... ignored.
#'
#' @details Unlike \code{setParams.ddsimca}, this re-runs \code{predict} so the
#' data-driven full-distance limits (which \code{setParams.ddsimca} does not
#' know about) are recomputed correctly.
#' @return the model with updated alpha/gamma and rebuilt calibration result.
#' @export
setParams.ddsimca3w <- function(obj, alpha = obj$alpha, gamma = obj$gamma, ...) {
   obj$alpha <- alpha
   obj$gamma <- gamma
   if (!is.null(obj$res) && !is.null(obj$res$cal)) {
      obj$res$cal <- predict(obj, obj$Xcal, obj$res$cal$simca$c.ref, alpha = alpha, gamma = gamma)
      obj$calres <- obj$res$cal
   }
   obj
}

#' @export
summary.ddsimca3w <- function(object, ncomp = object$ncomp.selected, res = object$res, ...) {
   title <- if (object$type == "parafac") "DD-SIMCA-PARAFAC" else "DD-SIMCA-TUCKER"
   fprintf("\n%s model for class '%s'%s\n\n", title, object$classname,
      sprintf(" (%d x %d)", object$dim[1], object$dim[2]))
   if (object$type == "tucker") {
      fprintf("Tucker ranks: (%dx%dx%d-%dx%dx%d) total, (%dx%dx%d) optimal\n",
         1, object$B.ncomp, object$C.ncomp, object$ncomp, object$B.ncomp, object$C.ncomp,
         object$ncomp.selected, object$B.ncomp, object$C.ncomp)
   }
   fprintf("Number of components: %d\n", object$ncomp)
   fprintf("Number of selected components: %d\n", ncomp)
   fprintf("Alpha: %s\n", object$alpha)
   fprintf("Gamma: %s\n\n", object$gamma)
   if (!is.null(object$res) && !is.null(object$res$cal)) {
      sum_data <- do.call(rbind, lapply(res, function(z) as.matrix(z)[ncomp, , drop = FALSE]))
      rownames(sum_data) <- capitalize(names(res))
      print(sum_data, 4)
      cat("\n")
   }
   return(invisible(object))
}

#' @export
print.ddsimca3w <- function(x, ...) {
   title <- if (x$type == "parafac") "DD-SIMCA-PARAFAC" else "DD-SIMCA-TUCKER"
   fprintf("\n%s model (class 'ddsimca3w')\n", title)
   cat("\nCall:\n"); print(x$call)
   cat("\nMain fields:\n")
   cat(" $ncomp - number of components\n")
   cat(" $ncomp.selected - number of selected components\n")
   cat(" $dim - 3-way dimensions (J x K)\n")
   cat(" $limParams - distance distribution parameters\n")
   cat(" $res - list with results (calibration)\n")
   return(invisible(x))
}

#' Plot a factor (loading) curve for a 3-way model.
#'
#' @param obj a 3-way DD-SIMCA model (class \code{ddsimca3w}).
#' @param comp 1-indexed component within the chosen sub-model.
#' @param mode which factor matrix to draw, "B" (mode-2) or "C" (mode-3).
#' @param ncomp which trained sub-model to read from (default selected).
#' @param type plot type.
#' @param ... arguments passed to \code{mdaplot}.
#' @export
plotFactors.ddsimca3w <- function(obj, comp = 1, mode = c("B", "C"),
   ncomp = obj$ncomp.selected, type = "l", ...) {
   mode <- match.arg(mode)
   if (ncomp < 1 || ncomp > obj$ncomp) stop("Wrong value for 'ncomp'.", call. = FALSE)
   M <- obj$models[[ncomp]][[mode]]
   if (comp < 1 || comp > ncol(M)) stop("Wrong value for 'comp'.", call. = FALSE)
   y <- M[, comp]
   pd <- matrix(y, nrow = 1)
   colnames(pd) <- seq_along(y)
   rownames(pd) <- sprintf("Comp %d", comp)
   attr(pd, "name") <- sprintf("%s factors (%s)",
      if (obj$type == "parafac") "PARAFAC" else "Tucker", mode)
   attr(pd, "xaxis.name") <- sprintf("Mode-%d index", if (mode == "B") 2 else 3)
   mdaplot(pd, type = type, ylab = "Factor values", ...)
}

#' @export
plotScores.ddsimca3w <- function(obj, comp = if (obj$ncomp.selected > 1) c(1, 2) else 1, ...) {
   plotScores(obj$res$cal, comp = comp, ...)
}

#' @export
plotLoadings.ddsimca3w <- function(obj, ...) {
   plotFactors(obj, ...)
}

# DoF methods read Nh/Nq/Nf directly from limParams (no $T2lim/$Qlim machinery).
#' @export
plotT2DoF.ddsimca3w <- function(obj, type = "b", labels = "values",
   xticks = seq_len(obj$ncomp), ylab = "Nh", est = "moments", ...) {
   y <- obj$limParams$T2[[est]]$Nu
   pd <- matrix(y, nrow = 1); colnames(pd) <- xticks; rownames(pd) <- "Nh"
   attr(pd, "name") <- "Degrees of freedom"
   mdaplot(pd, type = type, labels = labels, xticks = xticks, ylab = ylab, ...)
}

#' @export
plotQDoF.ddsimca3w <- function(obj, type = "b", labels = "values",
   xticks = seq_len(obj$ncomp), ylab = "Nq", est = "moments", ...) {
   y <- obj$limParams$Q[[est]]$Nu
   pd <- matrix(y, nrow = 1); colnames(pd) <- xticks; rownames(pd) <- "Nq"
   attr(pd, "name") <- "Degrees of freedom"
   mdaplot(pd, type = type, labels = labels, xticks = xticks, ylab = ylab, ...)
}

#' @export
plotDistDoF.ddsimca3w <- function(obj, type = "b", labels = "values",
   xticks = seq_len(obj$ncomp), est = "moments", ...) {
   pd <- rbind(
      Nh = obj$limParams$T2[[est]]$Nu,
      Nq = obj$limParams$Q[[est]]$Nu,
      Nf = obj$limParams$F[[est]]$Nu
   )
   colnames(pd) <- xticks
   attr(pd, "name") <- "Degrees of freedom"
   mdaplotg(pd, type = type, labels = labels, xticks = xticks, ...)
}
```

> `plotDistDoF.ddsimca3w` shows **Nh, Nq AND Nf** (3 series) because for `nf.as.sum = FALSE`, `Nf != Nh + Nq` — the informative case unique to 3-way. The `.pca` version shows only Nh/Nq. `selectCompNum.ddsimca3w` and `setParams.ddsimca3w` both rebuild the calibration result via `predict(obj, obj$Xcal, ...)`; `obj$Xcal` is stored on the model by `makeModel3w` (Task 8).

- [ ] **Step 5: Regenerate docs and run the suite**

Run: `Rscript -e 'roxygen2::roxygenize()'`
Run: `Rscript -e 'devtools::test(filter = "ddsimca3w")'`
Expected: PASS (all S3-glue tests). The "no applicable method" test confirms dropping `"pca"` works as intended.

- [ ] **Step 6: Full test suite + lint**

Run: `Rscript -e 'devtools::test()'`
Expected: PASS — no regressions in the ~10,000 existing tests.

Run: `Rscript -e 'lintr::lint("R/ddsimca3w.R")'`
Expected: no errors (warnings about line length must be fixed if any exceed 120).

- [ ] **Step 7: Commit**

```bash
git add R/ddsimca3w.R R/defaults.R NAMESPACE man/ tests/testthat/test-ddsimca3w.R
git commit -m "Add ddsimca3w S3 methods: selectCompNum, summary, print, plotFactors, DoF plots"
```

---

## Self-Review

**Spec coverage check (every spec section → task):**
- Public API (`ddsimca.parafac`/`ddsimca.tucker`/`predict.ddsimca3w`/`plotFactors`) → Tasks 8, 9, 10. ✓
- Classes without `"pca"` → Tasks 8, 9 (constructors) + Task 10 ("no applicable method" test). ✓
- Reused unchanged (`ddsimcares`, `classify`, `processMembers/Strangers`, limit fitters, result methods) → consumed in Task 8 `predict.ddsimca3w`. ✓
- Training (decompose → sign-fix → project → A_mean/S_pinv/SVD → H/Q → limParams → classify) → Tasks 2/3 (decompose), Task 8 (`buildSubModel3w`, `project3w`, `finalize3w`, `makeModel3w`). ✓
- Projection helpers (pinv3way, lstsq, Khatri-Rao order, Tucker A/B, T/U) → Task 1. ✓
- Sign-fix (Tucker only) → Task 3. ✓
- F + the one shared change (`classify(fp=)`, `f0.model`, three result methods) → Tasks 4, 5. ✓
- Decompositions (parafac.als/.nn, tucker.hooi/.nn, tol/iter, warning) → Tasks 2, 3. ✓
- Small S3 details (selectCompNum, summary/print, plotScores/plotLoadings adapters, as.matrix ncol(T2), DoF plots) → Tasks 5, 6, 10. ✓
- Edits to existing files (ddsimca.R, ddsimcares.R, pca.R, defaults.R, NAMESPACE/man) → Tasks 4, 5, 6, 10. ✓
- Testing (1e-9 golden, invariants, end-to-end, PCA regression) → Tasks 1–3 (invariants), 4 (PCA regression), 7 (fixtures), 8/9 (1e-9 golden + end-to-end). ✓

**Type-consistency check:** `project3w`/`finalize3w` return shapes match what `buildPcaresShape3w` and `classify` consume; `fp = list(u0, Nu)` matches `lp$F$moments`/`lp$F$robust` produced by `finalize3w`; `models[[nc]]` field names (`B`,`C`,`G`,`A_mean`,`S_pinv`,`A_V`,`A_sigma`) are produced by `buildSubModel3w` and read by `project3w`. ✓

**Known risk flags for the implementer (not placeholders — real watch-points):**
1. The mode-unfolding reshapes in `parafac.als` (`T3`/`X2`/`X3u`) and `tucker.*` (`.unfoldB`/`.unfoldC`) are the single most error-prone area. The Task 1/2/3 invariant tests (rank-1 recovery, reconstruction error, fold/unfold round-trip) are the gate — do not advance past a failing reshape test.
2. `diag(Nu, ncol)` in `finalize3w` must pass explicit `ncol` to avoid the scalar-`diag` trap; the golden test catches it.
3. `tucker.nn` mode-2/mode-3 multiplicative updates use the verified partner reshapes (`B`-update: `kron(A,C)` partner + `G2 = aperm(G, c(2,3,1))`; `C`-update: `kron(B,A)` + `G3 = aperm(G, c(3,1,2))`) — these were numerically validated (CP rel 7.7e-8, Tucker rel 5.4e-5, all factors non-negative). Do not "simplify" the partner order; the reconstruction assertion in the Task 3 test is the gate.
4. If `m.models[nc]` key names in the Python oracle differ from this plan's assumptions, read the actual Python (`ddsimca_parafac.py:178-184`, `ddsimca_tucker.py:223-229`) before dumping — never guess.
