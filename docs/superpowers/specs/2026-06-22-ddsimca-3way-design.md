# DD-SIMCA for 3-way data (PARAFAC & Tucker) — R design (simple)

Status: draft for review
Date: 2026-06-22
Branch: `develop`

## Principle

Almost everything is already in `R/ddsimca.R` / `R/ddsimcares.R`. The only new
work is: **decompose** the 3-way data (PARAFAC or Tucker), **project** to get
scores `A` and residual `E`, **compute `H`, `Q`, `F`**, and hand them to the
existing `classify()` → `ddsimcares()` machinery. The `ddsimcares` result class
is **identical** across all three methods. Python (`ddsimca-py`) is a numerical
oracle only.

**Explicitly out of scope** (per the package author): excluded rows/columns,
preprocessing (`prep`), and `center`/`scale`. The constructors take a clean
unfolded matrix and decompose it as-is. Also deferred: PCV (`$res$pv`) and JSON.

## Public API (new file `R/ddsimca3w.R`)

```r
ddsimca.parafac(x, classname, dim, ncomp = min(dim[1], dim[2], nrow(x) - 1),
               non.negative = FALSE, nf.as.sum = FALSE, alpha = 0.05, gamma = 0.01)

ddsimca.tucker(x, classname, dim, ncomp,   # ncomp = c(R_A_max, R_B, R_C)
              non.negative = FALSE, nf.as.sum = FALSE, alpha = 0.05, gamma = 0.01)

predict.ddsimca3w(object, x, c.ref = NULL,
                    alpha = object$alpha, gamma = object$gamma, ...)

plotFactors(obj, ...)                                       # new generic
plotFactors.ddsimca3w(obj, comp = 1, mode = c("B","C"), ncomp = obj$ncomp.selected, ...)
```

- `x` — numeric matrix `n x (J*K)`, tensor unfolded **row-major** (`j*K + k`).
- `dim = c(J, K)`, both integer `>= 2`.
- `ncomp` — PARAFAC: integer `<= min(J, K, n-1)`. Tucker: `c(R_A_max, R_B, R_C)`,
  with `R_B <= J`, `R_C <= K`, `R_A_max <= n-1`, `R_A_max <= R_B*R_C`.
- `non.negative` — non-negative decomposition (and Tucker Branch B projection).
- `nf.as.sum` — full-distance DoF rule (see "F" below). Default FALSE for 3-way.

Classes: model `c("ddsimca.parafac","ddsimca3w","ddsimca")` /
`c("ddsimca.tucker","ddsimca3w","ddsimca")` (no `"pca"` — see "Why no `pca`
class" below); result unchanged `c("ddsimcares","pcares","ldecomp")`.

## What is reused unchanged

`ddsimcares()`, `classify()`, `processMembers()`, `processStrangers()`,
`ddmoments.param()`/`ddrobust.param()`/`ldecomp.getLimParams()`, and every result
method (`summary/print/plot.ddsimcares`, `plotAcceptance/Distances/Extremes/
Aliens/FoM/FoMs/SelectivityArea/Sensitivity.ddsimcares`, `as.matrix/
as.data.frame/writeCSV.ddsimcares`). `setParams.ddsimca` reused (+ `fp`, below).

## Training (per sub-model rank `a = 1 … R_A_max`; sub-models are NOT nested)

1. Decompose `x` with rank `a` (PARAFAC) or `(a, R_B, R_C)` (Tucker).
2. Tucker only: canonical sign-fix (below) before projecting.
3. Scores `A`: PARAFAC = library factor; Tucker = re-project via the helper.
4. `A_mean = colMeans(A)`, `A_c = A - A_mean`, `S = crossprod(A_c)/(n-1)`
   (full `a x a`), `S_pinv = pinv3way(S)`; `svd(A_c)` → `V`, `sigma`.
5. `H[, a] = rowSums((A_c %*% S_pinv) * A_c)`   (Mahalanobis; NOT `T^2/eigenvalue`).
   `Q[, a] = rowSums(E_a^2)`   (`E_a` = reconstruction residual, below).

Then:
- `limParams$Q  <- ldecomp.getLimParams(Q)`, `limParams$T2 <- ldecomp.getLimParams(H)`
  (gives the `moments`/`robust` `(u0, Nu)` shape `classify()` expects; reuses the
  existing `clamp.dof` policy).
- `limParams$F` per `nf.as.sum` (below).
- Build a small `pcares`-shaped object (`scores`, `T2 = H`, `Q`, `residuals = E`,
  `expvar`, `ncomp`, `ncomp.selected`, dimnames), store the per-sub-model factors
  (`B`, `C`, Tucker `G`, `A_mean`, `S_pinv`, `V`, `sigma`) and per-sub-model
  scores (`scores.list`) on the model, plus metadata (`type`, `dim`,
  `non.negative`, `nf.as.sum`, `nrows`, `nclasses=1`, `classname`, `alpha`,
  `gamma`, `limType="moments"`, `call`; Tucker also `B.ncomp`, `C.ncomp`).
- `model$res$cal <- predict(model, x, rep(classname, n))`; sync `$calres`.

A fixed seed is set before each decomposition (deterministic results);
`.Random.seed` is saved/restored (`on.exit`) so the caller's RNG is untouched.

## Projection helpers (single source of truth at predict)

`pinv3way(X, rcond = 1e-15)` is a private SVD pseudo-inverse (NumPy semantics);
the exported `pinv()` is untouched. `lstsq(A, x) := pinv3way(A) %*% x`. R is
column-major, so all row-major reshapes (`j*K+k`, mode-1 unfolding `G_(1)`) are
forced via `aperm` + `matrix(byrow=TRUE)` in one tensor-utils spot.

- **PARAFAC:** `KR[, r] = kron(B[, r], C[, r])` (order `(B, C)`); `a_hat = lstsq(KR, x)`.
- **Tucker A (orthogonal):** `M = kron(B, C) %*% t(G_(1))`; `a_hat = lstsq(M, x)`.
- **Tucker B (non-negative):** `a_hat = x %*% S_kron %*% pinv3way(t(S_kron)%*%S_kron)
  %*% t(G_mat) %*% pinv3way(G_mat%*%t(G_mat))`, `S_kron = kron(B,C)`,
  `G_mat = G_(1)`. Intentionally **not** equal to Branch A for non-orthogonal
  B/C — implement literally.
- **Residual:** `x_hat = a_hat %*% t(KR)` (PARAFAC) or `%*% t(M)` (Tucker; Branch B
  materializes `M` only here). `E = x - x_hat`.
- **T / U** (selected sub-model only): `T = A_c %*% V` (R's `svd()$v` is `V`);
  `U = T / sigma` with a `1e-12` guard (degenerate columns → 0 in `U` only).

`_project()` returns `(H, Q, T, U, E)`: `H`/`Q` are `n x ncomp` (per-sub-model,
column-stacked — never broadcast one sub-model); `T`/`U` are `n x ncomp.selected`;
`E` is `n x (J*K)` from the final sub-model.

## Sign-fix (Tucker only)

Compute `sA`, `sB`, `sC` from the **unmodified** factors (each column's max-abs
entry's sign; exactly-zero → `+1`; first-index tie-break), THEN flip the factor
columns and apply the single core update `G[a,b,c] *= sA[a]*sB[b]*sC[c]`. PARAFAC
needs none (its ALS returns unnormalized factors; projection from fixed `B`/`C`
is stable).

## F and the one shared change

The full distance is the scaled `F = (H/h0)*Nh + (Q/q0)*Nq` (NOT `h + q`).

`classify()` gets a trailing optional `fp = NULL`:
- `fp = NULL` → existing sum rule (`Nf = Nh+Nq`, scale `f0.model = Nf`). PCA /
  scalar `ddsimca` behavior is **bit-identical** (regression-tested).
- `fp` supplied → per sub-model `Nf = fp$Nu`, `f0.model = fp$u0`,
  `fce = f0.model/Nf * qchisq(1-alpha, Nf)`, `fco` likewise with `gamma`.

`classify()` stores `f0.model` per outcome (distinct from `res$f0`, the alien fit
in `processStrangers`). Three result methods that assume `f0.model == Nf` need it
threaded (all reduce to current code when `f0.model == Nf`):
- `plotExtremes.ddsimcares`: `pchisq(f * Nf / f0.model, Nf)`.
- `plotDistances.ddsimcares`/`.ddsimca` (`f` branch): normalize by `f0.model`.
- `plotSelectivityArea.ddsimcares`: `fcrit <- f0.model/Nf * qchisq(1-alpha, Nf)`
  (keep alien `z <- fcrit/f0`).

Finalize:
```
if (nf.as.sum) {                 # f0.model = Nf = Nh+Nq; do NOT re-clamp the sum
  F$moments <- list(u0 = Nh_c+Nq_c, Nu = Nh_c+Nq_c)
  F$robust  <- list(u0 = Nh_r+Nq_r, Nu = Nh_r+Nq_r)
} else {                         # data-driven (paper eq.7)
  F$moments <- ddmoments.param((H/h0_c)*Nh_c + (Q/q0_c)*Nq_c)
  F$robust  <- ddrobust.param ((H/h0_r)*Nh_r + (Q/q0_r)*Nq_r)
  # ddmoments.param/ddrobust.param do NOT clamp (only ldecomp.getLimParams does);
  # clamp F's DoF to the package-wide integer [1,250] range explicitly:
  F$moments$Nu <- clamp.dof(F$moments$Nu)
  F$robust$Nu  <- clamp.dof(F$robust$Nu)
}
```

`setParams.ddsimca` makes two `classify()` calls; each gets its estimator
explicitly (`F[["moments"]]` / `F[["robust"]]`, both `NULL` for PCA).
`predict.ddsimca3w` forwards the same way.

## Decompositions (idiomatic R, zero dependencies)

`parafac.als` (CP-ALS), `parafac.als.nn` (multiplicative updates), `tucker.hooi`
(HOSVD init + HOOI), `tucker.nn` (multiplicative updates). `DECOMP_TOL = 1e-10`,
`max iter = 500`; non-convergence → `warning()`. These match tensorly's
algorithms (validated against goldens at ~1e-6 + invariants; byte-identity not
expected).

## Small S3 details

- `selectCompNum.ddsimca3w`: set `ncomp.selected` on model and each `$res`,
  re-point `$scores <- scores.list[[ncomp]]`, refresh aliases.
- `summary`/`print.ddsimca3w`: title `"DD-SIMCA-PARAFAC"`/`"DD-SIMCA-TUCKER"`;
  append ` (J x K)`; Tucker appends the rank-tuple descriptors
  (`(1xR_BxR_C-R_AmaxxR_BxR_C)` total, `(R_AselxR_BxR_C)` optimal, re-rendered each
  call). `invisible()` returns.
- `plotScores.ddsimca3w` → result's `plotScores.ldecomp` with a
  score-width-aware `comp` default; `plotLoadings.ddsimca3w` → `plotFactors`.
  (No `"pca"` in the model class, so PCA-only methods that need
  `$loadings`/`$eigenvals` — `plotBiplot`, `plotEigenvalues`, `plotResiduals`,
  `plotCumVariance`, `getProbabilities`, `asvector`, `getCalibrationData`,
  `asjson`/`writeJSON` — simply don't dispatch and yield R's native "no
  applicable method" error. No stub methods needed.)
- `as.matrix.ddsimcares`: key the component count/labels off `ncol(x$T2)` (equals
  `ncol(scores)` for PCA; correct for 3-way where `scores` is the selected
  sub-model). Backward-compatible.
- **DoF plots:** R's DoF plots are `plotT2DoF`/`plotQDoF`/`plotDistDoF` (currently
  plain functions in `pca.R` reading `$T2lim`/`$Qlim` row 4). Convert the three to
  S3 generics — the existing bodies become `.pca` methods (scalar `pca`/`ddsimca`
  unchanged) — and add `.ddsimca3w` methods that read DoF directly from
  `limParams`: `Nh = limParams$T2$<est>$Nu`, `Nq = limParams$Q$<est>$Nu`,
  `Nf = limParams$F$<est>$Nu` (no `$T2lim`/`$Qlim` machinery on the 3-way model).
  `plotDistDoF.ddsimca3w` shows **Nh, Nq, and Nf** — for `nf.as.sum = FALSE`,
  `Nf != Nh + Nq`, which is the informative case unique to 3-way.

## Edits to existing files

- `R/ddsimca.R`: `classify(..., fp = NULL)` storing `f0.model`; `setParams.ddsimca`
  + `plotDistances.ddsimca` (`f` branch) use `f0.model`.
- `R/ddsimcares.R`: `f0.model` in `plotExtremes`/`plotDistances`/
  `plotSelectivityArea`; `ncol(T2)` for the component count/labels in
  `as.matrix.ddsimcares` (and the `ncomp` bounds in `as.data.frame`/
  `plotAcceptance`/`plotDistances`); `plotScores.ddsimcares` adapter.
- `R/pca.R`: convert `plotT2DoF`/`plotQDoF`/`plotDistDoF` from plain functions to
  S3 generics (`UseMethod`) with their current bodies as `*.pca` methods (no
  behavior change for `pca`/`ddsimca`); add `*.ddsimca3w` methods reading
  `limParams` (Nh/Nq/Nf).
- `R/defaults.R`: `plotFactors` generic (and the `plotT2DoF`/`plotQDoF`/
  `plotDistDoF` generic stubs if kept here rather than in `pca.R`).
- `NAMESPACE`/`man` via roxygen (escape `%` as `\%`).

## Testing (Python = oracle)

1. **Fixed-factor 1e-9:** feed factors exported from `ddsimca-py/.venv` into the R
   projection → `H/Q/F` → limits → `classify`; exact `(u0, Nu)` → exact
   roles/decisions. (Export script needed — goldens currently lack factors; dump
   as CSV.)
2. **Invariants:** PARAFAC-on-`(B,C)` ≡ Tucker δ-core; Branch A ≡ Branch B under
   orthonormal `B`/`C` and ≠ `lstsq(M,x)` otherwise; sign-fix canonical
   invariance; a tiny analytic `(A,B,C)→X` fixture for Khatri-Rao order / reshape.
3. **End-to-end ~1e-6** on `demo/simdata.mat` (aggregate FoMs with tolerance;
   skip borderline role/decision flips).
4. **PCA regression:** `classify(fp=NULL)` bit-identical; existing `test-ddsimca.R`
   stays green.

## Why no `pca` class (resolved)

The model class is `c("ddsimca.<type>","ddsimca3w","ddsimca")` — `"pca"` is
deliberately **omitted**. Keeping it would silently route every PCA-only method
(`plotEigenvalues`, `plotBiplot`, `plotResiduals`, …) to a `*.pca` body that
reads `$loadings`/`$eigenvals` the 3-way model doesn't have, crashing with a
confusing subscript error — and forcing a list of `stop()` stubs to pre-empt it.
Dropping `"pca"` makes those calls fail cleanly with R's native "no applicable
method", needs zero stubs, and costs nothing: everything reused from
`ddsimca.R`/`ddsimcares.R` dispatches on `ddsimca3w`/`ddsimca` (model) or
`ddsimcares`/`pcares`/`ldecomp` (result), never on `pca` (model). The DoF plots
still dispatch correctly: `plotT2DoF(model)` finds `.ddsimca3w`; scalar
`pca`/`simca`/`ddsimca` objects still fall through to `.pca`.

Function naming is dotted — `ddsimca.parafac`/`ddsimca.tucker` — consistent with
the existing `ddsimca.fromjson`/`ddsimca.readJSON` plain helpers. `ddsimca` is not
an S3 generic, so these are inert plain constructors (no dispatch collision).
