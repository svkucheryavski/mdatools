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

test_that("DoF plot generics still dispatch for pca/ddsimca objects", {
   data(people)
   x <- people[, -4]
   m <- pca(x, ncomp = 3, center = TRUE, scale = TRUE, lim.type = "ddmoments")
   pdf(tempfile()); on.exit(dev.off())
   # plotT2DoF.pca draws (no show.plot arg); assert it dispatches without error
   expect_error(plotT2DoF(m), NA)
   expect_true(is.function(plotT2DoF))
})

context("ddsimca3w: PARAFAC golden")

gpath <- function(f) test_path("ddsimca3w", f)
rdcsv <- function(f) as.matrix(read.csv(gpath(f), header = FALSE))

test_that("PARAFAC project3w + finalize3w match the Python oracle (1e-9) on fixed factors", {
   skip_if_not(file.exists(gpath("parafac_X_train.csv")))
   X <- rdcsv("parafac_X_train.csv")
   ncomp <- 3; J <- 5; K <- 4
   # Build sub-models exactly as production does: PARAFAC's H-statistics come
   # from the ALS A factor (Python's factors[0], dumped as parafac_A_*.csv),
   # NOT a re-projection of X onto (B, C). Re-projecting here puts H ~4e-5 off
   # the oracle; feeding A_train reproduces it at machine precision.
   models <- lapply(seq_len(ncomp), function(nc) {
      A <- rdcsv(sprintf("parafac_A_%d.csv", nc))
      B <- rdcsv(sprintf("parafac_B_%d.csv", nc))
      C <- rdcsv(sprintf("parafac_C_%d.csv", nc))
      buildSubModel3w(A, B, C, NULL, X, FALSE, "parafac")
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
      # buildSubModel3w signature is (A, B, C, G, x, non.negative, type). For
      # Tucker the A arg is IGNORED (stats come from a re-projection, matching
      # the Python oracle's A_for_stats), so pass NULL here.
      buildSubModel3w(NULL, B, C, G, X, FALSE, "tucker")
   })
   model <- list(type = "tucker", models = models, non.negative = FALSE,
      dim = c(J, K), ncomp = ncomp, ncomp.selected = ncomp)
   pj <- project3w(model, X)
   expect_equal(unname(pj$H), unname(rdcsv("tucker_H.csv")), tolerance = 1e-9)
   expect_equal(unname(pj$Q), unname(rdcsv("tucker_Q.csv")), tolerance = 1e-9)

   lp <- finalize3w(pj$H, pj$Q, nf.as.sum = FALSE)
   # classic (moments)
   tf_c <- rdcsv("tucker_fParams_classic.csv")
   expect_equal(unname(lp$F$moments$u0), unname(tf_c[1, ]), tolerance = 1e-9)
   expect_equal(unname(lp$F$moments$Nu), unname(tf_c[2, ]), tolerance = 1e-9)
   th_c <- rdcsv("tucker_hParams_classic.csv")
   expect_equal(unname(lp$T2$moments$u0), unname(th_c[1, ]), tolerance = 1e-9)
   expect_equal(unname(lp$T2$moments$Nu), unname(th_c[2, ]), tolerance = 1e-9)
   tq_c <- rdcsv("tucker_qParams_classic.csv")
   expect_equal(unname(lp$Q$moments$u0), unname(tq_c[1, ]), tolerance = 1e-9)
   expect_equal(unname(lp$Q$moments$Nu), unname(tq_c[2, ]), tolerance = 1e-9)
   # robust
   tf_r <- rdcsv("tucker_fParams_robust.csv")
   expect_equal(unname(lp$F$robust$u0), unname(tf_r[1, ]), tolerance = 1e-9)
   expect_equal(unname(lp$F$robust$Nu), unname(tf_r[2, ]), tolerance = 1e-9)
   th_r <- rdcsv("tucker_hParams_robust.csv")
   expect_equal(unname(lp$T2$robust$u0), unname(th_r[1, ]), tolerance = 1e-9)
   expect_equal(unname(lp$T2$robust$Nu), unname(th_r[2, ]), tolerance = 1e-9)
   tq_r <- rdcsv("tucker_qParams_robust.csv")
   expect_equal(unname(lp$Q$robust$u0), unname(tq_r[1, ]), tolerance = 1e-9)
   expect_equal(unname(lp$Q$robust$Nu), unname(tq_r[2, ]), tolerance = 1e-9)
})

test_that("ddsimca.tucker trains and predicts end-to-end", {
   set.seed(21)
   J <- 5; K <- 4; n <- 30; R <- 3
   A0 <- matrix(rnorm(n * R), n, R)
   B0 <- matrix(rnorm(J * R), J, R)
   C0 <- matrix(rnorm(K * R), K, R)
   x <- A0 %*% t(khatriRao3w(B0, C0)) + 0.01 * matrix(rnorm(n * J * K), n, J * K)
   m <- ddsimca.tucker(x, "T", dim = c(J, K), ncomp = c(3, 5, 4))
   expect_s3_class(m, "ddsimca.tucker")
   expect_s3_class(m, "ddsimca3w")
   expect_s3_class(m$res$cal, "ddsimcares")
   expect_equal(ncol(m$res$cal$T2), 3L)
   expect_equal(m$B.ncomp, 5)
})

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
   expect_error(plotQDoF(m), NA)
   expect_error(plotDistDoF(m), NA)
})

test_that("a PCA-only method has no applicable method on the 3-way model", {
   m <- make_pf()
   expect_error(plotEigenvalues(m), "no applicable method|not.*found")
})

context("ddsimca3w: 3-way array input")

# A 3-way array and the equivalent COLUMN-MAJOR (R-native dim()<-) unfolded matrix
# must give identical models: an array is unfolded for the user, and matrix input
# is interpreted column-major and refolded through `dim`. So `dim(X) <- c(n, J*K)`
# (the natural R idiom) is correct.

# low-rank tensor (+ tiny noise) so PARAFAC/Tucker converge cleanly. `mat` is the
# column-major (R-native dim()<-) unfold of the array, the layout the matrix API
# now expects, so `array` and `mat` are equivalent inputs.
make3w <- function(seed, J, K, n, R) {
   set.seed(seed)
   A0 <- matrix(rnorm(n * R), n, R)
   B0 <- matrix(rnorm(J * R), J, R)
   C0 <- matrix(rnorm(K * R), K, R)
   xmat <- A0 %*% t(khatriRao3w(B0, C0)) + 0.01 * matrix(rnorm(n * J * K), n, J * K)
   A3 <- aperm(array(t(xmat), dim = c(K, J, n)), c(3, 2, 1))   # fold to [n, J, K]
   list(array = A3, mat = matrix(A3, nrow = n))                # column-major unfold
}

test_that("ddsimca.parafac accepts a 3-way array identical to the unfolded matrix", {
   J <- 5; K <- 4; n <- 30
   d <- make3w(123, J, K, n, 3)

   # PARAFAC has random ALS init; seed identically so both fits are comparable
   set.seed(3); m_mat <- ddsimca.parafac(d$mat,   "T", dim = c(J, K), ncomp = 3)
   set.seed(3); m_arr <- ddsimca.parafac(d$array, "T", ncomp = 3)   # dim inferred

   expect_s3_class(m_arr, "ddsimca.parafac")
   expect_equal(m_arr$dim, c(J, K))
   expect_equal(m_arr$res$cal$T2, m_mat$res$cal$T2)
   expect_equal(m_arr$res$cal$Q,  m_mat$res$cal$Q)
})

test_that("ddsimca.tucker accepts a 3-way array identical to the unfolded matrix", {
   J <- 5; K <- 4; n <- 30
   d <- make3w(321, J, K, n, 3)

   m_mat <- ddsimca.tucker(d$mat,   "T", dim = c(J, K), ncomp = c(3, 5, 4))
   m_arr <- ddsimca.tucker(d$array, "T", ncomp = c(3, 5, 4))   # dim inferred

   expect_s3_class(m_arr, "ddsimca.tucker")
   expect_equal(m_arr$dim, c(J, K))
   expect_equal(m_arr$res$cal$T2, m_mat$res$cal$T2)
   expect_equal(m_arr$res$cal$Q,  m_mat$res$cal$Q)
})

test_that("predict.ddsimca3w accepts a 3-way array identical to the unfolded matrix", {
   J <- 5; K <- 4
   tr <- make3w(123, J, K, 30, 3)
   set.seed(3); m <- ddsimca.parafac(tr$array, "T", ncomp = 3)

   nw <- make3w(55, J, K, 12, 3)   # new data, array + its row-major unfold
   r_arr <- predict(m, nw$array, rep("T", 12))
   r_mat <- predict(m, nw$mat,   rep("T", 12))

   expect_s3_class(r_arr, "ddsimcares")
   expect_equal(r_arr$T2, r_mat$T2)
   expect_equal(r_arr$Q,  r_mat$Q)
})

test_that("matrix input is interpreted column-major (R-native dim()<- idiom)", {
   J <- 5; K <- 4; n <- 30
   set.seed(8)
   A3 <- array(rnorm(n * J * K), dim = c(n, J, K))
   Mcol <- A3; dim(Mcol) <- c(n, J * K)            # the exact dim()<- idiom users reach for

   set.seed(4); m_arr <- ddsimca.parafac(A3,   "T", ncomp = 3)
   set.seed(4); m_mat <- ddsimca.parafac(Mcol, "T", dim = c(J, K), ncomp = 3)
   expect_equal(m_mat$res$cal$T2, m_arr$res$cal$T2)
   expect_equal(m_mat$res$cal$Q,  m_arr$res$cal$Q)
})

test_that("ddsimca.parafac errors when a matrix is given without 'dim'", {
   x <- matrix(rnorm(30 * 20), 30, 20)
   expect_error(ddsimca.parafac(x, "T", ncomp = 3), "dim")
})

test_that("ddsimca.parafac errors when 'dim' conflicts with array dimensions", {
   A3 <- array(rnorm(30 * 5 * 4), dim = c(30, 5, 4))
   expect_error(ddsimca.parafac(A3, "T", dim = c(4, 5), ncomp = 3), "match")
})

context("ddsimca3w: factor normalization")

# Build a positive low-rank tensor (returned as a 3-way array, so it is
# convention-independent) — true B/C loadings are positive bumps.
posTensor <- function(seed, J, K, n, R) {
   set.seed(seed)
   A0 <- matrix(abs(rnorm(n * R)) + 0.5, n, R)
   B0 <- matrix(abs(rnorm(J * R)) + 0.5, J, R)
   C0 <- matrix(abs(rnorm(K * R)) + 0.5, K, R)
   xmat <- A0 %*% t(khatriRao3w(B0, C0)) + 0.01 * matrix(rnorm(n * J * K), n, J * K)
   aperm(array(t(xmat), dim = c(K, J, n)), c(3, 2, 1))   # [n, J, K] array
}

test_that("ddsimca.parafac factors follow tensorly convention: unit-norm columns, positive mean", {
   J <- 5; K <- 4
   x <- posTensor(11, J, K, 30, 2)
   set.seed(1)
   m <- suppressWarnings(ddsimca.parafac(x, "T", ncomp = 2))
   B <- m$models[[2]]$B; C <- m$models[[2]]$C
   expect_equal(sqrt(colSums(B^2)), rep(1, 2), tolerance = 1e-8)   # unit-norm
   expect_equal(sqrt(colSums(C^2)), rep(1, 2), tolerance = 1e-8)
   expect_true(all(colMeans(B) >= 0))                             # positive-mean sign convention
   expect_true(all(colMeans(C) >= 0))
})

test_that("ddsimca.parafac calibration Q is consistent with the stored factors", {
   # Guards against the cal pipeline transforming the input twice: the
   # constructor unfolds the array, then the (former) predict() refolded the
   # already-unfolded matrix column-major and scrambled it. The rank-nc
   # sub-model's reported Q must equal the residual of re-projecting its own
   # stored B/C onto the same internal unfold used to fit them.
   J <- 5; K <- 4; n <- 30
   x <- posTensor(12, J, K, n, 2)
   set.seed(2)
   m <- suppressWarnings(ddsimca.parafac(x, "T", ncomp = 2))
   xi <- matrix(aperm(x, c(1, 3, 2)), nrow = n)   # internal row-major unfold
   for (nc in seq_len(m$ncomp)) {
      pj <- mdatools:::projectParafac(m$models[[nc]]$B, m$models[[nc]]$C, xi)
      expect_equal(unname(rowSums(pj$E^2)), unname(m$res$cal$Q[, nc]), tolerance = 1e-8)
   }
})

context("ddsimca3w: convergence reporting")

test_that("ddsimca.parafac is silent by default even when a sub-model hits max.iter", {
   # posTensor is ~rank-2; ncomp=3 over-factors so the rank-3 CP-ALS sub-model
   # swamps and runs to max.iter. The constructor must stay silent by default.
   x <- posTensor(12, 5, 4, 30, 2)
   set.seed(2)
   expect_silent(m <- ddsimca.parafac(x, "T", ncomp = 3))
})

# capture every message emitted while evaluating expr
captureMessages3w <- function(expr) {
   msgs <- character(0)
   withCallingHandlers(force(expr),
      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
   msgs
}

test_that("ddsimca.parafac verbose = TRUE names the over-factored sub-model by ncomp", {
   x <- posTensor(12, 5, 4, 30, 2)   # ~rank-2; the rank-3 sub-model over-factors
   set.seed(2)
   msgs <- captureMessages3w(ddsimca.parafac(x, "T", ncomp = 3, verbose = TRUE))
   expect_true(any(grepl("ncomp = 3", msgs, fixed = TRUE)))
   expect_true(any(grepl("reached max iterations", msgs, fixed = TRUE)))
   expect_true(all(nchar(trimws(msgs)) <= 70))   # message must fit 70-column width
})

test_that("ddsimca.tucker is silent by default and names the sub-model by ncomp when verbose", {
   x <- posTensor(7, 6, 5, 40, 2)
   set.seed(3)
   expect_silent(m <- ddsimca.tucker(x, "T", ncomp = c(2, 3, 3)))
   set.seed(3)
   msgs <- captureMessages3w(ddsimca.tucker(x, "T", ncomp = c(2, 3, 3), verbose = TRUE))
   expect_true(any(grepl("ncomp = 2", msgs, fixed = TRUE)))
   expect_true(any(grepl("reached max iterations", msgs, fixed = TRUE)))
   expect_true(all(nchar(trimws(msgs)) <= 70))   # message must fit 70-column width
})
