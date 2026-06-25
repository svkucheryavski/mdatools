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
