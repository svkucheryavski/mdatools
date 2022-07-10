#####################################
# Tests for crossval class method   #
#####################################

nobj <- 24

# tests for performance statistics
context(sprintf("crossval: cross-validation"))

test_that("leave one out works correctly", {
   cv.ind <- matrix(seq_len(nobj), ncol = 1)
   expect_equivalent(cv.ind, crossval(cv = "loo", nobj))
})


test_that("random cross-validation with one repetition works correctly", {

   # random full cross-validation
   nseg = nobj
   set.seed(42)
   cv.ind <- matrix(rep(seq_len(nseg), length.out = nobj)[sample(nobj)], ncol = 1)
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = 1, nobj))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("rand", 24), nobj))


   # number of objets is a multuply of number segments
   nseg = 4
   set.seed(42)
   cv.ind <- matrix(rep(seq_len(nseg), length.out = nobj)[sample(nobj)], ncol = 1)
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = 4, nobj))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("rand", 4), nobj))

   # number of objets is not a multuply of number segments
   nseg = 5
   set.seed(42)
   cv.ind <- matrix(rep(seq_len(nseg), length.out = nobj)[sample(nobj)], ncol = 1)
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = 5, nobj))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("rand", 5), nobj))
})


test_that("random cross-validation with several repetition works correctly", {
   # number of objets is a multuply of number segments
   set.seed(42)
   nseg <- 4
   nrep <- 3
   cv.ind <- sapply(seq_len(nrep), function(i) rep(seq_len(nseg), length.out = nobj)[sample(nobj)])
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("rand", 4, 3), nobj))

   # number of objets is not a multuply of number segments
   set.seed(42)
   nseg <- 5
   nrep <- 4
   cv.ind <- sapply(seq_len(nrep), function(i) rep(seq_len(nseg), length.out = nobj)[sample(nobj)])

   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("rand", 5, 4), nobj))
})


test_that("systematic cross-validation works correctly", {
   # number of objets is a multuply of number segments
   nseg <- 4
   set.seed(42)
   resp <- rnorm(nobj)
   cv.ind <- matrix(rep(seq_len(nseg), length.out = nobj)[order(resp)], ncol = 1)
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("ven", 4), resp = resp))

   # number of objets is not a multuply of number segments
   nseg <- 5
   set.seed(42)
   resp <- rnorm(nobj)
   cv.ind <- matrix(rep(seq_len(nseg), length.out = nobj)[order(resp)], ncol = 1)
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("ven", 5), resp = resp))

   # it also works well without any response
   expect_equivalent(crossval(list("ven", 4), 10), matrix(c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2), ncol = 1))
})

