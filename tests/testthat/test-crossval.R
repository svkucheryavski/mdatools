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
   cv.ind <- matrix(rep(seq_len(nseg), length.out = nobj)[order(order(resp))], ncol = 1)
   cv.computed <- crossval(cv = list("ven", 4), resp = resp);
   expect_equivalent(cv.ind, cv.computed)

   # check segment wise
   # sort resp values and split sorted values into segments using systematic indices: 1, 2, 3, 4, 1, 2, 3, 4, 1, 2...
   resp.sorted = sort(resp)
   cv.ind.sorted = rep(seq_len(nseg), length.out = nobj)

   # check that both splits are identicall
   for (i in 1:4) {
      expect_equivalent(sort(resp[cv.computed == i]), resp.sorted[cv.ind.sorted == i])
   }

   # number of objets is not a multuply of number segments
   nseg <- 5
   set.seed(42)
   resp <- rnorm(nobj)
   cv.ind <- matrix(rep(seq_len(nseg), length.out = nobj)[order(order(resp))], ncol = 1)
   expect_equivalent(cv.ind, crossval(cv = list("ven", 5), resp = resp))

   # it also works well without any response
   expect_equivalent(crossval(list("ven", 4), 10), matrix(c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2), ncol = 1))

   # to catch a possible bug
   #     1  2  3  4   5   6  7  8  9  10 11  12  - normal order
   #     2  4  8  7  11   3 10  9  1   6 12   5  - order of indices for sorted y-vales
   #     1  2  3  4   1   2  3  4  1   2  3   4  - normal order of segments
   #     1  1  2  2   4   2  4  3  4   3  1   3
   y = c(9, 1, 6, 2, 12, 10, 4, 3, 8,  7, 5, 11)
   expect_equivalent(crossval(list("ven", 4), 12, y), matrix(c(1, 1, 2, 2, 4, 2, 4, 3, 4, 3, 1, 3), ncol = 1))

})

