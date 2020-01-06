#####################################
# Tests for crossval class method   #
#####################################

nobj <- 24

# tests for performance statistics
context(sprintf("crossval: cross-validation"))

test_that("leave one out works correctly", {
   set.seed(42)
   cv.ind <- array(sample(nobj), dim = c(nobj, 1, 1))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = 1, nobj))
})


test_that("random cross-validation with one repetition works correctly", {
   # number of objets is a multuply of number segments
   set.seed(42)
   cv.ind <- array(matrix(sample(nobj), nrow = 4, byrow = TRUE), dim = c(4, nobj/4, 1))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = 4, nobj))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("rand", 4), nobj))

   # number of objets is not a multuply of number segments
   set.seed(42)
   cv.ind <- array(matrix(c(sample(nobj), NA), nrow = 5, byrow = TRUE), dim = c(5, ceiling(nobj/5), 1))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = 5, nobj))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("rand", 5), nobj))
})


test_that("random cross-validation with several repetition works correctly", {
   # number of objets is a multuply of number segments
   set.seed(42)
   cv.ind <- array(
      cbind(
         matrix(sample(nobj), nrow = 4, byrow = TRUE),
         matrix(sample(nobj), nrow = 4, byrow = TRUE),
         matrix(sample(nobj), nrow = 4, byrow = TRUE)
      ), dim = c(4, nobj/4, 3)
   )
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("rand", 4, 3), nobj))
   expect_equivalent(as.numeric(table(crossval(cv = list("rand", 4, 3), nobj))), rep(3, nobj))

   # number of objets is not a multuply of number segments
   set.seed(42)
   cv.ind <- array(
      cbind(
         matrix(c(sample(nobj), NA), nrow = 5, byrow = TRUE),
         matrix(c(sample(nobj), NA), nrow = 5, byrow = TRUE),
         matrix(c(sample(nobj), NA), nrow = 5, byrow = TRUE),
         matrix(c(sample(nobj), NA), nrow = 5, byrow = TRUE)
      ), dim = c(5, ceiling(nobj/5), 4)
   )

   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("rand", 5, 4), nobj))
   expect_equivalent(as.numeric(table(crossval(cv = list("rand", 5, 4), nobj), useNA = "ifany")), rep(4, nobj + 1))
})


test_that("random cross-validation with several repetition works correctly", {
   # number of objets is a multuply of number segments
   set.seed(42)
   resp <- rnorm(nobj)
   cv.ind <- array(matrix(order(resp), nrow = 4, byrow = FALSE), dim = c(4, nobj/4, 1))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("ven", 4), resp = resp))

   # number of objets is not a multuply of number segments
   set.seed(42)
   resp <- rnorm(nobj)
   cv.ind <- array(matrix(c(order(resp), NA), nrow = 5, byrow = FALSE), dim = c(5, ceiling(nobj/5), 1))
   set.seed(42)
   expect_equivalent(cv.ind, crossval(cv = list("ven", 5), resp = resp))
})

