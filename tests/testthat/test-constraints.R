data(carbs)
D <- carbs$D
S <- carbs$S
C <- carbs$C

context("mcrals: test implemented constrains")

for (cl in getImplementedConstraints()) {
   test_that(sprintf("%s works correctly with default parameters", cl$name), {
      expect_silent(cn <- constraint(cl$name, cl$params))
      expect_silent(Sr <- employ(cn, S, D))
      expect_true(is.matrix(Sr))
      expect_equal(dim(Sr), dim(S))
   })
}

#!############################
#! Non-negativity constraint #
#!############################

# add one of the existing constraints with default parameters
test_that("Non-negativity constraint works correctly", {
   S[sample(1:length(S))[1:10]] <- -10
   cn <- constraint("nonneg")
   Sr <- employ(cn, S, D)
   expect_equal(sum(Sr < 0), 0)
})

test_that("Non-negativity constraint with user defined parameters (wrong)", {
   expect_error(cn <- constraint("non-negativity", params = list(type = "xxx")))
})

#!###########################
#! Normalization constraint #
#!###########################

# add one of the existing constraints with default parameters
test_that("Normalization constraint works correctly", {
   cn <- constraint("norm")
   Sr <- employ(cn, S, D)

   # all columns should have unit length
   expect_equivalent(colSums(Sr^2), rep(1, ncol(Sr)))
})

test_that("Normalization constraint with user defined parameters (correct)", {
   cn <- constraint("norm", params = list(type = "area"))
   Sr <- employ(cn, S, D)

   # all columns should have unit area
   expect_equivalent(colSums(abs(Sr)), rep(1, ncol(Sr)))

   cn <- constraint("norm", params = list(type = "sum"))
   Sr <- employ(cn, S, D)

   # all columns should have unit sum
   expect_equivalent(colSums(Sr), rep(1, ncol(Sr)))
})

test_that("Normalization constraint with user defined parameters (wrong)", {
   expect_error(cn <- constraint("norm", params = list(method = "xxx")))
})


#!###################
#! Angle constraint #
#!###################

# add one of the existing constraints with default parameters
test_that("Angle constraint works correctly", {
   cn <- constraint("angle")
   Sr <- employ(cn, S, D)

   # compute the expected values
   m <- apply(D, 2, mean)
   m <- m / sqrt(sum(m^2))
   S <- t(prep.norm(t(S), "length"))
   Se <- S * 0.95 + 0.05 * matrix(m, nrow(S), ncol(S))

   # compare the expected and constrained values
   expect_equivalent(Sr, Se)
})

test_that("Angle constraint with user defined parameters (correct)", {
   cn <- constraint("angle", params = list(weight = 0))
   Sr <- employ(cn, S, D)

   # all columns should have unit length
   expect_equivalent(Sr, t(prep.norm(t(S), "length")))

   cn <- constraint("angle", params = list(weight = 1))
   Sr <- employ(cn, S, D)

   # compute the expected values
   m <- apply(D, 2, mean)
   m <- m / sqrt(sum(m^2))

   expect_equivalent(Sr, matrix(m, nrow(S), ncol(S)))
})

test_that("Angle constraint with user defined parameters (wrong)", {
   expect_error(cn <- constraint("norm", params = list(method = "xxx")))
})


#!###############
#! Unimodality  #
#!###############

# generate data with extra peaks
x  <- 1:500
y1 <- dnorm(x, m = 100, s = 20) * 0.8  + dnorm(x, m = 200, s = 10) * 0.2
y2 <- dnorm(x, m = 100, s = 10) * 0.2  + dnorm(x, m = 200, s = 20) * 0.8
y3 <- dnorm(x, m = 250, s = 20)
y <- cbind(y1, y2, y3)
y <- y + matrix(rnorm(length(y), 0, max(y) * 0.05), nrow(y), ncol(y))

check_unimodality <- function(y, tol = 0) {
   n <- length(y)
   im <- which.max(y)
   dl <- round(diff(y[im:1]) / abs(y[im:2]), 3)
   dr <- round(diff(y[im:n]) / abs(y[im:(n-1)]), 3)
   return(sum(dl > tol) + sum(dr > tol))
}

test_that("Unimodality constraint works correctly", {
   cn1 <- constraint("unimod")
   y.new1 <- employ(cn1, y, NULL)

   cn2 <- constraint("unimod", params = list(tol = 0.2))
   y.new2 <- employ(cn2, y, NULL)

   expect_true(all(apply(y, 2, check_unimodality, tol = 0) > 0))
   expect_true(all(apply(y.new1, 2, check_unimodality, tol = 0) < 0.00000001))
   expect_true(all(apply(y.new2, 2, check_unimodality, tol = 0.2) < 0.20))
})


#!###############
#! Closure      #
#!###############

# generate data with violation of the closure
x  <- 1:50
y1 <- x/10
y2 <- 1 / (1 + exp(-(20 - x))) - (x - 20) / 100 + 0.5
y3 <- max(y2) - y2 - (x - 20) / 100 + 0.5
y <- cbind(y1, y2, y3)
y <- y + matrix(rnorm(length(y), 0, max(y) * 0.01), nrow(y), ncol(y))

test_that("Closure constraint works correctly", {
   cn1 <- constraint("closure")
   y.new1 <- employ(cn1, y, NULL)

   cn2 <- constraint("closure", params = list(sum = 10))
   y.new2 <- employ(cn2, y, NULL)

   expect_true(length(unique(rowSums(y))) > 1)
   expect_true(all(abs(rowSums(y.new1) - 1) < 10^-8))
   expect_true(all(abs(rowSums(y.new2) - 10) < 10^-8))
})

