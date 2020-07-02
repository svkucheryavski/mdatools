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
   cn <- constraint("non-negativity")
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


# add one of the existing constraints with non-default parameters

# check all available existing constraints with default parameters

# add user defined constraint correctly

# add user defined constraint with wrong method

# add user defined constraint with wrong parameters

