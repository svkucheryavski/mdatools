data(simdata)
D <- simdata$spectra.c[1:5, ]
D <- t(prep.snv(D))
cn <- constraint("norm")
Dr <- employ(cn, D)


context("mcrals: test implemented constrains")

for (cl in getImplementedConstraints()) {
   test_that(sprintf("%s works correctly with default parameters", cl$name), {
      expect_silent(cn <- constraint(cl$name, cl$params))
      expect_silent(Dr <- employ(cn, D))
      expect_true(is.matrix(D))
      expect_equal(dim(Dr), dim(D))
   })
}

# add one of the existing constraints with default parameters
test_that("Non-negativity constraint", {
   cn <- constraint("non-negativity")
   Dr <- employ(cn, D)
   expect_equal(sum(Dr < 0), 0)
})

test_that("Non-negativity constraint with user defined parameters (wrong)", {
   expect_error(cn <- constraint("non-negativity", params = list(type = "xxx")))
})

# add one of the existing constraints with default parameters
test_that("Normalization constraint", {
   cn <- constraint("norm")
   Dr <- employ(cn, D)

   # all columns should have unit area
   expect_equivalent(colSums(abs(Dr)), rep(1, ncol(Dr)))
})

test_that("Non-negativity constraint with user defined parameters (correct)", {
   cn <- constraint("norm", params = list(type = "length"))
   Dr <- employ(cn, D)

   # all columns should have unit length
   expect_equivalent(sqrt(colSums(Dr^2)), rep(1, ncol(Dr)))
})

test_that("Non-negativity constraint with user defined parameters (wrong)", {
   expect_error(cn <- constraint("norm", params = list(method = "xxx")))
})

# add one of the existing constraints with non-default parameters

# add one of the existing constraints with wrong parameters
#cn <- constraint("non-negativity", thema = "dark")
#Dr <- employ(cn, D)

# check all available existing constraints with default parameters

# add user defined constraint correctly

# add user defined constraint with wrong method

# add user defined constraint with wrong parameters

