########################################################
# Tests for basic functionality of misc mda.* methods  #
########################################################

# Prepare datasets

factor_labels <- c("F1", "F2")
row_names <- paste0("O", 1:5)
x <- factor(c(1, 1, 2, 2, 1), labels = factor_labels)
d <- data.frame(Height = c(180, 190, 178, 167, 199), Gender = x, Weight = c(80, 90, 78, 67, 88))
rownames(d) <- row_names
d <- mda.exclrows(d, 2)
attr(d, "yaxis.name") <- "Objects"
attr(d, "yaxis.values") <- seq_len(5) * 10

context("misc: small utilities")

test_that("mda.df2mat() and mda.t() work correctly", {
   m <- mda.df2mat(d)

   expect_true(is.matrix(m))
   expect_equal(ncol(m), ncol(d))
   expect_equal(nrow(m), nrow(d))
   expect_equal(rownames(m), row_names)
   expect_equal(colnames(m), c(colnames(d)[c(1, 3)], factor_labels[1]))
   expect_equal(attr(m, "exclrows"), attr(d, "exclrows"))
   expect_equal(attr(m, "yaxis.name"), attr(d, "yaxis.name"))
   expect_equal(attr(m, "yaxis.values"), attr(d, "yaxis.values"))

   m <- mda.df2mat(d, full = TRUE)

   expect_true(is.matrix(m))
   expect_equal(ncol(m), ncol(d) + 1)
   expect_equal(nrow(m), nrow(d))
   expect_equal(rownames(m), row_names)
   expect_equal(colnames(m), c(colnames(d)[c(1, 3)], factor_labels))
   expect_equal(attr(m, "exclrows"), attr(d, "exclrows"))
   expect_equal(attr(m, "yaxis.name"), attr(d, "yaxis.name"))
   expect_equal(attr(m, "yaxis.values"), attr(d, "yaxis.values"))


   m <- mda.exclcols(m, 2)
   attr(m, "xaxis.name") <- "Variables"
   attr(m, "xaxis.values") <- seq_len(ncol(m)) * 100
   tm <- mda.t(m)

   expect_equal(attr(tm, "exclcols"), attr(m, "exclrows"))
   expect_equal(attr(tm, "xaxis.name"), attr(m, "yaxis.name"))
   expect_equal(attr(tm, "xaxis.values"), attr(m, "yaxis.values"))

   expect_equal(attr(m, "exclcols"), attr(tm, "exclrows"))
   expect_equal(attr(m, "xaxis.name"), attr(tm, "yaxis.name"))
   expect_equal(attr(m, "xaxis.values"), attr(tm, "yaxis.values"))

})

context("test purge methods")

data(people)
attr(people, "xaxis.values") <- seq(10, 120, by = 10)
attr(people, "xaxis.name") <- "My variables"
attr(people, "yaxis.values") <- seq(0.01, 0.32, by = 0.01)
attr(people, "yaxis.name") <- "My objects"
exclrows <- c(1, 5, 10)
exclcols <- c(4, 10)

x1 <- x2 <- x3 <- x4 <- people
x2 <- mda.exclrows(x2, exclrows)
x3 <- mda.exclcols(x3, exclcols)
x4 <- mda.exclrows(x4, exclrows)
x4 <- mda.exclcols(x4, exclcols)

test_that("purge returns original data if nothing is excluded", {
   x1pr <- mda.purgeRows(x1)
   x1pc <- mda.purgeCols(x1)
   x1p <- mda.purge(x1)

   expect_equal(x1, x1pr)
   expect_equal(x1, x1pc)
   expect_equal(x1, x1pc)

   expect_equal(attributes(x1), attributes(x1pr))
   expect_equal(attributes(x1), attributes(x1pc))
   expect_equal(attributes(x1), attributes(x1pc))

})

test_that("purge rows works correctly", {
   x2pr <- mda.purgeRows(x2)
   x3pr <- mda.purgeRows(x3)
   x4pr <- mda.purgeRows(x4)

   expect_equivalent(x2pr, x1[-exclrows, , drop = FALSE])
   expect_equivalent(x3pr, x1)
   expect_equivalent(x4pr, x1[-exclrows, , drop = FALSE])

   expect_equal(attr(x2pr, "xaxis.name"), attr(x1, "xaxis.name"))
   expect_equal(attr(x3pr, "xaxis.name"), attr(x1, "xaxis.name"))
   expect_equal(attr(x4pr, "xaxis.name"), attr(x1, "xaxis.name"))

   expect_equal(attr(x2pr, "xaxis.values"), attr(x1, "xaxis.values"))
   expect_equal(attr(x3pr, "xaxis.values"), attr(x1, "xaxis.values"))
   expect_equal(attr(x4pr, "xaxis.values"), attr(x1, "xaxis.values"))

   expect_equal(attr(x2pr, "yaxis.name"), attr(x1, "yaxis.name"))
   expect_equal(attr(x3pr, "yaxis.name"), attr(x1, "yaxis.name"))
   expect_equal(attr(x4pr, "yaxis.name"), attr(x1, "yaxis.name"))

   expect_equal(attr(x2pr, "yaxis.values"), attr(x1, "yaxis.values")[-exclrows])
   expect_equal(attr(x3pr, "yaxis.values"), attr(x1, "yaxis.values"))
   expect_equal(attr(x4pr, "yaxis.values"), attr(x1, "yaxis.values")[-exclrows])

})

test_that("purge cols works correctly", {
   x2pc <- mda.purgeCols(x2)
   x3pc <- mda.purgeCols(x3)
   x4pc <- mda.purgeCols(x4)

   expect_equivalent(x2pc, x1)
   expect_equivalent(x3pc, x1[, -exclcols, drop = FALSE])
   expect_equivalent(x4pc, x1[, -exclcols, drop = FALSE])

   expect_equal(attr(x2pc, "xaxis.name"), attr(x1, "xaxis.name"))
   expect_equal(attr(x3pc, "xaxis.name"), attr(x1, "xaxis.name"))
   expect_equal(attr(x4pc, "xaxis.name"), attr(x1, "xaxis.name"))

   expect_equal(attr(x2pc, "xaxis.values"), attr(x1, "xaxis.values"))
   expect_equal(attr(x3pc, "xaxis.values"), attr(x1, "xaxis.values")[-exclcols])
   expect_equal(attr(x4pc, "xaxis.values"), attr(x1, "xaxis.values")[-exclcols])

   expect_equal(attr(x2pc, "yaxis.name"), attr(x1, "yaxis.name"))
   expect_equal(attr(x3pc, "yaxis.name"), attr(x1, "yaxis.name"))
   expect_equal(attr(x4pc, "yaxis.name"), attr(x1, "yaxis.name"))

   expect_equal(attr(x2pc, "yaxis.values"), attr(x1, "yaxis.values"))
   expect_equal(attr(x3pc, "yaxis.values"), attr(x1, "yaxis.values"))
   expect_equal(attr(x4pc, "yaxis.values"), attr(x1, "yaxis.values"))

})

test_that("purge both works correctly", {
   x2pc <- mda.purge(x2)
   x3pc <- mda.purge(x3)
   x4pc <- mda.purge(x4)

   expect_equivalent(x2pc, x1[-exclrows, , drop = FALSE])
   expect_equivalent(x3pc, x1[, -exclcols, drop = FALSE])
   expect_equivalent(x4pc, x1[-exclrows, -exclcols, drop = FALSE])

   expect_equal(attr(x2pc, "xaxis.name"), attr(x1, "xaxis.name"))
   expect_equal(attr(x3pc, "xaxis.name"), attr(x1, "xaxis.name"))
   expect_equal(attr(x4pc, "xaxis.name"), attr(x1, "xaxis.name"))

   expect_equal(attr(x2pc, "xaxis.values"), attr(x1, "xaxis.values"))
   expect_equal(attr(x3pc, "xaxis.values"), attr(x1, "xaxis.values")[-exclcols])
   expect_equal(attr(x4pc, "xaxis.values"), attr(x1, "xaxis.values")[-exclcols])

   expect_equal(attr(x2pc, "yaxis.name"), attr(x1, "yaxis.name"))
   expect_equal(attr(x3pc, "yaxis.name"), attr(x1, "yaxis.name"))
   expect_equal(attr(x4pc, "yaxis.name"), attr(x1, "yaxis.name"))

   expect_equal(attr(x2pc, "yaxis.values"), attr(x1, "yaxis.values")[-exclrows])
   expect_equal(attr(x3pc, "yaxis.values"), attr(x1, "yaxis.values"))
   expect_equal(attr(x4pc, "yaxis.values"), attr(x1, "yaxis.values")[-exclrows])

})