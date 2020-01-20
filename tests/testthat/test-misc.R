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