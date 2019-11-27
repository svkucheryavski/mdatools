##########################################################
# Tests for basic functionality of plotseries() methods  #
##########################################################

# types of plots
all_plots <- c("p", "d", "l", "b", "h", "e")
nonscatter_plots <- c("l", "b", "h", "e")
scatter_plots <- c("p", "d")

main_plots <- c("p", "l", "b", "h")
linebar_plots <- c("l", "b", "h")


#################################################
# Block 1. preparePlotData()                    #
#################################################

context("plotseries: prepare data for plot")

## empty data always raises an error
test_that("if data is empty it raises an error for all plots", {
   expect_error(preparePlotData(NULL))
})


## vectors should be fine for main plots
test_that("vector without dimension handles correctly", {
   x <- 1:4
   pd <- preparePlotData(x)
   expect_false(is.null(dim(pd$visible)))
   expect_equal(dim(pd$visible), c(1, 4))
   expect_equal(attr(pd$visible, "xaxis.values"), 1:4)
   expect_equal(attr(pd$visible, "yaxis.values"), 1)
   expect_null(pd$excluded)
   expect_null(pd$excluded_rows)
   expect_null(pd$excluded_cols)
})


## data frame should be converted to matrix
test_that("vector without dimension handles correctly", {
   x <- data.frame(x = 1:4)
   pd <- preparePlotData(x)
   expect_false(is.list(dim(pd$visible)))
   expect_false(is.null(dim(pd$visible)))
   expect_equal(dim(pd$visible), c(4, 1))
})


## test names for rows and columns - same expectations for all plots
x <- matrix(1:9, ncol = 3)
rownames <- paste0("O", 1:3)
colnames <- paste0("X", 1:3)

test_that("if no rownames are not specified they are added", {
   expect_equal(rownames(preparePlotData(x)$visible), rownames)
})

test_that("if no colnames are not specified they are added", {
   expect_equal(colnames(preparePlotData(x)$visible), colnames)
})

x <- matrix(1:9, ncol = 3)
rownames(x) <- paste0("Obj", 1:3)
colnames(x) <- paste0("Var", 1:3)

test_that("if rownames are specified they are kept", {
   expect_equal(rownames(preparePlotData(x)$visible), rownames(x))
})

test_that("if colnames are specified they are kept", {
   expect_equal(colnames(preparePlotData(x)$visible), colnames(x))
})

## remove one column of three - similar behaviour for all plots
xn1c <- mda.exclcols(x, c(FALSE, TRUE, FALSE))
test_that("missing columns are handled correctly (more than 1 is left)", {
   for (p in all_plots) {
      pd <- preparePlotData(xn1c)
      expect_equal(dim(pd$visible), c(3, 2))
      expect_equal(colnames(pd$visible), c("Var1", "Var3"))
      expect_equal(rownames(pd$visible), rownames(x))
      expect_equal(attr(pd$visible, "xaxis.values"), c(1, 3))
      expect_equal(attr(pd$visible, "yaxis.values"), 1:3)

      expect_equal(pd$excluded_cols, 2)
      expect_null(pd$excluded_rows)
      expect_null(pd$excluded)
   }
})

## remove one row of three - similar behaviour for all plots
xn1r <- mda.exclrows(x, c(FALSE, TRUE, FALSE))
test_that("missing rows are handled correctly (more than 1 is left)", {
   pd <- preparePlotData(xn1r)

   # visible part of data
   expect_equal(dim(pd$visible), c(2, 3))
   expect_equal(colnames(pd$visible), colnames(x))
   expect_equal(rownames(pd$visible), c("Obj1", "Obj3"))
   expect_equal(attr(pd$visible, "xaxis.values"), 1:3)
   expect_equal(attr(pd$visible, "yaxis.values"), c(1, 3))

   # excluded values
   expect_equal(dim(pd$excluded), c(1, 3))
   expect_equal(colnames(pd$excluded), colnames(x))
   expect_equal(rownames(pd$excluded), c("Obj2"))
   expect_equal(attr(pd$excluded, "xaxis.values"), 1:3)
   expect_equal(attr(pd$excluded, "yaxis.values"), 2)
   expect_equal(pd$excluded_rows, 2)
   expect_null(pd$excluded_cols)
})


## remove two columns of three
xn2c <- mda.exclcols(x, c(TRUE, FALSE, TRUE))
test_that("missing columns are handled correctly (only one is left)", {
   pd <- preparePlotData(xn2c)
   expect_equal(dim(pd$visible), c(3, 1))
   expect_equal(colnames(pd$visible), c("Var2"))
   expect_equal(rownames(pd$visible), rownames(x))
   expect_equal(attr(pd$visible, "xaxis.values"), 2)
   expect_equal(attr(pd$visible, "yaxis.values"), 1:3)

   expect_equal(pd$excluded_cols, c(1, 3))
   expect_null(pd$excluded_rows)
   expect_null(pd$excluded)
})

## remove two rows of three
xn2r <- mda.exclrows(x, c(TRUE, FALSE, TRUE))
test_that("missing rows are handled correctly (only one is left)", {
   pd <- preparePlotData(xn2r)

   # visible part of data
   expect_equal(dim(pd$visible), c(1, 3))
   expect_equal(colnames(pd$visible), colnames(x))
   expect_equal(rownames(pd$visible), c("Obj2"))
   expect_equal(attr(pd$visible, "xaxis.values"), 1:3)
   expect_equal(attr(pd$visible, "yaxis.values"), 2)

   # excluded values
   expect_equal(dim(pd$excluded), c(2, 3))
   expect_equal(colnames(pd$excluded), colnames(x))
   expect_equal(rownames(pd$excluded), c("Obj1", "Obj3"))
   expect_equal(attr(pd$excluded, "xaxis.values"), 1:3)
   expect_equal(attr(pd$excluded, "yaxis.values"), c(1, 3))
   expect_equal(pd$excluded_rows, c(1, 3))
   expect_null(pd$excluded_cols)
})

#################################################
# Block 2. splitPlotData()                      #
#################################################

context("plotseries: split plot data")

x <- matrix(1:9, ncol = 3)

test_that("data is split correctly for scatter plots", {
   res <- splitPlotData(preparePlotData(x)$visible, "p")
   expect_equivalent(res$x_values, x[, 1])
   expect_equivalent(res$y_values, x[, 2])
   expect_equal(attr(res$x_values, "name"), "X1")
   expect_equal(attr(res$y_values, "name"), "X2")
   expect_equal(names(res$x_values), c("O1", "O2", "O3"))
})

test_that("data is split correctly for line plots", {
   res <- splitPlotData(preparePlotData(x)$visible, "l")
   expect_equivalent(res$x_values, 1:3)
   expect_equivalent(res$y_values, x)
   expect_equal(attr(res$x_values, "name"), "Variables")
   expect_equal(attr(res$y_values, "name"), "")
   expect_equal(names(res$x_values), c("X1", "X2", "X3"))
})


test_that("data is split correctly for bar plots", {
   res <- splitPlotData(preparePlotData(x)$visible, "h")
   expect_equivalent(res$x_values, 1:3)
   expect_equivalent(res$y_values, x)
   expect_equal(attr(res$x_values, "name"), "Variables")
   expect_equal(attr(res$y_values, "name"), "")
   expect_equal(names(res$x_values), c("X1", "X2", "X3"))
})


test_that("data is split correctly for errorbar plots", {
   res <- splitPlotData(preparePlotData(x)$visible, "e")
   errorbar_data <- rbind(x[1, ], x[1, ] - x[2, ], x[1, ] + x[3, ])
   expect_equivalent(res$x_values, 1:3)
   expect_equivalent(res$y_values, errorbar_data)
   expect_equal(attr(res$x_values, "name"), "Variables")
   expect_equal(attr(res$y_values, "name"), "")
   expect_equal(names(res$x_values), c("X1", "X2", "X3"))
})

test_that("xaxis.name argument is handled correctly", {
   xaxis.name <- "My x variable"
   attr(x, "xaxis.name") <- xaxis.name
   for (type in nonscatter_plots) {
      res <- splitPlotData(preparePlotData(x)$visible, type)
      expect_equal(attr(res$x_values, "name"), xaxis.name)
   }

   for (type in scatter_plots) {
      res <- splitPlotData(preparePlotData(x)$visible, type)
      expect_equal(attr(res$x_values, "name"), "X1")
   }
})

test_that("yaxis.name argument is handled correctly", {
   yaxis.name <- "My y variable"
   attr(x, "yaxis.name") <- yaxis.name
   for (type in nonscatter_plots) {
      res <- splitPlotData(preparePlotData(x)$visible, type)
      expect_equal(attr(res$y_values, "name"), yaxis.name)
   }

   for (type in scatter_plots) {
      res <- splitPlotData(preparePlotData(x)$visible, type)
      expect_equal(attr(res$y_values, "name"), "X2")
   }
})

###########################################
# Block 3. Test excluded rows processing  #
###########################################

context("plotseries: excluded rows processing")

x1 <- matrix(1:12, ncol = 2)
x2 <- matrix(1:12, ncol = 2)
x1 <- mda.exclrows(x1, c(TRUE, FALSE, TRUE, FALSE))

test_that("excluded rows are processed correctly for scatter plot", {
   x1r <- plotseries(x1, type = "p")
   expect_equivalent(x1r$x_values_excluded, x1[c(1, 3), 1])
   expect_equivalent(x1r$y_values_excluded, x1[c(1, 3), 2])

   x2r <- plotseries(x2, type = "p")
   expect_equal(x2r$x_values_excluded, NULL)
   expect_equal(x2r$y_values_excluded, NULL)
})

test_that("excluded rows are processed correctly for line plots", {
   x1r <- plotseries(x1, type = "l")
   expect_equivalent(x1r$x_values_excluded, x1r$x_values)
   expect_equivalent(x1r$y_values_excluded, x1[c(1, 3), ])

   x2r <- plotseries(x2, type = "l")
   expect_equal(x2r$x_values_excluded, NULL)
   expect_equal(x2r$y_values_excluded, NULL)
})

test_that("trying to call this function for other plot gives error (if excluded rows exist)", {
   expect_error(plotseries(x1, type = "e"))
   expect_silent(plotseries(x2, type = "e"))
})


###########################################
# Block 4. Test data labels preparation   #
###########################################

context("plotseries: prepare data labels")

excluded_rows <- c(TRUE, FALSE, TRUE, FALSE)
excluded_cols <- c(FALSE, TRUE, FALSE)
x1 <- x2 <- x3 <- matrix(1:12, ncol = 3)
x2 <- mda.exclrows(x2, excluded_rows)
x3 <- mda.exclcols(x3, excluded_cols)

row_labels <- c("A", "B", "C", "D")
row_labels_wrong <- c("A", "C")
col_labels <- c("X", "Y", "Z")
col_labels_wrong <- c("X", "Z")

tf <- function(data, type) {
   ps <- list()
   ps$type <- type
   plot_data <- preparePlotData(data)
   ps <- c(ps, splitPlotData(plot_data$visible, type))

   ps$excluded_cols <- plot_data$excluded_cols
   ps$excluded_rows <- plot_data$excluded_rows

   if (!is.null(plot_data$excluded)) {
      ps <- c(ps, splitExcludedData(plot_data$excluded, type))
   }

   return(ps)
}

test_that("processing of labels provided by user", {

   # data without hidden rows and cols (scatter plots)
   x1r <- tf(x1, "p")
   expect_error(getDataLabels(x1r, row_labels_wrong))

   x1r <- getDataLabels(x1r, row_labels)
   expect_equal(x1r$labels, row_labels)
   expect_equal(x1r$labels_excluded, NULL)

   # data with hidden rows (scatter plots)
   x2r <- tf(x2, "p")
   expect_error(getDataLabels(x2r, row_labels_wrong))

   x2r <- getDataLabels(x2r, row_labels)
   expect_equal(x2r$labels, row_labels[!excluded_rows])
   expect_equal(x2r$labels_excluded, row_labels[excluded_rows])

   # data with hidden cols (scatter plots)
   x3r <- tf(x3, "p")
   expect_error(getDataLabels(x3r, row_labels_wrong))

   x3r <- getDataLabels(x3r, row_labels)
   expect_equal(x3r$labels, row_labels)
   expect_equal(x3r$labels_excluded, NULL)


   for (type in c("l", "b")) {
      # data with hidden rows (non scatter plots)
      x2r <- tf(x2, type)
      expect_error(getDataLabels(x2r, col_labels_wrong))

      x2r <- getDataLabels(x2r, col_labels)
      expect_equivalent(x2r$labels, col_labels)
      expect_equivalent(x2r$labels_excluded, NULL)
   }

   for (type in nonscatter_plots) {
      # data without hidden rows and cols (non scatter plots)
      x1r <- tf(x1, type)
      expect_error(getDataLabels(x1r, col_labels_wrong))

      x1r <- getDataLabels(x1r, col_labels)
      expect_equivalent(x1r$labels, col_labels)
      expect_equal(x1r$labels_excluded, NULL)

      # data with hidden cols (scatter plots)
      x3r <- tf(x3, type)
      expect_error(getDataLabels(x3r, col_labels_wrong))

      x3r <- getDataLabels(x3r, col_labels)
      expect_equivalent(x3r$labels, col_labels[!excluded_cols])
      expect_equivalent(x3r$labels_excluded, col_labels[excluded_cols])
   }
})


test_that("processing of labels which were not specified (scatter)", {
   x1r <- tf(x1, "p")
   x1r <- getDataLabels(x1r, NULL)
   expect_equal(x1r$labels, c("O1", "O2", "O3", "O4"))
   expect_equal(x1r$labels_excluded, NULL)

   x2r <- tf(x2, "p")
   x2r <- getDataLabels(x2r, NULL)
   expect_equal(x2r$labels, c("O2", "O4"))
   expect_equal(x2r$labels_excluded, c("O1", "O3"))

   x3r <- tf(x3, "p")
   x3r <- getDataLabels(x3r, NULL)
   expect_equal(x3r$labels, c("O1", "O2", "O3", "O4"))
   expect_equal(x3r$labels_excluded, NULL)
})

test_that("processing of labels which were not specified (other plots)", {

   for (type in c("l", "b")) {
      x2r <- tf(x2, type)
      x2r <- getDataLabels(x2r, NULL)
      expect_equal(x2r$labels, c("X1", "X2", "X3"))
      expect_equal(x2r$labels_excluded, NULL)
   }

   for (type in nonscatter_plots) {
      x1r <- tf(x1, type)
      x1r <- getDataLabels(x1r, NULL)
      expect_equal(x1r$labels, c("X1", "X2", "X3"))
      expect_equal(x1r$labels_excluded, NULL)

      x3r <- tf(x3, type)
      x3r <- getDataLabels(x3r, NULL)
      expect_equal(x3r$labels, c("X1", "X3"))
      expect_equal(x3r$labels_excluded, NULL)
   }
})

test_that("processing of labels specified as 'values' (scatter)", {
   x1r <- tf(x1, "p")
   x1r <- getDataLabels(x1r, "values")
   expect_equivalent(x1r$labels, x1[, 2])
   expect_equivalent(x1r$labels_excluded, NULL)

   x2r <- tf(x2, "p")
   x2r <- getDataLabels(x2r, "values")
   expect_equivalent(x2r$labels, x2[!excluded_rows, 2])
   expect_equivalent(x2r$labels_excluded, x2[excluded_rows, 2])

   x3r <- tf(x3, "p")
   x3r <- getDataLabels(x3r, "values")
   expect_equivalent(x3r$labels, x3[, 3])
   expect_equivalent(x3r$labels_excluded, NULL)
})

test_that("processing of labels specified as 'values' (other plots)", {

   for (type in c("l", "b")) {
      x2r <- tf(x2, type)
      x2r <- getDataLabels(x2r, "values")
      expect_equivalent(x2r$labels, apply(x2[!excluded_rows, , drop = F], 2, max))
      expect_equivalent(x2r$labels_excluded, apply(x2[excluded_rows, , drop = F], 2, max))
   }

   for (type in c("h", "e")) {
      x1r <- tf(x1, type)
      x1r <- getDataLabels(x1r, "values")
      expect_equivalent(x1r$labels, x1[1, ])
      expect_equivalent(x1r$labels_excluded, NULL)

      x3r <- tf(x3, type)
      x3r <- getDataLabels(x3r, "values")
      expect_equivalent(x3r$labels, x3[1, !excluded_cols, drop = F])
      expect_equivalent(x3r$labels_excluded, NULL)
   }

})

test_that("processing of labels specified as 'indices' (scatter)", {
   x1r <- tf(x1, "p")
   x1r <- getDataLabels(x1r, "indices")
   expect_equivalent(x1r$labels, 1:4)
   expect_equivalent(x1r$labels_excluded, NULL)

   x2r <- tf(x2, "p")
   x2r <- getDataLabels(x2r, "indices")
   expect_equivalent(x2r$labels, c(2, 4))
   expect_equivalent(x2r$labels_excluded, c(1, 3))

   x3r <- tf(x3, "p")
   x3r <- getDataLabels(x3r, "indices")
   expect_equivalent(x3r$labels, 1:4)
   expect_equivalent(x3r$labels_excluded, NULL)
})


test_that("processing of labels specified as 'indices' (other plots)", {

   for (type in c("l", "b")) {
      x2r <- tf(x2, type)
      x2r <- getDataLabels(x2r, "indices")
      expect_equivalent(x2r$labels, 1:3)
      expect_equivalent(x2r$labels_excluded, NULL)
   }

   for (type in nonscatter_plots) {
      x1r <- tf(x1, type)
      x1r <- getDataLabels(x1r, "indices")
      expect_equivalent(x1r$labels, 1:3)
      expect_equivalent(x1r$labels_excluded, NULL)

      x3r <- tf(x3, type)
      x3r <- getDataLabels(x3r, "indices")
      expect_equivalent(x3r$labels, c(1, 3))
      expect_equivalent(x3r$labels_excluded, NULL)
   }
})
