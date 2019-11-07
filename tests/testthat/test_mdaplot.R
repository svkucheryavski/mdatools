# prepare dataset
all_plots <- c("p", "l", "b", "h", "e")
nonscatter_plots <- c("l", "b", "h", "e")
main_plots <- c("p", "l", "b", "h")
linebar_plots <- c("l", "b", "h")

par(mfrow = c(2, 2))

data("people")
people <- people[, -6]
people <- as.data.frame(people)

m <- apply(people, 2, mean)
s <- apply(people, 2, sd)
stat1 <- rbind(m, s)
stat2 <- rbind(m, s, 1.5 * s)
colnames(stat1) <- colnames(people)
colnames(stat2) <- colnames(people)

####################################################
# Block 1. Test prepare data for plot function     #
####################################################

context("mdaplot: prepare data for plots")

## empty data always raises an error
test_that("if data is empty it raises an error for all plots", {
   for (p in all_plots) {
      expect_error(prepare_data_for_plot(NULL, !!p))
   }
})

## vectors should be fine for main plots
test_that("vector without dimension handles correctly", {
   x <- 1:4
   for (p in linebar_plots) {
      expect_equal(dim(prepare_data_for_plot(x, !!p)$data), c(1, 4))
   }

   # for scatter plot it is converted to column plus extra column is added
   expect_equal(dim(prepare_data_for_plot(x, "p")$data), c(4, 2))

   # for error bar it raises an error
   expect_error(prepare_data_for_plot(x, "e"))
})

## for errorbar it should give an error


## test names for rows and columns - same expectations for all plots

x <- matrix(1:9, ncol = 3)
rownames <- paste0("O", 1:3)
colnames <- paste0("X", 1:3)

test_that("if no rownames are not specified they are added", {
   for (p in all_plots) {
      expect_equal(rownames(prepare_data_for_plot(x, !!p)$data), rownames)
   }
})

test_that("if no colnames are not specified they are added", {
   for (p in all_plots) {
      expect_equal(colnames(prepare_data_for_plot(x, !!p)$data), colnames)
   }
})

x <- matrix(1:9, ncol = 3)
rownames(x) <- paste0("Obj", 1:3)
colnames(x) <- paste0("Var", 1:3)

test_that("if rownames are specified they are kept", {
   for (p in all_plots) {
      expect_equal(rownames(prepare_data_for_plot(x, !!p)$data), rownames(x))
   }
})

test_that("if colnames are specified they are kept", {
   for (p in all_plots) {
      expect_equal(colnames(prepare_data_for_plot(x, !!p)$data), colnames(x))
   }
})

## remove one column of three - similar behaviour for all plots
xn1c <- mda.exclcols(x, c(FALSE, TRUE, FALSE))
test_that("missing columns are handled correctly (more than 1 is left)", {
   for (p in all_plots) {
      expect_equal(dim(prepare_data_for_plot(xn1c, p)$data), c(3, 2))
      expect_equal(colnames(prepare_data_for_plot(xn1c, p)$data), c("Var1", "Var3"))
      expect_equal(rownames(prepare_data_for_plot(xn1c, p)$data), rownames(x))
      expect_equal(prepare_data_for_plot(xn1c, p)$excluded_cols, 2)
   }
})

## remove one row of three - similar behaviour for all plots
xn1r <- mda.exclrows(x, c(FALSE, TRUE, FALSE))
test_that("missing rows are handled correctly (more than 1 is left)", {
   for (p in all_plots) {
      expect_equal(dim(prepare_data_for_plot(xn1r, p)$data), c(2, 3))
      expect_equal(colnames(prepare_data_for_plot(xn1r, p)$data), colnames(x))
      expect_equal(rownames(prepare_data_for_plot(xn1r, p)$data), c("Obj1", "Obj3"))
      expect_equal(prepare_data_for_plot(xn1r, p)$excluded_rows, 2)
   }
})

## remove two columns of three 
xn2c <- mda.exclcols(x, c(TRUE, FALSE, TRUE))
test_that("missing columns are handled correctly (only one is left)", {
   for (p in nonscatter_plots) {
      expect_equal(dim(prepare_data_for_plot(xn2c, p)$data), c(3, 1))
      expect_equal(colnames(prepare_data_for_plot(xn2c, p)$data), c("Var2"))
      expect_equal(rownames(prepare_data_for_plot(xn2c, p)$data), rownames(x))
      expect_equal(prepare_data_for_plot(xn2c, p)$excluded_cols, c(1, 3))
   }

   # for scatter is a special one (since a column will be added if less than 2)
   expect_equal(dim(prepare_data_for_plot(xn2c, "p")$data), c(3, 2))
   expect_equal(colnames(prepare_data_for_plot(xn2c, "p")$data), c("Objects", "Var2"))
   expect_equal(rownames(prepare_data_for_plot(xn2c, "p")$data), rownames(x))
   expect_equal(prepare_data_for_plot(xn2c, "p")$excluded_cols, c(1, 3))
})

## remove two rows of three 
xn2r <- mda.exclrows(x, c(TRUE, FALSE, TRUE))
test_that("missing rows are handled correctly (only one is left)", {
   for (p in main_plots) {
      expect_equal(dim(prepare_data_for_plot(xn2r, p)$data), c(1, 3))
      expect_equal(colnames(prepare_data_for_plot(xn2r, p)$data), colnames(x))
      expect_equal(rownames(prepare_data_for_plot(xn2r, p)$data), c("Obj2"))
      expect_equal(prepare_data_for_plot(xn2r, p)$excluded_rows, c(1, 3))
   }

   # for error it raises an error
   expect_error(prepare_data_for_plot(xn2r, "e"))
})


####################################################
# Block 2. Test split plot data function     #
####################################################

context("mdaplot: split plot data")

x <- matrix(1:9, ncol = 3)

test_that("data is split correctly for scatter plots", {
   res <- split_plot_data(prepare_data_for_plot(x, "p"), "p")
   expect_equivalent(res$x_values, x[, 1])
   expect_equivalent(res$y_values, x[, 2])
   expect_equal(attr(res$x_values, "name"), "X1")
   expect_equal(attr(res$y_values, "name"), "X2")
   expect_equal(names(res$x_values), c("O1", "O2", "O3"))
})

test_that("data is split correctly for line plots", {
   res <- split_plot_data(prepare_data_for_plot(x, "l"), "l")
   expect_equivalent(res$x_values, 1:3)   
   expect_equivalent(res$y_values, x)
   expect_equal(attr(res$x_values, "name"), "Variables")
   expect_equal(attr(res$y_values, "name"), "")
   expect_equal(names(res$x_values), c("X1", "X2", "X3"))
})

test_that("data is split correctly for bar plots", {
   res <- split_plot_data(prepare_data_for_plot(x, "h"), "h")
   expect_equivalent(res$x_values, 1:3)
   expect_equivalent(res$y_values, x[1, ])
   expect_equal(attr(res$x_values, "name"), "Variables")
   expect_equal(attr(res$y_values, "name"), "")
   expect_equal(names(res$x_values), c("X1", "X2", "X3"))
})

test_that("data is split correctly for errorbar plots", {
   res <- split_plot_data(prepare_data_for_plot(x, "e"), "e")
   expect_equivalent(res$x_values, 1:3)
   expect_equivalent(res$y_values, x[1, ])
   expect_equivalent(res$lower, x[2, ])
   expect_equivalent(res$upper, x[3, ])
   expect_equal(attr(res$x_values, "name"), "Variables")
   expect_equal(attr(res$y_values, "name"), "O1")
   expect_equal(names(res$x_values), c("X1", "X2", "X3"))
})

test_that("data is split correctly for errorbar plots (two rows case)", {
   res <- split_plot_data(prepare_data_for_plot(x[1:2, ], "e"), "e")
   expect_equivalent(res$x_values, 1:3)
   expect_equivalent(res$y_values, x[1, ])
   expect_equivalent(res$lower, x[2, ])
   expect_equivalent(res$upper, x[2, ])
   expect_equal(attr(res$x_values, "name"), "Variables")
   expect_equal(attr(res$y_values, "name"), "O1")
   expect_equal(names(res$x_values), c("X1", "X2", "X3"))
})

###########################################
# Block 3. Test excluded rows processing  #
###########################################

context("mdaplot: excluded rows processing")

x1 <- matrix(1:12, ncol = 2)
x2 <- matrix(1:12, ncol = 2)
x1 <- mda.exclrows(x1, c(TRUE, FALSE, TRUE, FALSE))

test_that("excluded rows are processed correctly for scatter plot", {
   x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, "p"), "p"), "p")
   expect_equivalent(x1r$x_values_excluded, x1[c(1, 3), 1])
   expect_equivalent(x1r$y_values_excluded, x1[c(1, 3), 2])

   x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, "p"), "p"), "p")
   expect_equal(x2r$x_values_excluded, NULL)
   expect_equal(x2r$y_values_excluded, NULL)
})

test_that("excluded rows are processed correctly for line plots", {
   x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, "l"), "l"), "l")
   expect_equivalent(x1r$x_values_excluded, x1r$x_values)
   expect_equivalent(x1r$y_values_excluded, x1[c(1, 3), ])

   x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, "l"), "l"), "l")
   expect_equal(x2r$x_values_excluded, NULL)
   expect_equal(x2r$y_values_excluded, NULL)
})

test_that("trying to call this function for other plot gives error (if excluded rows exist)", {
   expect_error(process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, "h"), "h"), "h"))
   expect_error(process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, "e"), "e"), "e"))
   expect_silent(process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, "h"), "h"), "h"))
   expect_silent(process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, "e"), "e"), "e"))
})

###########################################
# Block 4. Test data labels preparation   #
###########################################

context("mdaplot: prepare data labels")

excluded_rows <- c(TRUE, FALSE, TRUE, FALSE)
excluded_cols <- c(FALSE, TRUE, FALSE)
x1 <- x2 <- x3 <- matrix(1:12, ncol = 3)
x2 <- mda.exclrows(x2, excluded_rows)
x3 <- mda.exclcols(x3, excluded_cols)

row_labels <- c("A", "B", "C", "D")
row_labels_wrong <- c("A", "C")
col_labels <- c("X", "Y", "Z")
col_labels_wrong <- c("X", "Z")

test_that("processing of labels provided by user", {

   # data without hidden rows and cols (scatter plots)
   x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, "p"), "p"), "p")
   expect_error(prepare_plot_data_labels(x1r, "p", row_labels_wrong))
   x1r <- prepare_plot_data_labels(x1r, "p", row_labels)
   expect_equal(x1r$labels, row_labels)
   expect_equal(x1r$labels_excluded, NULL)

   # data with hidden rows (scatter plots)
   x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, "p"), "p"), "p")
   expect_error(prepare_plot_data_labels(x2r, "p", row_labels_wrong))
   x2r <- prepare_plot_data_labels(x2r, "p", row_labels)
   expect_equal(x2r$labels, row_labels[!excluded_rows])
   expect_equal(x2r$labels_excluded, row_labels[excluded_rows])

   # data with hidden cols (scatter plots)
   x3r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x3, "p"), "p"), "p")
   expect_error(prepare_plot_data_labels(x3r, "p", row_labels_wrong))
   x3r <- prepare_plot_data_labels(x3r, "p", row_labels)
   expect_equal(x3r$labels, row_labels)
   expect_equal(x3r$labels_excluded, NULL)


   for (p in c("l", "b")) {
      # data with hidden rows (non scatter plots)
      x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, p), p), p)
      expect_error(prepare_plot_data_labels(x2r, p, col_labels_wrong))
      x2r <- prepare_plot_data_labels(x2r, p, col_labels)
      expect_equivalent(x2r$labels, col_labels)
      expect_equivalent(x2r$labels_excluded, NULL)
   }

   for (p in nonscatter_plots) {
      # data without hidden rows and cols (non scatter plots)
      x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, p), p), p)   
      expect_error(prepare_plot_data_labels(x1r, p, col_labels_wrong))
      x1r <- prepare_plot_data_labels(x1r, p, col_labels)
      expect_equivalent(x1r$labels, col_labels)
      expect_equal(x1r$labels_excluded, NULL)

      # data with hidden cols (scatter plots)
      x3r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x3, p), p), p)
      expect_error(prepare_plot_data_labels(x3r, p, col_labels_wrong))
      x3r <- prepare_plot_data_labels(x3r, p, col_labels)
      expect_equivalent(x3r$labels, col_labels[!excluded_cols])
      expect_equivalent(x3r$labels_excluded, col_labels[excluded_cols])
   }
})

test_that("processing of labels which were not specified (scatter)", {
   x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, "p"), "p"), "p")   
   x1r <- prepare_plot_data_labels(x1r, "p", NULL)
   expect_equal(x1r$labels, c("O1", "O2", "O3", "O4"))
   expect_equal(x1r$labels_excluded, NULL)

   x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, "p"), "p"), "p")
   x2r <- prepare_plot_data_labels(x2r, "p", NULL)
   expect_equal(x2r$labels, c("O2", "O4"))
   expect_equal(x2r$labels_excluded, c("O1", "O3"))

   x3r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x3, "p"), "p"), "p")   
   x3r <- prepare_plot_data_labels(x3r, "p", NULL)
   expect_equal(x3r$labels, c("O1", "O2", "O3", "O4"))
   expect_equal(x3r$labels_excluded, NULL)
})

test_that("processing of labels which were not specified (other plots)", {

   for (p in c("l", "b")) {
      x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, p), p), p)
      x2r <- prepare_plot_data_labels(x2r, p, NULL)
      expect_equal(x2r$labels, c("X1", "X2", "X3"))
      expect_equal(x2r$labels_excluded, NULL)
   }

   for (p in nonscatter_plots) {
      x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, p), p), p)
      x1r <- prepare_plot_data_labels(x1r, p, NULL)
      expect_equal(x1r$labels, c("X1", "X2", "X3"))
      expect_equal(x1r$labels_excluded, NULL)

      x3r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x3, p), p), p)   
      x3r <- prepare_plot_data_labels(x3r, p, NULL)
      expect_equal(x3r$labels, c("X1", "X3"))
      expect_equal(x3r$labels_excluded, NULL)
   }
})

test_that("processing of labels specified as 'values' (scatter)", {
   x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, "p"), "p"), "p")   
   x1r = prepare_plot_data_labels(x1r, "p", "values")
   expect_equivalent(x1r$labels, x1[, 2])
   expect_equivalent(x1r$labels_excluded, NULL)

   x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, "p"), "p"), "p")
   x2r <- prepare_plot_data_labels(x2r, "p", "values")
   expect_equivalent(x2r$labels, x2[!excluded_rows, 2])
   expect_equivalent(x2r$labels_excluded, x2[excluded_rows, 2])

   x3r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x3, "p"), "p"), "p")   
   x3r = prepare_plot_data_labels(x3r, "p", "values")
   expect_equivalent(x3r$labels, x3[, 3])
   expect_equivalent(x3r$labels_excluded, NULL)
})

test_that("processing of labels specified as 'values' (other plots)", {

   for (p in c("l", "b")) {
      x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, p), p), p)
      x2r <- prepare_plot_data_labels(x2r, p, "values")
      expect_equivalent(x2r$labels, apply(x2[!excluded_rows, , drop = F], 2, max))
      expect_equivalent(x2r$labels_excluded, apply(x2[excluded_rows, , drop = F], 2, max))
   }

   for (p in c("h", "e")) {
      x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, p), p), p)
      x1r <- prepare_plot_data_labels(x1r, p, "values")
      expect_equivalent(x1r$labels, x1[1, ])
      expect_equivalent(x1r$labels_excluded, NULL)

      x3r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x3, p), p), p)   
      x3r <- prepare_plot_data_labels(x3r, p, "values")
      expect_equivalent(x3r$labels, x3[1, !excluded_cols, drop = F])
      expect_equivalent(x3r$labels_excluded, NULL)
   }

})

test_that("processing of labels specified as 'indices' (scatter)", {
   x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, "p"), "p"), "p")   
   x1r <- prepare_plot_data_labels(x1r, "p", "indices")
   expect_equivalent(x1r$labels, 1:4)
   expect_equivalent(x1r$labels_excluded, NULL)

   x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, "p"), "p"), "p")
   x2r <- prepare_plot_data_labels(x2r, "p", "indices")
   expect_equivalent(x2r$labels, c(2, 4))
   expect_equivalent(x2r$labels_excluded, c(1, 3))

   x3r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x3, "p"), "p"), "p")
   x3r <- prepare_plot_data_labels(x3r, "p", "indices")
   expect_equivalent(x3r$labels, 1:4)
   expect_equivalent(x3r$labels_excluded, NULL)
})

test_that("processing of labels specified as 'indices' (other plots)", {

   for (p in c("l", "b")) {
      x2r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x2, p), p), p)
      x2r <- prepare_plot_data_labels(x2r, p, "indices")
      expect_equivalent(x2r$labels, 1:3)
      expect_equivalent(x2r$labels_excluded, NULL)
   }

   for (p in nonscatter_plots) {
      x1r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x1, p), p), p)
      x1r <- prepare_plot_data_labels(x1r, p, "indices")
      expect_equivalent(x1r$labels, 1:3)
      expect_equivalent(x1r$labels_excluded, NULL)

      x3r <- process_excluded_rows(split_plot_data(prepare_data_for_plot(x3, p), p), p)
      x3r <- prepare_plot_data_labels(x3r, p, "indices")
      expect_equivalent(x3r$labels, c(1, 3))
      expect_equivalent(x3r$labels_excluded, NULL)
   }
})

####################################################
# Block 5. General tests for mdaplot functionality #
####################################################

context("mdaplot: main functionality")

## test function (shortcut)
tf <- function(...) mdaplot(people, ...)

test_that("can create four main plots without extra parameters", {
   expect_silent(tf(type = "p"))
   expect_silent(tf(type = "l"))
   expect_silent(tf(type = "b"))
   expect_silent(tf(type = "h"))
})

test_that("can create error bar plots as well", {
   expect_silent(mdaplot(stat1, type = "e"))
   expect_silent(mdaplot(stat2, type = "e"))
   expect_silent(mdaplot(stat1, type = "e"))
   expect_silent(mdaplot(stat2, type = "e"))
})

test_that("can work when matrix has only one column", {
   expect_silent(mdaplot(people[, 1, drop = F], type = "p"))
   expect_silent(mdaplot(people[, 1, drop = F], type = "l"))
   expect_silent(mdaplot(people[, 1, drop = F], type = "b"))
   expect_silent(mdaplot(people[, 1, drop = F], type = "h"))
})

test_that("can make errorbar plot when matrix has only one column", {
   expect_silent(mdaplot(stat1[, 1, drop = F], type = "e"))
   expect_silent(mdaplot(stat2[, 1, drop = F], type = "e"))
   expect_silent(mdaplot(stat1[, 1, drop = F], type = "e"))
   expect_silent(mdaplot(stat2[, 1, drop = F], type = "e"))
})

test_that("can work when matrix has only one row", {
   expect_silent(mdaplot(people[1, , drop = F], type = "p"))
   expect_silent(mdaplot(people[1, , drop = F], type = "l"))
   expect_silent(mdaplot(people[1, , drop = F], type = "b"))
   expect_silent(mdaplot(people[1, , drop = F], type = "h"))
})

test_that("can work with vectors as data source (part 1)", {
   expect_silent(mdaplot(people[1, ], type = "p"))
   expect_silent(mdaplot(people[1, ], type = "l"))
   expect_silent(mdaplot(people[1, ], type = "b"))
   expect_silent(mdaplot(people[1, ], type = "h"))
})

test_that("can work with vectors as data source (part 2)", {
   expect_silent(mdaplot(people[, 1], type = "p"))
   expect_silent(mdaplot(people[, 1], type = "l"))
   expect_silent(mdaplot(people[, 1], type = "b"))
   expect_silent(mdaplot(people[, 1], type = "h"))
})

par(mfrow = c(1, 1))
test_that("several plots can be combined together using 'show.axes' parameter", {
   m = apply(people, 2, mean)
   expect_silent(mdaplot(m, type = 'h', col = 'pink', bwd = 0.25, show.labels = T))
   expect_silent(mdaplot(stat1, type = 'e', col = 'red', show.axes = F, show.grid = F))
})
par(mfrow = c(2, 2))


######################################
# Block 6. Handling data attributes  #
######################################

data("simdata")
spectra <- simdata$spectra.c

context("mdaplot: handling data attributes")

attr(people, "name") <- "People"
test_that("handle 'name' attribute correctly", {
   expect_silent(tf(type = "p"))
   expect_silent(tf(type = "l"))
   expect_silent(tf(type = "b"))
   expect_silent(tf(type = "h"))
})

attr(stat1, "name") <- "People statistics"
attr(stat2, "name") <- "People statistics"
test_that("handle 'name' attribute correctly (errorbar)", {
   expect_silent(mdaplot(stat1, type = "e"))
   expect_silent(mdaplot(stat2, type = "e"))
   expect_silent(mdaplot(stat1, type = "e"))
   expect_silent(mdaplot(stat2, type = "e"))
})


attr(spectra, "xaxis.values") <- simdata$wavelength
attr(spectra, "xaxis.name")   <- "Wavelength, nm"
test_that("handle 'xaxis.values' and 'xaxis.names' attribute correctly", {
   expect_silent(mdaplot(spectra, type = "p"))
   expect_silent(mdaplot(spectra, type = "l"))
   expect_silent(mdaplot(spectra, type = "b"))
   expect_silent(mdaplot(spectra, type = "h"))
})

attr(spectra, "yaxis.name") <- "Absorbance"
attr(spectra, "name") <- "UV/Vis spectra"
test_that("handle and 'yaxis.names' attribute correctly", {
   expect_silent(mdaplot(spectra, type = "p"))
   expect_silent(mdaplot(spectra, type = "l"))
   expect_silent(mdaplot(spectra, type = "b"))
   expect_silent(mdaplot(spectra, type = "h"))
})

attr(spectra, "xaxis.values") = 10^7 / simdata$wavelength
attr(spectra, "xaxis.name") = expression("Wavenumber, cm"^-1)
test_that("attribute 'xaxis.values' can be in reverse order", {
   expect_silent(mdaplot(spectra, type = "p"))
   expect_silent(mdaplot(spectra, type = "l"))
   expect_silent(mdaplot(spectra, type = "b"))
   expect_silent(mdaplot(spectra, type = "h"))
})


#########################################
# Block 7. Manual ticks and ticklabels  #
#########################################


context("mdaplot: manual ticks and ticklabels")

xticks <- seq(1, 11, by = 3)
xticklabels <- colnames(people)[xticks]

test_that("manual xticks work correctly", {
   expect_silent(tf(type = "p", xticks = c(165, 180, 195)))
   expect_silent(tf(type = "l", xticks = xticks))
   expect_silent(tf(type = "b", xticks = xticks))
   expect_silent(tf(type = "h", xticks = xticks))
})

test_that("'xticklabels' can not be specified without 'xticks'", {
   expect_error(tf(type = "l", xticklabels = xticklabels))
   expect_error(tf(type = "b", xticklabels = xticklabels))
   expect_error(tf(type = "h", xticklabels = xticklabels))
})


test_that("'xticks' and 'xticklabels' work correctly together", {
   expect_silent(tf(type = "p", xticks = c(165, 180, 195), xticklabels = c("L", "M", "H")))
   expect_silent(tf(type = "l", xticks = xticks, xticklabels = xticklabels))
   expect_silent(tf(type = "b", xticks = xticks, xticklabels = xticklabels))
   expect_silent(tf(type = "h", xticks = xticks, xticklabels = xticklabels))
})

test_that("'xlas' parameter is handling correctly", {
   expect_silent(tf(type = "p", xlas = 2, xticks = c(165, 180, 195), xticklabels=c("L", "M", "H")))
   expect_silent(tf(type = "l", xlas = 2, xticks = xticks, xticklabels = xticklabels))
   expect_silent(tf(type = "b", xlas = 2, xticks = xticks, xticklabels = xticklabels))
   expect_silent(tf(type = "h", xlas = 2, xticks = xticks, xticklabels = xticklabels))
})


yticks <- c(100, 200, 300, 400, 500)
yticklabels <- c("XS", "S", "M", "L", "XL")

test_that("manual yticks work correctly", {
   expect_silent(tf(type = "p", yticks = c(50, 70, 100)))
   expect_silent(tf(type = "l", yticks = yticks))
   expect_silent(tf(type = "b", yticks = yticks))
   expect_silent(tf(type = "h", yticks = yticks))
})

test_that("'yticklabels' can not be specified without 'xticks'", {
   expect_error(tf(type = "l", yticklabels = yticklabels))
   expect_error(tf(type = "b", yticklabels = yticklabels))
   expect_error(tf(type = "h", yticklabels = yticklabels))
})


test_that("'yticks' and 'yticklabels' work correctly together", {
   expect_silent(tf(type = "p", yticks = c(50, 70, 100), yticklabels = c("L", "M", "H")))
   expect_silent(tf(type = "l", yticks = yticks, yticklabels = yticklabels))
   expect_silent(tf(type = "b", yticks = yticks, yticklabels = yticklabels))
   expect_silent(tf(type = "h", yticks = yticks, yticklabels = yticklabels))
})

test_that("'ylas' parameter is handling correctly", {
   expect_silent(tf(type = "p", ylas = 2, yticks = c(50, 70, 100), yticklabels = c("L", "M", "H")))
   expect_silent(tf(type = "l", ylas = 2, yticks = yticks, yticklabels = yticklabels))
   expect_silent(tf(type = "b", ylas = 2, yticks = yticks, yticklabels = yticklabels))
   expect_silent(tf(type = "h", ylas = 2, yticks = yticks, yticklabels = yticklabels))
})



#####################################
# Block 8. Color groups and labels  #
#####################################

context("mdaplot: color groups and labels")


test_that("row names are used as labels by default", {
   expect_silent(tf(type = "p", show.labels = T))
   expect_silent(tf(type = "l", show.labels = T))
   expect_silent(tf(type = "b", show.labels = T))
   expect_silent(tf(type = "h", show.labels = T))
})

test_that("in case of vector, numbers or names are used as labels", {
   expect_silent(mdaplot(people[1, ], type = "p", show.labels = T))
   expect_silent(mdaplot(people[1, ], type = "l", show.labels = T))
   expect_silent(mdaplot(people[1, ], type = "b", show.labels = T))
   expect_silent(mdaplot(people[1, ], type = "h", show.labels = T))
})

test_that("user can specify the labels manually", {
   expect_silent(tf(type = "p", show.labels = T, labels = people[, 1]))
   expect_silent(tf(type = "l", show.labels = T, labels = people[1, ]))
   expect_silent(tf(type = "b", show.labels = T, labels = people[1, ]))
   expect_silent(tf(type = "h", show.labels = T, labels = people[1, ]))
})

test_that("user can specify the labels as 'indices'", {
   expect_silent(tf(type = "p", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "l", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "b", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "h", show.labels = T, labels = "indices"))
})

test_that("user can specify the labels as 'values'", {
   expect_silent(tf(type = "p", show.labels = T, labels = "values"))
   expect_silent(tf(type = "l", show.labels = T, labels = "values"))
   expect_silent(tf(type = "b", show.labels = T, labels = "values"))
   expect_silent(tf(type = "h", show.labels = T, labels = "values"))
})

test_that("user can specify the labels as 'names'", {
   expect_silent(tf(type = "p", show.labels = T, labels = "names"))
   expect_silent(tf(type = "l", show.labels = T, labels = "names"))
   expect_silent(tf(type = "b", show.labels = T, labels = "names"))
   expect_silent(tf(type = "h", show.labels = T, labels = "names"))
})


test_that("user can use color grouping by numeric values", {
   expect_silent(tf(type = "p", cgroup = people[, 1]))
   expect_silent(tf(type = "l", cgroup = people[, 6]))
   expect_silent(tf(type = "b", cgroup = people[, 3]))
   expect_silent(tf(type = "h", cgroup = people[1, ]))
})

cgroup <- factor(people[, 8], labels = c("Males", "Females"))
test_that("user can use color grouping by factors", {
   expect_silent(tf(type = "p", cgroup = cgroup))
   expect_silent(tf(type = "l", cgroup = cgroup))
   expect_silent(tf(type = "b", cgroup = cgroup))
   expect_silent(tf(type = "h", cgroup = factor(rep(c(1, 2), 6), labels = c("A", "B"))))
})

test_that("user if both color grouping and color are specified first wins", {
   expect_silent(tf(type = "p", col = "red", cgroup = people[, 1]))
   expect_silent(tf(type = "l", col = "red", cgroup = people[, 6]))
   expect_silent(tf(type = "b", col = "red", cgroup = people[, 3]))
   expect_silent(tf(type = "h", col = "red", cgroup = people[1, ]))
})

test_that("user can hide the colorbar", {
   expect_silent(tf(type = "p", show.colorbar = FALSE, cgroup = people[, 1]))
   expect_silent(tf(type = "l", show.colorbar = FALSE, cgroup = people[, 6]))
   expect_silent(tf(type = "b", show.colorbar = FALSE, cgroup = people[, 3]))
   expect_silent(tf(type = "h", show.colorbar = FALSE, cgroup = people[1, ]))
})

test_that("different color maps can be used", {
   expect_silent(tf(type = "p", cgroup = people[, 1], colmap = "old"))
   expect_silent(tf(type = "l", cgroup = people[, 6], colmap = "gray"))
   expect_silent(tf(type = "b", cgroup = people[, 7], colmap = "jet"))
   expect_silent(tf(type = "h", cgroup = people[1, ], colmap = c("green", "red")))
})


#######################################
# Block 9. Excluding rows and columns  #
#######################################

## only rows
context("mdaplot: can handle excluded rows")

attr(people, "exclrows") <- NULL
attr(people, "exclcols") <- NULL
people <- mda.exclrows(people, people[, "Beer"] > 300)

test_that("excluded values are hidden by default", {
   expect_silent(tf(type = "p", show.labels = T))
   expect_silent(tf(type = "l", show.labels = T))
   expect_silent(tf(type = "b", show.labels = T))
   expect_silent(tf(type = "h", show.labels = T))
})

test_that("hidden excluded values and labels (indices) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "l", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "b", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "h", show.labels = T, labels = "indices"))
})

test_that("hidden excluded values and labels (indices) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "values"))
   expect_silent(tf(type = "l", show.labels = T, labels = "values"))
   expect_silent(tf(type = "b", show.labels = T, labels = "values"))
   expect_silent(tf(type = "h", show.labels = T, labels = "values"))
})

test_that("excluded values can be shown on all plots except bar and errorbar", {
   tf(type = "p", show.excluded = T)
   expect_silent(tf(type = "p", show.excluded = T))
   expect_silent(tf(type = "l", show.excluded = T))
   expect_silent(tf(type = "b", show.excluded = T))
   expect_error(tf(type = "h", show.excluded = T))
   expect_error(tf(type = "e", show.excluded = T))
})

par(mfrow = c(2, 2))
test_that("excluded values and labels (names) work find", {
   expect_silent(tf(type = "p", show.labels = T, show.excluded = T))
   expect_silent(tf(type = "l", show.labels = T, show.excluded = T))
   expect_silent(tf(type = "b", show.labels = T, show.excluded = T))
})

par(mfrow = c(2, 2))
test_that("excluded values and labels (values) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "values", show.excluded = T))
   expect_silent(tf(type = "l", show.labels = T, labels = "values", show.excluded = T))
   expect_silent(tf(type = "b", show.labels = T, labels = "values", show.excluded = T))
})

par(mfrow = c(2, 2))
test_that("excluded values and labels (indices) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "indices", show.excluded = T))
   expect_silent(tf(type = "l", show.labels = T, labels = "indices", show.excluded = T))
   expect_silent(tf(type = "b", show.labels = T, labels = "indices", show.excluded = T))
})

## only columns
context("mdaplot: can handle excluded columns")

attr(people, "exclrows") <- NULL
attr(people, "exclcols") <- NULL
people <- mda.exclcols(people, c("Height", "Beer"))

test_that("excluded columns are always hidden", {
   expect_silent(tf(type = "p", show.labels = T))
   expect_silent(tf(type = "l", show.labels = T))
   expect_silent(tf(type = "b", show.labels = T))
   expect_silent(tf(type = "h", show.labels = T))
})

test_that("excluded columns and labels (indices) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "l", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "b", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "h", show.labels = T, labels = "indices"))
})

test_that("excluded columns and labels (indices) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "values"))
   expect_silent(tf(type = "l", show.labels = T, labels = "values"))
   expect_silent(tf(type = "b", show.labels = T, labels = "values"))
   expect_silent(tf(type = "h", show.labels = T, labels = "values"))
})

test_that("show excluded does not change anything", {
   tf(type = "p", show.excluded = T)
   expect_silent(tf(type = "p", show.excluded = T))
   expect_silent(tf(type = "l", show.excluded = T))
   expect_silent(tf(type = "b", show.excluded = T))
   expect_silent(tf(type = "h", show.excluded = T))
})

## both rows and columns
context("mdaplot: can handle both excluded rows and columns")

attr(people, "exclrows") <- NULL
attr(people, "exclcols") <- NULL
people <- mda.exclrows(people, people[, "Beer"] > 300)
people <- mda.exclcols(people, c("Height", "Beer"))

test_that("excluded values are hidden by default", {
   expect_silent(tf(type = "p", show.labels = T))
   expect_silent(tf(type = "l", show.labels = T))
   expect_silent(tf(type = "b", show.labels = T))
   expect_silent(tf(type = "h", show.labels = T))
})

test_that("hidden excluded values and labels (indices) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "l", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "b", show.labels = T, labels = "indices"))
   expect_silent(tf(type = "h", show.labels = T, labels = "indices"))
})

test_that("hidden excluded values and labels (indices) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "values"))
   expect_silent(tf(type = "l", show.labels = T, labels = "values"))
   expect_silent(tf(type = "b", show.labels = T, labels = "values"))
   expect_silent(tf(type = "h", show.labels = T, labels = "values"))
})

test_that("excluded values can be shown on all plots except bar and errorbar", {
   tf(type = "p", show.excluded = T)
   expect_silent(tf(type = "p", show.excluded = T))
   expect_silent(tf(type = "l", show.excluded = T))
   expect_silent(tf(type = "b", show.excluded = T))
   expect_error(tf(type = "h", show.excluded = T))
   expect_error(tf(type = "e", show.excluded = T))
})

par(mfrow = c(2, 2))
test_that("excluded values and labels (names) work find", {
   expect_silent(tf(type = "p", show.labels = T, show.excluded = T))
   expect_silent(tf(type = "l", show.labels = T, show.excluded = T))
   expect_silent(tf(type = "b", show.labels = T, show.excluded = T))
})

par(mfrow = c(2, 2))
test_that("excluded values and labels (values) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "values", show.excluded = T))
   expect_silent(tf(type = "l", show.labels = T, labels = "values", show.excluded = T))
   expect_silent(tf(type = "b", show.labels = T, labels = "values", show.excluded = T))
})

par(mfrow = c(2, 2))
test_that("excluded values and labels (indices) work fine", {
   expect_silent(tf(type = "p", show.labels = T, labels = "indices", show.excluded = T))
   expect_silent(tf(type = "l", show.labels = T, labels = "indices", show.excluded = T))
   expect_silent(tf(type = "b", show.labels = T, labels = "indices", show.excluded = T))
})


################################
# Block 10. New functionality  #
################################

## color grouping with excluded values
context("mdaplot: color grouping works fine with excluded values")

attr(people, "exclrows") <- NULL
attr(people, "exclcols") <- NULL
people <- mda.exclrows(people, people[, "Beer"] > 300)
people <- mda.exclcols(people, c("Height", "Beer"))

par(mfrow = c(2, 2))
test_that("excluded values and color grouping work fine together", {
   expect_silent(tf(type = "p", cgroup = people$Wine, show.excluded = T, show.labels = T))
   expect_silent(tf(type = "l", cgroup = people$Wine, show.excluded = T))
   expect_silent(tf(type = "b", cgroup = people$Wine, show.excluded = T))
   expect_silent(tf(type = "h", cgroup = rep(c(1, 2), 5), show.excluded = F))
})


## compare different colormaps
context("mdaplot: colormaps, opacity and markers")

x <- rnorm(10000)
y <- rnorm(10000)
d <- sqrt(x^2 + y^2)
d <- 2.5 - d
d[d < 0] <- 0
pd <- cbind(x, y)

par(mfrow = c(3, 2))
test_that("excluded values and color grouping work fine together", {
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: default"))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: old", colmap = "old"))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: gray", colmap = "gray"))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: jet", colmap = "jet"))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: user", colmap = c("red", "blue")))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: user", colmap = c("red", "green", "blue")))
})

par(mfrow = c(2, 2))
test_that("opacity parameter works well", {
   expect_silent(tf(type = "p", opacity = 0.5))
   expect_silent(tf(type = "l", opacity = 0.5))
   expect_silent(tf(type = "b", opacity = 0.5))
   expect_silent(tf(type = "h", opacity = 0.5))
})

par(mfrow = c(2, 2))
test_that("opacity parameter works well together with color grouping", {
   expect_silent(tf(type = "p", cgroup = people[, "Wine"], opacity = 0.5))
   expect_silent(tf(type = "l", cgroup = people[, "Wine"], opacity = 0.5))
   expect_silent(tf(type = "b", cgroup = people[, "Wine"], opacity = 0.5))
   expect_silent(tf(type = "h", cgroup = people[, "Wine"], opacity = 0.5))
})

par(mfrow = c(4, 2))
test_that("different pch values works fine", {
   expect_silent(tf(type = "p", pch = 16))
   expect_silent(tf(type = "p", pch = 21))
   expect_silent(tf(type = "p", pch = 16, col = "red"))
   expect_silent(tf(type = "p", pch = 21, col = "red"))
   expect_silent(tf(type = "p", pch = 16, col = "red"))
   expect_silent(tf(type = "p", pch = 21, bg = "red", col = "blue"))
   expect_silent(tf(type = "p", pch = 16, cgroup = people[, "Wine"]))
   expect_silent(tf(type = "p", pch = 21, cgroup = people[, "Wine"]))
})

par(mfrow = c(2, 2))
test_that("special parameters for pch values work fine", {
   expect_silent(tf(type = "p", pch = 21))
   expect_silent(tf(type = "p", pch = 21, cgroup = people[, "Height"]))
   expect_silent(tf(type = "p", pch = 21, pch.colinv = TRUE))
   expect_silent(tf(type = "p", pch = 21, cgroup = people[, "Height"], 
      bg = "white", lwd = 0.5, cex = 1.2, pch.colinv = TRUE))
})


g1 <- people[, "Sex"]
f1 <- factor(g1, labels = c("M", "F"))
g2 <- people[, "Region"]
f2 <- factor(g2, labels = c("S", "M"))

par(mfrow = c(2, 2))
test_that("color grouping with factors is shown with discrete colorbar", {
   expect_silent(tf(type = "p", cgroup = g1))
   expect_silent(tf(type = "p", cgroup = f1))
   expect_silent(tf(type = "p", cgroup = interaction(f1, f2)))
})

## add convex hull
context("mdaplot: add convex hull")

data(people)
attr(people, "exclrows") <- NULL
attr(people, "exclcols") <- NULL
gw <- people[, "Sex"]
g1 <- factor(people[, "Sex"], labels = c("Male", "Female"))
g2 <- factor(people[, "Region"], labels = c("Scan", "Med"))
g <- interaction(g1, g2)

tf <- function(type, cgroup, ...) {
   p <- mdaplot(people, type = type, cgroup = cgroup)
   add_convex_hull(p, ...)
}

par(mfrow = c(2, 2))
test_that("add_convex_hull returns error if wrong type of plot is used", {
   expect_error(tf(type = "l", cgroup = g1))
   expect_error(tf(type = "b", cgroup = g1))
   expect_error(tf(type = "h", cgroup = g1))
   expect_error(tf(type = "e", cgroup = g1))
})

test_that("add_convex_hull returns error if cgroup is not a factor", {
   expect_error(tf(type = "p", cgroup = gw))
})

par(mfrow = c(2, 2))
test_that("add_convex_hull works as expected with default styles", {
   expect_silent(tf(type = "p", cgroup = g1))
   expect_silent(tf(type = "p", cgroup = g2))
   expect_silent(tf(type = "p", cgroup = g))
})

par(mfrow = c(2, 2))
test_that("add_convex_hull works as expected with different opacity", {
   expect_silent(tf(type = "p", opacity = 0.8, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.4, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, cgroup = g2))
   expect_silent(tf(type = "p", opacity = 0.1, cgroup = g))
})

par(mfrow = c(2, 2))
test_that("add_convex_hull works as expected with different line parameters", {
   expect_silent(tf(type = "p", opacity = 0.4, lwd = 1, lty = 1, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 2, lty = 2, cgroup = g2))
   expect_silent(tf(type = "p", opacity = 0.1, lwd = 3, lty = 3, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.1, lwd = 3, lty = 4, cgroup = g))
})

## add convex hull for excluded data
context("mdaplot: add convex hull for excluded data")

data(people)
attr(people, "exclrows") <- NULL
attr(people, "exclcols") <- NULL
people <- mda.exclrows(people, people[, "Beer"] > 300)
people <- mda.exclcols(people, c("Height", "Beer"))

gw <- people[, "Sex"]
g1 <- factor(people[, "Sex"], labels = c("Male", "Female"))
g2 <- factor(people[, "Region"], labels = c("Scan", "Med"))
g <- interaction(g1, g2)

p <- mdaplot(people, type = "p", cgroup = g1)
add_convex_hull(p)

tf <- function(type, cgroup, ...) {
   p <- mdaplot(people, type = type, cgroup = cgroup)
   add_convex_hull(p, ...)
}

par(mfrow = c(2, 2))
test_that("add_convex_hull returns error if wrong type of plot is used", {
   expect_error(tf(type = "l", cgroup = g1))
   expect_error(tf(type = "b", cgroup = g1))
   expect_error(tf(type = "h", cgroup = g1))
   expect_error(tf(type = "e", cgroup = g1))
})

test_that("add_convex_hull returns error if cgroup is not a factor", {
   expect_error(tf(type = "p", cgroup = gw))
})

par(mfrow = c(2, 2))
test_that("add_convex_hull works as expected with default styles", {
   expect_silent(tf(type = "p", cgroup = g1))
   expect_silent(tf(type = "p", cgroup = g2))
   expect_silent(tf(type = "p", cgroup = g))
})

par(mfrow = c(2, 2))
test_that("add_convex_hull works as expected with different opacity", {
   expect_silent(tf(type = "p", opacity = 0.8, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.4, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, cgroup = g2))
   expect_silent(tf(type = "p", opacity = 0.1, cgroup = g))
})

par(mfrow = c(2, 2))
test_that("add_convex_hull works as expected with different line parameters", {
   expect_silent(tf(type = "p", opacity = 0.4, lwd = 1, lty = 1, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 2, lty = 2, cgroup = g2))
   expect_silent(tf(type = "p", opacity = 0.1, lwd = 3, lty = 3, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.1, lwd = 3, lty = 4, cgroup = g))
})

## add confidence ellipse
context("mdaplot: add confidence ellipse")

data(people)
attr(people, "exclrows") <- NULL
attr(people, "exclcols") <- NULL
gw <- people[, "Sex"]
g1 <- factor(people[, "Sex"], labels = c("Male", "Female"))
g2 <- factor(people[, "Region"], labels = c("Scan", "Med"))
g <- interaction(g1, g2)

tf <- function(type, cgroup, ...) {
   p <- mdaplot(people, type = type, cgroup = cgroup)
   add_confidence_ellipse(p, ...)
}

par(mfrow = c(2, 2))
test_that("add_confidence_ellipse returns error if wrong type of plot is used", {
   expect_error(tf(type = "l", cgroup = g1))
   expect_error(tf(type = "b", cgroup = g1))
   expect_error(tf(type = "h", cgroup = g1))
   expect_error(tf(type = "e", cgroup = g1))
})

test_that("add_confidence_ellipse returns error if cgroup is not a factor", {
   expect_error(tf(type = "p", cgroup = gw))
})

par(mfrow = c(2, 2))
test_that("add_confidence_ellipse works as expected with default styles", {
   expect_silent(tf(type = "p", cgroup = g1))
   expect_silent(tf(type = "p", cgroup = g2))
   expect_silent(tf(type = "p", cgroup = g))
})

par(mfrow = c(2, 2))
test_that("add_confidence_ellipse works as expected with different opacity", {
   expect_silent(tf(type = "p", opacity = 0.8, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.4, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.2, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.1, cgroup = g))
})

par(mfrow = c(2, 2))
test_that("add_confidence_ellipse works as expected with different confidence level", {
   expect_silent(tf(type = "p", opacity = 0.8, conf.level = 0.80, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.4, conf.level = 0.90, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.2, conf.level = 0.95, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.1, conf.level = 0.99, cgroup = g))
})

par(mfrow = c(2, 2))
test_that("add_confidence_ellipse works as expected with different line parameters", {
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 1, lty = 1, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 2, lty = 2, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 3, lty = 3, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 3, lty = 4, cgroup = g1))
})

