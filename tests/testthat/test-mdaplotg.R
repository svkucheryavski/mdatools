################################################################
# Block 1: using matrix as data source plus groupby parameter  #
################################################################

setup({
   pdf(file = tempfile("mdatools-test-mdaplotg-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-mdaplotg-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

context("mdaplotg: plots with matrix as data and groupby parameter")
par(mfrow = c(2, 2))
data(people)
data <- people[, -6]
groupby <- as.factor(people[, "Sex"])
groupby_wrong <- people[, "Sex"]

groupby_df <- as.data.frame(people[, c("Sex", "Region")])
groupby_df$Sex <- factor(groupby_df$Sex, labels = c("M", "F"))
groupby_df$Region <- factor(groupby_df$Region, labels = c("S", "M"))

## shortcut for plotting function
tf <- function(groupby, ...) mdaplotg(data, groupby = groupby, ...)

# simple plots with groupby parameter

test_that("providing correct 'groupy' parameter works for any plot", {
   expect_silent(tf(groupby, type = "p"))
   expect_silent(tf(groupby, type = "l"))
   expect_silent(tf(groupby, type = "b"))
   expect_silent(tf(groupby, type = "h"))
})

test_that("providing 'groupy' as data frame parameter works for any plot", {
   expect_silent(tf(groupby_df, type = "p"))
   expect_silent(tf(groupby_df, type = "l"))
   expect_silent(tf(groupby_df, type = "b"))
   expect_silent(tf(groupby_df, type = "h"))
})

test_that("providing wrong 'groupy' parameter raises error for any plot", {
   expect_error(tf(groupby_wrong, type = "p"))
   expect_error(tf(groupby_wrong, type = "l"))
   expect_error(tf(groupby_wrong, type = "b"))
   expect_error(tf(groupby_wrong, type = "h"))
})

## checking labels

tfl <- function(groupby, ...) {
   mdaplotg(data, groupby = groupby, showl.labels = TRUE, ...)
}

par(mfrow = c(2, 2))
test_that("parameter 'show.labels' is accepted", {
   expect_silent(tf(groupby, type = "p", show.labels = TRUE))
   expect_silent(tf(groupby, type = "l", show.labels = TRUE))
   expect_silent(tf(groupby, type = "b", show.labels = TRUE))
   expect_silent(tf(groupby, type = "h", show.labels = TRUE))
})

labels <- "values"
test_that("values can be used as labels", {
   expect_silent(tf(groupby, type = "p", show.labels = TRUE, labels = labels))
   expect_silent(tf(groupby, type = "l", show.labels = TRUE, labels = labels))
   expect_silent(tf(groupby, type = "b", show.labels = TRUE, labels = labels))
   expect_silent(tf(groupby, type = "h", show.labels = TRUE, labels = labels))
})

labels <- "indices"
test_that("indices can be used as labels", {
   expect_silent(tf(groupby, type = "p", show.labels = TRUE, labels = labels))
   expect_silent(tf(groupby, type = "l", show.labels = TRUE, labels = labels))
   expect_silent(tf(groupby, type = "b", show.labels = TRUE, labels = labels))
   expect_silent(tf(groupby, type = "h", show.labels = TRUE, labels = labels))
})

## checking legend

test_that("legend can be hidden or have another position", {
   expect_silent(tf(groupby, type = "p", show.legend = FALSE))
   expect_silent(tf(groupby, type = "l", legend.position = "topleft"))
   expect_silent(tf(groupby, type = "b", legend.position = "top"))
   expect_silent(tf(groupby, type = "h", legend.position = "bottomleft"))
})

legend <- c("a", "b")
test_that("legend values can be modified", {
   expect_silent(tf(groupby, type = "p", legend = legend))
   expect_silent(tf(groupby, type = "l", legend = legend))
   expect_silent(tf(groupby, type = "b", legend = legend))
   expect_silent(tf(groupby, type = "h", legend = legend))
})

wrong_legend <- c("a", "b", "c", "d")
test_that("if number of provided legend values is wrong - error raises", {
   expect_error(tf(groupby, type = "p", legend = wrong_legend))
   expect_error(tf(groupby, type = "l", legend = wrong_legend))
   expect_error(tf(groupby, type = "b", legend = wrong_legend))
   expect_error(tf(groupby, type = "h", legend = wrong_legend))
})

## checking paremeters for plot series

test_that("'cex' can be specified as a single value or individual for each group", {
   expect_silent(tf(groupby, cex = 2, type = "p"))
   expect_silent(tf(groupby, cex = 2, type = "b"))
   expect_silent(tf(groupby, cex = 1:2, type = "p"))
   expect_silent(tf(groupby, cex = 1:2, type = "b"))
   expect_error( tf(groupby, cex = 1:3, type = "b"))
})

test_that("'pch' can be specified as a single value or individual for each group", {
   expect_silent(tf(groupby, pch = 2, type = "p"))
   expect_silent(tf(groupby, pch = 2, type = "b"))
   expect_silent(tf(groupby, pch = 1:2, type = "p"))
   expect_silent(tf(groupby, pch = 1:2, type = "b"))
   expect_error( tf(groupby, pch = 1:3, type = "b"))
})

test_that("'lty' can be specified as a single value or individual for each group", {
   expect_silent(tf(groupby, lty = 1, type = "l"))
   expect_silent(tf(groupby, lty = 1, type = "b"))
   expect_silent(tf(groupby, lty = 1:2, type = "l"))
   expect_silent(tf(groupby, lty = 1:2, type = "b"))
   expect_error( tf(groupby, lty = 1:3, type = "l"))
})

test_that("'lwd' can be specified as a single value or individual for each group", {
   expect_silent(tf(groupby, lwd = 1, type = "l"))
   expect_silent(tf(groupby, lwd = 1, type = "b"))
   expect_silent(tf(groupby, lwd = 1:2, type = "l"))
   expect_silent(tf(groupby, lwd = 1:2, type = "b"))
   expect_error( tf(groupby, lwd = 1:3, type = "l"))
})

test_that("'col' can be specified individual for each group", {
   expect_silent(tf(groupby, col = c("red", "green"), type = "p"))
   expect_silent(tf(groupby, col = c("red", "green"), type = "h"))
   expect_silent(tf(groupby, col = c("red", "green"), type = "l"))
   expect_silent(tf(groupby, col = c("red", "green"), type = "b"))
   expect_error( tf(groupby, col = c("red", "green", "blue"), type = "p"))
})

test_that("'col' can be specified as a single value", {
   expect_silent(tf(groupby, col = c("red"), type = "p", pch = 1:2))
   expect_silent(tf(groupby, col = c("red"), type = "h", border = c("green", "blue")))
   expect_silent(tf(groupby, col = c("red"), type = "l", lty = 1:2))
   expect_silent(tf(groupby, col = c("red"), type = "b", lty = 1:2, pch = 1:2))
})


## checking ticks and ticklabels for x-axis

xticks   <- c(2, 4, 6, 8, 10, 12)
xticks_p <- c(150, 170, 190)
xticklabels   <- colnames(people)[xticks]
xticklabels_p <- c("Small", "Normal", "Large")

par(mfrow = c(2, 2))

test_that("'xticks' can be user defined", {
   expect_silent(tf(groupby, type = "l", xticks = xticks))
   expect_silent(tf(groupby, type = "b", xticks = xticks))
   expect_silent(tf(groupby, type = "h", xticks = xticks))
   expect_silent(tf(groupby, type = "p", xticks = xticks_p))
})

test_that("'xticklabels' can be user defined", {
   expect_silent(tf(groupby, type = "l", xticks = xticks, xticklabels = xticklabels))
   expect_silent(tf(groupby, type = "b", xticks = xticks, xticklabels = xticklabels))
   expect_silent(tf(groupby, type = "h", xticks = xticks, xticklabels = xticklabels))
   expect_silent(tf(groupby, type = "p", xticks = xticks_p, xticklabels = xticklabels_p))
})

test_that("'xticklabels' can not be used without 'xticks' or if they have different length", {
   expect_error(tf(groupby, type = "l", xticklabels = xticklabels))
   expect_error(tf(groupby, type = "p", xticklabels = xticklabels_p))
   expect_error(tf(groupby, type = "l", xticks = 1:10, xticklabels = xticklabels))
   expect_error(tf(groupby, type = "p", xticks = 1:10, xticklabels = xticklabels_p))
})

## checking ticks and ticklabels for y-axis

yticks   <- c(100, 300)
yticks_p <- c(50, 70, 90)
yticklabels   <- c("Mini", "Mega")
yticklabels_p <- c("Small", "Normal", "Large")

test_that("'yticks' can be user defined", {
   expect_silent(tf(groupby, type = "l", yticks = yticks))
   expect_silent(tf(groupby, type = "b", yticks = yticks))
   expect_silent(tf(groupby, type = "h", yticks = yticks))
   expect_silent(tf(groupby, type = "p", yticks = yticks_p))
})

test_that("'yticklabels' can be user defined", {
   expect_silent(tf(groupby, type = "b", yticks = yticks,   yticklabels = yticklabels))
   expect_silent(tf(groupby, type = "l", yticks = yticks,   yticklabels = yticklabels))
   expect_silent(tf(groupby, type = "h", yticks = yticks,   yticklabels = yticklabels))
   expect_silent(tf(groupby, type = "p", yticks = yticks_p, yticklabels = yticklabels_p))
})

test_that("'yticklabels' can not be used without 'yticks' or if they have different length", {
   expect_error(tf(groupby, type = "l", yticklabels = yticklabels))
   expect_error(tf(groupby, type = "p", yticklabels = yticklabels_p))
   expect_error(tf(groupby, type = "l", yticks = 1:10, yticklabels = yticklabels))
   expect_error(tf(groupby, type = "p", yticks = 1:10, yticklabels = yticklabels_p))
})


## using different colormaps and opacity colors

test_that("different colormaps can be used", {
   expect_silent(tf(groupby_df, type = "p", ))
   expect_silent(tf(groupby_df, type = "p", colmap = "old"))
   expect_silent(tf(groupby_df, type = "p", colmap = "gray"))
   expect_silent(tf(groupby_df, type = "p", colmap = "jet"))
})

test_that("different colormaps can be used with single opacity value", {
   expect_silent(tf(groupby_df, opacity = 0.5, type = "p", ))
   expect_silent(tf(groupby_df, opacity = 0.5, type = "p", colmap = "old"))
   expect_silent(tf(groupby_df, opacity = 0.5, type = "p", colmap = "gray"))
   expect_silent(tf(groupby_df, opacity = 0.5, type = "p", colmap = "jet"))
})

test_that("different colormaps can be used with individual opacity values", {
   expect_silent(tf(groupby_df, opacity = c(0.2, 0.4, 0.6, 0.8), type = "p", ))
   expect_silent(tf(groupby_df, opacity = c(0.2, 0.4, 0.6, 0.8), type = "p", colmap = "old"))
   expect_silent(tf(groupby_df, opacity = c(0.2, 0.4, 0.6, 0.8), type = "p", colmap = "gray"))
   expect_silent(tf(groupby_df, opacity = c(0.2, 0.4, 0.6, 0.8), type = "p", colmap = "jet"))
})

## test grid parameters
test_that("grid can be shown", {
   expect_silent(tf(groupby, type = "p", show.grid = TRUE))
   expect_silent(tf(groupby, type = "l", show.grid = TRUE))
   expect_silent(tf(groupby, type = "h", show.grid = TRUE))
   expect_silent(tf(groupby, type = "b", show.grid = TRUE))
})

test_that("grid can be hidden", {
   expect_silent(tf(groupby, type = "p", show.grid = FALSE))
   expect_silent(tf(groupby, type = "l", show.grid = FALSE))
   expect_silent(tf(groupby, type = "h", show.grid = FALSE))
   expect_silent(tf(groupby, type = "b", show.grid = FALSE))
})

test_that("grid color and thickness can be changed", {
   expect_silent(tf(groupby, type = "p", show.grid = TRUE, grid.col = "red", grid.lwd = "2"))
   expect_silent(tf(groupby, type = "l", show.grid = TRUE, grid.col = "red", grid.lwd = "2"))
   expect_silent(tf(groupby, type = "h", show.grid = TRUE, grid.col = "red", grid.lwd = "2"))
   expect_silent(tf(groupby, type = "b", show.grid = TRUE, grid.col = "red", grid.lwd = "2"))
})


## handling hidden columns

data(people)
data <- people[, -6]
data <- mda.exclcols(data, "Height")

test_that("hidden data is not shown by default", {
   expect_silent(tf(groupby, type = "p"))
   expect_silent(tf(groupby, type = "l"))
   expect_silent(tf(groupby, type = "b"))
   expect_silent(tf(groupby, type = "h"))
})

test_that("hidden data is shown as gray if needed and labels are produced correctly", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T))
})

test_that("hidden data can be used with labels as values", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T, labels = "values"))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T, labels = "values"))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T, labels = "values"))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T, labels = "values"))
})

test_that("hidden data can be used with labels as names", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T, labels = "names"))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T, labels = "names"))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T, labels = "names"))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T, labels = "names"))
})


test_that("hidden data can be used with labels as indices", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T, labels = "indices"))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T, labels = "indices"))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T, labels = "indices"))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T, labels = "indices"))
})


## handling hidden rows

data(people)
data <- people[, -6]
data <- mda.exclrows(data, data[, "Beer"] > 300)

test_that("hidden data is not shown by default", {
   expect_silent(tf(groupby, type = "p"))
   expect_silent(tf(groupby, type = "l"))
   expect_silent(tf(groupby, type = "b"))
   expect_silent(tf(groupby, type = "h"))
})

test_that("hidden data is shown as gray if needed and labels are produced correctly", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T))
})

test_that("hidden data can be used with labels as values", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T, labels = "values"))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T, labels = "values"))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T, labels = "values"))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T, labels = "values"))
})

test_that("hidden data can be used with labels as names", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T, labels = "names"))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T, labels = "names"))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T, labels = "names"))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T, labels = "names"))
})

test_that("hidden data can be used with labels as indices", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T, labels = "indices"))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T, labels = "indices"))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T, labels = "indices"))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T, labels = "indices"))
})

## handling both hidden columns and rows

data(people)
data <- people[, -6]
data <- mda.exclrows(data, data[, "Beer"] > 300)
data <- mda.exclcols(data, "Height")
tf(groupby, type = "h")

test_that("hidden data is not shown by default", {
   expect_silent(tf(groupby, type = "p"))
   expect_silent(tf(groupby, type = "l"))
   expect_silent(tf(groupby, type = "b"))
   expect_silent(tf(groupby, type = "h"))
})

test_that("hidden data is shown as gray if needed and labels are produced correctly", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T))
})

test_that("hidden data can be used with labels as values", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T, labels = "values"))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T, labels = "values"))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T, labels = "values"))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T, labels = "values"))
})

test_that("hidden data can be used with labels as names", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T, labels = "names"))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T, labels = "names"))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T, labels = "names"))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T, labels = "names"))
})

test_that("hidden data can be used with labels as indices", {
   expect_silent(tf(groupby, type = "p", show.excluded = TRUE, show.labels = T, labels = "indices"))
   expect_silent(tf(groupby, type = "l", show.excluded = TRUE, show.labels = T, labels = "indices"))
   expect_silent(tf(groupby, type = "b", show.excluded = TRUE, show.labels = T, labels = "indices"))
   expect_silent(tf(groupby, type = "h", show.excluded = TRUE, show.labels = T, labels = "indices"))
})

#######################################
# Block 2: using list as data source  #
#######################################

# TODO: extend the tests for block 2

context("mdaplotg: plots with list as data")

par(mfrow = c(2, 2))
data(people)
males <- people[people[, "Sex"] == -1, -6]
females <- people[people[, "Sex"] == 1, -6]
data <- list("males" = males, "females" = females)

## shortcut for plotting function
tf <- function(...) mdaplotg(data, ...)

# simple plots with groupby parameter

test_that("providing data as list works for any plot", {
   expect_silent(tf(type = "p"))
   expect_silent(tf(type = "l"))
   expect_silent(tf(type = "b"))
   expect_silent(tf(type = "h"))
})


#######################""################
# Block 3: using matrix as data source  #
#########################""##############

# TODO: extend the tests for block 3

context("mdaplotg: plots with matrix as data without groupby")

## shortcut for plotting function
tf <- function(...) mdaplotg(people[1:4, -6], ...)

test_that("providing data as matrix works for line and bar plots only", {
   expect_silent(tf(type = "l"))
   expect_silent(tf(type = "b"))
   expect_silent(tf(type = "h"))
   expect_error(tf(type = "p"))
})

