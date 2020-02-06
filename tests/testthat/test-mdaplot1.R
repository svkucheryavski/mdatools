###################################################################
# Tests for basic functionality of mdaplot() and related methods  #
###################################################################

setup({
   pdf(file = tempfile("mdatools-test-mdaplot1-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-mdaplot1-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# prepare dataset
all_plots <- c("p", "l", "b", "h", "e")
nonscatter_plots <- c("l", "b", "h", "e")
main_plots <- c("p", "l", "b", "h")
linebar_plots <- c("l", "b", "h")

par(mfrow = c(2, 2))

data("people")
people <- people[, -6]

m <- apply(people, 2, mean)
s <- apply(people, 2, sd)
stat1 <- rbind(m, s)
stat2 <- rbind(m, s, 1.5 * s)
colnames(stat1) <- colnames(people)
colnames(stat2) <- colnames(people)

####################################################
# Block 1. General tests for mdaplot functionality #
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
   m <- apply(people, 2, mean)
   expect_silent(mdaplot(m, type = "h", col = "pink", bwd = 0.25, show.labels = T))
   expect_silent(mdaplot(stat1, type = "e", col = "red", show.axes = F, show.grid = F))
})
par(mfrow = c(2, 2))

######################################
# Block 2. Handling data attributes  #
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

attr(spectra, "xaxis.values") <- 10^7 / simdata$wavelength
attr(spectra, "xaxis.name") <- expression("Wavenumber, cm"^-1)
test_that("attribute 'xaxis.values' can be in reverse order", {
   expect_silent(mdaplot(spectra, type = "p"))
   expect_silent(mdaplot(spectra, type = "l"))
   expect_silent(mdaplot(spectra, type = "b"))
   expect_silent(mdaplot(spectra, type = "h"))
})


#########################################
# Block 3. Manual ticks and ticklabels  #
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
   xticklabels2 <- c("L", "M", "H")
   expect_silent(tf(type = "p", xlas = 2, xticks = c(165, 180, 195), xticklabels = xticklabels2))
   expect_silent(tf(type = "l", xlas = 2, xticks = xticks, xticklabels = xticklabels))
   expect_silent(tf(type = "b", xlas = 2, xticks = xticks, xticklabels = xticklabels))
   expect_silent(tf(type = "h", xlas = 2, xticks = xticks, xticklabels = xticklabels))
})


yticks <- c(100, 200, 300, 400, 500)
yticklabels <- c("XS", "S", "M", "L", "XL")

test_that("manual yticks work correctly", {
   expect_silent(tf(type = "p", yticks = c(50, 70, 90)))
   expect_silent(tf(type = "l", yticks = yticks))
   expect_silent(tf(type = "b", yticks = yticks))
   expect_silent(tf(type = "h", yticks = yticks))
})

test_that("'yticklabels' can not be specified without 'yticks'", {
   expect_error(tf(type = "l", yticklabels = yticklabels))
   expect_error(tf(type = "b", yticklabels = yticklabels))
   expect_error(tf(type = "h", yticklabels = yticklabels))
})


test_that("'yticks' and 'yticklabels' work correctly together", {
   expect_silent(tf(type = "p", yticks = c(50, 70, 90), yticklabels = c("L", "M", "H")))
   expect_silent(tf(type = "l", yticks = yticks, yticklabels = yticklabels))
   expect_silent(tf(type = "b", yticks = yticks, yticklabels = yticklabels))
   expect_silent(tf(type = "h", yticks = yticks, yticklabels = yticklabels))
})

test_that("'ylas' parameter is handling correctly", {
   expect_silent(tf(type = "p", ylas = 2, yticks = c(50, 70, 90), yticklabels = c("L", "M", "H")))
   expect_silent(tf(type = "l", ylas = 2, yticks = yticks, yticklabels = yticklabels))
   expect_silent(tf(type = "b", ylas = 2, yticks = yticks, yticklabels = yticklabels))
   expect_silent(tf(type = "h", ylas = 2, yticks = yticks, yticklabels = yticklabels))
})


#####################################
# Block 4. Color groups and labels  #
#####################################

context("mdaplot: using labels")
tf <- function(...) mdaplot(people, ...)

test_that("row names are used as labels by default", {
   expect_silent(tf(type = "p", show.labels = T))
   expect_silent(tf(type = "l", show.labels = T))
   expect_silent(tf(type = "b", show.labels = T))
   expect_silent(tf(type = "h", show.labels = T))
})

test_that("in case of vector, numbers or names are used as labels", {
   expect_silent(mdaplot(people[, 1, drop = F], type = "p", show.labels = T))
   expect_silent(mdaplot(people[1, ,drop = F], type = "l", show.labels = T))
   expect_silent(mdaplot(people[1, ,drop = F], type = "b", show.labels = T))
   expect_silent(mdaplot(people[1, ,drop = F], type = "h", show.labels = T))
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


context("mdaplot: using color grouping")

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
   expect_silent(tf(type = "h", cgroup = factor(rep(c(1, 2), 6), labels = c("A", "B"))[1:11]))
})


test_that("user if both color grouping and color are specified color wins", {
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
# Block 5. Excluding rows and columns  #
#######################################

## only rows
context("mdaplot: can handle excluded rows")

attr(people, "exclrows") <- NULL
attr(people, "exclcols") <- NULL
people <- mda.exclrows(people, people[, "Beer"] > 300)
tf <- function(...) mdaplot(people, ...)

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

par(mfrow = c(2, 2))
test_that("excluded values can be shown on all plots except errorbar", {
   expect_silent(tf(type = "p", show.excluded = T))
   expect_silent(tf(type = "l", show.excluded = T))
   expect_silent(tf(type = "b", show.excluded = T))
   expect_silent(tf(type = "h", show.excluded = T))
   expect_error(tf(type = "e", show.excluded = T))
})

par(mfrow = c(2, 2))
test_that("excluded values and labels (names) work fine", {
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

par(mfrow = c(2, 2))
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
   expect_silent(tf(type = "p", show.excluded = T))
   expect_silent(tf(type = "l", show.excluded = T))
   expect_silent(tf(type = "b", show.excluded = T))
   expect_silent(tf(type = "h", show.excluded = T))
   expect_error(tf(type = "e", show.excluded = T))
})

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
