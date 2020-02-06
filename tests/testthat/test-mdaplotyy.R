
#####################################################
# Tests for basic functionality of mdaplotyy()      #
#####################################################

setup({
   pdf(file = tempfile("mdatools-test-mdaplotyy-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-mdaplotyy-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

data(people)

context("mdaplotyy: basic functionality")

plot_data <- t(mda.subset(people, select = c("Height", "Income")))
attr(plot_data, "name") <- "People"
attr(plot_data, "xaxis.name") <- "Persons"

wrong_plot_data <- t(mda.subset(people, select = c("Height", "IQ", "Income")))

test_that("mdaplotyy returns error if some parameters are wrong", {
   expect_error(mdaplotyy(wrong_plot_data))
   expect_error(mdaplotyy(plot_data, type = "p"))
   expect_error(mdaplotyy(plot_data, type = "h"))
   expect_error(mdaplotyy(plot_data, col = "red"))
   expect_error(mdaplotyy(plot_data, lty = 1))
   expect_error(mdaplotyy(plot_data, lwd = 1))
   expect_error(mdaplotyy(plot_data, pch = 1))
   expect_error(mdaplotyy(plot_data, ylab = "xxx"))
})

test_that("mdaplotyy works fine with basic settings", {
   par(mfrow = c(2, 2))
   expect_silent(mdaplotyy(plot_data))
   expect_silent(mdaplotyy(plot_data, col = c("red", "green")))
   expect_silent(mdaplotyy(plot_data, type = "b", pch = c(1, 2)))
   expect_silent(mdaplotyy(plot_data, type = "b", lty = c(2, 3), lwd = c(1, 2)))
})

test_that("labels are handled correctly", {
   par(mfrow = c(2, 2))
   expect_silent(mdaplotyy(plot_data, show.labels = T))
   expect_silent(mdaplotyy(plot_data, show.labels = T, labels = "values"))
   expect_silent(mdaplotyy(plot_data, show.labels = T, labels = "indices", type = "b"))
   expect_silent(mdaplotyy(plot_data, show.labels = T, labels = "names", type = "b"))
})

test_that("labels, titles and ticks can be changed", {
   xticks <- seq(1, 32, by = 7)
   xticklabels <- rownames(people)[xticks]
   par(mfrow = c(2, 2))
   expect_silent(mdaplotyy(plot_data, main = "Data", xlab = "Rows", ylab = c("Y1", "Y2")))
   expect_silent(mdaplotyy(plot_data, xticks = xticks))
   expect_silent(mdaplotyy(plot_data, xticks = xticks, xticklabels = xticklabels))
   expect_silent(mdaplotyy(plot_data, xticks = xticks, xticklabels = xticklabels, ylas = 2,
      xlas = 2, xlab = "", ylab = c("", ""), legend = c("Y1", "Y2")))
})

test_that("limits can be changed", {
   par(mfrow = c(2, 2))
   expect_silent(mdaplotyy(plot_data, xlim = c(0, 40)))
   expect_silent(mdaplotyy(plot_data, ylim = list(c(0, 200), NULL)))
   expect_silent(mdaplotyy(plot_data, ylim = list(c(0, 200), c(0, 60000))))
})

test_that("legend and its position can be changed", {
   par(mfrow = c(2, 2))
   expect_silent(mdaplotyy(plot_data, show.legend = FALSE))
   expect_silent(mdaplotyy(plot_data, legend = c("11", "22")))
   expect_silent(mdaplotyy(plot_data, legend.position = "top"))
   expect_silent(mdaplotyy(plot_data, legend.position = "bottom"))
})
