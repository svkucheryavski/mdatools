#############################################################
# Tests for new features for mdaplot() and related methods  #
#############################################################

setup({
   pdf(file = tempfile("mdatools-test-mdaplot2-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-mdaplot2-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

par(mfrow = c(2, 2))

data("people")
people <- people[, -6]
people <- as.data.frame(people)
tf <- function(...) mdaplot(people, ...)

#############################################################
#  Block 1: colormaps and color groupins                    #
#############################################################

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
   expect_silent(tf(type = "h", cgroup = rep(c(1, 2), 6)[1:11], show.excluded = F))
})

## compare different colormaps
context("mdaplot: colormaps, opacity and markers")

x <- rnorm(10000)
y <- rnorm(10000)
d <- sqrt(x^2 + y^2)
d <- 2.5 - d
d[d < 0] <- 0
pd <- cbind(x, y)

usr_cmap1 <- c("red", "blue")
usr_cmap2 <- c("red", "green", "blue")

par(mfrow = c(3, 2))
test_that("excluded values and color grouping work fine together", {
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: default"))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: old", colmap = "old"))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: gray", colmap = "gray"))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: jet", colmap = "jet"))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: user", colmap = usr_cmap1))
   expect_silent(mdaplot(pd, cgroup = d, main = "Colormap: user", colmap = usr_cmap2))
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
   expect_silent(tf(type = "h", cgroup = people[1, ], opacity = 0.5))
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

#############################################################
#  Block 2: convex hull and confidence ellipses             #
#############################################################

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
   plotConvexHull(p, ...)
}

par(mfrow = c(2, 2))
test_that("plotConvexHull() returns error if wrong type of plot is used", {
   expect_error(tf(type = "l", cgroup = g1))
   expect_error(tf(type = "b", cgroup = g1))
   expect_error(tf(type = "h", cgroup = g1))
   expect_error(tf(type = "e", cgroup = g1))
})

test_that("plotConvexHull() returns error if cgroup is not a factor", {
   expect_error(tf(type = "p", cgroup = gw))
})

par(mfrow = c(2, 2))
test_that("plotConvexHull() works as expected with default styles", {
   expect_silent(tf(type = "p", cgroup = g1))
   expect_silent(tf(type = "p", cgroup = g2))
   expect_silent(tf(type = "p", cgroup = g))
})

par(mfrow = c(2, 2))
test_that("plotConvexHull works as expected with different opacity", {
   expect_silent(tf(type = "p", opacity = 0.8, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.4, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, cgroup = g2))
   expect_silent(tf(type = "p", opacity = 0.1, cgroup = g))
})

par(mfrow = c(2, 2))
test_that("plotConvexHull works as expected with different line parameters", {
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

tf <- function(type, cgroup, ...) {
   p <- mdaplot(people, type = type, cgroup = cgroup)
   plotConvexHull(p, ...)
}

par(mfrow = c(2, 2))
test_that("plotConvexHull returns error if wrong type of plot is used", {
   expect_error(tf(type = "l", cgroup = g1))
   expect_error(tf(type = "b", cgroup = g1))
   expect_error(tf(type = "h", cgroup = g1))
   expect_error(tf(type = "e", cgroup = g1))
})

test_that("plotConvexHull returns error if cgroup is not a factor", {
   expect_error(tf(type = "p", cgroup = gw))
})

par(mfrow = c(2, 2))
test_that("plotConvexHull works as expected with default styles", {
   expect_silent(tf(type = "p", cgroup = g1))
   expect_silent(tf(type = "p", cgroup = g2))
   expect_silent(tf(type = "p", cgroup = g))
})

par(mfrow = c(2, 2))
test_that("plotConvexHull works as expected with different opacity", {
   expect_silent(tf(type = "p", opacity = 0.8, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.4, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, cgroup = g2))
   expect_silent(tf(type = "p", opacity = 0.1, cgroup = g))
})

tf(type = "p", opacity = 0.1, lwd = 3, lty = 3, cgroup = g)

par(mfrow = c(2, 2))
test_that("plotConvexHull works as expected with different line parameters", {
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
   plotConfidenceEllipse(p, ...)
}

par(mfrow = c(2, 2))
test_that("plotConfidenceEllipse returns error if wrong type of plot is used", {
   expect_error(tf(type = "l", cgroup = g1))
   expect_error(tf(type = "b", cgroup = g1))
   expect_error(tf(type = "h", cgroup = g1))
   expect_error(tf(type = "e", cgroup = g1))
})

test_that("plotConfidenceEllipse returns error if cgroup is not a factor", {
   expect_error(tf(type = "p", cgroup = gw))
})

par(mfrow = c(2, 2))
test_that("plotConfidenceEllipse works as expected with default styles", {
   expect_silent(tf(type = "p", cgroup = g1))
   expect_silent(tf(type = "p", cgroup = g2))
   expect_silent(tf(type = "p", cgroup = g))
})

par(mfrow = c(2, 2))
test_that("plotConfidenceEllipse works as expected with different opacity", {
   expect_silent(tf(type = "p", opacity = 0.8, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.4, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.2, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.1, cgroup = g))
})

par(mfrow = c(2, 2))
test_that("plotConfidenceEllipse works as expected with different confidence level", {
   expect_silent(tf(type = "p", opacity = 0.8, conf.level = 0.80, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.4, conf.level = 0.90, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.2, conf.level = 0.95, cgroup = g))
   expect_silent(tf(type = "p", opacity = 0.1, conf.level = 0.99, cgroup = g))
})

par(mfrow = c(2, 2))
test_that("plotConfidenceEllipse works as expected with different line parameters", {
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 1, lty = 1, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 2, lty = 2, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 3, lty = 3, cgroup = g1))
   expect_silent(tf(type = "p", opacity = 0.2, lwd = 3, lty = 4, cgroup = g1))
})


#############################################################
#  Block 3: density plot                                    #
#############################################################

## test density plots
context("mdaplot: density plot (hexbin)")

N <- 100000
x <- rnorm(N, 0, 1)
y <- rnorm(N, 1, 5)
d <- cbind(x, y)
colnames(d) <- c("Var 1", "Var 2")

tf <- function(x, ...) mdaplot(x, type = "d", ...)
par(mfrow = c(2, 2))
test_that("density plot works as expected with main parameters", {
   expect_silent(tf(d))
   expect_silent(tf(d, nbins = 10, colmap = "gray"))
   expect_silent(tf(d, nbins = 100, colmap = "jet"))
   expect_silent(tf(d, nbins = 50, colmap = c("blue", "green", "red")))
})
