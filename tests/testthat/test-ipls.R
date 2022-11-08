#############################
# Tests for ipls() methods  #
#############################

setup({
   pdf(file = tempfile("mdatools-test-ipls-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-ipls-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# Prepare datasets
datasets <- list()

## small data no names no extra arguments
data(people)
xc <- people[, -4, drop = F]
yc <- people[, 4, drop = F]
datasets[["people"]] <- list(xc = xc, yc = yc, center = T, scale = T, ncomp = 8)

## spectral data with extra attributes, three responses, hidden values and test set
data(simdata)
xc <- simdata$spectra.c
yc <- simdata$conc.c[, 1, drop = F]
xc <- mda.exclrows(xc, c(1, 10, 20, 30, 40))
xc <- mda.exclcols(xc, c(1, 100:110))
yc <- mda.exclrows(yc, c(1, 10, 20, 30, 40))
attr(xc, "name") <- "Spectra, cal"
attr(xc, "xaxis.name") <- "Wavelength, nm"
attr(xc, "xaxis.values") <- simdata$wavelength

xt <- simdata$spectra.t
yt <- simdata$conc.t[, 1, drop = F]
xt <- mda.exclrows(xt, c(15, 35))
yt <- mda.exclrows(yt, c(15, 35))
attr(xt, "name") <- "Spectra, val"
attr(xt, "xaxis.name") <- "Wavelength, nm"
attr(xt, "xaxis.values") <- simdata$wavelength

datasets[["spectra"]] <- list(xc = xc, yc = yc, xt = xt, yt = yt, center = T, scale = F, ncomp = 6)

context("ipls: test forward method")

for (i in seq_along(datasets)) {
   d <- datasets[[i]]
   name <- names(datasets)[i]

   expect_error(ipls(d$xc, d$yc, glob.ncomp = d$ncomp, int.num = 1))
   expect_error(ipls(d$xc, d$yc, glob.ncomp = d$ncomp, int.num = 10000))
   expect_error(ipls(d$xc, d$yc, glob.ncomp = d$ncomp, int.width = 0))
   expect_error(ipls(d$xc, d$yc, glob.ncomp = d$ncomp, int.width = 10000))

   expect_silent(m1 <- ipls(d$xc, d$yc, glob.ncomp = d$ncomp, center = d$center,
      scale = d$scale, int.num = 4, silent = TRUE))

   expect_output(m1 <- ipls(d$xc, d$yc, glob.ncomp = d$ncomp, center = d$center,
      scale = d$scale, int.num = 10))

   if (ncol(d$xc) < 20) {
      expect_output(m2 <- ipls(d$xc, d$yc, glob.ncomp = d$ncomp,
         center = d$center, scale = d$scale, int.width = 1))
   } else {
      expect_output(m2 <- ipls(d$xc, d$yc, glob.ncomp = d$ncomp,
         center = d$center, scale = d$scale, int.width = 20))
   }
   #expect_output(m3 <- ipls(d$xc, d$yc, glob.ncomp = d$ncomp, center = d$center,
   #   scale = d$scale, int.width = 2))

   expect_output(summary(m1))
   expect_output(summary(m2))
   #expect_output(summary(m3))

   expect_silent(plot(m1))
   expect_silent(plot(m2))
   #expect_silent(plot(m3))

   expect_silent(plotRMSE(m1))
   expect_silent(plotRMSE(m2))
   #expect_silent(plotRMSE(m3))

   cat("\nOutput for forward iPLS\n\n")
   summary(m1)
   summary(m2)
   #summary(m3)

   print(m1)
   print(m2)
   #print(m3)
}

context("ipls: test backward method")

for (i in seq_along(datasets)) {
   d <- datasets[[i]]
   name <- names(datasets)[i]

   expect_error(ipls(d$xc, d$yc, method = "backward", glob.ncomp = d$ncomp, int.num = 1))
   expect_error(ipls(d$xc, d$yc, method = "backward", glob.ncomp = d$ncomp, int.num = 10000))
   expect_error(ipls(d$xc, d$yc, method = "backward", glob.ncomp = d$ncomp, int.width = 0))
   expect_error(ipls(d$xc, d$yc, method = "backward", glob.ncomp = d$ncomp, int.width = 10000))

   expect_silent(m1 <- ipls(d$xc, d$yc, method = "backward", glob.ncomp = d$ncomp, center = d$center,
      scale = d$scale, int.num = 4, silent = TRUE))

   expect_output(m1 <- ipls(d$xc, d$yc, method = "backward", glob.ncomp = d$ncomp,
      center = d$center, scale = d$scale, int.num = 10))

   if (ncol(d$xc) < 20) {
      expect_output(m2 <- ipls(d$xc, d$yc, method = "backward", glob.ncomp = d$ncomp,
         center = d$center, scale = d$scale, int.width = 1))
   } else {
      expect_output(m2 <- ipls(d$xc, d$yc, method = "backward", glob.ncomp = d$ncomp,
         center = d$center, scale = d$scale, int.width = 20))
   }

   #expect_output(m3 <- ipls(d$xc, d$yc, method = "backward", glob.ncomp = d$ncomp,
   #   center = d$center, scale = d$scale, int.width = 2))

   expect_output(summary(m1))
   expect_output(summary(m2))
   #expect_output(summary(m3))

   expect_silent(plot(m1))
   expect_silent(plot(m2))
   #expect_silent(plot(m3))

   expect_silent(plotRMSE(m1))
   expect_silent(plotRMSE(m2))
   #expect_silent(plotRMSE(m3))

   cat("\nOutput for backward iPLS\n\n")
   summary(m1)
   summary(m2)
   #summary(m3)

   print(m1)
   print(m2)
   #print(m3)
}

context("ipls: check if works with test set")

d <- datasets[["spectra"]]

# test that if test set is bad it raises an error
expect_error(ipls(d$xc, d$yc, int.width = 10, glob.ncomp = d$ncomp, x.test = d$xt[, -1], y.test = d$yt))
expect_error(ipls(d$xc, d$yc, int.width = 10, glob.ncomp = d$ncomp, x.test = d$xt[-1, ], y.test = d$yt))
expect_error(ipls(d$xc, d$yc, int.width = 10, glob.ncomp = d$ncomp, x.test = d$xt, y.test = d$yt[-1, ]))

# test forward method with test set
expect_silent(m <- ipls(d$xc, d$yc, method = "forward", int.width = 10, glob.ncomp = d$ncomp,
   x.test = d$xt, y.test = d$yt, silent = TRUE))

expect_output(summary(m))
expect_silent(plot(m))
expect_equivalent(m$int.selected, c(6, 11, 4, 9))

expect_silent(m <- ipls(d$xc, d$yc, method = "backward", int.width = 10, glob.ncomp = d$ncomp,
   x.test = d$xt, y.test = d$yt, silent = TRUE))

expect_output(summary(m))
expect_silent(plot(m))
expect_equivalent(m$int.selected, c(3, 4, 8, 13))


### ! added 2022-11-07 ###
context("ipls: check that parameters 'full' and 'int.niter' work correctly.")

data(simdata)
X = simdata$spectra.c
y = simdata$conc.c[, 2, drop = FALSE]

# check forward method

## - default settings (full = FALSE)
m1 = ipls(X, y, glob.ncomp = 4, int.num = 15)
expect_equivalent(m1$int.selected, c(6, 1, 11, 3, 12))

## - with full = TRUE
m2 = ipls(X, y, glob.ncomp = 4, int.num = 15, full = TRUE)
expect_equal(length(m2$int.selected), 15)
expect_equivalent(m2$int.selected[1:5], c(6, 1, 11, 3, 12))

## - default settings (full = FALSE) + different int.num
m3 = ipls(X, y, glob.ncomp = 4, int.num = 50)
expect_equal(length(m3$int.selected), 20)

## - with full = TRUE and different int.num
m4 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE)
expect_equal(length(m4$int.selected), 30)

## - with full = TRUE and different int.num + larger max.niter
m5 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE, int.niter = 40)
expect_equal(length(m5$int.selected), 40)

## - with full = TRUE and different int.num + very large max.niter
m6 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE, int.niter = 100)
expect_equal(length(m6$int.selected), 50)

# check backward method

## - default settings (full = FALSE)
m1 = ipls(X, y, glob.ncomp = 4, int.num = 15, method = "backward")
expect_equivalent(m1$int.selected, c(1, 3, 5, 9, 11))

## - with full = TRUE
m2 = ipls(X, y, glob.ncomp = 4, int.num = 15, full = TRUE, method = "backward")
expect_equal(length(m2$int.selected), 2)
expect_equivalent(m2$int.selected, c(1, 5))

## - default settings (full = FALSE) + different int.num
m3 = ipls(X, y, glob.ncomp = 4, int.num = 50, method = "backward")
expect_equal(length(m3$int.selected), 20)

## - with full = TRUE and different int.num
m4 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE, method = "backward")
expect_equal(length(m4$int.selected), 20)

## - with full = TRUE and different int.num + larger max.niter
m5 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE, int.niter = 40, method = "backward")
expect_equal(length(m5$int.selected), 10)

## - with full = TRUE and different int.num + very large max.niter
m6 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE, int.niter = 100, method = "backward")
expect_equal(length(m6$int.selected), 2)

# compare with the paper: Anderson, Bro, JChem, 2010

d <- read.csv2("Beer.csv")
y <- d[1:40, 1, drop = FALSE]
X <- d[1:40, 2:ncol(d)]
w <- seq(400, 2250, by = 2)
attr(X, "xaxis.values") <- w
attr(X, "xaxis.name") <- "Wavelength, nm"

# 20 intervals: 1240-1330
mb <- ipls(X, y, cv = 1, int.num = 20)
expect_equal(length(mb$int.selected), 1)
expect_equal(max(w[mb$var.selected]), 1330)
expect_equal(max(w[mb$var.selected]), 1330)
expect_equal(round(mb$int.stat$RMSE[2], 2), 0.14)
expect_equivalent(mb$glob.stat$nComp[2:21], c(2, 1, 1, 1, 1, 3, 3, 6, 5, 7, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1))