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


###############################
# 1. Prepare datasets         #
###############################
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


###############################
# 2. Test forward method      #
###############################
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


###############################
# 3. Test backward method     #
###############################
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

###############################
# 4. Test ipls with test set  #
###############################
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



###############################
# 5. Test new parameters      #  - added 2022-11-07
###############################
context("ipls: check that parameters 'full' and 'int.niter' work correctly.")

data(simdata)
X = simdata$spectra.c
y = simdata$conc.c[, 2, drop = FALSE]

# check forward method

## - default settings (full = FALSE)
m1 = ipls(X, y, glob.ncomp = 4, int.num = 15)
expect_equal(nrow(m1$int.stat), 7)
expect_equal(length(m1$int.selected), 6)
expect_equivalent(m1$int.selected, c(6, 11, 1, 2, 13, 12))

## - with full = TRUE
m2 = ipls(X, y, glob.ncomp = 4, int.num = 15, full = TRUE)
expect_equal(nrow(m2$int.stat), 16)
expect_equal(length(m2$int.selected), 6)
expect_equivalent(m2$int.selected, c(6, 11, 1, 2, 13, 12))

## - default settings (full = FALSE) + different int.num
m3 = ipls(X, y, glob.ncomp = 4, int.num = 50)
expect_equal(nrow(m3$int.stat), 15)
expect_equal(length(m3$int.selected), 14)

## - with full = TRUE and different int.num
m4 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE)
expect_equal(nrow(m4$int.stat), 31)
expect_equal(length(m4$int.selected), 27)

## - with full = TRUE and different int.num + larger max.niter
m5 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE, int.niter = 40)
expect_equal(nrow(m5$int.stat), 41)
expect_equal(length(m5$int.selected), 27)

## - with full = TRUE and different int.num + very large max.niter
m6 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE, int.niter = 100)
expect_equal(nrow(m6$int.stat), 51)
expect_equal(length(m6$int.selected), 27)

# check backward method - in this case number of iterations is less than number of intervals - 1

## - default settings (full = FALSE)
m1 = ipls(X, y, glob.ncomp = 4, int.num = 15, method = "backward")
expect_equal(nrow(m1$int.stat), 10)
expect_equivalent(m1$int.selected, c(1, 3, 5, 9, 11, 13))

## - with full = TRUE
m2 = ipls(X, y, glob.ncomp = 4, int.num = 15, full = TRUE, method = "backward")
expect_equal(nrow(m2$int.stat), 15)
expect_equivalent(m2$int.selected, c(1, 3, 5, 9, 11, 13))

## - default settings (full = FALSE) + different int.num
m3 = ipls(X, y, glob.ncomp = 4, int.num = 50, method = "backward")
expect_equal(nrow(m3$int.stat), 31)
expect_equal(length(m3$int.selected), 20) # 50 in total minus 30 excluded

## - with full = TRUE and different int.num + large max.niter
m4 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE, int.niter = 40, method = "backward")
expect_equal(nrow(m4$int.stat), 41)
expect_equal(length(m4$int.selected), 13) # 50 in total minus 37 excluded

## - with full = TRUE and different int.num + very large max.niter
m5 = ipls(X, y, glob.ncomp = 4, int.num = 50, full = TRUE, int.niter = 100, method = "backward")
expect_equal(nrow(m5$int.stat), 50) # not 51 because last interval can not be excluded, so general model + 49 iterations
expect_equal(length(m5$int.selected), 13) # 50 in total minus 37 excluded


###############################
# 6. Test using Beer data     #  - added 2022-11-09
###############################

# load data used in the paper: Anderson, Bro, JChem, 2010
d <- read.csv2("Beer.csv")
# d <- read.csv2("tests/testthat/Beer.csv")
y <- d[1:40, 1, drop = FALSE]
X <- d[1:40, 2:ncol(d)]
w <- seq(400, 2250, by = 2)
attr(X, "xaxis.values") <- w
attr(X, "xaxis.name") <- "Wavelength, nm"

## 20 intervals: 1240-1330 must be selected
mb <- ipls(X, y, cv = 1, int.num = 20)
expect_equal(length(mb$int.selected), 1)
expect_equal(max(w[mb$var.selected]), 1330)
expect_equal(max(w[mb$var.selected]), 1330)
expect_equal(round(mb$int.stat$RMSE[2], 2), 0.14)
expect_equivalent(mb$glob.stat$nComp[2:21], c(2, 1, 1, 1, 1, 3, 3, 6, 5, 7, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1))

## test for global minimum
mb <- ipls(X, y, cv = 1, int.num = 20, full = TRUE)

int.selected <- c(10, 9, 7, 5, 6, 4, 3, 8)
var.selected <- do.call(c, lapply(int.selected, function(x) mb$int.limits[x, 1]:mb$int.limits[x, 2]))
expect_equal(nrow(mb$int.stat), 21)
expect_equal(length(mb$int.selected), 8)
expect_equivalent(mb$int.selected, int.selected)
expect_equivalent(mb$var.selected, var.selected)
expect_equivalent(mb$int.stat$selected, c(FALSE, rep(TRUE, 8), rep(FALSE, 12)))

# ! add a test for manual interval boundaries, so they have different length

###############################
# 7. Test using manual splits #  - added 2022-11-09
###############################

int.limits <- matrix(c(1, 50, 51, 100, 101, 200, 201, 220, 221, 240, 241, 260, 261, 280, 281, 300, 301, 320, 321, 340, 341, 360, 361, 380, 381, 400, 401, 420, 421, 440, 441, 460, 461, 480, 481, 500, 501, 750, 751, 926), ncol = 2, byrow = TRUE)

## full = FALSE
mb <- ipls(X, y, cv = 1, int.limits = int.limits)
int.selected <- c(14, 12, 8, 9)
var.selected <- do.call(c, lapply(int.selected, function(x) int.limits[x, 1]:int.limits[x, 2]))

expect_equal(nrow(mb$int.stat), length(int.selected) + 1)
expect_equal(length(mb$int.selected), length(int.selected))
expect_equivalent(mb$int.selected, int.selected)
expect_equivalent(mb$var.selected, var.selected)
expect_equivalent(mb$int.stat$selected, c(FALSE, rep(TRUE, length(int.selected))))

## full = TRUE
mb <- ipls(X, y, cv = 1, int.limits = int.limits, full = TRUE, int.niter = 100)
int.selected <- c(14, 12, 8, 9, 5, 10, 4, 6)
var.selected <- do.call(c, lapply(int.selected, function(x) int.limits[x, 1]:int.limits[x, 2]))

expect_equal(nrow(mb$int.stat), nrow(int.limits) + 1)
expect_equal(length(mb$int.selected), length(int.selected))
expect_equivalent(mb$int.selected, int.selected)
expect_equivalent(mb$var.selected, var.selected)
expect_equivalent(mb$int.stat$selected, c(FALSE, rep(TRUE, length(int.selected)), rep(FALSE, nrow(int.limits) - length(int.selected))))

