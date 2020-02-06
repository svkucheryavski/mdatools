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
yt <- simdata$conc.t
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
