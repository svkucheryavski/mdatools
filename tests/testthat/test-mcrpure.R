########################################################
# Block 1: testing methods implementing pca algorithms #
########################################################

setup({
   #pdf(file = "mdatools-test-mcrpure.pdf")
   pdf(file = tempfile("mdatools-test-mcrpure-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-mcrpure-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

library(mdatools)
data(carbs)

params_ok <- list(
   list(ncomp = 1),
   list(ncomp = 1, offset = 0),
   list(ncomp = 1, offset = 0.05),
   list(ncomp = 4),
   list(ncomp = 4, offset = 0),
   list(ncomp = 4, offset = 0.05),
   list(ncomp = 3),
   list(ncomp = 3, offset = 0),
   list(ncomp = 3, offset = 0.05),
   list(ncomp = 3, offset = 0.25),

   # with excluded columns
   list(ncomp = 3, exclcols = 1:9),

   # with excluded rows
   list(ncomp = 3, exclrows = 1:3)
)

params_err <- list(
   list(ncomp = 0),
   list(ncomp = 1, offset = -1),
   list(ncomp = 1, offset = 1),
   list(ncomp = 30)
)

C <- carbs$C
S <- carbs$S
D <- carbs$D


for (p in params_err) {
   context("mcrpure: testing error cases")
   expect_error(m <- do.call(mcrpure, c(list(x = D), p)))
}

n <- 1
for (p in params_ok) {

   context(sprintf("mcrpure: testing case %d", n))

   par(mfrow = c(1, 1))
   plot.new()
   plot.window(xlim = c(0, 1), ylim = c(0, 1))
   text(0.25, 1, paste0("mcrpure - case - ", n), pos = 4, font = 2)
   text(0.25, 0.5, paste0(capture.output(str(p)), collapse="\n"), pos = 4)

   expect_silent(m <- do.call(mcrpure, c(list(x = D), p)))
   summary(m)

   # check dimensions and names
   expect_equal(nrow(m$resspec), ncol(D))
   expect_equal(ncol(m$resspec), p$ncomp)
   expect_equal(nrow(m$purityspec), ncol(D))
   expect_equal(ncol(m$purityspec), p$ncomp)
   expect_equal(nrow(m$rescont), nrow(D))
   expect_equal(colnames(m$rescont), colnames(m$resspec))

   # excluded rows and columns
   expect_equal(attr(m$rescont, "exclrows"), p$exclrows)
   expect_equal(attr(m$resspec, "exclrows"), p$exclcols)
   expect_equal(attr(m$purityspec, "exclrows"), p$exclcols)

   # variance plots
   par(mfrow = c(2, 2))
   expect_silent(plotVariance(m))
   expect_silent(plotCumVariance(m))
   expect_silent(plotVariance(m, show.labels = TRUE))
   expect_silent(plotCumVariance(m, type = "h", col = "red"))

   # purity plots
   par(mfrow = c(2, 2))
   expect_silent(plotPurity(m))
   expect_silent(plotPurity(m, type = "h"))
   expect_silent(plotPurity(m, type = "h", col = "red", show.labels = TRUE))
   expect_silent(plotPurity(m, type = "h", col = "red", show.labels = TRUE, labels = "names"))

   # purity spectra
   par(mfrow = c(2, 2))
   expect_silent(plotPuritySpectra(m))
   expect_silent(plotPuritySpectra(m, comp = 1, type = "b", col = "black"))
   expect_silent(if (m$ncomp > 1) plotPuritySpectra(m, comp = 2, type = "h"))
   expect_silent(if (m$ncomp > 2) plotPuritySpectra(m, comp = 3))

   # resolved contributions
   par(mfrow = c(2, 2))
   expect_silent(plotContributions(m))
   expect_silent(plotContributions(m, comp = 1, type = "b", col = "black"))
   expect_silent(if (m$ncomp > 1) plotContributions(m, comp = 2, type = "h"))
   expect_silent(if (m$ncomp > 2) plotContributions(m, comp = 3))

   # resolved spectra
   par(mfrow = c(2, 2))
   expect_silent(plotSpectra(m))
   expect_silent(plotSpectra(m, comp = 1, type = "b", col = "black"))
   expect_silent(if (m$ncomp > 1) plotSpectra(m, comp = 2, type = "h"))
   expect_silent(if (m$ncomp > 2) plotSpectra(m, comp = 3))

   # resolved spectra vs new
   par(mfrow = c(2, 2))
   for (i in 1:min(c(m$ncomp, ncol(S)))) {
      mdaplotg(
         list(
            original = prep.norm(mda.subset(mda.t(S), i), "area"),
            resolved = prep.norm(mda.subset(mda.t(m$resspec), i), "area")
         ), type = "l", col = c("gray", "red"), lwd = c(2, 1), opacity = c(0.5, 1),
         , xlim = c(1600, 200), xticks = seq(1600, 200, by = -200)
      )
   }

   # resolved contributions vs new
   par(mfrow = c(2, 2))
   for (i in 1:min(c(m$ncomp, ncol(C)))) {
      mdaplotg(
         list(
            original = prep.norm(mda.subset(mda.t(C), i), "area"),
            resolved = prep.norm(mda.subset(mda.t(m$rescont), i), "area")
         ), type = "l", col = c("gray", "red"), lwd = c(2, 1), opacity = c(0.5, 1)
      )
   }


   # check that predictions work correctly
   res <- predict(m, D)
   expect_equivalent(res, m$rescont)

   n <- n + 1
}

context("mcrpure: testing for simdata")

data(simdata)

simdata$spectra.c <- simdata$spectra.c[order(simdata$conc.c[, 1]), ]
attr(simdata$spectra.c, "yaxis.name") <- "Time, s"
attr(simdata$spectra.c, "yaxis.values") <- seq(0, 10, length.out = nrow(simdata$spectra.c))


expect_silent(m <- mcrpure(simdata$spectra.c, 3, offset = 0.0004))
summary(m)

par(mfrow = c(2, 1))
expect_silent(plotSpectra(m))
expect_silent(plotContributions(m))

expect_silent(m <- mcrpure(simdata$spectra.c, 3, offset = 0.0004, exclcols = 141:150))
summary(m)

par(mfrow = c(2, 1))
expect_silent(plotSpectra(m))
expect_silent(plotContributions(m))

expect_silent(m <- mcrpure(simdata$spectra.c, 3, offset = 0.0004, exclrows = 1:10))
summary(m)

par(mfrow = c(2, 1))
expect_silent(plotSpectra(m))
expect_silent(plotContributions(m, type = "l"))
