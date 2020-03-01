#####################################################
# Tests for basic functionality of plot() class  #
#####################################################

setup({
   pdf(file = tempfile("mdatools-test-pcaplots-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-pcaplots-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# prepare cases

## 1. full data no test set
data(people)
x1 <- people
m1 <- pca(x1, 10, scale = T)

## 2. calibration and test sets
ind <- seq(1, 32, by = 4)
x2 <- people[-ind, ]
x2.test <- people[ind, ]
m2 <- pca(x2, 10, scale = T, x.test = x2.test)

## 3. calibration and test sets with excluded data
x3 <- people[-ind, ]
x3 <- mda.exclrows(x3, c(1, 10, 20))
x3 <- mda.exclcols(x3, c(3, 12))
x3.test <- people[ind, ]
m3 <- pca(x3, 10, scale = T, x.test = x3.test)

## combine all together
x <- list(x1, x2, x3)
x.test <- list(NULL, x2.test, x3.test)
m <- list("full" = m1, "full + test" = m2, "excluded + test" = m3)

#########################################
# Block 1: Variance plot                #
#########################################

context("pca: variance plots")

tf <- function(model, name) {

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("PCA - variance plots - ", name), pos = 4)

   # basic variance plots
   par(mfrow = c(2, 2))
   test_name <- "variance plot works fine with different settings"
   test_that(test_name, {
      expect_silent(plotVariance(model))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotVariance(model, type = "h", show.labels = T))
      expect_silent(plotVariance(model, show.labels = T, col = c("red", "green")[seq_along(model$res)]))
      expect_silent(plotVariance(model, res = list("cal" = model$res[["cal"]])))
   })

   # basic cumulative plots
   par(mfrow = c(2, 2))
   test_name <- "cumulative variance plot works fine with different settings"
   test_that(test_name, {
      expect_silent(plotCumVariance(model))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotCumVariance(model, type = "h", show.labels = T))
      expect_silent(plotCumVariance(model, show.labels = T, col = c("red", "green")[seq_along(model$res)]))
      expect_silent(plotCumVariance(model, res = list("cal" = model$res[["cal"]])))
   })

}

for (i in seq_len(length(m))) {
   tf(m[[i]], names(m)[[i]])
}

#########################################
# Block 2: Scores plot                  #
#########################################

context("pca: scores plots")

tf <- function(model, x.cal, name) {
   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("PCA - score plots - ", name), pos = 4)

   col = c("red", "green")[seq_along(model$res)]

   # basic plots (scatter)
   par(mfrow = c(2, 2))
   test_name <- "scores plot works fine with basic settings"
   test_that(test_name, {
      expect_silent(plotScores(model))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotScores(model, c(1, 3)))
      expect_silent(plotScores(model, c(1, 3), show.labels = T, col = col))
      expect_silent(plotScores(model, c(1, 3), res = list("cal" = model$res[["cal"]])))
   })

   # basic plots (scatter, one component)
   par(mfrow = c(2, 2))
   test_name <- "scatter plot for single component"
   test_that(test_name, {
      expect_silent(plotScores(model, 1))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotScores(model, 2))
      expect_silent(plotScores(model, 2, show.labels = T, col = col))
      expect_silent(plotScores(model, 2, res = list("cal" = model$res[["cal"]])))
   })

   # line and bar plots
   if (length(model$res) > 1) {
      test_name <- "line and bar plot can nt be made if more than one result object is available"
      test_that(test_name, {
            expect_error(plotScores(model, c(1, 3), type = "l"))
            expect_error(plotScores(model, c(1, 3), type = "h"))
            expect_error(plotScores(model, c(1, 3), type = "b"))
      })
   }

   if (length(model$res) == 1) {
      par(mfrow = c(2, 2))
      test_name <- "line and bar plot for one component"
      test_that(test_name, {
         expect_silent(plotScores(model, 1, type = "l"))
         mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
         expect_silent(plotScores(model, 2, type = "l"))
         expect_silent(plotScores(model, 1, type = "b"))
         expect_silent(plotScores(model, 1, type = "h"))
      })

      par(mfrow = c(2, 2))
      test_name <- "line and bar plot for two or three components"
      test_that(test_name, {
         expect_silent(plotScores(model, c(1, 2), type = "l"))
         mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
         expect_silent(plotScores(model, c(1, 2), type = "h"))
         expect_silent(plotScores(model, c(1, 2, 3), type = "h"))
         expect_silent(plotScores(model, c(1, 2), type = "b"))
      })
   }

   # playing with legend
   if (length(model$res) > 1) {
      par(mfrow = c(2, 2))
      test_name <- "legend position can be changed or hidden"
      test_that(test_name, {
         expect_silent(plotScores(model, show.legend = FALSE))
         mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
         expect_silent(plotScores(model, c(1, 3), legend.position = "top"))
         expect_silent(plotScores(model, c(1, 3), show.labels = T, legend.position = "topleft"))
         expect_silent(plotScores(model, c(1, 3), legend.position = "bottom"))
      })
   }

   # show hidden values
   par(mfrow = c(2, 2))
   test_name <- "hidden values can be shown"
   test_that(test_name, {
      expect_silent(plotScores(model, show.excluded = TRUE))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotScores(model, c(1, 3), show.excluded = TRUE))
      expect_silent(plotScores(model, c(1, 3), show.excluded = TRUE, show.labels = T, col = "red"))
      expect_silent(plotScores(model, c(1, 3), show.excluded = TRUE, res = list("cal" = model$res[["cal"]])))
   })

   # add Hotelling ellipse
   test_name <- "Hotelling ellipse can be added"
   if (is.null(model$res$test)) {
      par(mfrow = c(2, 2))
      test_that(test_name, {
         expect_silent({p <- plotScores(model); plotHotellingEllipse(p)})
         mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
         expect_silent({p <- plotScores(model, c(2, 3)); plotHotellingEllipse(p)})
         expect_silent({p <- plotScores(model); plotHotellingEllipse(p)})
         expect_silent({p <- plotScores(model, c(2, 3)); plotHotellingEllipse(p)})
      })
   } else {
      expect_error({p <- plotScores(model); plotHotellingEllipse(p)})
      expect_error({p <- plotScores(model, c(2, 3)); plotHotellingEllipse(p)})
      par(mfrow = c(2, 2))
      test_that(test_name, {
         expect_silent({p <- plotScores(model); plotHotellingEllipse(p[[1]])})
         mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
         expect_silent({p <- plotScores(model, c(2, 3)); plotHotellingEllipse(p[[1]])})
         expect_silent({p <- plotScores(model); plotHotellingEllipse(p[[1]])})
         expect_silent({p <- plotScores(model, c(2, 3)); plotHotellingEllipse(p[[1]])})
      })
   }

}

for (i in seq_along(m)) {
   tf(m[[i]], x[[i]], names(m)[[i]])
}

#########################################
# Block 3: Residual distance plot       #
#########################################

context("pca: residual distance plots")

tf <- function(x.cal, x.test, name) {

   m1 <- pca(x.cal, 4, scale = T, x.test = x.test)
   m2 <- pca(x.cal, 4, scale = T, x.test = x.test, lim.type = "chisq")
   m3 <- pca(x.cal, 4, scale = T, x.test = x.test, lim.type = "ddmoments", alpha = 0.1, gamma = 0.05)
   m4 <- pca(x.cal, 4, scale = T, x.test = x.test, lim.type = "ddrobust")

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("PCA - residual distance plots - ", name), pos = 4)

   # basic plots (scatter)
   par(mfrow = c(2, 2))
   test_name <- "residual distance plot with default settings work fine"
   test_that(test_name, {
      expect_silent(plotResiduals(m1))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2))
      expect_silent(plotResiduals(m3))
      expect_silent(plotResiduals(m4))
   })

   par(mfrow = c(2, 2))
   test_name <- "residual distance plot work fine with labels and ncomp value"
   test_that(test_name, {
      expect_silent(plotResiduals(m1, show.labels = T))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, ncomp = 2, show.labels = T))
      expect_silent(plotResiduals(m3, ncomp = 2, show.labels = T))
      expect_silent(plotResiduals(m4, ncomp = 2, show.labels = T))
   })

   par(mfrow = c(2, 2))
   test_name <- "residual distance plot can show color groups by categories"
   test_that(test_name, {
      expect_silent(plotResiduals(m1, show.labels = T, cgroup = "categories"))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, show.labels = T, cgroup = "categories"))
      expect_silent(plotResiduals(m3, show.labels = T, cgroup = "categories"))
      expect_silent(plotResiduals(m4, show.labels = T, cgroup = "categories"))
   })

   par(mfrow = c(2, 2))
   test_name <- "residual distance plot can work with log transform"
   test_that(test_name, {
      expect_silent(plotResiduals(m1, show.labels = T, log = T, cgroup = "categories"))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, show.labels = T, log = T, cgroup = "categories"))
      expect_silent(plotResiduals(m3, show.labels = T, log = T, cgroup = "categories"))
      expect_silent(plotResiduals(m4, show.labels = T, log = T, cgroup = "categories"))
   })

   par(mfrow = c(2, 2))
   test_name <- "residual distance plot can work without normaliszation"
   test_that(test_name, {
      expect_silent(plotResiduals(m1, show.labels = T, norm = F, cgroup = "categories"))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, show.labels = T, norm = F, cgroup = "categories"))
      expect_silent(plotResiduals(m3, show.labels = T, norm = F, cgroup = "categories"))
      expect_silent(plotResiduals(m4, show.labels = T, norm = F, cgroup = "categories"))
   })

   par(mfrow = c(2, 2))
   test_name <- "residual distance plot can hide limits"
   test_that("residual distance plot can hide limits", {
      expect_silent(plotResiduals(m1, show.labels = T, cgroup = "categories", show.limits = F))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, show.labels = T, cgroup = "categories", show.limits = F))
      expect_silent(plotResiduals(m3, show.labels = T, cgroup = "categories", show.limits = F))
      expect_silent(plotResiduals(m4, show.labels = T, cgroup = "categories", show.limits = F))
   })

   par(mfrow = c(2, 2))
   test_name <- "residual distance plot can be used with manual title and axis labels"
   test_that(test_name, {
      expect_silent(plotResiduals(m1, main = "Distances", xlab = "T2", ylab = "Q"))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, main = "Distances", xlab = "T2", ylab = "Q"))
      expect_silent(plotResiduals(m3, main = "Distances", xlab = "T2", ylab = "Q"))
      expect_silent(plotResiduals(m4, main = "Distances", xlab = "T2", ylab = "Q"))
   })

   par(mfrow = c(2, 2))
   test_name <- "residual distance plot can show excluded values"
   test_that(test_name, {
      expect_silent(plotResiduals(m1, show.excluded = T))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, show.excluded = T))
      expect_silent(plotResiduals(m3, show.excluded = T))
      expect_silent(plotResiduals(m4, show.excluded = T))
   })


   test_name <- "residual distance plot understands two values for show.limits"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotResiduals(m1, show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, show.limits = c(FALSE, TRUE)))
      expect_silent(plotResiduals(m3, show.limits = c(TRUE, FALSE)))
      expect_silent(plotResiduals(m4, show.limits = c(TRUE, TRUE)))
   })

   test_name <- "show.limits works fine with log transform"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotResiduals(m1, log = TRUE, show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, log = TRUE, show.limits = c(FALSE, TRUE)))
      expect_silent(plotResiduals(m3, log = TRUE, show.limits = c(TRUE, FALSE)))
      expect_silent(plotResiduals(m4, log = TRUE, show.limits = c(TRUE, TRUE)))
   })

   test_name <- "show.limits works fine with none normalized values"
   par(mfrow = c(2, 2))
   test_that(test_name, {
      expect_silent(plotResiduals(m1, norm = FALSE, show.limits = c(FALSE, FALSE)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotResiduals(m2, norm = FALSE, show.limits = c(FALSE, TRUE)))
      expect_silent(plotResiduals(m3, norm = FALSE, show.limits = c(TRUE, FALSE)))
      expect_silent(plotResiduals(m4, norm = FALSE, show.limits = c(TRUE, TRUE)))
   })
}

for (i in seq_len(length(m))) {
   tf(x[[i]], x.test[[i]], names(m)[[i]])
}


#########################################
# Block 4: Loadings plot                #
#########################################

context("pca: loadings plots")

tf <- function(model, x.cal, name) {
   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("PCA - loadings plots - ", name), pos = 4)

   # basic plots (scatter)
   par(mfrow = c(2, 2))
   test_name <- "loadings plot works fine with basic settings"
   test_that(test_name, {
      expect_silent(plotLoadings(model))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotLoadings(model, c(1, 3)))
      expect_silent(plotLoadings(model, c(1, 3), show.labels = T, col = c("red", "green")))
      expect_silent(plotLoadings(model, 1))
   })


   par(mfrow = c(2, 2))
   test_name <- "line and bar plot for one component"
   test_that(test_name, {
      expect_silent(plotLoadings(model, 1, type = "l"))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotLoadings(model, 2, type = "l"))
      expect_silent(plotLoadings(model, 1, type = "b"))
      expect_silent(plotLoadings(model, 1, type = "h"))
   })

   par(mfrow = c(2, 2))
   test_name <- "line and bar plot for one component and labels"
   test_that(test_name, {
      expect_silent(plotLoadings(model, 1, type = "l", show.labels = T))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotLoadings(model, 2, type = "l", show.labels = T))
      expect_silent(plotLoadings(model, 1, type = "b", show.labels = T))
      expect_silent(plotLoadings(model, 1, type = "h", show.labels = T))
   })

   par(mfrow = c(2, 2))
   test_name <- "line and bar plot for one component and xticklabels"
   test_that(test_name, {
      expect_silent(plotLoadings(model, 1, type = "l", xticks = 1:ncol(x.cal), xlas = 2, xticklabels = colnames(x.cal)))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotLoadings(model, 2, type = "l", xticks = 1:ncol(x.cal), xlas = 2, xticklabels = colnames(x.cal)))
      expect_silent(plotLoadings(model, 1, type = "b", xticks = 1:ncol(x.cal), xlas = 2, xticklabels = colnames(x.cal)))
      expect_silent(plotLoadings(model, 1, type = "h", xticks = 1:ncol(x.cal), xlas = 2, xticklabels = colnames(x.cal)))
   })

   par(mfrow = c(2, 2))
   test_name <- "line and bar plot for two or three components"
   test_that(test_name, {
      expect_silent(plotLoadings(model, c(1, 2), type = "l"))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotLoadings(model, c(1, 2), type = "h"))
      expect_silent(plotLoadings(model, c(1, 2, 3), type = "h"))
      expect_silent(plotLoadings(model, c(1, 2), type = "b"))
   })

   # playing with legend
   par(mfrow = c(2, 2))
   test_name <- "legend position can be changed or hidden"
   test_that(test_name, {
      expect_silent(plotLoadings(model, type = "l", show.legend = FALSE))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotLoadings(model, c(1, 2), type = "l", legend.position = "top"))
      expect_silent(plotLoadings(model, c(1, 2), type = "l", show.labels = T, legend.position = "topleft"))
      expect_silent(plotLoadings(model, c(1, 2), type = "h", legend.position = "bottom"))
   })

   # show hidden values
   par(mfrow = c(2, 2))
   test_name <- "hidden values can be shown"
   test_that(test_name, {
      expect_silent(plotLoadings(model, show.excluded = TRUE))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotLoadings(model, c(1, 3), show.excluded = TRUE))
      expect_silent(plotLoadings(model, c(1, 3), show.excluded = TRUE, show.labels = T, col = "red"))
      expect_silent(plotLoadings(model, c(1, 3), type = "l", show.excluded = TRUE))
   })
}

for (i in seq_along(m)) {
   tf(m[[i]], x[[i]], names(m)[[i]])
}

#########################################
# Block 5: Biplot plot                #
#########################################

context("pca: biplots")

tf <- function(model, x.cal, name) {
   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("PCA - biplots - ", name), pos = 4)

   # basic plots (scatter)
   par(mfrow = c(2, 2))
   test_name <- "biplot works fine with basic settings"
   test_that(test_name, {
      expect_silent(plotBiplot(model))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotBiplot(model, c(1, 3), show.labels = T))
      expect_silent(plotBiplot(model, c(1, 3), lty = 2, lwd = 2))
      expect_silent(plotBiplot(model, c(1, 3), show.labels = T, col = c("red", "green")))
      expect_error(plotBiplot(model, 1))
   })
}

for (i in seq_along(m)) {
   tf(m[[i]], x[[i]], names(m)[[i]])
}


#########################################
# Block 6: PCA for spectral data        #
#########################################

data(simdata)
x.cal <- simdata$spectra.c
x.test <- simdata$spectra.t
attr(x.cal, "xaxis.name") <- attr(x.test, "xaxis.name") <- "Wavenumbers, cm-1"
attr(x.cal, "xaxis.values") <- attr(x.test, "xaxis.values") <- 10^7 / simdata$wavelength

ms <- pca(x.cal, 4, x.test = x.test)

par(mfrow = c(1, 1))
plot.new()
text(0, 0, "PCA - spectral data (simdata)", pos = 4)

par(mfrow = c(2, 2))
expect_silent(plot(ms))

par(mfrow = c(2, 2))
expect_silent(plotLoadings(ms, type = "l"))
expect_silent(plotLoadings(ms, type = "l", ncomp = 1))
expect_silent(plotScores(ms))
expect_silent(plotScores(ms, type = "l", ncomp = 1, res = list(ms$calres)))

par(mfrow = c(2, 2))
expect_silent(plotResiduals(ms))
expect_silent(plotResiduals(ms, ncomp = 2))
expect_silent(plotResiduals(ms, ncomp = 2, log = T))
expect_silent(plotResiduals(ms, ncomp = 4, log = T, cgroup = "categories",
   res = list("cal" = ms$calres)))


#########################################
# Block 7: DoF plots                    #
#########################################

m1 <- pca(x1, 10, scale = TRUE, lim.type = "ddmoments")
m2 <- pca(simdata$spectra.c, 10, lim.type = "ddmoments")

context("pca: DoF plots for T2")

par(mfrow = c(1, 1))
plot.new()
text(0, 0, "PCA - DoF plot for T2 distance", pos = 4)

# basic variance plots
par(mfrow = c(2, 2))
test_that("DoF plot for T2 distance with People data", {
   expect_silent(plotT2DoF(m1))
   expect_silent(plotT2DoF(m1, type = "h", show.labels = T))
   expect_silent(plotT2DoF(m1, type = "b", col = "red"))
   expect_silent(plotT2DoF(m1, type = "l", lty = 2, lwd = 2, col = "green"))
})

par(mfrow = c(2, 2))
test_that("DoF plot for T2 distance with Simdata data", {
   expect_silent(plotT2DoF(m2))
   expect_silent(plotT2DoF(m2, type = "h", show.labels = T))
   expect_silent(plotT2DoF(m2, type = "b", col = "red"))
   expect_silent(plotT2DoF(m2, type = "l", lty = 2, lwd = 2, col = "green"))
})

context("pca: DoF plots for Q")

par(mfrow = c(1, 1))
plot.new()
text(0, 0, "PCA - DoF plot for Q distance", pos = 4)

# basic variance plots
par(mfrow = c(2, 2))
test_that("DoF plot for Q distance with People data", {
   expect_silent(plotQDoF(m1))
   expect_silent(plotQDoF(m1, type = "h", show.labels = T))
   expect_silent(plotQDoF(m1, type = "b", col = "red"))
   expect_silent(plotQDoF(m1, type = "l", lty = 2, lwd = 2, col = "green"))
})

par(mfrow = c(2, 2))
test_that("DoF plot for Q distance with Simdata data", {
   expect_silent(plotQDoF(m2))
   expect_silent(plotQDoF(m2, type = "h", show.labels = T))
   expect_silent(plotQDoF(m2, type = "b", col = "red"))
   expect_silent(plotQDoF(m2, type = "l", lty = 2, lwd = 2, col = "green"))
})

context("pca: DoF plots for both T2 and Q")

par(mfrow = c(1, 1))
plot.new()
text(0, 0, "PCA - DoF plot for both distances", pos = 4)

# basic variance plots
par(mfrow = c(2, 2))
test_that("DoF plot for both distances with People data", {
   expect_silent(plotDistDoF(m1))
   expect_silent(plotDistDoF(m1, type = "l", show.labels = T))
   expect_silent(plotDistDoF(m1, type = "b", col = c("red", "blue")))
   expect_silent(plotDistDoF(m1, type = "l", lty = c(1, 2), lwd = c(1, 2), col = c("red", "red")))
})

par(mfrow = c(2, 2))
test_that("DoF plot for both distances with Simdata data", {
   expect_silent(plotDistDoF(m2))
   expect_silent(plotDistDoF(m2, type = "l", show.labels = T))
   expect_silent(plotDistDoF(m2, type = "b", col = c("red", "blue")))
   expect_silent(plotDistDoF(m2, type = "l", lty = c(1, 2), lwd = c(1, 2), col = c("red", "red")))
})

#########################################
# Block 8: Extreme plots                    #
#########################################

context("pca: extreme plots")

tf <- function(model, name) {

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("PCA - extreme plots - ", name), pos = 4)

   # basic extreme plots
   par(mfrow = c(2, 2))
   test_name <- "extreme plot works fine with default settings"
   test_that(test_name, {
      expect_silent(plotExtreme(model))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotExtreme(model, comp = 1:3, pch = 21:23))
      expect_silent(plotExtreme(model, comp = 1:2, bg = c("red", "green")))
      expect_silent(plotExtreme(model, pch = 1, lwd = 1, col = "orange"))
   })

   # extreme plots for test set
   if (!is.null(model$res[["test"]])) {
      test_name <- "extreme plot works fine with manual result object"
      test_that(test_name, {
         expect_silent(plotExtreme(model, res = model$res[["test"]]))
         mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
         expect_silent(plotExtreme(model, res = model$res[["test"]], comp = 1:3, pch = 21:23))
         expect_silent(plotExtreme(model, res = model$res[["test"]], comp = 1:2, bg = c("red", "green")))
         expect_silent(plotExtreme(model, res = model$res[["test"]], pch = 1, lwd = 1, col = "orange"))
      })
   }

   # ddrobust extreme plots
   model <- setDistanceLimits(model, lim.type = "ddrobust")
   par(mfrow = c(2, 2))
   test_name <- "extreme plot works fine with ddrobust lim type"
   test_that(test_name, {
      expect_silent(plotExtreme(model))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotExtreme(model, comp = 1:3, pch = 21:23))
      expect_silent(plotExtreme(model, comp = 1:2, bg = c("red", "green")))
      expect_silent(plotExtreme(model, pch = 1, lwd = 1, col = "orange"))
   })

   # chisq extreme plots
   model <- setDistanceLimits(model, lim.type = "chisq")
   par(mfrow = c(2, 2))
   test_name <- "extreme plot works fine with chisq lim type"
   test_that(test_name, {
      expect_silent(plotExtreme(model))
      mtext(test_name, side = 3, line = -1, outer = TRUE, cex = 0.75, col = "gray")
      expect_silent(plotExtreme(model, comp = 1:3, pch = 21:23))
      expect_silent(plotExtreme(model, comp = 1:2, bg = c("red", "green")))
      expect_silent(plotExtreme(model, pch = 1, lwd = 1, col = "orange"))
   })

}


tf(ms, "spectral data")


# just output to check in txt file
print(ms)
print(summary(ms))
