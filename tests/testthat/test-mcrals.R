setup({
   #pdf(file = "mdatools-test-mcrpure.pdf")
   pdf(file = tempfile("mdatools-test-mcrpure-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-mcrpure-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})


#!###################################
#! Block 1. Testing solver functins #
#!###################################

context("mcrals: testing solvers")

# simple test from Bro article about fast nnls
C <- matrix(c(73, 87, 72, 80, 71, 74, 2, 89, 52, 46, 7, 71), ncol = 3)
D <- matrix(c(49, 67, 68, 20), ncol = 1)
S.ols <- c(1.123, 0.917, -2.068)
S.nnls <- c(0.650, 0, 0)

test_that("mcrals.ols works correctly", {
   res <- mcrals.ols(D, C)
   expect_equal(dim(res), c(ncol(D), ncol(C)))
   expect_equal(round(as.numeric(res), 3), S.ols)
})

test_that("mcrals.nnls works correctly", {
   res <- mcrals.nnls(D, C)
   expect_equal(dim(res), c(ncol(D), ncol(C)))
   expect_equal(round(as.numeric(res), 3), S.nnls)
})

test_that("mcrals.fcnnls works correctly", {
   res <- mcrals.fcnnls(D, C)
   expect_equal(dim(res), c(ncol(D), ncol(C)))
   expect_equal(round(as.numeric(res), 3), S.nnls)
})

# simple test from Benthem and Keenan article about fast nnls
C <- matrix(c(95, 23, 61, 42, 89, 76, 46, 2, 82, 44, 62, 79), ncol = 3)
D <- matrix(c(92, 74, 18, 41, 99, 19, 41, 61, 80, 43, 51, 39), ncol = 3)
S.ols  <- c(-0.468, 0.833, 0.284, 0.748, -0.150, 0.286, 0.661, .272, 0.335)
S.nnls  <- c(0, 0.687, 0.284, 0.627, 0, 0.286, 0.352, 0.287, 0.335)

test_that("mcrals.ols works correctly", {
   res <- mcrals.ols(D, C)
   expect_equal(dim(res), c(ncol(D), ncol(C)))
   expect_equal(round(as.numeric(res), 3), S.ols)
})

test_that("mcrals.nnls works correctly", {
   res <- mcrals.nnls(D, C)
   expect_equal(dim(res), c(ncol(D), ncol(C)))
   expect_equal(round(as.numeric(res), 3), S.nnls)
})

test_that("mcrals.fcnnls works correctly", {
   res <- mcrals.fcnnls(D, C)
   expect_equal(dim(res), c(ncol(D), ncol(C)))
   expect_equal(round(as.numeric(res), 3), S.nnls)
})

#!####################################
#! Block 2. Simple tests for SIMDATA #
#!####################################


# concentration constrains
cc <- list(
   constraint("non-negativity")
)

# spectral constrains
sc <- list(
   constraint("non-negativity"),
   constraint("norm", params = list(type = "area"))
)

# solvers
solvers <- list(mcrals.ols, mcrals.nnls, mcrals.fcnnls)


context("mcrals: testing for simdata")

data(simdata)

D <- simdata$spectra.c[order(simdata$conc.c[, 1]), ]
attr(D, "yaxis.name") <- "Time, s"
attr(D, "yaxis.values") <- seq(0, 10, length.out = nrow(simdata$spectra.c))


n <- 1
for (s in solvers) {

   set.seed(1)
   par(mfrow = c(1, 1))
   plot.new()
   plot.window(xlim = c(0, 1), ylim = c(0, 1))
   text(0.25, 1, paste0("mcrals - solver - ", n), pos = 4, font = 2)

   expect_silent(m <- mcrals(D, 3, spec.constraints = sc, cont.constraints = cc,
      spec.solver = s, cont.solver = s))
   summary(m)

   par(mfrow = c(2, 1))
   expect_silent(plotSpectra(m))
   expect_silent(plotContributions(m))

   par(mfrow = c(2, 1))
   expect_silent(plotSpectra(m))
   expect_silent(plotContributions(m))

   expect_silent(m <- mcrals(D, 3, spec.constraints = sc, cont.constraints = cc,
      exclcols = 141:150,
      spec.solver = s, cont.solver = s))
   summary(m)

   par(mfrow = c(2, 1))
   expect_silent(plotSpectra(m))
   expect_silent(plotContributions(m))

   expect_silent(mcrals(D, 3, spec.constraints = sc, cont.constraints = cc,
      exclrows = 1:10,
      spec.solver = s, cont.solver = s))
   summary(m)

   par(mfrow = c(2, 1))
   expect_silent(plotSpectra(m))
   expect_silent(plotContributions(m, type = "l"))

   n <- n + 1
}

context("mcrals: compare nnls and fcnnls")

test_that("nnls and fcnnls give identical results", {
   set.seed(6)
   m1 <- mcrals(D, 3, spec.solver = mcrals.nnls, cont.solver = mcrals.nnls)
   set.seed(6)
   m2 <- mcrals(D, 3, spec.solver = mcrals.fcnnls, cont.solver = mcrals.fcnnls)

   expect_equal(m1$resspec, m2$resspec)
   expect_equal(m1$rescont, m2$rescont)
   expect_equal(m1$expvar, m2$expvar)
})

#!#########################################
#! Block 3. Thorough tests based on Carbs #
#!#########################################

# prepare data for other tests
data(carbs)
C <- carbs$C
S <- carbs$S
D <- carbs$D
D.ini <- pca(D, 3)$loadings
D.ini <- abs(D.ini)

params_ok <- list(
   list(ncomp = 1),
   list(ncomp = 3),

   # initial spectra by PCA
   list(ncomp = 1 , spec.ini = D.ini[, 1, drop = FALSE]),
   list(ncomp = 3 , spec.ini = D.ini),

   # with constraints
   list(ncomp = 1, spec.constraints = sc, cont.constraints = cc),
   list(ncomp = 3, spec.constraints = sc, cont.constraints = cc),

   # with constraints and manual estimated initial spectra
   list(ncomp = 1, spec.ini = D.ini[, 1, drop = F], spec.constraints = sc, cont.constraints = cc),
   list(ncomp = 3, spec.ini = D.ini, spec.constraints = sc, cont.constraints = cc),

   # with excluded columns
   list(ncomp = 3, exclcols = 1:9),
   list(ncomp = 3, exclcols = 1:9, spec.ini = D.ini),
   list(ncomp = 3, spec.ini = D.ini, spec.constraints = sc, cont.constraints = cc),

   # with excluded rows
   list(ncomp = 3, exclrows = 1:3),
   list(ncomp = 3, exclrows = 1:3, spec.ini = D.ini),
   list(ncomp = 3, exclrows = 1:3, spec.ini = D.ini, spec.constraints = sc, cont.constraints = cc)
)

params_err <- list(
   list(ncomp = 0),
   list(ncomp = 30)
)

context("mcrals: testing error cases")
for (p in params_err) {
   expect_error(m <- do.call(mcrals, c(list(x = D), p)))
}

n <- 1
for (s in solvers) {
   context(sprintf("mcrals: testing solver %d", n))

   par(mfrow = c(1, 1))
   plot.new()
   plot.window(xlim = c(0, 1), ylim = c(0, 1))
   text(0.25, 1, paste0("mcrals - solver - ", n), pos = 4, font = 2)

   ncase <- 1
   for (p in params_ok) {


      par(mfrow = c(1, 1))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      text(0.25, 1, paste0("mcrals - case - ", ncase), pos = 4, font = 2)
      text(0.25, 0.5,
         paste0(capture.output(str(p, max.level = 1, give.attr = FALSE)), collapse="\n"), pos = 4)

      expect_output(m <- do.call(mcrals, c(list(x = D, verbose = TRUE), p)))

      set.seed(6)
      expect_silent(m <- do.call(mcrals, c(list(x = D), p)))

      # check dimensions and names
      expect_equal(nrow(m$resspec), ncol(D))
      expect_equal(ncol(m$resspec), p$ncomp)
      expect_equal(nrow(m$rescont), nrow(D))
      expect_equal(colnames(m$rescont), colnames(m$resspec))

      # excluded rows and columns
      expect_equal(attr(m$rescont, "exclrows"), p$exclrows)
      expect_equal(attr(m$resspec, "exclrows"), p$exclcols)

      # variance plots
      par(mfrow = c(2, 2))
      expect_silent(plotVariance(m))
      expect_silent(plotCumVariance(m))
      expect_silent(plotVariance(m, show.labels = TRUE))
      expect_silent(plotCumVariance(m, type = "h", col = "red"))

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
            xlim = c(1600, 200), xticks = seq(1600, 200, by = -200)
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
      if (is.null(p$exclrows)) {
         res <- predict(m, D)
         expect_equivalent(res, m$rescont, tolerance = 0.01)
      }

      ncase <- ncase + 1
   }

   n <- n + 1
}

