################################
# Tests for pca class methods  #
################################

## create several datasets with different layout based on People data
data(people)

x1 <- people
x2 <- x1[1:2, , drop = F]
x3 <- x1[, 1:2, drop = F]

## same but scaled
x4 <- scale(people, center = TRUE, scale = TRUE)
x5 <- x4[1:2, , drop = F]
x6 <- x5[, 1:2, drop = F]

## large random data 1000 x 10
x7 <- cbind(
   rnorm(1000, 0, 50),
   rnorm(1000, 0, 10),
   rnorm(1000, 0, 7),
   rnorm(1000, 0, 5),
   rnorm(1000, 0, 3),
   rnorm(1000, 0, 2),
   rnorm(1000, 0, 1),
   rnorm(1000, 0, 0.1),
   rnorm(1000, 0, 0.01),
   rnorm(1000, 0, 0.001)
)

## combine all to a list
x <- list(x1, x2, x3, x4, x5, x6, x7)


########################################################
# Block 1: testing methods implementing pca algorithms #
########################################################

context('pca: testing different algorithms')

## testing function for testing results of one model
tf <- function(x, f) {

   # set of correct and incorrect number of components
   ncomp_correct <- c(1, min(ncol(x), nrow(x) - 1), round(min(ncol(x), nrow(x) - 1)/2) + 1)
   ncomp_wrong <- c(0, 9999)

   for (ncomp in ncomp_correct) {
      res <- f(x, ncomp)

      # check size of scores, loadings and eigenvalues
      expect_equal(dim(res$loadings), c(ncol(x), ncomp))
      expect_equal(dim(res$scores), c(nrow(x), ncomp))
      expect_equal(length(res$eigenvals), ncomp)

      # check computed values
      expect_equivalent(res$scores, x %*% res$loadings)
      expect_equivalent(res$eigenvals, colSums(res$scores^2)/(nrow(x) - 1))
   }

   for (ncomp in ncomp_wrong) {
      expect_error(f(x, ncomp))
   }
}

test_that("pca.svd works correctly", {
   for (x.c in x) {
      tf(x.c, pca.svd)
   }
})

test_that("pca.nipals works correctly", {
   for (x.c in x) {
      tf(x.c, pca.nipals)
   }
})

# testing function for direct comparing two methods
tf <- function(x, f1, f2) {

   ncomp_correct <- c(1, min(ncol(x), nrow(x) - 1), round(min(ncol(x), nrow(x) - 1)/2) + 1)

   for (ncomp in ncomp_correct) {
      res1 <- f1(x, ncomp)
      res2 <- f2(x, ncomp)

      expect_equivalent(abs(res1$scores), abs(res2$scores), tolerance = 10^-3)
      expect_equivalent(abs(res1$loadings), abs(res2$loadings), tolerance = 10^-3)
      expect_equivalent(res1$eigenvals, res2$eigenvals, tolerance = 10^-3)
   }
}

test_that("SVD and NIPALS gives similar results works correctly", {
   for (x.c in x) {
      tf(x.c, pca.nipals, pca.svd)
   }
})

######################################
# Block 2: testing pca.run() method  #
######################################

context("pca: testing pca.run()")

## testing function for comparing pca.run() results with algorithms
tf <- function(x, f1, method, rand, tol = 10^-6) {

   ncomp_correct <- c(1, min(ncol(x), nrow(x) - 1), round(min(ncol(x), nrow(x) - 1)/2) + 1)

   for (ncomp in ncomp_correct) {
      res1 <- f1(x, ncomp)
      res2 <- pca.run(x, ncomp, method, rand)

      expect_equivalent(abs(res1$scores), abs(res2$scores), tolerance = tol)
      expect_equivalent(abs(res1$loadings), abs(res2$loadings), tolerance = tol)
      expect_equivalent(res1$eigenvals, res2$eigenvals, tolerance = tol)

      expect_equal(dim(res2$residuals), dim(x))
      expect_equivalent(res2$residuals, x - tcrossprod(res2$scores, res2$loadings))
   }
}

test_that("pca.run works correctly without rand", {
   for (x.c in x) {
      tf(x.c, pca.svd, method = "svd", rand = NULL)
      tf(x.c, pca.nipals, method = "nipals", rand = NULL)
   }
})

test_that("pca.run works correctly with rand (svd)", {
   for (x.c in x[1:7]) {
      ## the bigger q or p in rand = c(p, q) the more precise results are
      tf(x.c, pca.svd, method = "svd", rand = c(15, 0), tol = 10^-2)
      tf(x.c, pca.svd, method = "svd", rand = c(15, 1), tol = 10^-3)
      tf(x.c, pca.svd, method = "svd", rand = c(15, 2), tol = 10^-4)
   }
})

test_that("pca.run works correctly with rand (nipals)", {
   for (x.c in x) {
      ## the bigger q or p in rand = c(p, q) the more precise results are
      tf(x.c, pca.nipals, method = "svd", rand = c(15, 0), tol = 10^-2)
      tf(x.c, pca.nipals, method = "svd", rand = c(15, 1), tol = 10^-3)
      tf(x.c, pca.nipals, method = "svd", rand = c(15, 2), tol = 10^-4)
   }
})

######################################
# Block 3: testing pca.cal() method  #
######################################

context("pca: testing pca.cal()")


## setup two datasets - one without attributes and with
x1 <- people
x2 <- people
attr(x2, "yaxis.name") <- "Persons"
attr(x2, "yaxis.values") <- 1:32
attr(x2, "xaxis.name") <- "Measurements"
attr(x2, "xaxis.values") <- (1:12)/10

## number of components to test
ncomp_correct <- c(1, 5, 12)

tf <- function(x, ncomp, center, scale) {
   m <- pca.cal(x, ncomp, center = center, scale = scale, method = "svd")
   expect_equal(m$ncomp, ncomp)

   # excluded rows to compute scale and center correctly
   attrs <- attributes(x)
   if (length(attrs$exclrows) > 0) x <- x[-attrs$exclrows, , drop = FALSE]
   expect_equivalent(m$scale,
      (if (is.logical(scale) && scale == TRUE) apply(x, 2, sd) else scale))
   expect_equivalent(m$center,
      (if (is.logical(center) && center == TRUE) apply(x, 2, mean) else center))

   # excluded columns to compute other thinngs correctly
   xt <- prep.autoscale(x, center = center, scale = scale)
   if (length(attrs$exclcols) > 0) xt <- xt[, -attrs$exclcols, drop = F]
   mt <- pca.svd(xt, ncomp)

   expect_equivalent(m$eigenvals, mt$eigenvals)
   expect_equivalent(length(m$eigenvals), ncomp)
   expect_equivalent(names(m$eigenvals), colnames(m$loadings))

   expect_equivalent(dim(m$loadings), c(ncol(x), ncomp))
   expect_equivalent(rownames(m$loadings), colnames(x))

   if (is.null(attr(x, "xaxis.name"))) {
      expect_equivalent(attr(m$loadings, "yaxis.name"), "Variables")
   } else {
      expect_equivalent(attr(m$loadings, "yaxis.name"), attr(x, "xaxis.name"))
   }

   expect_equivalent(attr(m$loadings, "yaxis.values"), attr(x, "xaxis.values"))

   if (length(attrs$exclcols) > 0) {
      expect_equivalent(m$loadings[-attrs$exclcols, , drop = FALSE], mt$loadings)
      expect_equivalent(m$loadings[attrs$exclcols, , drop = FALSE],
      matrix(0, length(attrs$exclcols), ncomp))
      expect_equivalent(attr(m$loadings, "exclrows"), attrs$exclcols)
   } else {
      expect_equivalent(m$loadings, mt$loadings)
   }
}

test_that("pca.cal works fine with full data", {

   for (x in list(x1, x2)) {
      for (ncomp in ncomp_correct) {
         tf(x, ncomp, center = FALSE, scale = FALSE)
         tf(x, ncomp, center = TRUE, scale = FALSE)
         tf(x, ncomp, center = TRUE, scale = TRUE)
         tf(x, ncomp, center = apply(x, 2, median), scale = TRUE)
         tf(x, ncomp, center = apply(x, 2, median), scale = apply(x, 2, IQR))
      }
   }
})

test_that("pca() calculates limit parameters correctly (full data)", {
   for (x in list(x1, x2)) {
      m <- pca(x, 11, center = TRUE, scale = TRUE, method = "svd")
      # parameters computed using DDSimca Toolbox
      # using people data, autoscaled

      momentsParamsQ <- list(
         u0 = c(5.39623603, 3.22376502, 1.65661907, 0.68981822, 0.38111636, 0.22105099, 0.12476396,
            0.07148856, 0.04489746, 0.02151256, 0.00678829),
         Nu = c(10.00000000, 4.00000000, 4.00000000, 3.00000000, 3.00000000, 3.00000000, 5.0000000,
            4.00000000, 2.00000000, 2.00000000, 1.00000000),
         nobj = nrow(x)
      )

      momentsParamsT2 <- list(
         u0 = c(0.96875000, 1.93750000, 2.90625000, 3.87500000, 4.84375000, 5.81250000, 6.78125000,
            7.75000000, 8.71875000, 9.68750000, 10.65625000),
         Nu = c(2.00000000, 10.00000000, 13.00000000, 14.00000000, 13.00000000, 12.00000000, 8.000,
            9.00000000, 10.00000000, 12.00000000, 15.00000000),
         nobj = nrow(x)
      )

      robustParamsQ <- list(
         u0 = c(5.59398247, 3.14412554, 1.68373781, 0.63464055, 0.33473045, 0.21529782, 0.13614951,
            0.06862920, 0.04666029, 0.02054372, 0.00374660),
         Nu = c(18.00000000, 10.00000000, 6.00000000, 3.00000000, 4.00000000, 3.00000000, 5.0000000,
            3.00000000, 4.00000000, 2.00000000, 1.00000000),
         nobj = nrow(x)
      )

      robustParamsT2 <- list(
         u0 = c(1.23877706, 1.84445603, 2.67812507, 3.91673370, 4.51759976, 5.66332668, 5.96116787,
            7.19886896, 8.06726598, 8.95166844, 10.01296469),
         Nu = c(3.00000000, 22.00000000, 17.00000000, 11.00000000, 15.00000000, 18.000000, 16.00000,
            18.00000000, 18.00000000, 10.00000000, 14.00000000),
         nobj = nrow(x)
      )

      # check parameters for Q residuals (moments)
      expect_equivalent(m$limParams$Q$moments$u0, momentsParamsQ$u0, tolerance = 10^-5)
      expect_equivalent(round(m$limParams$Q$moments$Nu), momentsParamsQ$Nu)
      expect_equivalent(m$limParams$Q$moments$nobj, momentsParamsQ$nobj)

      # check parameters for Q residuals (robust)
      expect_equivalent(m$limParams$Q$robust$u0, robustParamsQ$u0, tolerance = 10^-5)
      expect_equivalent(round(m$limParams$Q$robust$Nu), robustParamsQ$Nu)
      expect_equivalent(m$limParams$Q$robust$nobj, robustParamsQ$nobj)

      # check parameters for T2 residuals (moments)
      expect_equivalent(m$limParams$T2$moments$u0, momentsParamsT2$u0, tolerance = 10^-5)
      expect_equivalent(round(m$limParams$T2$moments$Nu), momentsParamsT2$Nu)
      expect_equivalent(m$limParams$T2$moments$nobj, momentsParamsT2$nobj)

      # check parameters for T2 residuals (robust)
      expect_equivalent(m$limParams$T2$robust$u0, robustParamsT2$u0, tolerance = 10^-5)
      expect_equivalent(round(m$limParams$T2$robust$Nu), robustParamsT2$Nu)
      expect_equivalent(m$limParams$T2$robust$nobj, robustParamsT2$nobj)

   }
})

## exclude several columns and rows and redefine components
x1 <- mda.exclrows(x1, c(1, 10, 20))
x2 <- mda.exclcols(x1, c(3, 12))
ncomp_correct <- c(1, 5, 9)

test_that("pca.cal works fine with excluded data", {

   for (x in list(x1, x2)) {
      for (ncomp in ncomp_correct) {
         tf(x, ncomp, center = FALSE, scale = FALSE)
         tf(x, ncomp, center = TRUE, scale = FALSE)
         tf(x, ncomp, center = TRUE, scale = TRUE)
         tf(x, ncomp, center = apply(x, 2, median), scale = TRUE)
         tf(x, ncomp, center = apply(x, 2, median), scale = apply(x, 2, IQR))
      }
   }
})

test_that("pca() calculates limit parameters correctly (excluded data)", {
   for (x in list(x2)) {
      m <- pca(x, 9, center = TRUE, scale = TRUE, method = "svd")

      # parameters computed using DDSimca Toolbox
      # using people data, autoscaled

      momentsParamsQ <- list(
         u0 = c(4.43617757, 2.19049660, 0.44027710, 0.25151771, 0.14625433, 0.08237345, 0.05117557,
            0.02675517, 0.00723308),
         Nu = c(7.00000000, 3.00000000, 5.00000000, 3.00000000, 4.00000000, 3.00000000, 3.00000000,
            2.00000000, 1.00000000),
         nobj = nrow(x) - length(attr(x, "exclrows"))
      )

      momentsParamsT2 <- list(
         u0 = c(0.96551724, 1.93103448, 2.89655172, 3.86206897, 4.82758621, 5.79310345, 6.75862069,
            7.72413793, 8.68965517),
         Nu = c(3.00000000, 13.00000000, 15.00000000, 15.00000000, 14.00000000, 14.000000, 20.00000,
            14.00000000, 17.00000000),
         nobj = nrow(x) - length(attr(x, "exclrows"))
      )

      robustParamsQ <- list(
         u0 = c(4.62537076, 2.40195599, 0.43916799, 0.23975306, 0.13645627, 0.07560735, 0.04351790,
            0.02600730, 0.00787551),
         Nu = c(16.00000000, 3.00000000, 4.00000000, 3.00000000, 5.00000000, 4.00000000, 2.0000000,
            2.00000000, 1.00000000),
         nobj = nrow(x) - length(attr(x, "exclrows"))
      )

      robustParamsT2 <- list(
         u0 = c(1.39125723, 2.06314350, 2.61209799, 3.71494609, 5.03000803, 5.67830732, 6.11739113,
            7.31069105, 8.12330118),
         Nu = c(3.00000000, 19.00000000, 25.00000000, 12.00000000, 9.00000000, 12.0000, 37.0000000,
            21.00000000, 15.00000000),
         nobj = nrow(x) - length(attr(x, "exclrows"))
      )

      # check parameters for Q residuals (moments)
      expect_equivalent(m$limParams$Q$moments$u0, momentsParamsQ$u0, tolerance = 10^-5)
      expect_equivalent(round(m$limParams$Q$moments$Nu), momentsParamsQ$Nu)
      expect_equivalent(m$limParams$Q$moments$nobj, momentsParamsQ$nobj)

      # check parameters for Q residuals (robust)
      expect_equivalent(m$limParams$Q$robust$u0, robustParamsQ$u0, tolerance = 10^-5)
      expect_equivalent(round(m$limParams$Q$robust$Nu), robustParamsQ$Nu)
      expect_equivalent(m$limParams$Q$robust$nobj, robustParamsQ$nobj)

      # check parameters for T2 residuals (moments)
      expect_equivalent(m$limParams$T2$moments$u0, momentsParamsT2$u0, tolerance = 10^-5)
      expect_equivalent(round(m$limParams$T2$moments$Nu), momentsParamsT2$Nu)
      expect_equivalent(m$limParams$T2$moments$nobj, momentsParamsT2$nobj)

      # check parameters for T2 residuals (robust)
      expect_equivalent(m$limParams$T2$robust$u0, robustParamsT2$u0, tolerance = 10^-5)
      expect_equivalent(round(m$limParams$T2$robust$Nu), robustParamsT2$Nu)
      expect_equivalent(m$limParams$T2$robust$nobj, robustParamsT2$nobj)

   }
})


############################################################
# Block 4: testing getQlimits() and getT2Limits() methods  #
############################################################

x <- people
m <- pca(x, 11, center = TRUE, scale = TRUE, method = "svd")

context("pca: testing getQlimits()")

test_that("getQLimits() works fine (jm)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim1 <- rbind(
      c(13.98208375, 8.91523789, 4.86682122, 1.81125667, 0.98188720, 0.57560902, 0.32312846,
         0.18598893, 0.12950361, 0.06958773, 0.02625457),
      c(37.92391547, 27.28432368, 18.46348557, 5.51382242, 2.82348737, 1.68333989, 0.89628978,
         0.47263837, 0.37503680, 0.22537117, 0.09561530)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim2 <- rbind(
      c(11.08714606, 6.88643872, 3.62417784, 1.41387882, 0.77404080, 0.45201939, 0.25555116,
         0.14863431, 0.09998945, 0.05180038, 0.01849205),
      c(29.65164966, 20.68354357, 13.16120605, 4.16057052, 2.16549979, 1.28564550, 0.69582345,
         0.37714904, 0.29044220, 0.17075224, 0.07133111)
   )

   Qlim1 <- ldecomp.getQLimits("jm", 0.05, 0.01, m$limParams, m$res$cal$residuals, m$eigenvals)
   Qlim2 <- ldecomp.getQLimits("jm", 0.10, 0.05, m$limParams, m$res$cal$residuals, m$eigenvals)

   expect_equivalent(Qlim1[1:2, ], expQlim1, tolerance = 10^-5)
   expect_equivalent(Qlim2[1:2, ], expQlim2, tolerance = 10^-5)
})

test_that("getQLimits() works fine (chisq)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim1 <- rbind(
      c(9.25512512, 7.37431769, 3.64879767, 1.55968461, 0.82955551, 0.52538528, 0.23914330,
         0.15607468, 0.11587272, 0.04727562, 0.04666972),
      c(16.87605393, 16.33648881, 8.08326204, 4.19946801, 2.23358734, 1.25780278, 0.52977943,
         0.37365182, 0.31198858, 0.15980027, 0.15775224)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim2 <- rbind(
      c(8.03234581, 6.04655385, 2.99182277, 1.19880757, 0.63761443, 0.42028176, 0.19608497,
         0.12485188, 0.08906230, 0.03329627, 0.03286953),
      c(14.57482653, 13.53946609, 6.69930078, 3.35119154, 1.78241124, 1.02644677, 0.43907419,
         0.30492357, 0.24896808, 0.12254836, 0.12097776)
   )

   Qlim1 <- ldecomp.getQLimits("chisq", 0.05, 0.01, m$limParams, m$res[["cal"]]$residuals, m$eigenvals)
   Qlim2 <- ldecomp.getQLimits("chisq", 0.10, 0.05, m$limParams, m$res[["cal"]]$residuals, m$eigenvals)

   expect_equivalent(Qlim1[1:2, ], expQlim1, tolerance = 10^-5)
   expect_equivalent(Qlim2[1:2, ], expQlim2, tolerance = 10^-5)
})

test_that("getQLimits() works fine (ddmoments)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expQlim1 <- rbind(
      c(11.34617157, 19.08855942, 11.42533747, 6.34336609, 3.34064053, 1.84178155, 0.55799538,
         0.39965752, 0.47200886, 0.25476038, 0.17850637),
      c(19.47199649, 31.77017452, 18.32858290, 10.17605927, 5.41986135, 3.02485080, 0.94234321,
         0.67494207, 0.81004899, 0.42401219, 0.28960906)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expQlim2 <- rbind(
      c(10.00967132, 16.97646011, 10.25820924, 5.69537459, 2.99072482, 1.64367017, 0.49436318,
         0.35408172, 0.41640949, 0.22657181, 0.15980870),
      c(17.04539107, 28.00480337, 16.29329975, 9.04606673, 4.80551079, 2.67446799, 0.82790983,
         0.59298052, 0.70910047, 0.37375866, 0.25678138)
   )

   Qlim1 <- ldecomp.getQLimits("ddmoments", 0.05, 0.01, m$limParams, m$res[["cal"]]$residuals, m$eigenvals)
   Qlim2 <- ldecomp.getQLimits("ddmoments", 0.10, 0.05, m$limParams, m$res[["cal"]]$residuals, m$eigenvals)
   expect_equivalent(Qlim1[1:2, ], expQlim1, tolerance = 10^-5)
   expect_equivalent(Qlim2[1:2, ], expQlim2, tolerance = 10^-5)
})

test_that("getQLimits() works fine (ddrobust)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expQlim1 <- rbind(
      c(10.15325434, 14.52405375, 9.87019680, 5.01044538, 2.52248917, 2.34463394, 0.88961632,
         0.74738500, 0.39573108, 0.21597701, 0.09364917),
      c(15.68388669, 20.95236870, 15.01145008, 8.33916906, 3.96573316, 3.62179177, 1.37420388,
         1.15449701, 0.60641832, 0.37065397, 0.15380476)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expQlim2 <- rbind(
      c(9.20367955, 13.38917413, 8.98186792, 4.45605266, 2.27646590, 2.12535397, 0.80641568,
         0.67748643, 0.35943919, 0.19053642, 0.08357579),
      c(14.06529482, 19.09677672, 13.51144603, 7.35081860, 3.54188909, 3.24801945, 1.23238475,
         1.03535183, 0.54485629, 0.32446297, 0.13598882)
   )

   Qlim1 <- ldecomp.getQLimits("ddrobust", 0.05, 0.01, m$limParams, m$res[["cal"]]$residuals, m$eigenvals)
   Qlim2 <- ldecomp.getQLimits("ddrobust", 0.10, 0.05, m$limParams, m$res[["cal"]]$residuals, m$eigenvals)
   expect_equivalent(Qlim1[1:2, ], expQlim1, tolerance = 10^-5)
   expect_equivalent(Qlim2[1:2, ], expQlim2, tolerance = 10^-5)
})

context('pca: testing getT2limits()')

test_that("getT2Limits() works fine (chisq and jm)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "tsqlim" function from PLS_Toolbox
   expT2lim1 <- rbind(
      c(4.15961510, 6.85271430, 9.40913034, 12.01947857, 14.76453307, 17.69939358, 20.87303997,
         24.33584211, 28.14388532, 32.36253393, 37.07020868),
      c(16.43763399, 22.07592005, 27.49293159, 33.11628022, 39.14366436, 45.72590957, 53.01070123,
         61.16173704, 70.37250826, 80.88041709, 92.98424615)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "tsqlim" function from PLS_Toolbox
   expT2lim2 <- rbind(
      c(2.87478394, 5.14334644, 7.32156708, 9.55303105, 11.90044572, 14.40708180, 17.11153965,
         20.05346964, 23.27685391, 26.83267841, 30.78172763),
      c(11.96070234, 16.61285114, 21.04123968, 25.60304588, 30.45666797, 35.71776166, 41.49572457,
         47.90883728, 55.09430867, 63.21788003, 72.48520544)
   )

   T2lim11 <- ldecomp.getT2Limits("chisq", 0.05, 0.01, m$limParams)
   T2lim21 <- ldecomp.getT2Limits("chisq", 0.10, 0.05, m$limParams)
   T2lim12 <- ldecomp.getT2Limits("jm", 0.05, 0.01, m$limParams)
   T2lim22 <- ldecomp.getT2Limits("jm", 0.10, 0.05, m$limParams)

   expect_equivalent(T2lim11[1:2, ], expT2lim1, tolerance = 10^-5)
   expect_equivalent(T2lim21[1:2, ], expT2lim2, tolerance = 10^-5)
   expect_equivalent(T2lim12[1:2, ], expT2lim1, tolerance = 10^-5)
   expect_equivalent(T2lim22[1:2, ], expT2lim2, tolerance = 10^-5)
})

test_that("getT2Limits() works fine (ddmoments)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expT2lim1 <- rbind(
      c(10.18450977, 4.58893048, 6.16731299, 7.63572085, 9.79787223, 12.10733682, 18.95532441,
         19.25620258, 18.33211759, 19.12054367, 18.68127639),
      c(17.47838355, 7.63761785, 9.89363402, 12.24926117, 15.89608595, 19.88449041, 32.01177265,
         32.51989603, 31.46109039, 31.82340773, 30.30853722)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expT2lim2 <- rbind(
      c(8.98484522, 4.08117729, 5.53730578, 6.85571191, 8.77159317, 10.80501016, 16.79371331,
         17.06028019, 16.17272140, 17.00490539, 16.72450430),
      c(15.30022271, 6.73241458, 8.79500317, 10.88905154, 14.09423741, 17.58117559, 28.12442515,
         28.57084459, 27.54040088, 28.05172743, 26.87301266)
   )

   T2lim1 <- ldecomp.getT2Limits("ddmoments", 0.05, 0.01, m$limParams)
   T2lim2 <- ldecomp.getT2Limits("ddmoments", 0.10, 0.05, m$limParams)
   expect_equivalent(T2lim1[1:2, ], expT2lim1, tolerance = 10^-5)
   expect_equivalent(T2lim2[1:2, ], expT2lim2, tolerance = 10^-5)
})

test_that("getT2Limits() works fine (ddrobust)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expT2lim1 <- rbind(
      c(13.49051625, 3.87287602, 5.54095383, 8.43336946, 9.07842620, 10.27911626, 12.17217082,
         13.06617385, 15.20430507, 18.82185388, 17.87728459),
      c(20.83900601, 5.58700262, 8.42716245, 14.03613618, 14.27265426, 15.87830752, 18.80253774,
         20.18351784, 23.29907765, 32.30156273, 29.36076688)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expT2lim2 <- rbind(
      c(12.22882677, 3.57025747, 5.04226172, 7.50023912, 8.19298966, 9.31777034, 11.03377852,
         11.84417065, 13.80994172, 16.60477015, 15.95431303),
      c(18.68840097, 5.09220428, 7.58508672, 12.37258654, 12.74724151, 14.23965123, 16.86209814,
         18.10055980, 20.93381524, 28.27613333, 25.95977002)
   )

   T2lim1 <- ldecomp.getT2Limits("ddrobust", 0.05, 0.01, m$limParams)
   T2lim2 <- ldecomp.getT2Limits("ddrobust", 0.10, 0.05, m$limParams)
   expect_equivalent(T2lim1[1:2, ], expT2lim1, tolerance = 10^-5)
   expect_equivalent(T2lim2[1:2, ], expT2lim2, tolerance = 10^-5)
})

#########################################
# Block 5: testing categorize() method  #
#########################################

tf <- function(m, outliers, extremes) {

   expect_equivalent(which(categorize(m, m$calres, 1) == "outlier"), outliers[[1]])
   expect_equivalent(which(categorize(m, m$calres, 3) == "outlier"), outliers[[2]])
   expect_equivalent(which(categorize(m, m$calres, 5) == "outlier"), outliers[[3]])

   expect_equivalent(which(categorize(m, m$calres, 1) == "extreme"), extremes[[1]])
   expect_equivalent(which(categorize(m, m$calres, 3) == "extreme"), extremes[[2]])
   expect_equivalent(which(categorize(m, m$calres, 5) == "extreme"), extremes[[3]])
}

context("pca: testing categorize()")

x <- people
m <- pca(x, 11, scale = TRUE)

test_that("categorize() works well (chisq)", {

   ## outliers and extremes computed using PLS toolbox

   m <- setDistanceLimits.pca(m, "chisq", 0.05, 0.01)
   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(1, 18, 27, 28), c(7), c(4, 20, 28, 29, 31))
   tf(m, outliers, extremes)

   m <- setDistanceLimits.pca(m, "chisq", 0.10, 0.05)
   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(1, 17, 18, 27, 28), c(4, 7, 19), c(4, 19, 20, 28, 29, 31))
   tf(m, outliers, extremes)

})

test_that("categorize() works well (ddmoments)", {

   ## outliers and extremes computed using DD SIMCA toolbox

   m <- setDistanceLimits.pca(m, "ddmoments", 0.05, 0.01)

   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(28), c(1, 28), c(1, 4, 28))
   tf(m, outliers, extremes)

   m <- setDistanceLimits.pca(m, "ddmoments", 0.10, 0.05)
   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(18, 28), c(1, 18, 28), c(1, 4, 19, 20, 27, 28))
   tf(m, outliers, extremes)

})

test_that("categorize() works well (ddrobust)", {

   ## outliers and extremes computed using DD SIMCA toolbox

   m <- setDistanceLimits.pca(m, "ddrobust", 0.05, 0.01)

   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(28), c(1, 18, 28), c(1, 4, 19, 20, 28))
   tf(m, outliers, extremes)

   m <- setDistanceLimits.pca(m, "ddrobust", 0.10, 0.05)
   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(17, 18, 27, 28), c(1, 7, 17, 18, 19, 27, 28), c(1, 4, 7, 19, 20, 27, 28))
   tf(m, outliers, extremes)

})

x <- people
x <- mda.exclrows(x, c(1, 10, 20))
x <- mda.exclcols(x, c(3, 12))
m <- pca(x, 11, scale = TRUE)
test_that("categorize() works well for excluded data (chisq)", {

   ## outliers and extremes computed using PLS toolbox

   m <- setDistanceLimits.pca(m, "chisq", 0.05, 0.01)
   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(16, 25), c(9, 18, 28), c(17, 19, 28))
   tf(m, outliers, extremes)

   m <- setDistanceLimits.pca(m, "chisq", 0.10, 0.05)
   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(15, 16, 25), c(9, 11, 18, 26, 28), c(17, 19, 28))
   tf(m, outliers, extremes)

})

test_that("categorize() works well for excluded data (ddmoments)", {

   ## outliers and extremes computed using DD SIMCA toolbox

   m <- setDistanceLimits.pca(m, "ddmoments", 0.05, 0.01)

   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(16), c(16), c(28))
   tf(m, outliers, extremes)

   m <- setDistanceLimits.pca(m, "ddmoments", 0.10, 0.05)
   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(15, 16, 25), c(16, 25), c(9, 16, 17, 28))
   tf(m, outliers, extremes)

})

test_that("categorize() works well for excluded data (ddrobust)", {

   ## outliers and extremes computed using DD SIMCA toolbox

   m <- setDistanceLimits.pca(m, "ddrobust", 0.05, 0.01)

   outliers <- list(integer(0), integer(0), integer(0))
   extremes <- list(c(15, 16, 25), c(15, 16, 22, 25), c(17, 28))
   tf(m, outliers, extremes)

   m <- setDistanceLimits.pca(m, "ddrobust", 0.10, 0.05)
   outliers <- list(integer(0), c(16), integer(0))
   extremes <- list(c(15, 16, 25), c(15, 22, 23, 25), c(17, 28))
   tf(m, outliers, extremes)

})


##################################
# Block 6: testing pca()         #
##################################

tf <- function(x, ncomp) {
   m <- pca(x, ncomp = ncomp, scale = TRUE)
   expect_equal(m$calres, predict(m, x))
   expect_equal(dim(m$Qlim), c(4, ncomp))
   expect_equal(dim(m$T2lim), c(4, ncomp))
   expect_equal(m$exclcols, attr(x, "exclcols"))
   expect_equal(m$exclrows, attr(x, "exclrows"))
   expect_output(summary(m))
   expect_output(print(m))
}

context("pca: testing pca()")

# full data
x <- people

tf(x, 1)
tf(x, 5)
tf(x, 12)

# excluded data
x <- people
x <- mda.exclrows(x, c(1, 10, 20))
x <- mda.exclcols(x, c(3, 12))

tf(x, 1)
tf(x, 5)
tf(x, 10)


#########################################
# Block 7: testing pca() with test set  #
#########################################


tf <- function(x, x.test, ncomp) {
   m <- pca(x, ncomp = ncomp, scale = TRUE, x.test = x.test)

   # results for calibration and test set are correct
   expect_equal(m$calres, predict(m, x))
   expect_equal(m$testres, predict(m, x.test))

   # model$res returns correct list
   expect_equivalent(m$res , list("cal" = m$calres, "test" = m$testres))

   # dimension of limits is correct
   expect_equal(dim(m$Qlim), c(4, ncomp))
   expect_equal(dim(m$T2lim), c(4, ncomp))

   # excluded data is stored in model
   expect_equal(m$exclcols, attr(x, "exclcols"))
   expect_equal(m$exclrows, attr(x, "exclrows"))

   # loadings matrix has zeros for excluded variables
   if (length(attr(x, "exclcols")) > 0) {
      nexcols <- length(attr(x, "exclcols"))
      expect_equivalent(m$loadings[attr(x, "exclcols"), ], matrix(0, nrow = nexcols, ncol = ncomp))
   }

   # print and summery produce some output
   expect_output(summary(m))
   expect_output(print(m))
}

context("pca: testing pca() for test set")

# full data
ind <- seq(1, 32, by = 4)
x <- people[-ind, ]
x.test <- people[ind, ]

tf(x, x.test, 1)
tf(x, x.test, 5)
tf(x, x.test, 12)

# excluded data
x <- people[-ind, ]
x <- mda.exclrows(x, c(1, 10, 20))
x <- mda.exclcols(x, c(3, 12))
x.test <- people[ind, ]

tf(x, x.test, 1)
tf(x, x.test, 5)
tf(x, x.test, 10)


#########################################
# Block 8: testing pca.mvreplace()      #
#########################################


data(people)
context("pca: test pca.mvreplace() function")

test_that("for simple dataset it can reconstruct values in full", {
   data <- odata <- matrix(c(1:10, 2 * (1:10), 3 * (1:10)), ncol = 3)
   data[3, 1] <- data[7, 3] <- NA
   expect_equivalent(odata, pca.mvreplace(data), tolerance = 10^-2)
   expect_equivalent(odata, pca.mvreplace(data, scale = TRUE), tolerance = 10^-2)
})

test_that("if too many values are missing it returns error", {
   data <- people
   data[seq(1, 32, by = 4), ] <- NA
   expect_error(pca.mvreplace(data))
   expect_error(pca.mvreplace(data, scale = T))
})

test_that("for people data it can partially reconstruct the values", {
   data <- odata <- people
   data[c(3, 10), 1] <- NA
   data[c(25, 32), 2] <- NA
   data[c(5, 10, 15), 11] <- NA
   data[c(14), 12] <- NA

   expect_equivalent(
      odata,
      pca.mvreplace(data, scale = T),
      tolerance = 0.1
   )
})

#########################################
# Block 9: testing bug fixes            #
#########################################

data(people)
context("pca: test recent bug fixes")

test_that("row and column names are kept when data frame is used as data source", {
   m1 <- pca(people)
   m2 <- pca(as.data.frame(people))
   m3 <- pca(people, 5, scale = TRUE)
   m4 <- pca(as.data.frame(people), 5, scale = TRUE)

   expect_equal(rownames(m1$res$cal$scores), rownames(people))
   expect_equal(rownames(m2$res$cal$scores), rownames(people))
   expect_equal(rownames(m3$res$cal$scores), rownames(people))
   expect_equal(rownames(m4$res$cal$scores), rownames(people))

   expect_equal(rownames(m1$loadings), colnames(people))
   expect_equal(rownames(m2$loadings), colnames(people))
   expect_equal(rownames(m3$loadings), colnames(people))
   expect_equal(rownames(m4$loadings), colnames(people))

})

###########################################
# Block 10: testing preprocessing in PCA  #
###########################################

context("pca: testing preprocessing")

test_that("preprocessing is applied to training set and becomes part of a model", {
   data(simdata)
   xc <- simdata$spectra.c
   xt <- simdata$spectra.t

   xc <- mda.exclrows(xc, c(20, 30))
   xc <- mda.exclcols(xc, c(20, 100))
   xt <- mda.exclrows(xc, c(10, 40))

   p <- list(
      prep("savgol", width = 7, porder = 2, dorder = 2),
      prep("norm", type = "snv")
   )

   pm <- prep.fit(p, xc)
   xcp <- prep.apply(pm, xc)
   xtp <- prep.apply(pm, xt)

   m1 <- pca(xc, 5, prep = p)
   m2 <- pca(xcp, 5)

   r1 <- predict(m1, xt)
   r2 <- predict(m2, xtp)

   expect_equivalent(m1$loadings, m2$loadings)
   expect_equivalent(m1$eigenvals, m2$eigenvals)
   expect_equivalent(m1$T2lim, m2$T2lim)
   expect_equivalent(m1$Qlim, m2$Qlim)
   expect_equivalent(m1$calres$scores, m2$calres$scores)
   expect_equivalent(m1$calres$residuals, m2$calres$residuals)

   expect_equivalent(r1$scores, r2$scores)
   expect_equivalent(r1$residuals, r2$residuals)

})


###########################################
# Block 11: testing JSON methods in PCA   #
###########################################

context("pca: testing JSON methods")

# expected vector for case 1: People, A = 5, center, scale
v.exp1 <- c(12,5,173.125,64.46875,0,39.90625,34.4375,27437.5,249.5,131.59375,0,81.5,0,115.125,10.057095072152922,15.191220671291395,1.016001016001524,3.8967262934996616,9.517174628177514,8929.608234732723,90.59694504569623,49.496731319253016,1.016001016001524,7.3176763454510745,1.016001016001524,12.164862143392416,0.37528576400253966,0.3811357231221713,-0.3377734617313739,0.37769704201578835,0.1429453571041515,0.19046640122719039,0.3246657899569954,-0.12414873126517677,-0.3517825945127897,0.3649037252959106,-0.14412124582660282,-0.0440672793224172,-0.1354585410268543,-0.11144717582485403,0.15016339845720608,-0.15080609754202087,0.061462758699004895,0.286893202417765,0.3082852282854548,-0.5542002831368638,0.23167088884816056,-0.11242522260275464,-0.5952589880490593,-0.12260413937914583,0.07237117350560711,0.06810859929961677,-0.07901752221168583,0.0012780596082086185,-0.7202286133731134,-0.5859268287185029,0.18756640487343607,-0.21158206832514032,0.052092420014533486,0.13537492864570158,-0.1302454774560239,0.06234357833906442,-0.07431852638691923,0.03313762191848652,-0.11431601174307382,-0.06569226677135552,0.05531501931059969,0.08459843779894861,0.039697043799153016,-0.1251825244452304,-0.05139188333241813,-0.08102086608087655,-0.022163848756660727,0.9689350446949387,-0.18590711265802426,-0.10026968601336853,-0.6601722295988279,-0.1517388171557198,0.028543074436581533,-0.06321381310217518,-0.23112668047508714,-0.41519416601127884,-0.3128402579814745,-0.33560912070844556,0.15112219981533864,-0.18048721502797108,6.429691844266282,2.24255071844727,1.6176990433911904,0.9979879761227338,0.31865997784265254,5.396236025867017,3.223765017371203,1.6566190690859912,0.6898182172170945,0.3811163636820247,10,4,4,3,3,5.593988468173895,3.144122518950645,1.6837411648431408,0.6346374607315604,0.33472986127762133,18,10,6,3,4,0.9687499999999998,1.937500000000004,2.9062500000000098,3.8750000000000124,4.843750000000016,2,10,13,14,13,1.238771043821127,1.844456172046646,2.6781248187392355,3.916736065678244,4.517604926884154,3,22,17,11,15,14.389962735645376,8.596706712989873,4.417650850895974,1.839515245912252,1.0163103031520657,3,1,1,3,6,19.35930713762732,6.874566199215227,2.0421381833736065,1.6530767247400835,1.0809996425225796,1,1,2,3,7,2.583333333333339,5.166666666666698,7.750000000000045,10.33333333333338,12.916666666666723,4,5,5,4,6,4.07527851608531,5.016654209223919,6.764958860491086,8.009701813600438,12.231079713226494,4,36,3,3,6
)

# expected vector for case 2: People, A = 5, center, scale, excluded outliers: 1, 18
v.exp2 <- c(12,5,171.83333333333334,62.9,0.06666666666666667,39.43333333333333,34.233333333333334,27216.666666666668,242.16666666666666,130.6,0.06666666666666667,80.5,0,115.5,8.855480801734318,14.312642108916544,1.01483252680985,3.530002767736859,9.409032099927774,8391.994680208258,87.35740433388588,50.339947803186945,1.01483252680985,6.366669675632677,1.0170952554312156,12.232658379701558,0.37568343086065487,0.38244367204802776,-0.3416012902152422,0.37809405743410734,0.12729431443575767,0.18120811956046676,0.31484036483660816,-0.14369855357005445,-0.3576133200431884,0.36652901430302803,-0.14765864153375458,-0.010217906573563303,0.13303530050003334,0.09781266535931817,-0.13631806788514908,0.1496782844412064,0.08778201168818146,-0.16341275403830627,-0.36398298129971224,0.5684730553390754,-0.24458534891904193,0.07420134431772461,0.6063596941674132,0.09323915837774206,0.11886304300404579,0.07853451050580801,-0.08563870926226136,0.018754049996189154,-0.7130509547168177,-0.6480398830659602,0.10624365853950413,-0.08047851879007488,0.02015996188108685,0.15923367799966606,-0.007183246790827171,0.006554114131705765,0.04975918527599104,-0.03302166586801164,0.0791404837935771,0.07082952894635605,0.014755301663878815,-0.0358463269637331,-0.05085708507926133,0.09281922027731736,0.001578124828469866,0.05520129203849313,0.010643197030011164,-0.9845232730183615,-0.1330842544720288,-0.14597746945919735,-0.7239560227552181,-0.16486433205986334,0.025374679003891735,-0.044731467580512975,-0.2083392264235022,-0.4330271626719794,-0.18888319543535773,-0.31027341130705777,0.1754458195606194,-0.11574267348124469,6.380594383603337,2.2231523206093438,1.6904669204512461,1.006851916275124,0.2835530869022018,5.432092095850103,3.2830448525943887,1.6489268294915231,0.6756366437589066,0.40153532642011075,9,5,5,2,4,5.908391935903243,3.0181380565508498,1.628287190896415,0.5971479418174565,0.4167720898695836,13,9,8,5,3,0.9666666666666675,1.9333333333333422,2.9000000000000177,3.866666666666687,4.833333333333358,3,26,18,16,14,1.263363484309913,1.9911712311324234,2.72147012577123,3.600621577806794,4.50171453388795,2,37,24,24,15,13.580230239625253,8.207612131485973,4.122317073728808,1.6890916093972663,1.0038383160502768,3,1,1,2,6,18.600359269261418,6.730239685206099,1.8324476846250017,1.2621648370919205,0.9687040662656421,1,1,2,7,11,2.4166666666666705,4.8333333333333615,7.250000000000038,9.66666666666671,12.08333333333339,4,4,5,4,6,3.782853508758194,4.798912943211668,6.397395588204569,7.804847670649446,11.765882254884279,4,14,4,4,5
)

# expected vector for case 3: People, A = 5, no center, no scale, excluded outliers: 1, 18, preprocessing: median centering, pareto scaling
v.exp3 <- c(12,5,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,-0.00651284588370324,-0.01133568001657788,-0.0005258754388312916,-0.004888331165770438,-0.02607419591820334,-0.9987088622635989,-0.034950646871687675,0.02074170228900641,0.00026811522396751107,-0.005839053029597615,0.003962036529544309,0.0020068691408697266,-0.16946760520513476,-0.2205059352515631,0.05718517081742274,-0.09857317901518337,0.07616210389958658,0.04125447465512403,-0.7909011125601083,0.4969225515458174,0.046780658426530564,-0.15094456721061905,0.06824424603185021,0.009049231984647747,0.35733402285338434,0.42575411619329456,-0.12820647737424773,0.24708341425499508,0.1276387722052461,-0.004013492178614407,0.15362354305472292,0.6856397836600245,-0.15724896188348214,0.26802234808334774,0.06968827884208464,0.046142959529780314,0.03954654968818249,-0.03361519315178758,0.048674776361779915,0.021874699613568464,-0.02123158831272278,-0.0011469222680858743,0.016841350856262226,0.059074252375181754,0.034574925506410996,0.04835118288136612,-0.015241292237606115,-0.9932106017083606,0.28824077970675455,0.37022880995132984,-0.1764118008830507,0.17555257800444155,0.22096435708705747,-0.005187744525341536,-0.5497009949800385,-0.5245504222019238,-0.19894453758334568,0.19859836172820636,0.10111068696467598,-0.04987937065554586,9370.944068259547,116.77179392960575,33.968389458233844,12.444978752322351,6.560543623349179,168.3830326831894,55.50363188457031,22.667522074944248,10.637375947699304,4.295517111795095,5,5,5,3,3,171.33560717914833,55.96266085441812,25.100294181156578,9.64910833936484,4.587567519411123,7,5,6,3,2,0.9666666666666669,1.933333333333336,2.9000000000000044,3.8666666666666734,4.833333333333341,1,4,7,9,11,1.0478335426191867,2.1923300357709428,2.8698480275381346,3.6645625598420137,4.4492593777205425,1,4,10,20,9,420.95758170797353,138.75907971142576,56.66880518736062,26.59343986924827,10.738792779487738,1,2,1,3,1,279.86272970390064,110.18604438748238,34.49128942336145,33.18149311946683,12.709877996016669,1,2,3,3,1,2.4166666666666687,4.83333333333334,7.2500000000000115,9.66666666666668,12.083333333333352,1,1,1,1,2,0.0069141664319162365,1.6354308963160635,6.35238263602421,11.675049099206504,16.59596894436647,1,1,1,1,1)

# expected vector for case 4: People, A = 5, no center, no scale, excluded outliers: 1, 18, preprocessing: median centering, pareto scaling, excluded variables: 2, 7, 8
v.exp4 <- c(9,5,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0.006439085455613531,0.0005530714130669823,0.004846691510345359,0.026137205755654785,0.9995992896687314,-0.00024725031505873005,0.005772797111865385,-0.003929750974699295,-0.00200478902348955,-0.6527918842836696,0.23984775742079564,-0.41390007631831655,-0.03846214833252279,0.010305320547422654,0.2520730372064941,-0.5284607909968222,0.024537576783166218,0.0048252491609501276,0.012122140012400818,0.06131374462793386,-0.0016493747635598697,-0.05744481103670594,-0.0009114675549412477,0.05308929111056793,0.03291133717853047,-0.03546303120375522,-0.993796239371183,0.08785146554955976,0.22853210456310086,-0.12010213987393817,-0.8184788302250368,0.018858811613371513,0.3022716245881942,0.2815538405649365,-0.27411474246227646,0.09791746870315325,0.21112100979685539,0.5795189821621545,-0.10450746628509879,0.2555651854887526,-0.0070941745210474794,0.3896342451827627,0.27749605424463336,0.5597686164592998,0.03376621932484878,9354.453615457029,19.554594181724475,12.497407958413431,4.200395556954394,1.2315830709943258,36.982336873284765,18.07956249761766,5.998734804484688,1.9383524327620996,0.7478221308009176,6,3,2,2,4,34.45177900466982,18.975968943597355,5.047895565922706,1.6436459267347594,0.7980768840069301,6,4,2,4,4,0.9666666666666702,1.9333333333333396,2.9000000000000106,3.86666666666668,4.833333333333347,1,3,6,7,12,1.051793299956409,1.743635843683613,2.9287103846480305,3.782628404324857,4.430159278242252,1,2,7,5,10,123.27445624428252,60.265208325392194,19.995782681615623,6.461174775873667,2.4927404360030585,2,1,1,3,4,118.10893361526027,22.20382607266186,19.102541478833622,5.589133549014062,3.1942254535888455,2,4,5,11,5,3.222222222222223,6.444444444444465,9.666666666666696,12.888888888888934,16.11111111111116,1,1,1,3,7,0.0009855834490419266,5.314329537810869,9.793316124306301,14.355841046777572,14.622772945567407,2,1,1,2,7)

compareResults <- function(rr, rjs) {
   expect_equal(abs(rjs$scores), abs(rr$scores))
   expect_equal(rjs$Q, rr$Q)
   expect_equal(rjs$T2, rr$T2)
   expect_equal(rjs$expvar, rr$expvar)
   expect_equal(rjs$cumexpvar, rr$cumexpvar)

   rr.res <- mda.purgeCols(rr$residuals)
   rjs.res <- mda.purgeCols(rjs$residuals)

   attr(rr.res, "exclcols") <- NULL
   attr(rjs.res, "exclcols") <- NULL
   attr(rr.res, "prep:center") <- NULL
   attr(rjs.res, "prep:center") <- NULL
   attr(rr.res, "prep:scale") <- NULL
   attr(rjs.res, "prep:scale") <- NULL

   expect_equal(rr.res, rjs.res)
}

compareModels <- function(mr, mjs) {

   if (is.null(attr(mr$loadings, "yaxis.values"))) {
      attr(mjs$loadings, "yaxis.values") <- NULL
   }

   expect_equal(abs(mjs$loadings), abs(mda.purge(mr$loadings)))
   expect_equal(mr$eigenvals, mjs$eigenvals)

   # we need to round DoF in R as they are rounded in web version
   mr$limParams$Q$moments$Nu <- round(mr$limParams$Q$moments$Nu)
   mr$limParams$Q$robust$Nu <- round(mr$limParams$Q$robust$Nu)
   mr$limParams$T2$moments$Nu <- round(mr$limParams$T2$moments$Nu)
   mr$limParams$T2$robust$Nu <- round(mr$limParams$T2$robust$Nu)
   testList(mr$limParams, mjs$limParams)

   # same for Qlim/T2lim
   expect_equal(mr$Qlim, mjs$Qlim)
   expect_equal(mr$T2lim, mjs$T2lim)

   # prep
   ncr <- ncol(mr$loadings)
   ncjs <- ncol(mjs$loadings)
   if (!is.null(mr$prep) && ncr == ncjs) {
      for (n in seq_along(mr$prep)) {
         if (!is.list(mr$prep[[n]])) next
         testList(mr$prep[[n]]$params, mjs$prep[[n]]$params)
      }
   }
}

test_that("asvector.pca works correctly", {
   data(people)

   # case 1
   m <- pca(people, 5, scale = TRUE)
   v <- asvector(m)
   expect_equivalent(abs(v), abs(v.exp1), tolerance = 0.0001)

   writeJSON(m, "./jsonfiles/pca-model-1-r.json")
   m1 <- readJSON("./jsonfiles/pca-model-1-r.json")
   m2 <- readJSON("./jsonfiles/pca-model-1.json")
   compareModels(m, m1)

   r <- predict(m, people)
   r1 <- predict(m1, people)
   r2 <- predict(m2, people)
   compareResults(r, r1)
   compareResults(r, r2)


   # case 2
   m <- pca(people, 5, scale = TRUE, exclrows = c(1, 18))
   v <- asvector(m)
   expect_equivalent(abs(v), abs(v.exp2), tolerance = 0.0001)

   writeJSON(m, "./jsonfiles/pca-model-2-r.json")
   m1 <- readJSON("./jsonfiles/pca-model-2-r.json")
   m2 <- readJSON("./jsonfiles/pca-model-2.json")
   compareModels(m, m1)

   r <- predict(m, people)
   r1 <- predict(m1, people)
   r2 <- predict(m2, people)
   compareResults(r, r1)
   compareResults(r, r2)


   # case 3
   p <- list(
      prep("center", type = "median"),
      prep("scale", type = "pareto")
   )

   m <- pca(people, 5, center = FALSE, scale = FALSE, exclrows = c(1, 18), prep = p)
   v <- asvector(m)
   expect_equivalent(abs(v), abs(v.exp3), tolerance = 0.0001)

   writeJSON(m, "./jsonfiles/pca-model-3-r.json")
   m1 <- readJSON("./jsonfiles/pca-model-3-r.json")
   m2 <- readJSON("./jsonfiles/pca-model-3.json")
   compareModels(m, m1)

   r <- predict(m, people)
   r1 <- predict(m1, people)
   r2 <- predict(m2, people)
   compareResults(r, r1)
   compareResults(r, r2)

   # case 4: without comparing with web-app
   m <- pca(people, 5, center = FALSE, scale = FALSE, exclrows = c(1, 18), prep = p, exclcols = c(2, 7, 8))
   v <- asvector(m)
   expect_equivalent(abs(v), abs(v.exp4), tolerance = 0.0001)

   writeJSON(m, "./jsonfiles/pca-model-4-r.json")
   m1 <- pca.readJSON("./jsonfiles/pca-model-4-r.json")
   compareModels(m, m1)

   r <- predict(m, people)
   r1 <- predict(m1, people)
   compareResults(r, r1)

   # case 4: with comparing with web-app (to do this we remove columns)
   people_exclvars <- people[, -c(2, 7, 8)]
   m <- pca(people_exclvars, 5, center = FALSE, scale = FALSE, exclrows = c(1, 18), prep = p)
   v <- asvector(m)
   expect_equivalent(abs(v), abs(v.exp4), tolerance = 0.0001)

   writeJSON(m, "./jsonfiles/pca-model-4-r.json")
   m1 <- pca.readJSON("./jsonfiles/pca-model-4-r.json")
   m2 <- pca.readJSON("./jsonfiles/pca-model-4.json")
   compareModels(m, m1)
   compareModels(m, m2)

   r <- predict(m, people_exclvars)
   r1 <- predict(m1, people_exclvars)
   r2 <- predict(m2, people_exclvars)
   compareResults(r, r1)
   compareResults(r, r2)

})


