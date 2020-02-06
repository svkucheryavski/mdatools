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
      tf(x.c, pca.svd, method = "svd", rand = c(5, 0), tol = 10^-2)
      tf(x.c, pca.svd, method = "svd", rand = c(5, 1), tol = 10^-3)
      tf(x.c, pca.svd, method = "svd", rand = c(5, 2), tol = 10^-4)
   }
})

test_that("pca.run works correctly with rand (nipals)", {
   for (x.c in x) {
      ## the bigger q or p in rand = c(p, q) the more precise results are
      tf(x.c, pca.nipals, method = "svd", rand = c(5, 0), tol = 10^-2)
      tf(x.c, pca.nipals, method = "svd", rand = c(5, 1), tol = 10^-3)
      tf(x.c, pca.nipals, method = "svd", rand = c(5, 2), tol = 10^-4)
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
