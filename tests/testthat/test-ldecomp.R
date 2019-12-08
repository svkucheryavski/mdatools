#####################################################
# Tests for basic functionality of ldecomp() class  #
#####################################################

# Prepare 3 datasets

data(people)

## plain data, no names no attributes
X1 <- list()
X1$data <- people
rownames(X1$data) <- NULL
colnames(X1$data) <- NULL
X1$exp_exclrows <- NULL
X1$exp_exclcols <- NULL

## same data but with names and attributes
X2 <- list()
X2$data <- people
X2$exp_exclrows <- NULL
X2$exp_exclcols <- NULL
attr(X2$data, "name") <- "People"
attr(X2$data, "xaxis.name") <- "Properties"
attr(X2$data, "xaxis.values") <- 12:1
attr(X2$data, "yaxis.name") <- "Values"
attr(X2$data, "yaxis.values") <- 32:1

## for full PCA model of People data we know some results - use them for double check
Q_exp <- NULL
T2_exp <- NULL
expvar_exp <- NULL

## data with names and attributes and excluded rows
X3 <- X2
X3$exp_exclrows <- c(1, 10, 20)
X3$data <- mda.exclrows(X3$data, X3$exp_exclrows)

## data with names and attributes, excluded rows and columns
X4 <- X3
X4$data <- mda.exclcols(X4$data, c("Hairleng", "IQ"))

# function to get scores, loadings and residuals from data
getPCARes <- function(X, ncomp) {
   rows_excluded <- attr(X, "exclrows")
   cols_excluded <- attr(X, "exclcols")

   if (is.null(rownames(X))) {
      rownames(X) <- paste0('O', 1:nrow(X))
   }

   if (is.null(colnames(X))) {
      colnames(X) <- paste0('X', 1:ncol(X))
   }

   rownames <- rownames(X)
   colnames <- colnames(X)

   X_cal <- X

   # remove excluded rows from calibration data
   if (length(rows_excluded) > 0) {
      X_cal <- X_cal[-rows_excluded, , drop = FALSE]
   }

   # get mean and center and do autoscaling of the calibration data
   m <- apply(X_cal, 2, mean)
   s <- apply(X_cal, 2, sd)
   X_cal <- scale(X_cal, center = m, scale = s)

   # remove excluded columns
   if (length(cols_excluded) > 0) {
      X_cal <- X_cal[, -cols_excluded, drop = FALSE]
   }

   # find loadings
   loadings_visible <- svd(X_cal)$v[, 1:ncomp, drop = F]
   loadings <- matrix(0, nrow = ncol(X), ncol = ncomp)
   if (length(cols_excluded) > 0) {
      loadings[-cols_excluded, ] <- loadings_visible
   } else {
      loadings <- loadings_visible
   }

   # eigenvalues using only visible rows
   tnorm <- sqrt(colSums((X_cal %*% loadings_visible)^2) / (nrow(X_cal) - 1))

   X <- scale(X, center = m, scale = s)
   scores <- X %*% loadings
   residuals <- X - tcrossprod(scores, loadings)

   scores <- mda.setattr(scores, mda.getattr(X), type = 'row')
   residuals <- mda.setattr(residuals, mda.getattr(X))

   attr(loadings, "exclrows") <- attr(X, "exclcols")
   attr(loadings, "yaxis.name") <- attr(X, "xaxis.name")
   attr(loadings, "yaxis.values") <- attr(X, "xaxis.values")

   rownames(scores) <- rownames(residuals) <- rownames
   rownames(loadings) <- colnames(residuals) <- colnames
   colnames(scores) <- colnames(loadings) <- paste("Comp", 1:ncomp)

   return(list(scores = scores, loadings = loadings, residuals = residuals,
      tnorm = tnorm, totvar = sum(X_cal^2)))
}

##################################################
# Block 1. Computing distances and tnorm values  #
##################################################

context("ldecomp: computing distances and tnorm")

## shortcut for testing function
tf <- function(X) {

   for (A in c(1, 5, ncol(X$data) - length(attr(X$data, "exclcols")))) {
      m <- getPCARes(X$data, A)
      res <- ldecomp.getDistances(m$scores, m$loadings, m$residuals)

      # check dimensions
      expect_equal(dim(res$Q), c(32, A))
      expect_equal(dim(res$T2), c(32, A))

      # check values
      for (a in 1:A) {
         E <- m$residuals
         if (a < A) {
            E <- E + m$scores[, (a + 1):A, drop = F] %*% t(m$loadings[, (a + 1):A, drop = F])
         }
         U <- scale(m$scores[, 1:a], center = FALSE, scale = m$tnorm[1:a])

         expect_equal(res$Q[, a], rowSums(E^2))
         expect_equal(res$T2[, a], rowSums(U^2))
         expect_equivalent(res$tnorm[1:a], m$tnorm[1:a])
      }


      # check names
      expect_equal(rownames(res$Q), rownames(m$scores))
      expect_equal(colnames(res$Q), colnames(m$loadings))
      expect_equal(rownames(res$T2), rownames(m$scores))
      expect_equal(colnames(res$T2), colnames(m$loadings))
      expect_equal(names(res$tnorm), colnames(m$loadings))

      # check attributes
      expect_equal(attr(res$Q, "exclrows"), X$exp_exclrows)
      expect_equal(attr(res$Q, "yaxis.name"), attr(X$data, "yaxis.name"))
      expect_equal(attr(res$Q, "yaxis.values"), attr(X$data, "yaxis.values"))
      expect_equal(attr(res$T2, "exclrows"), X$exp_exclrows)
      expect_equal(attr(res$T2, "yaxis.name"), attr(X$data, "yaxis.name"))
      expect_equal(attr(res$T2, "yaxis.values"), attr(X$data, "yaxis.values"))
   }
}


test_that("tnorm, Q and T2 are correct for plain data", { tf(X1) })
test_that("tnorm, Q and T2 are correct for data attributes", { tf(X2) })
test_that("tnorm, Q and T2 are correct for data with excluded rows", { tf(X3) })
test_that("tnorm, Q and T2 are correct for data with excluded columns", { tf(X3) })


#######################################
# Block 2. Computing variance values  #
#######################################

context("ldecomp: computing variance")

## shortcut for testing function
tf <- function(X) {

   for (A in c(1, 5, ncol(X$data) - length(attr(X$data, "exclcols")))) {
      m <- getPCARes(X$data, A)
      dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals)
      res <- ldecomp.getVariances(m$scores, m$loadings, m$residuals, dist$Q)

      # check dimensions
      expect_equal(length(res$expvar), A)
      expect_equal(length(res$cumexpvar), A)

      # check values
      for (a in 1:A) {
         E <- m$residuals
         if (a < A) {
            E <- E + m$scores[, (a + 1):A, drop = F] %*% t(m$loadings[, (a + 1):A, drop = F])
         }

         # excluded rows should not be taken into account for calculation of variance
         rows_excluded <- attr(m$scores, "exclrows")
         if (length(rows_excluded) > 0) {
            E <- E[-rows_excluded, ]
         }

         expect_equivalent(res$cumexpvar[a], 100 * (1 - sum(E^2) / m$totvar))
      }

      # check names
      expect_equal(names(res$expvar), colnames(m$loadings))
      expect_equal(names(res$cumexpvar), colnames(m$loadings))
   }
}

test_that("explained variance is correct for plain data", { tf(X1) })
test_that("explained variance is correct for data attributes", { tf(X2) })
test_that("explained variance is correct for data with excluded rows", { tf(X3) })
test_that("explained variance is correct for data with excluded columns", { tf(X3) })


#########################################################
# Block 3. Computing critical limits (chisq/hotelling)  #
#########################################################

## chisq/hotelling for full data
context("ldecomp: computing limits (chisq/hotelling, full data)")

# in this case we use predefined values computed for full People data with autoscaling
# the values we got from PLS_Toolbox and DD-SIMCA Toolbox

m <- getPCARes(X1$data, 12)
dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals)

test_that("critical limits for Q are correct (chisq)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim1 <- rbind(
      c(9.25512512, 7.37431769, 3.64879767, 1.55968461, 0.82955551, 0.52538528, 0.23914330,
         0.15607468, 0.11587272, 0.04727562, 0.04666972, 0.00000000),
      c(16.87605393, 16.33648881, 8.08326204, 4.19946801, 2.23358734, 1.25780278, 0.52977943,
         0.37365182, 0.31198858, 0.15980027, 0.15775224, 0.00000000),
      apply(dist$Q, 2, mean),
      apply(dist$Q, 2, function(x) 2 * (mean(x) / sd(x))^2)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim2 <- rbind(
      c(8.03234581, 6.04655385, 2.99182277, 1.19880757, 0.63761443, 0.42028176, 0.19608497,
         0.12485188, 0.08906230, 0.03329627, 0.03286953, 0.00000000),
      c(14.57482653, 13.53946609, 6.69930078, 3.35119154, 1.78241124, 1.02644677, 0.43907419,
         0.30492357, 0.24896808, 0.12254836, 0.12097776, 0.00000000),
      apply(dist$Q, 2, mean),
      apply(dist$Q, 2, function(x)  2 * (mean(x) / sd(x))^2)
   )

   lim1 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.05, gamma = 0.01, lim.type = "chisq")
   lim2 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.10, gamma = 0.05, lim.type = "chisq")

   expect_equivalent(chisq.crit(dist$Q[, 3], 3, 0.05, 0.01), expQlim1[, 3], tolerance = 10^-5)
   expect_equivalent(chisq.crit(dist$Q[, 4], 4, 0.10, 0.05), expQlim2[, 4], tolerance = 10^-5)
   expect_equivalent(lim1$Qlim, expQlim1, tolerance = 10^-5)
   expect_equivalent(lim2$Qlim, expQlim2, tolerance = 10^-5)

})

test_that("critical limits for T2 are correct (hotelling)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "tsqlim" function from PLS_Toolbox
   expT2lim1 <- rbind(
      c(4.15961510, 6.85271430, 9.40913034, 12.01947857, 14.76453307, 17.69939358, 20.87303997,
         24.33584211, 28.14388532, 32.36253393, 37.07020868, 42.36299868),
      c(16.43763399, 22.07592005, 27.49293159, 33.11628022, 39.14366436, 45.72590957, 53.01070123,
         61.16173704, 70.37250826, 80.88041709, 92.98424615, 107.06768987),
      apply(dist$T2, 2, mean),
      nrow(dist$T2) - (1:ncol(dist$T2))
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "tsqlim" function from PLS_Toolbox
   expT2lim2 <- rbind(
      c(2.87478394, 5.14334644, 7.32156708, 9.55303105, 11.90044572, 14.40708180, 17.11153965,
         20.05346964, 23.27685391, 26.83267841, 30.78172763, 35.19795488),
      c(11.96070234, 16.61285114, 21.04123968, 25.60304588, 30.45666797, 35.71776166, 41.49572457,
         47.90883728, 55.09430867, 63.21788003, 72.48520544, 83.15679204),
      apply(dist$T2, 2, mean),
      nrow(dist$T2) - (1:ncol(dist$T2))
   )


   lim1 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.05, gamma = 0.01, lim.type = "chisq")
   lim2 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.10, gamma = 0.05, lim.type = "chisq")

   expect_equivalent(hotelling.crit(dist$T2[, 3], 3, 0.05, 0.01), expT2lim1[, 3], tolerance = 10^-5)
   expect_equivalent(hotelling.crit(dist$T2[, 4], 4, 0.10, 0.05), expT2lim2[, 4], tolerance = 10^-5)
   expect_equivalent(lim1$T2lim, expT2lim1, tolerance = 10^-5)
   expect_equivalent(lim2$T2lim, expT2lim2, tolerance = 10^-5)
})


## chisq/hotelling for excluded data data
context("ldecomp: computing limits (chisq/hotelling, excluded data)")

m <- getPCARes(X4$data, 10)
dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals)

rows_excluded <- X4$exp_exclrows
cols_excluded <- X4$exp_exclcols
test_that("critical limits for Q are correct when rows and columns are excluded (chisq)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim1 <- rbind(
      c(8.48240511, 4.44080558, 0.89616079, 0.56211531, 0.32748958, 0.18707463, 0.11620620,
         0.05214964, 0.02356848, 0.00000000),
      c(16.22421281, 11.81101671, 1.85705202, 1.49503355, 0.71805040, 0.49755423, 0.30906855,
         0.17377461, 0.07853561, 0.00000000),
      apply(dist$Q[-rows_excluded, , drop = F], 2, mean),
      apply(dist$Q[-rows_excluded, , drop = F], 2, function(x) 2 * (mean(x) / sd(x))^2)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim2 <- rbind(
      c(7.24620424, 3.41329992, 0.74768643, 0.43205407, 0.26852429, 0.14378964, 0.08931862, 0.03672904, 0.01659930, 0.00000000),
      c(13.81464099, 9.39586381, 1.55391044, 1.18932451, 0.59370604, 0.39581282, 0.24586927, 0.13272621, 0.05998421, 0.00000000),
      apply(dist$Q[-rows_excluded, , drop = F], 2, mean),
      apply(dist$Q[-rows_excluded, , drop = F], 2, function(x)  2 * (mean(x) / sd(x))^2)
   )

   lim1 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.05, gamma = 0.01, lim.type = "chisq")
   lim2 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.10, gamma = 0.05, lim.type = "chisq")

   expect_equivalent(chisq.crit(dist$Q[-rows_excluded, 3], 3, 0.05, 0.01), expQlim1[, 3], tolerance = 10^-5)
   expect_equivalent(chisq.crit(dist$Q[-rows_excluded, 4], 4, 0.10, 0.05), expQlim2[, 4], tolerance = 10^-5)
   expect_equivalent(lim1$Qlim, expQlim1, tolerance = 10^-5)
   expect_equivalent(lim2$Qlim, expQlim2, tolerance = 10^-5)

})

test_that("critical limits for T2 are correct when rows and columns are excluded (chisq/hotelling)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "residuallimit" function from PLS_Toolbox
   expT2lim1 <- rbind(
      c(4.19597182, 6.95671580, 9.61203588, 12.35902290, 15.28714920, 18.46287368, 21.94998503,
         25.81826344, 30.14945777, 35.04323329),
      c(16.57969094, 22.52147374, 28.33766821, 34.48651357, 41.20190951, 48.68181023, 57.13678256,
         66.81570665, 78.02909222, 91.17751461),
      apply(dist$T2[-rows_excluded, , drop = F], 2, mean),
      nrow(dist$T2) - (1:ncol(dist$T2)) - length(rows_excluded)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "residuallimit" function from PLS_Toolbox
   expT2lim2 <- rbind(
      c(2.89384646, 5.20718834, 7.45497071, 9.78540224, 12.26769498, 14.95365699, 17.89298964,
         21.13982013, 24.75715160, 28.82121257),
      c(11.94402185, 16.77785866, 21.45842445, 26.36104904, 31.66743997, 37.52399529, 44.08105163,
         51.51194381, 60.02878045, 69.90058746),
      apply(dist$T2[-rows_excluded, , drop = F], 2, mean),
      nrow(dist$T2) - (1:ncol(dist$T2)) - length(rows_excluded)
   )

   lim1 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.05, gamma = 0.01, lim.type = "chisq")
   lim2 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.10, gamma = 0.05, lim.type = "chisq")

   expect_equivalent(hotelling.crit(dist$T2[-rows_excluded, 3], 3, 0.05, 0.01), expT2lim1[, 3], tolerance = 10^-5)
   expect_equivalent(hotelling.crit(dist$T2[-rows_excluded, 4], 4, 0.10, 0.05), expT2lim2[, 4], tolerance = 10^-5)
   expect_equivalent(lim1$T2lim, expT2lim1, tolerance = 10^-5)
   expect_equivalent(lim2$T2lim, expT2lim2, tolerance = 10^-5)

})

###################################################
# Block 4. Computing critical limits (ddmoments)  #
###################################################


## ddmoments for full data data
context("ldecomp: computing limits (ddmoments, full data)")

### in this case we use 11 components as the 12th explain mostly noise for Q
m <- getPCARes(X1$data, 11)
dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals)

test_that("critical limits for Q are correct (ddmoments)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expQlim1 <- rbind(
      c(11.34617157, 19.08855942, 11.42533747, 6.34336609, 3.34064053, 1.84178155, 0.55799538,
         0.39965752, 0.47200886, 0.25476038, 0.17850637),
      c(19.47199649, 31.77017452, 18.32858290, 10.17605927, 5.41986135, 3.02485080, 0.94234321,
         0.67494207, 0.81004899, 0.42401219, 0.28960906),
      c(5.39623603, 3.22376502, 1.65661907, 0.68981822, 0.38111636, 0.22105099, 0.12476396,
         0.07148856, 0.04489746, 0.02151256, 0.00678829),
      c(10.00000000, 4.00000000, 4.00000000, 3.00000000, 3.00000000, 3.00000000, 5.00000000,
         4.00000000, 2.00000000, 2.00000000, 1.00000000)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expQlim2 <- rbind(
      c(10.00967132, 16.97646011, 10.25820924, 5.69537459, 2.99072482, 1.64367017, 0.49436318,
         0.35408172, 0.41640949, 0.22657181, 0.15980870),
      c(17.04539107, 28.00480337, 16.29329975, 9.04606673, 4.80551079, 2.67446799, 0.82790983,
         0.59298052, 0.70910047, 0.37375866, 0.25678138),
      c(5.39623603, 3.22376502, 1.65661907, 0.68981822, 0.38111636, 0.22105099, 0.12476396,
         0.07148856, 0.04489746, 0.02151256, 0.00678829),
      c(10.00000000, 4.00000000, 4.00000000, 3.00000000, 3.00000000, 3.00000000, 5.00000000,
         4.00000000, 2.00000000, 2.00000000, 1.00000000)
   )

   lim1 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.05, gamma = 0.01, lim.type = "ddmoments")
   lim2 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.10, gamma = 0.05, lim.type = "ddmoments")

   expect_equivalent(ddmoments.crit(dist$Q[, 3], 3, 0.05, 0.01)[3:4], expQlim1[3:4, 3], tolerance = 10^-4)
   expect_equivalent(ddmoments.crit(dist$Q[, 4], 4, 0.10, 0.05)[3:4], expQlim2[3:4, 4], tolerance = 10^-4)
   expect_equivalent(lim1$Qlim, expQlim1, tolerance = 10^-4)
   expect_equivalent(lim2$Qlim, expQlim2, tolerance = 10^-4)
})


test_that("critical limits for T2 are correct (ddmoments)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expT2lim1 <- rbind(
      c(10.18450977, 4.58893048, 6.16731299, 7.63572085, 9.79787223, 12.10733682, 18.95532441,
         19.25620258, 18.33211759, 19.12054367, 18.68127639),
      c(17.47838355, 7.63761785, 9.89363402, 12.24926117, 15.89608595, 19.88449041, 32.01177265,
         32.51989603, 31.46109039, 31.82340773, 30.30853722),
      c(0.96875000, 1.93750000, 2.90625000, 3.87500000, 4.84375000, 5.81250000, 6.78125000,
         7.75000000, 8.71875000, 9.68750000, 10.65625000),
      c(2.00000000, 10.00000000, 13.00000000, 14.00000000, 13.00000000, 12.00000000, 8.00000000,
         9.00000000, 10.00000000, 12.00000000, 15.00000000)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expT2lim2 <- rbind(
      c(8.98484522, 4.08117729, 5.53730578, 6.85571191, 8.77159317, 10.80501016, 16.79371331,
         17.06028019, 16.17272140, 17.00490539, 16.72450430),
      c(15.30022271, 6.73241458, 8.79500317, 10.88905154, 14.09423741, 17.58117559, 28.12442515,
         28.57084459, 27.54040088, 28.05172743, 26.87301266),
      c(0.96875000, 1.93750000, 2.90625000, 3.87500000, 4.84375000, 5.81250000, 6.78125000,
         7.75000000, 8.71875000, 9.68750000, 10.65625000),
      c(2.00000000, 10.00000000, 13.00000000, 14.00000000, 13.00000000, 12.00000000, 8.00000000,
         9.00000000, 10.00000000, 12.00000000, 15.00000000)
   )

   lim1 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.05, gamma = 0.01, lim.type = "ddmoments")
   lim2 <- ldecomp.getLimits(dist$Q, dist$T2, alpha = 0.10, gamma = 0.05, lim.type = "ddmoments")

   expect_equivalent(ddmoments.crit(dist$T2[, 3], 3, 0.05, 0.01)[3:4], expT2lim1[3:4, 3], tolerance = 10^-4)
   expect_equivalent(ddmoments.crit(dist$T2[, 4], 4, 0.10, 0.05)[3:4], expT2lim2[3:4, 4], tolerance = 10^-4)
   expect_equivalent(lim1$T2lim, expT2lim1, tolerance = 10^-4)
   expect_equivalent(lim2$T2lim, expT2lim2, tolerance = 10^-4)
})