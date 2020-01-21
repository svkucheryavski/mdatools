#############################################################
# Tests for basic functionality of critical limits methods  #
#############################################################

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
   eigenvals <- colSums((X_cal %*% loadings_visible)^2) / (nrow(X_cal) - 1)

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
      eigenvals = eigenvals, ncomp = ncomp, totvar = sum(X_cal^2)))
}


#########################################################
#region 1. Computing critical limits (chisq/hotelling)  #
#########################################################

## chisq/hotelling for full data
context("misc: computing limits (chisq/hotelling, full data)")

# in this case we use predefined values computed for full People data with autoscaling
# the values we got from PLS_Toolbox and DD-SIMCA Toolbox


test_that("critical limits for Q are correct (jm) full decomposition", {

   # test for both full and partial decomposition
   for (a in c(4, 12)) {
      m <- getPCARes(X1$data, a)

      # expected limits for alpha = 0.05 and gamma = 0.01
      # computed using "residuallimit" function from PLS_Toolbox
      expQlim1 <- rbind(
         c(13.98208375, 8.91523789, 4.86682122, 1.81125667, 0.98188720, 0.57560902, 0.32312846,
            0.18598893, 0.12950361, 0.06958773, 0.02625457, 0),
         c(37.92391547, 27.28432368, 18.46348557, 5.51382242, 2.82348737, 1.68333989, 0.89628978,
            0.47263837, 0.37503680, 0.22537117, 0.09561530, 0)
      )

      # expected limits for alpha = 0.10 and gamma = 0.05
      # computed using "residuallimit" function from PLS_Toolbox
      expQlim2 <- rbind(
         c(11.08714606, 6.88643872, 3.62417784, 1.41387882, 0.77404080, 0.45201939, 0.25555116,
            0.14863431, 0.09998945, 0.05180038, 0.01849205, 0),
         c(29.65164966, 20.68354357, 13.16120605, 4.16057052, 2.16549979, 1.28564550, 0.69582345,
            0.37714904, 0.29044220, 0.17075224, 0.07133111, 0)
      )

      expect_equivalent(jm.crit(m$residuals, m$eigenvals, 0.05, 0.01), expQlim1[, seq_len(a)], tolerance = 10^-5)
      expect_equivalent(jm.crit(m$residuals, m$eigenvals, 0.10, 0.05), expQlim2[, seq_len(a)], tolerance = 10^-5)
   }

})

m <- getPCARes(X1$data, 12)
dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals, m$eigenvals)

test_that("critical limits for Q are correct (chisq)", {

   # expected distribution parameters
   expParams <- list(
      u0 = apply(dist$Q, 2, mean),
      Nu = apply(dist$Q, 2, function(x) 2 * (mean(x) / sd(x))^2),
      nobj = nrow(dist$Q)
   )

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim1 <- rbind(
      c(9.25512512, 7.37431769, 3.64879767, 1.55968461, 0.82955551, 0.52538528, 0.23914330,
         0.15607468, 0.11587272, 0.04727562, 0.04666972, 0.00000000),
      c(16.87605393, 16.33648881, 8.08326204, 4.19946801, 2.23358734, 1.25780278, 0.52977943,
         0.37365182, 0.31198858, 0.15980027, 0.15775224, 0.00000000)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim2 <- rbind(
      c(8.03234581, 6.04655385, 2.99182277, 1.19880757, 0.63761443, 0.42028176, 0.19608497,
         0.12485188, 0.08906230, 0.03329627, 0.03286953, 0.00000000),
      c(14.57482653, 13.53946609, 6.69930078, 3.35119154, 1.78241124, 1.02644677, 0.43907419,
         0.30492357, 0.24896808, 0.12254836, 0.12097776, 0.00000000)
   )

   # parameters
   p <- ddmoments.param(dist$Q)
   expect_equivalent(p, expParams, tolerance = 10^-5)

   # limits
   expect_equivalent(chisq.crit(p, 0.05, 0.01), expQlim1, tolerance = 10^-5)
   expect_equivalent(chisq.crit(p, 0.10, 0.05), expQlim2, tolerance = 10^-5)
})

test_that("critical limits for T2 are correct (hotelling)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "tsqlim" function from PLS_Toolbox
   expT2lim1 <- rbind(
      c(4.15961510, 6.85271430, 9.40913034, 12.01947857, 14.76453307, 17.69939358, 20.87303997,
         24.33584211, 28.14388532, 32.36253393, 37.07020868, 42.36299868),
      c(16.43763399, 22.07592005, 27.49293159, 33.11628022, 39.14366436, 45.72590957, 53.01070123,
         61.16173704, 70.37250826, 80.88041709, 92.98424615, 107.06768987)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "tsqlim" function from PLS_Toolbox
   expT2lim2 <- rbind(
      c(2.87478394, 5.14334644, 7.32156708, 9.55303105, 11.90044572, 14.40708180, 17.11153965,
         20.05346964, 23.27685391, 26.83267841, 30.78172763, 35.19795488),
      c(11.96070234, 16.61285114, 21.04123968, 25.60304588, 30.45666797, 35.71776166, 41.49572457,
         47.90883728, 55.09430867, 63.21788003, 72.48520544, 83.15679204)
   )

   nobj <- nrow(dist$T2)
   ncomp <- seq_len(ncol(dist$T2))
   expect_equivalent(hotelling.crit(nobj, ncomp, 0.05, 0.01), expT2lim1, tolerance = 10^-5)
   expect_equivalent(hotelling.crit(nobj, ncomp, 0.10, 0.05), expT2lim2, tolerance = 10^-5)
})

## chisq/hotelling for excluded data data
context("misc: computing limits (chisq/hotelling, excluded data)")

m <- getPCARes(X4$data, 10)
dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals, m$eigenvals)
rows_excluded <- X4$exp_exclrows

test_that("critical limits for Q are correct for excluded data (chisq)", {

   # expected distribution parameters
   expParams <- list(
      u0 = apply(dist$Q[-rows_excluded, , drop = F], 2, mean),
      Nu = apply(dist$Q[-rows_excluded, , drop = F], 2, function(x) 2 * (mean(x) / sd(x))^2),
      nobj = nrow(dist$Q) - length(rows_excluded)
   )

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim1 <- rbind(
      c(8.48240511, 4.44080558, 0.89616079, 0.56211531, 0.32748958, 0.18707463, 0.11620620,
         0.05214964, 0.02356848, 0.00000000),
      c(16.22421281, 11.81101671, 1.85705202, 1.49503355, 0.71805040, 0.49755423, 0.30906855,
         0.17377461, 0.07853561, 0.00000000)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "residuallimit" function from PLS_Toolbox
   expQlim2 <- rbind(
      c(7.24620424, 3.41329992, 0.74768643, 0.43205407, 0.26852429, 0.14378964, 0.08931862,
         0.03672904, 0.01659930, 0.00000000),
      c(13.81464099, 9.39586381, 1.55391044, 1.18932451, 0.59370604, 0.39581282, 0.24586927,
         0.13272621, 0.05998421, 0.00000000)
   )

   # parameters
   p <- ddmoments.param(dist$Q[-rows_excluded, ])
   expect_equivalent(p, expParams, tolerance = 10^-5)

   # limits
   expect_equivalent(chisq.crit(p, 0.05, 0.01), expQlim1, tolerance = 10^-5)
   expect_equivalent(chisq.crit(p, 0.10, 0.05), expQlim2, tolerance = 10^-5)
})

test_that("critical limits for T2 are correct for excluded data (chisq/hotelling)", {

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using "residuallimit" function from PLS_Toolbox
   expT2lim1 <- rbind(
      c(4.19597182, 6.95671580, 9.61203588, 12.35902290, 15.28714920, 18.46287368, 21.94998503,
         25.81826344, 30.14945777, 35.04323329),
      c(16.57969094, 22.52147374, 28.33766821, 34.48651357, 41.20190951, 48.68181023, 57.13678256,
         66.81570665, 78.02909222, 91.17751461)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using "residuallimit" function from PLS_Toolbox
   expT2lim2 <- rbind(
      c(2.89384646, 5.20718834, 7.45497071, 9.78540224, 12.26769498, 14.95365699, 17.89298964,
         21.13982013, 24.75715160, 28.82121257),
      c(11.94402185, 16.77785866, 21.45842445, 26.36104904, 31.66743997, 37.52399529, 44.08105163,
         51.51194381, 60.02878045, 69.90058746)
   )

   # limits
   nobj <- nrow(dist$T2) - length(rows_excluded)
   ncomp <- 1:ncol(dist$T2)
   expect_equivalent(hotelling.crit(nobj, ncomp, 0.05, 0.01), expT2lim1, tolerance = 10^-5)
   expect_equivalent(hotelling.crit(nobj, ncomp, 0.10, 0.05), expT2lim2, tolerance = 10^-5)
})

###################################################
# Block 4. Computing critical limits (ddmoments)  #
###################################################

## ddmoments for full data data
context("misc: computing limits (ddmoments, full data)")

### in this case we use 11 components as the 12th explain mostly noise for Q
m <- getPCARes(X1$data, 11)
dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals, m$eigenvals)
pQ <- ddmoments.param(dist$Q)
pT2 <- ddmoments.param(dist$T2)

test_that("critical limits for Q are correct (ddmoments)", {

   # expected distribution parameters
   expParams <- list(
      u0 = c(5.39623603, 3.22376502, 1.65661907, 0.68981822, 0.38111636, 0.22105099, 0.12476396,
         0.07148856, 0.04489746, 0.02151256, 0.00678829),
      Nu = c(10.00000000, 4.00000000, 4.00000000, 3.00000000, 3.00000000, 3.00000000, 5.00000000,
         4.00000000, 2.00000000, 2.00000000, 1.00000000),
      nobj = nrow(dist$Q)
   )

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

   # check parameters (round DoF first)
   pQ$Nu <- round(pQ$Nu)
   expect_equivalent(pQ, expParams, tolerance = 10^-5)

   # limits
   lim1 <- scale(dd.crit(pQ, pT2, 0.05, 0.01), center = FALSE, scale = pQ$Nu / pQ$u0)
   lim2 <- scale(dd.crit(pQ, pT2, 0.10, 0.05), center = FALSE, scale = pQ$Nu / pQ$u0)

   expect_equivalent(lim1, expQlim1, tolerance = 10^-5)
   expect_equivalent(lim2, expQlim2, tolerance = 10^-5)
})

test_that("critical limits for T2 are correct (ddmoments)", {

   # expected distribution parameters
   expParams <- list(
      u0 = c(0.96875000, 1.93750000, 2.90625000, 3.87500000, 4.84375000, 5.81250000, 6.78125000,
         7.75000000, 8.71875000, 9.68750000, 10.65625000),
      Nu = c(2.00000000, 10.00000000, 13.00000000, 14.00000000, 13.00000000, 12.00000000, 8.00000000,
         9.00000000, 10.00000000, 12.00000000, 15.00000000),
      nobj = nrow(dist$T2)
   )

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

   # check parameters (round DoF first)
   pT2$Nu <- round(pT2$Nu)
   expect_equivalent(pT2, expParams, tolerance = 10^-5)

   # limits
   lim1 <- scale(dd.crit(pQ, pT2, 0.05, 0.01), center = FALSE, scale = pT2$Nu / pT2$u0)
   lim2 <- scale(dd.crit(pQ, pT2, 0.10, 0.05), center = FALSE, scale = pT2$Nu / pT2$u0)

   expect_equivalent(lim1, expT2lim1, tolerance = 10^-5)
   expect_equivalent(lim2, expT2lim2, tolerance = 10^-5)
})

## ddmoments for excluded data
context("misc: computing limits (ddmoments, excluded data)")

m <- getPCARes(X4$data, 9)
dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals, m$eigenvals)
rows_excluded <- X4$exp_exclrows
pQ <- ddmoments.param(dist$Q[-rows_excluded, ])
pT2 <- ddmoments.param(dist$T2[-rows_excluded, ])

test_that("critical limits for Q are correct for excluded data (ddmoments)", {

   # expected distribution parameters
   expParams <- list(
      u0 = c(4.43617757, 2.19049660, 0.44027710, 0.25151771, 0.14625433, 0.08237345, 0.05117557,
         0.02675517, 0.00723308),
      Nu = c(7.00000000, 3.00000000, 5.00000000, 3.00000000, 4.00000000, 3.00000000, 3.00000000,
         2.00000000, 1.00000000),
      nobj = nrow(dist$Q) - length(rows_excluded)
   )

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expQlim1 <- rbind(
      c(11.60190138, 19.20059700, 2.76585832, 2.42038062, 1.05556531, 0.75748212, 0.59999005,
         0.35177998, 0.20881386),
      c(20.51821091, 30.94285188, 4.28239312, 3.81776382, 1.66498567, 1.20720511, 0.90716035,
         0.56691340, 0.32937052)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expQlim2 <- rbind(
      c(10.13171606, 17.18942863, 2.50182801, 2.17893453, 0.95026695, 0.68010334, 0.54599027,
         0.31493275, 0.18798355),
      c(17.78135669, 27.40039229, 3.82840757, 3.39793186, 1.48189048, 1.07182017, 0.81563549,
         0.50201092, 0.29315029)
   )

   # check parameters (round DoF first)
   pQ$Nu <- round(pQ$Nu)
   expect_equivalent(pQ, expParams, tolerance = 10^-5)

   # limits
   lim1 <- scale(dd.crit(pQ, pT2, 0.05, 0.01), center = FALSE, scale = pQ$Nu / pQ$u0)
   lim2 <- scale(dd.crit(pQ, pT2, 0.10, 0.05), center = FALSE, scale = pQ$Nu / pQ$u0)

   expect_equivalent(lim1, expQlim1, tolerance = 10^-5)
   expect_equivalent(lim2, expQlim2, tolerance = 10^-5)
})

test_that("critical limits for T2 are correct for excluded data (ddmoments)", {

   # expected distribution parameters
   expParams <- list(
      u0 = c(0.96551724, 1.93103448, 2.89655172, 3.86206897, 4.82758621, 5.79310345, 6.75862069,
         7.72413793, 8.68965517),
      Nu = c(3.00000000, 13.00000000, 15.00000000, 15.00000000, 14.00000000, 14.00000000, 20.00000000,
         14.00000000, 17.00000000),
      nobj = nrow(dist$T2) - length(rows_excluded)
   )

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expT2lim1 <- rbind(
      c(5.89192305, 3.90607053, 6.06546174, 7.43301692, 9.95493338, 11.41536020, 11.88586181,
         14.50826197, 14.75672477),
      c(10.41999203, 6.29485437, 9.39118662, 11.72439691, 15.70231730, 18.19274770, 17.97093555,
         23.38088765, 23.27637623)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expT2lim2 <- rbind(
      c(5.14530244, 3.49692880, 5.48644953, 6.69153318, 8.96187480, 10.24925118, 10.81612080,
         12.98859268, 13.28466146),
      c(9.03010480, 5.57419464, 8.39560708, 10.43508810, 13.97556442, 16.15247798, 16.15781902,
         20.70415152, 20.71671903)
   )

   # check parameters (round DoF first)
   pT2$Nu <- round(pT2$Nu)
   expect_equivalent(pT2, expParams, tolerance = 10^-5)

   # limits
   lim1 <- scale(dd.crit(pQ, pT2, 0.05, 0.01), center = FALSE, scale = pT2$Nu / pT2$u0)
   lim2 <- scale(dd.crit(pQ, pT2, 0.10, 0.05), center = FALSE, scale = pT2$Nu / pT2$u0)

   expect_equivalent(lim1, expT2lim1, tolerance = 10^-5)
   expect_equivalent(lim2, expT2lim2, tolerance = 10^-5)
})


###################################################
# Block 5. Computing critical limits (ddrobust)   #
###################################################

## ddmoments for full data data
context("misc: computing limits (ddrobust, full data)")

### in this case we use 11 components as the 12th explain mostly noise for Q
m <- getPCARes(X1$data, 11)
dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals, m$eigenvals)
pQ <- ddrobust.param(dist$Q)
pT2 <- ddrobust.param(dist$T2)

test_that("critical limits for Q are correct (ddrobust)", {

   # expected distribution parameters
   expParams <- list(
      u0 = c(5.59398247, 3.14412554, 1.68373781, 0.63464055, 0.33473045, 0.21529782, 0.13614951,
         0.06862920, 0.04666029, 0.02054372, 0.00374660),
      Nu = c(18.00000000, 10.00000000, 6.00000000, 3.00000000, 4.00000000, 3.00000000, 5.00000000,
         3.00000000, 4.00000000, 2.00000000, 1.00000000),
      nobj = nrow(dist$Q)
   )

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

   # check parameters (round DoF first)
   pQ$Nu <- round(pQ$Nu)
   expect_equivalent(pQ, expParams, tolerance = 10^-5)

   # limits
   lim1 <- scale(dd.crit(pQ, pT2, 0.05, 0.01), center = FALSE, scale = pQ$Nu / pQ$u0)
   lim2 <- scale(dd.crit(pQ, pT2, 0.10, 0.05), center = FALSE, scale = pQ$Nu / pQ$u0)

   expect_equivalent(lim1, expQlim1, tolerance = 10^-5)
   expect_equivalent(lim2, expQlim2, tolerance = 10^-5)
})

test_that("critical limits for T2 are correct (ddrobust)", {

   # expected distribution parameters
   expParams <- list(
      u0 = c(1.23877706, 1.84445603, 2.67812507, 3.91673370, 4.51759976, 5.66332668, 5.96116787,
         7.19886896, 8.06726598, 8.95166844, 10.01296469),
      Nu = c(3.00000000, 22.00000000, 17.00000000, 11.00000000, 15.00000000, 18.00000000, 16.00000000,
         18.00000000, 18.00000000, 10.00000000, 14.00000000),
      nobj = nrow(dist$T2)
   )

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

   # check parameters (round DoF first)
   pT2$Nu <- round(pT2$Nu)
   expect_equivalent(pT2, expParams, tolerance = 10^-5)

   # limits
   lim1 <- scale(dd.crit(pQ, pT2, 0.05, 0.01), center = FALSE, scale = pT2$Nu / pT2$u0)
   lim2 <- scale(dd.crit(pQ, pT2, 0.10, 0.05), center = FALSE, scale = pT2$Nu / pT2$u0)

   expect_equivalent(lim1, expT2lim1, tolerance = 10^-5)
   expect_equivalent(lim2, expT2lim2, tolerance = 10^-5)
})


## ddmoments for excluded data
context("misc: computing limits (ddrobust, excluded data)")

m <- getPCARes(X4$data, 9)
dist <- ldecomp.getDistances(m$scores, m$loadings, m$residuals, m$eigenvals)
rows_excluded <- X4$exp_exclrows
pQ <- ddrobust.param(dist$Q[-rows_excluded, ])
pT2 <- ddrobust.param(dist$T2[-rows_excluded, ])

test_that("critical limits for Q are correct for excluded data (ddrobust)", {

   # expected distribution parameters
   expParams <- list(
      u0 = c(4.62537076, 2.40195599, 0.43916799, 0.23975306, 0.13645627, 0.07560735, 0.04351790,
         0.02600730, 0.00787551),
      Nu = c(16.00000000, 3.00000000, 4.00000000, 3.00000000, 5.00000000, 4.00000000, 2.00000000,
         2.00000000, 1.00000000),
      nobj = nrow(dist$Q) - length(rows_excluded)
   )

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expQlim1 <- rbind(
      c(8.71406195, 27.16167173, 4.67241416, 1.99760592, 0.64638798, 0.49704699, 1.18743408,
         0.45737021, 0.20709631),
      c(13.61362649, 41.37423787, 6.80661271, 3.25834093, 1.06828931, 0.80101944, 1.65599362,
         0.69152500, 0.33374746)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expQlim2 <- rbind(
      c(7.86416256, 24.67071673, 4.29149155, 1.78273328, 0.57486684, 0.44498376, 1.10230325,
         0.41620638, 0.18540399),
      c(12.14433627, 37.13321925, 6.17541759, 2.87712373, 0.94039845, 0.70931558, 1.51864614,
         0.62175594, 0.29553873)
   )

   # check parameters (round DoF first)
   pQ$Nu <- round(pQ$Nu)
   expect_equivalent(pQ, expParams, tolerance = 10^-5)

   # limits
   lim1 <- scale(dd.crit(pQ, pT2, 0.05, 0.01), center = FALSE, scale = pQ$Nu / pQ$u0)
   lim2 <- scale(dd.crit(pQ, pT2, 0.10, 0.05), center = FALSE, scale = pQ$Nu / pQ$u0)

   expect_equivalent(lim1, expQlim1, tolerance = 10^-5)
   expect_equivalent(lim2, expQlim2, tolerance = 10^-5)
})

test_that("critical limits for T2 are correct for excluded data (ddrobust)", {

   # expected distribution parameters
   expParams <- list(
      u0 = c(1.39125723, 2.06314350, 2.61209799, 3.71494609, 5.03000803, 5.67830732, 6.11739113,
         7.31069105, 8.12330118),
      Nu = c(3.00000000, 19.00000000, 25.00000000, 12.00000000, 9.00000000, 12.00000000, 37.0000000,
         21.00000000, 15.00000000),
      nobj = nrow(dist$Q) - length(rows_excluded)
   )

   # expected limits for alpha = 0.05 and gamma = 0.01
   # computed using DDSimca toolbox
   expT2lim1 <- rbind(
      c(13.97913356, 3.68373634, 4.44651850, 7.73816836, 13.23719406, 12.44317048, 9.02269130,
         12.24451916, 14.24081027),
      c(21.83903490, 5.61128141, 6.47753566, 12.62190422, 21.87719039, 20.05287550, 12.58303051,
         18.51321075, 22.94987405)
   )

   # expected limits for alpha = 0.10 and gamma = 0.05
   # computed using DDSimca toolbox
   expT2lim2 <- rbind(
      c(12.61572152, 3.34590657, 4.08401223, 6.90581165, 11.77253307, 11.13980940, 8.37582660,
         11.14249858, 12.74915524),
      c(19.48199356, 5.03610347, 5.87685672, 11.14517509, 19.25815009, 17.75714331, 11.53939875,
         16.64538321, 20.32248206)
   )

   # check parameters (round DoF first)
   pT2$Nu <- round(pT2$Nu)
   expect_equivalent(pT2, expParams, tolerance = 10^-5)

   # limits
   lim1 <- scale(dd.crit(pQ, pT2, 0.05, 0.01), center = FALSE, scale = pT2$Nu / pT2$u0)
   lim2 <- scale(dd.crit(pQ, pT2, 0.10, 0.05), center = FALSE, scale = pT2$Nu / pT2$u0)

   expect_equivalent(lim1, expT2lim1, tolerance = 10^-5)
   expect_equivalent(lim2, expT2lim2, tolerance = 10^-5)
})
