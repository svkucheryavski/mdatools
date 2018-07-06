data(iris)
data = iris[1:50, 1:4]

# limits for T2-residuals by Hotteling method

context('PCA: residual limits for standardized data (Iris, Setosa)')

# limits for T2-residuals by Hotelling's T2 distribution
T2lim.hotelling = matrix(
   c(
      4.4036, 7.6057, 1.9026, 1.0000,
      1.9659, 3.3955, 0.9009, 1.0000,
      0.5671, 0.9794, 0.2464, 1.0000,
      0.0000, 0.0000, 0.0000, 1.0000
   ),
   ncol = 4
)

# limits for Q-residuals by chisq method
Qlim.chisq = matrix(
   c(
      4.4036, 15.8443, 1.9026, 1.6598,
      1.9659,  7.0734, 0.9009, 1.7604,
      0.5671,  2.0403, 0.2464, 1.6694,
      0.0000,  0.0000, 0.0000, 1.0000
      ),
   ncol = 4
)

m = pca(data, ncomp = 4, scale = T, lim.type = 'chisq')
m$Qlim
test_that("Chisq limits are correct (standardized)", {
   expect_equal(sum(round(m$Qlim, 4) - Qlim.chisq), 0)
})
