#####################################################
# Tests for plotting methods of pcares() class      #
#####################################################

setup({
   pdf(file = tempfile("mdatools-test-pcaresplots-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-pcaresplots-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

data(people)
x <- people
m <- pca(x, scale = TRUE, ncomp = 10)

context("pcares: residual distance plots")

test_that("residuals plot works fine in general.", {
   expect_silent({
      par(mfrow = c(2, 2))
      plotResiduals(m$calres)
      plotResiduals(m$calres, 1, show.labels = T)
      plotResiduals(m$calres, 2, show.labels = T, colmap = c("red", "green"))
      plotResiduals(m$calres, 3, cgroup = x[, 1], show.labels = T, colmap = c("red", "green"), pch = 17)
   })
})

test_that("residuals plot work find with different norm and log settings.", {
   expect_silent({
      par(mfrow = c(2, 2))
      plotResiduals(m$calres, 1, norm = F, log = F, show.labels = T)
      plotResiduals(m$calres, 1, norm = F, log = T, show.labels = T)
      plotResiduals(m$calres, 1, norm = T, log = F, show.labels = T)
      plotResiduals(m$calres, 1, norm = T, log = T, show.labels = T)
   })
})

test_that("residuals plot works well with categorize coloring.", {
   m1 <- m
   m1 <- setDistanceLimits(m1, alpha = 0.10, gamma = 0.05)
   cgroup1 <- categorize(m1, m1$calres, 1)

   m2 <- m
   m2 <- setDistanceLimits(m2, lim.type="ddrobust", alpha = 0.40, gamma = 0.30)
   cgroup2 <- categorize(m2, m2$calres, 3)

   expect_silent({
      par(mfrow = c(2, 2))
      plotResiduals(m$calres, 1, norm = F, cgroup = cgroup1, show.labels = T)
      plotResiduals(m$calres, 1, cgroup = cgroup1, show.labels = T)
      plotResiduals(m$calres, 1, log = T, cgroup = cgroup1, show.labels = T)
      plotResiduals(m$calres, 3, cgroup = cgroup2, show.labels = T)
   })
})

test_that("residuals plot returns correct plot data.", {
   pd <- plotResiduals(m$calres, ncomp = 4, show.plot = FALSE)
   expect_true("matrix" %in% class(pd))
   expect_equivalent(pd[, 1], m$calres$T2[, 4])
   expect_equivalent(pd[, 2], m$calres$Q[, 4])

   pd <- plotResiduals(m$calres, ncomp = 4, norm = TRUE, show.plot = FALSE)
   expect_true("matrix" %in% class(pd))
   expect_equivalent(pd[, 1], m$calres$T2[, 4] / m$T2lim[3, 4])
   expect_equivalent(pd[, 2], m$calres$Q[, 4] / m$Qlim[3, 4])

   pd <- plotResiduals(m$calres, ncomp = 4, norm = TRUE, log = TRUE, show.plot = FALSE, col = "#00AAFFFF",
      labels = "values")
   expect_true("matrix" %in% class(pd))
   expect_equivalent(pd[, 1], log(1 + m$calres$T2[, 4]/m$T2lim[3, 4]))
   expect_equivalent(pd[, 2], log(1 + m$calres$Q[, 4]/m$Qlim[3, 4]))
})


context('pcares: scores plots')

test_that("scores plot can return plot data.", {
   pd <- plotScores(m$calres, c(1, 3), show.plot = FALSE)
   expect_true("matrix" %in% class(pd))
   expect_equivalent(pd[, 1], m$calres$scores[, 1])
   expect_equivalent(pd[, 2], m$calres$scores[, 3])
})

test_that("scores plot works fine with different attributes.", {
   expect_silent({
      par(mfrow = c(2, 2))
      plotScores(m$calres)
      plotScores(m$calres, c(1, 3), show.labels = T)
      plotScores(m$calres, c(1, 2), cgroup = x[, 1], show.labels = T, show.colorbar = F, show.axes = F)
      plotScores(m$calres, 1, show.labels = T, colmap = c("red", "green"), pch = 17)
   })

   expect_silent({
      par(mfrow = c(2, 2))
      plotScores(m$calres, type = "l")
      plotScores(m$calres, c(1, 3), type = "h", show.labels = T)
      plotScores(m$calres, c(1, 2), type = "h", cgroup = x[1, ], show.labels = T, show.axes = F)
      plotScores(m$calres, 1, type = "l", show.labels = T, colmap = c("red", "green"), pch = 17)
   })

})


context('pcares: variance plots')

test_that("varance plot can return plot data.", {
   pd <- plotVariance(m$calres, show.plot = FALSE)
   expect_equal(class(pd), "numeric")
   expect_equivalent(pd, m$calres$expvar)
})

test_that("variance plot works fine.", {
   expect_silent({
      par(mfrow = c(2, 2))
      plotVariance(m$calres)
      plotVariance(m$calres, type = "h")
      plotVariance(m$calres, type = "h", col = "red")
      plotVariance(m$calres, type = "h", col = "red", show.labels = TRUE)
   })
})

test_that("cumulative varance plot can return plot data.", {
   pd <- plotCumVariance(m$calres, show.plot = FALSE)
   expect_equal(class(pd), "numeric")
   expect_equivalent(pd, m$calres$cumexpvar)
})

test_that("variance plot works fine.", {
   expect_silent({
      par(mfrow = c(2, 2))
      plotCumVariance(m$calres)
      plotCumVariance(m$calres, type = "h")
      plotCumVariance(m$calres, type = "h", col = "red", show.labels = TRUE)
      plotCumVariance(m$calres, type = "h", col = "red", show.labels = TRUE, labels = "names")
   })
})

context('pcares: overall plot and summary')

test_that("print() and summary() produce output", {
   expect_output(print(m$calres))
   expect_output(summary(m$calres))
})

test_that("print() and summary() produce output", {
   expect_silent(plot(m$calres))
   expect_silent(plot(m$calres, comp = c(1, 3)))
   expect_silent(plot(m$calres, comp = c(1, 3), ncomp = 4, show.labels = T))
})

# just output to check in txt file
print(m$calres)
print(summary(m$calres))
