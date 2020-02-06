######################################
# Tests for regcoeffs class methods  #
######################################

setup({
   pdf(file = tempfile("mdatools-test-regcoeffs-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-regcoeffs-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# create several datasets

## number of predictors, components and y-variables
npred <- 5
ncomp <- 3
nresp <- 2
nseg <- 10

## basic coefficients for npred x nresp
b2 <- cbind(c(8, 6, 10, 4, 2), c(1, 9, 5, 7, 3))
b1 <- b2 + cbind(runif(npred, -0.5, 0.5), runif(npred, -0.5, 0.5))
b3 <- b2 + cbind(runif(npred, -0.75, 0.75), runif(npred, -0.75, 0.75))

## create array for all components
coeffs <- array(rbind(b1, b2, b3), dim = c(npred, ncomp, nresp))

## create ci.coeffs for each of 10 segments by adding some error
## smallest error is for ncomp = 2
ci.coeffs <- array(0, dim = c(dim(coeffs), nseg))
for (i in seq_len(nseg)) {
   noise1 <- matrix(rnorm(npred * nresp, 0, 0.25), npred, nresp)
   noise2 <- matrix(rnorm(npred * nresp, 0, 0.75), npred, nresp)
   noise3 <- matrix(rnorm(npred * nresp, 0, 1.00), npred, nresp)
   ci.coeffs[, , , i] <- coeffs + array(rbind(noise1, noise2, noise3), dim = dim(coeffs))
}

test_data <- list()
## 1. data with one component and one response
test_data[["one comp + one resp"]] = list(
   values = coeffs[, 1, 1, drop = F],
   ci.coeffs = ci.coeffs[, 1, 1, , drop = F]
)

## 2. data with one component and two responses
test_data[["one comp + two resp"]] = list(
   values = coeffs[, 1, , drop = F],
   ci.coeffs = ci.coeffs[, 1, , , drop = F]
)

## 3. data with three components and one response
test_data[["three comp + one resp"]] = list(
   values = coeffs[, , 1, drop = F],
   ci.coeffs = ci.coeffs[, , 1, , drop = F]
)

## 4. data with three components and two responses
test_data[["full"]] = list(
   values = coeffs,
   ci.coeffs = ci.coeffs
)

## 5. same as #4 but with all attributes and excluded variables
attr(coeffs, "name") <- "Regression coefficients"
attr(coeffs, "yaxis.name") <- "Wavelength, nm"
attr(coeffs, "yaxis.values") <- seq(100, 500, length.out = npred)
attr(coeffs, "exclrows") <- c(3)

dimnames(coeffs) <- list(
   paste0("X", seq_len(npred)),
   paste("Comp", seq_len(ncomp)),
   paste0("Y", seq_len(nresp))
)

test_data[["full + attrs"]] = list(
   values = coeffs,
   ci.coeffs = ci.coeffs[-3, , , , drop = FALSE]
)

# tests for performance statistics
context(sprintf("regcoeffs: general tests"))
test_that("if coeffs has wrong dimension it raises an error", {
   # for coeffs
   expect_error(regcoeffs(1:30))
   expect_error(regcoeffs(matrix(1:30, 10, 3)))
   expect_silent(regcoeffs(array(1:30, dim = c(10, 3, 1))))
   # for ci.coeffs
   expect_error(regcoeffs(coeffs, 1:30))
   expect_error(regcoeffs(coeffs, matrix(1:30, 10, 3)))
   expect_error(regcoeffs(coeffs, array(1:30, dim = c(10, 3, 1))))
   expect_error(regcoeffs(array(1:30, dim = c(10, 3, 1)), array(1:30, dim = c(10, 1, 1, 3))))
   expect_silent(regcoeffs(array(1:10, dim = c(10, 1, 1)), array(1:30, dim = c(10, 1, 1, 3))))
})


for (i in seq_along(test_data)) {

   ## define most useful parameters
   td <- test_data[[i]]
   case_name <- names(test_data)[i]
   ci.coeffs <- td$ci.coeffs

   m <- apply(ci.coeffs, 1:3, mean)
   err <- sweep(ci.coeffs, 1:3, m, "-")
   ssq <- apply(err^2, 1:3, sum)
   se <- sqrt((nseg - 1)/nseg * ssq)
   t <- m / se
   p <- 2 * (1 - pt(abs(t), nseg - 1))

   attrs <- mda.getattr(td$values)
   row_ind <- seq_len(dim(td$values)[1])
   if (length(attrs$exclrows) > 0) {
      row_ind <- row_ind[-attrs$exclrows]
   }

   # tests for performance statistics
   context(sprintf("regcoeffs: statistics for %s", case_name))

   rc <- regcoeffs(td$values, td$ci.coeffs)

   test_that("DoF is correct", {
      expect_equal(rc$DoF, dim(td$ci.coeffs)[4] - 1)
   })

   test_that("standard error is computed correctly", {
      expect_equivalent(rc$se[row_ind, , , drop = FALSE], se)
      expect_equal(dimnames(rc$se), dimnames(td$values))
      expect_equal(attr(rc$se, "exclrows"), attr(td$values, "exclrows"))
      expect_equal(attr(rc$se, "xaxis.values"), attr(td$values, "xaxis.values"))
      expect_equal(attr(rc$se, "xaxis.name"), attr(td$values, "xaxis.name"))
   })

   test_that("t-values are computed correctly", {
      expect_equivalent(rc$t.values[row_ind, , , drop = FALSE], t)
      expect_equal(dimnames(rc$t.values), dimnames(td$values))
      expect_equal(attr(rc$t.values, "exclrows"), attr(td$values, "exclrows"))
      expect_equal(attr(rc$t.values, "xaxis.values"), attr(td$values, "xaxis.values"))
      expect_equal(attr(rc$t.values, "xaxis.name"), attr(td$values, "xaxis.name"))
   })

   test_that("p-values are computed correctly", {
      expect_equivalent(rc$p.values[row_ind, , , drop = FALSE], p)
      expect_equal(dimnames(rc$p.values), dimnames(td$values))
      expect_equal(attr(rc$p.values, "exclrows"), attr(td$values, "exclrows"))
      expect_equal(attr(rc$p.values, "xaxis.values"), attr(td$values, "xaxis.values"))
      expect_equal(attr(rc$p.values, "xaxis.name"), attr(td$values, "xaxis.name"))
   })

   test_that("confidence interval is computed correctly", {
      ci <- confint(rc)
      expect_equal(rownames(ci), dimnames(td$values)[[1]])
      expect_equal(attr(ci, "exclrows"), attr(td$values, "exclrows"))
      expect_equal(attr(ci, "xaxis.values"), attr(td$values, "xaxis.values"))
      expect_equal(attr(ci, "xaxis.name"), attr(td$values, "xaxis.name"))
      expect_equal(nrow(ci), dim(td$values)[1])

      f <- qt(0.975, rc$DoF)
      expect_equal(ci[row_ind, 1], td$values[row_ind, 1, 1] + rc$se[row_ind, 1, 1] * -f)
      expect_equal(ci[row_ind, 2], td$values[row_ind, 1, 1] + rc$se[row_ind, 1, 1] *  f)

      f <- qt(0.995, rc$DoF)
      ci <- confint(rc, level = 0.99)
      expect_equal(ci[row_ind, 1], td$values[row_ind, 1, 1] + rc$se[row_ind, 1, 1] * -f)
      expect_equal(ci[row_ind, 2], td$values[row_ind, 1, 1] + rc$se[row_ind, 1, 1] *  f)
   })

   # tests for plots
   context(sprintf("regcoeffs: plots for %s", case_name))

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("Regcoeffs - ", case_name), pos = 4)

   test_that("regcoeffs plot raise an error if parameters are wrong", {
      expect_error(plot(rc, ncomp = 0))
      expect_error(plot(rc, ny = 0))
      expect_error(plot(rc, ncomp = 1110))
      expect_error(plot(rc, ny = 1110))
      expect_error(plot(rc, ncomp = 1:2))
      expect_error(plot(rc, ny = 1:2))
   })

   test_that("regcoeffs plot works find with default settings", {
      par(mfrow = c(2, 2))
      expect_silent(plot(rc))
      expect_silent(plot(rc, main = "Regcoeffs", ylab = "coeffs", , ncomp = dim(td$values)[2]))
      expect_silent(plot(rc, ny = dim(td$values)[3], show.ci = FALSE, show.labels = T, alpha = 0.1))
      expect_silent(plot(rc, col = c("red", "pink"), show.excluded = T, show.labels = T))
   })

   test_that("regcoeffs plot works find as scatter line plot", {
      par(mfrow = c(2, 2))
      expect_silent(plot(rc, type = "b"))
      expect_silent(plot(rc, type = "b", main = "Regcoeffs", ylab = "coeffs", ncomp = dim(td$values)[2]))
      expect_silent(plot(rc, type = "b", ny = dim(td$values)[3], show.ci = FALSE, show.labels = T, alpha = 0.1))
      expect_silent(plot(rc, type = "b", col = c("red", "pink"), show.excluded = T, show.labels = T))
   })

   test_that("regcoeffs plot works find as line plot", {
      par(mfrow = c(2, 2))
      expect_silent(plot(rc, type = "l"))
      expect_silent(plot(rc, type = "l", main = "Regcoeffs", ylab = "coeffs", , ncomp = dim(td$values)[2]))
      expect_silent(plot(rc, type = "l", ny = dim(td$values)[3], show.ci = FALSE, show.labels = T, alpha = 0.1))
      expect_silent(plot(rc, type = "l", col = c("red", "pink"), show.excluded = T, show.labels = T))
   })

   fprintf("\nOutput for regcoeffs: %s\n", case_name)
   print(rc)
   summary(rc)

   # tests for performance statistics
   context(sprintf("regcoeffs: no ci data for %s", case_name))

   rc <- regcoeffs(td$values)
   fprintf("\nOutput for regcoeffs no ci: %s\n", case_name)
   print(rc)
   summary(rc)

   test_that("regcoeffs plot works find with default settings", {
      par(mfrow = c(3, 2))
      expect_silent(plot(rc))
      expect_silent(plot(rc, col = c("red", "pink"), show.excluded = T, show.labels = T))
      expect_silent(plot(rc, type = "b"))
      expect_silent(plot(rc, type = "b", col = c("red", "pink"), show.excluded = T, show.labels = T))
      expect_silent(plot(rc, type = "l"))
      expect_silent(plot(rc, type = "l", col = c("red", "pink"), show.excluded = T, show.labels = T))
   })

}
