###################################
# Tests for regres class methods  #
###################################

setup({
   pdf(file = tempfile("mdatools-test-regres-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-regres-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# create several datasets

## number of objects, components and y-variables
nobj <- 30
ncomp <- 3
nresp <- 2

## create x, y and yp for full data, for ncomp = 2 predictions should be best
x <- matrix(runif(nobj, 0, 10), nrow = nobj, ncol = 1)
b <- matrix(c(2, 1, -1, -1), ncol = 2)

y.ref <- cbind(1, x) %*% b
y.ref <- y.ref + cbind(rnorm(nobj, 0, 0.1), rnorm(nobj, 0, 0.2))
y.err <- array(
   rbind(
      matrix(rnorm(nobj * nresp, 0, 0.50), nobj, nresp),
      matrix(rnorm(nobj * nresp, 0, 0.25), nobj, nresp),
      matrix(rnorm(nobj * nresp, 0, 0.75), nobj, nresp)
   ), dim = c(nobj, ncomp, nresp)
)

y.pred <- array(rbind(y.ref, y.ref, y.ref), dim = c(nobj, ncomp, nresp)) + y.err

test_data <- list()
## 1. data with one component and one response
test_data[["one comp + one resp"]] = list(
   y.ref = y.ref[, 1, drop = FALSE],
   y.pred = y.pred[, 1, 1, drop = FALSE],
   y.err = y.err[, 1, 1, drop = FALSE]
)

SStot <- sum((y.ref[, 1, drop = FALSE] - mean(y.ref[, 1, drop = FALSE]))^2)
SSerr <- sum((y.ref[, 1, drop = FALSE] - y.pred[, 1, 1])^2)

## 2. data with one component and three responses
test_data[["one comp + three resp"]] = list(
   y.ref = y.ref,
   y.pred = y.pred[, 1, , drop = FALSE],
   y.err = y.err[, 1, , drop = FALSE]
)

## 3. data with three components and one response
test_data[["three comp + one resp"]] = list(
   y.ref = y.ref[, 1, drop = FALSE],
   y.pred = y.pred[, , 1, drop = FALSE],
   y.err = y.err[, , 1, drop = FALSE]
)

## 4. data with three components and three responses
test_data[["full"]] = list(
   y.ref = y.ref,
   y.pred = y.pred,
   y.err = y.err
)

## 5. data from #4 with excluded rows and all attributes
rownames(y.ref) <- dimnames(y.pred)[[1]] <- paste0("O", seq_len(nobj))
colnames(y.ref) <- dimnames(y.pred)[[3]] <- paste0("Resp", seq_len(nresp))
dimnames(y.pred)[[2]] <- paste0("Comp", seq_len(ncomp))
attr(y.pred, "exclrows") <- attr(y.ref, "exclrows") <- c(5, 15, 25)
attr(y.pred, "yaxis.name") <- "My dear objects"
attr(y.pred, "yaxis.values") <- seq_len(nobj) / 100
test_data[["full + attrs"]] = list(
   y.ref = y.ref,
   y.pred = y.pred,
   y.err = y.err
)


for (i in seq_along(test_data)) {

   ## define most useful parameters
   td <- test_data[[i]]
   case_name <- names(test_data)[i]
   y.pred <- td$y.pred
   y.ref <- td$y.ref
   y.err <- if (!is.null(td$y.err)) -td$y.err
   ncomp <- dim(y.pred)[2]
   nresp <- dim(y.pred)[3]
   ncomp.selected <- if (ncomp > 1) 2 else 1

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("Regres - ", case_name), pos = 4)

   ## compute statistics manually - hardcore: no loops, no apply functions
   exclrows <- attr(y.pred, "exclrows")
   if (length(exclrows) > 0) {
      y.err <- y.err[-exclrows, , , drop = FALSE]
      y.pred <- y.pred[-exclrows, , , drop = FALSE]
      y.ref <- y.ref[-exclrows, , drop = FALSE]
   }

   ### total variance for each response
   SStot <- colSums(scale(y.ref, center = TRUE, scale = FALSE)^2)

   ### error variance for each response and component
   SSerr <- t(colSums(y.err^2))

   ### r2 values
   r2 <- 1 - sweep(SSerr, 1, SStot, "/")
   rmse <- sqrt(SSerr / nrow(y.ref))
   bias <- apply(y.err, c(3, 2), mean)


   # tests for performance statistics
   context(sprintf("regres: performance statistics for %s", case_name))

   ## create regres object
   test_that("regres object can be created", {
      expect_silent(regres(td$y.pred, td$y.ref, ncomp.selected))
   })

   ## create object for calculations
   res <- regres(td$y.pred, td$y.ref, ncomp.selected)

   test_that("regres.err works correctly", {
      err <- regres.err(y.pred, y.ref)
      expect_equal(dim(err), dim(y.err))
      expect_equivalent(err, y.err, tolerance = 10 * .Machine$double.eps)
   })

   test_that("R2 is computed correctly", {
      expect_equal(dim(res$r2), c(nresp, ncomp))
      expect_equal(dimnames(res$r2), dimnames(y.pred)[c(3, 2)])
      expect_equivalent(res$r2, r2, tolerance = 10 * .Machine$double.eps)
   })

   test_that("RMSE is computed correctly", {
      expect_equal(dim(res$rmse), c(nresp, ncomp))
      expect_equal(dimnames(res$rmse), dimnames(y.pred)[c(3, 2)])
      expect_equivalent(res$rmse, rmse, tolerance = 10 * .Machine$double.eps)
   })

   test_that("Bias is computed correctly", {
      expect_equal(dim(res$bias), c(nresp, ncomp))
      expect_equal(dimnames(res$bias), dimnames(y.pred)[c(3, 2)])
      expect_equivalent(res$bias, bias, tolerance = 10 * .Machine$double.eps)
   })

   # tests for plots
   context(sprintf("regres: plots for %s", case_name))

   test_that("Prediction plot works correctly", {
      pd <- plotPredictions(res, show.plot = FALSE)
      expect_equivalent(pd, cbind(td$y.ref[, 1], td$y.pred[, ncomp.selected, 1]))
      expect_equal(attr(pd, "exclrows"), exclrows)

      pd <- plotPredictions(res, ny = nresp, ncomp = ncomp, show.plot = FALSE)
      expect_equivalent(pd, cbind(td$y.ref[, nresp], td$y.pred[, ncomp, nresp]))
      expect_equal(attr(pd, "exclrows"), exclrows)

      par(mfrow = c(2, 2))
      expect_silent(plotPredictions(res, main = case_name))
      expect_silent(plotPredictions(res, show.labels = T, axes.equal = F))
      expect_silent(plotPredictions(res, show.labels = T, show.line = F, show.stat = T))
      expect_silent(plotPredictions(res, ncomp = ncomp, ny = nresp, col = "red", show.excluded = T))
      expect_error(plotPredictions(res, ny = 1:2))
   })

   test_that("Residuals plot works correctly", {
      pd <- plotResiduals(res, show.plot = FALSE)
      expect_equivalent(pd, cbind(td$y.ref[, 1], -td$y.err[, ncomp.selected, 1]))
      expect_equal(attr(pd, "exclrows"), exclrows)

      pd <- plotResiduals(res, ny = nresp, ncomp = ncomp, show.plot = FALSE)
      expect_equivalent(pd, cbind(td$y.ref[, nresp], -td$y.err[, ncomp, nresp]))
      expect_equal(attr(pd, "exclrows"), exclrows)

      par(mfrow = c(2, 2))
      expect_silent(plotResiduals(res, main = case_name))
      expect_silent(plotResiduals(res, show.labels = T))
      expect_silent(plotResiduals(res, show.labels = T, show.line = F))
      expect_silent(plotResiduals(res, ncomp = ncomp, ny = nresp, col = "red", show.excluded = T))
      expect_error(plotResiduals(res, ny = 1:2))
   })

   test_that("RMSE plot works correctly", {
      pd <- plotRMSE(res, show.plot = FALSE)
      expect_equivalent(pd, rmse[1, , drop = FALSE])

      pd <- plotRMSE(res, ny = nresp, show.plot = FALSE)
      expect_equivalent(pd, rmse[nresp, , drop = FALSE])

      par(mfrow = c(2, 2))
      expect_silent(plotRMSE(res, main = case_name))
      expect_silent(plotRMSE(res, ny = nresp, show.labels = T, labels = "values"))
      expect_silent(plotRMSE(res, ny = nresp, type = "b", col = "red"))
      expect_silent(plotRMSE(res, ny = nresp, type = "l", col = "red"))
      expect_error(plotRMSE(res, ny = 1:3))
   })

   ## test print and output methods
   #test_that("as.matrix method works correctly", {
   #})

   # tests if there is no reference

   context(sprintf("regres: no reference for case for %s", case_name))
   ## create object for calculations

   res <- regres(td$y.pred, y.ref = NULL, ncomp.selected)

   test_that("plot methods work well", {
      expect_error(plotResiduals(res))
      expect_error(plotRMSE(res))

      pd <- plotPredictions(res, show.plot = FALSE)
      expect_equivalent(pd, td$y.pred[, ncomp.selected, 1])
      expect_equal(attr(pd, "exclrows"), exclrows)

      pd <- plotPredictions(res, ny = nresp, ncomp = ncomp, show.plot = FALSE)
      expect_equivalent(pd, td$y.pred[, ncomp, nresp])
      expect_equal(attr(pd, "exclrows"), exclrows)

      par(mfrow = c(2, 2))
      expect_silent(plotPredictions(res, main = case_name))
      expect_silent(plotPredictions(res, show.labels = T, axes.equal = F))
      expect_silent(plotPredictions(res, show.labels = T, show.line = F, show.stat = T))
      expect_silent(plotPredictions(res, ncomp = ncomp, ny = nresp, col = "red", show.excluded = T))
      expect_error(plotPredictions(res, ny = 1:2))
   })
}
