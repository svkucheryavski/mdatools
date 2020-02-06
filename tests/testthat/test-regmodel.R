######################################
# Tests for regmodel class methods   #
######################################

setup({
   pdf(file = tempfile("mdatools-test-regmodel-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-regmodel-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

# create functions for cal and predict
test.cal <- function(x, y, center, scale, ncomp = 1, method = "mlr", cv = FALSE) {

   exclcols <- attr(x, "exclrows")
   x <- prep.autoscale(x, center = center, scale = scale)
   y <- prep.autoscale(y, center = center, scale = scale)

   fit <- lm(y ~ x)
   coeffs <- array(fit$coefficients[2:(ncol(x) + 1)], dim = c(ncol(x), 1, 1))

   m <- list(
      center = center,
      scale = scale,
      xcenter = attr(x, "prep:center"),
      xscale = attr(x, "prep:scale"),
      ycenter = attr(y, "prep:center"),
      yscale = attr(y, "prep:scale"),
      coeffs = regcoeffs(coeffs),
      method = method,
      ncomp = 1,
      ncomp.selected = 1
   )
   class(m) <- c("testmodel", "regmodel")

   if (cv) {
      return(m)
   }

   dimnames(coeffs) <- list(colnames(x), "Comp 1", "Y1")
   m$coeffs = regcoeffs(coeffs)
   return(m)
}

predict.testmodel <- function(obj, x, y = NULL, cv = FALSE) {

   x <- prep.autoscale(x, center = obj$xcenter, scale = obj$xscale)
   y.pred <- x %*% obj$coeffs$values
   y.pred <- sweep(y.pred, 2, obj$yscale, "*")
   y.pred <- sweep(y.pred, 2, obj$ycenter, "+")
   dim(y.pred) <- c(nrow(x), 1, 1)
   dimnames(y.pred) <- list(
      rownames(x),
      dimnames(obj$coeffs$values)[[2]],
      dimnames(obj$coeffs$values)[[3]]
   )
   return(regres(y.pred = y.pred, y.ref = y, ncomp.selected = 1))
}

assign("predict.testmodel", predict.testmodel, envir = .GlobalEnv)

# we remove sex several other variables to ret rid of collinearity
data(people)
x <- people[, -c(4, 1, 3, 9, 10, 11)]
y <- people[,  4, drop = FALSE]

m <- test.cal(x, y, center = TRUE, scale = TRUE)
r <- predict(m, x)


# tests for performance statistics
context(sprintf("regmodel: cross-validation"))
test_that("leave one out works correctly", {
   # leave one out
   cvres <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal)
   expect_equal(dim(cvres$y.pred), c(nrow(x), 1, 1))
   expect_equal(dim(cvres$jk.coeffs), c(ncol(x), 1, 1, nrow(x)))
   expect_gt(cor(y, cvres$y.pred[, 1, 1]), 0.9)
   expect_lt(cor(y, cvres$y.pred[, 1, 1]), cor(y, r$y.pred))

   # random without repetitions
   cvres <- crossval.regmodel(m, x, y, cv = 8, cal.fun = test.cal)
   expect_equal(dim(cvres$y.pred), c(nrow(x), 1, 1))
   expect_equal(dim(cvres$jk.coeffs), c(ncol(x), 1, 1, 8))
   expect_gt(cor(y, cvres$y.pred[, 1, 1]), 0.9)
   expect_lt(cor(y, cvres$y.pred[, 1, 1]), cor(y, r$y.pred))

   # random with repetitions
   cvres <- crossval.regmodel(m, x, y, cv = list("rand", 8, 6), cal.fun = test.cal)
   expect_equal(dim(cvres$y.pred), c(nrow(x), 1, 1))
   expect_equal(dim(cvres$jk.coeffs), c(ncol(x), 1, 1, 8))
   expect_gt(cor(y, cvres$y.pred[, 1, 1]), 0.9)
   expect_lt(cor(y, cvres$y.pred[, 1, 1]), cor(y, r$y.pred))

   # systematic
   cvres <- crossval.regmodel(m, x, y, cv = list("ven", 10), cal.fun = test.cal)
   expect_equal(dim(cvres$y.pred), c(nrow(x), 1, 1))
   expect_equal(dim(cvres$jk.coeffs), c(ncol(x), 1, 1, 10))
   expect_gt(cor(y, cvres$y.pred[, 1, 1]), 0.9)
   expect_lt(cor(y, cvres$y.pred[, 1, 1]), cor(y, r$y.pred))

})

context(sprintf("regmodel: getRegCoeffs"))

test_that("getRegcoeffs works well without inference", {
   mfull <- lm(y~x)
   coeffs <- getRegcoeffs.regmodel(m, ncomp = 1)
   expect_equivalent(mfull$coefficients, coeffs[, 1])
})

test_that("getRegcoeffs works well without inference", {
   mfull <- lm(y~x)
   cvres <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal)
   m$coeffs <- regcoeffs(m$coeffs$values, cvres$jk.coeffs)
   coeffs <- getRegcoeffs.regmodel(m, ncomp = 1, full = TRUE)

   # check standard error deviation (max 50% difference)
   se_full <- summary(mfull)$coefficients[2:ncol(x), 2]
   se <- coeffs[2:ncol(x), 2]
   expect_lt( max(abs((se_full - se)/se_full)), 0.5)

   # check t-value deviation (max 40% difference)
   t_full <- summary(mfull)$coefficients[2:ncol(x), 3]
   t <- coeffs[2:ncol(x), 3]
   expect_lt( max(abs((t_full - t)/t_full)), 0.4)

   # check p-value deviation (max 10% difference)
   # TODO: implement later when figure out correctness of JK
   #p_full <- summary(mfull)$coefficients[2:ncol(x), 4]
   #p <- coeffs[2:ncol(x), 4]
   #expect_equal( sum(p_full < 0.05), sum(p < 0.05))

})


test_that("main plots works fine", {
   m <- test.cal(x, y, center = TRUE, scale = TRUE)
   cvres <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal)
   m$res <- list()
   m$res[["cal"]] <- predict(m, x, y)
   m$res[["cv"]] <- regres(cvres$y.pred, y)
   m$coeffs <- regcoeffs(m$coeffs$values, cvres$jk.coeffs)

   # RMSE
   par(mfrow = c(2, 2))
   expect_silent(plotRMSE(m))
   expect_silent(plotRMSE(m, type = "b", legend.position = "topleft", res = list("xres" = m$res$cal)))
   expect_silent(plotRMSE(m, show.legend = F, main = "My RMSE", ylab = "RMSE value", xlab = "Complexity"))
   expect_silent(plotRMSE(m, col = c("green", "orange"), show.labels = T, labels = "values"))

   # Predictions
   par(mfrow = c(2, 2))
   expect_silent(plotPredictions(m))
   expect_silent(plotPredictions(m, legend.position = "bottom", res = list("xres" = m$res$cal)))
   expect_silent(plotPredictions(m, show.legend = F, main = "My predictions", ylab = "pred", xlab = "ref"))
   expect_silent(plotPredictions(m, col = c("green", "orange"), show.labels = T, labels = "names"))

   # Resdiauls
   par(mfrow = c(2, 2))
   expect_silent(plotYResiduals(m))
   expect_silent(plotYResiduals(m, legend.position = "bottom", res = list("xres" = m$res$cal)))
   expect_silent(plotYResiduals(m, show.legend = F, main = "My residuals", ylab = "resid", xlab = "ref"))
   expect_silent(plotYResiduals(m, col = c("green", "orange"), show.labels = T, labels = "names"))

   # Regression coefficients
   par(mfrow = c(2, 2))
   expect_silent(plotRegcoeffs(m, show.ci = F, show.labels = T, labels = "values"))
   expect_silent(plotRegcoeffs(m, col = c("red", "pink"), show.excluded = T, show.labels = T))
   expect_silent(plotRegcoeffs(m, type = "b", col = c("red", "pink"), show.excluded = T, show.labels = T))
   expect_silent(plotRegcoeffs(m, type = "l", col = c("red", "pink"), show.excluded = T, show.labels = T))

})

test_that("text outcome works fine", {
   m <- test.cal(x, y, center = TRUE, scale = TRUE)
   cvres <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal)
   m$res <- list()
   m$res[["cal"]] <- predict(m, x, y)
   m$res[["cv"]] <- regres(cvres$y.pred, y)
   m$coeffs <- regcoeffs(m$coeffs$values, cvres$jk.coeffs)

   # just to make prints to sink to check the output
   cat("\nOutput for getRegcoeffs():\n\n")
   expect_output(print(round(getRegcoeffs.regmodel(m, ncomp = 1, full = TRUE), 4)))
   expect_output(print(round(getRegcoeffs.regmodel(m, ncomp = 1, full = TRUE, alpha = 0.10), 4)))
   expect_output(print(round(getRegcoeffs.regmodel(m, ncomp = 1, full = FALSE), 4)))

   # print and summary
   expect_output(summary(m))
   expect_output(summary(m, ny = 1))
   expect_output(summary(m, ny = 1, ncomp = 1))
   expect_output(summary(m, ny = 1, ncomp = 1, res = list("xres" = m$res$cal)))
   expect_output(print(m))
})
