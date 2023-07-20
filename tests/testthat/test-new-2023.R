###########################################################
# Tests for bugs and new functionalities appeared in 2023 #
###########################################################

setup({
   #pdf(file = tempfile("mdatools-test-new-2023-", fileext = ".pdf"))
   #sink(tempfile("mdatools-test-test-new-2923-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   #dev.off()
   #sink()
})

######################################
# new bugs                           #
######################################

######################################
# bug fixes for v. 0.14.1            #
######################################
context("tests fo bug to be fixed in 0.14.1.")


test_that("bug #114 is fixed", {
   data(people)
   X <- people[, -4]
   Y <- people[, 4]
   cv <- rep(1:4, 8)

   m <- pls(X, Y, 10, cv = cv)
   expect_output(summary(m))
   expect_equivalent(crossval.str(cv), "user defined with 4 segments")
   expect_equivalent(crossval.str(rep(1:10, 4)), "user defined with 10 segments")
})

test_that("bug #112 is fixed", {
   data(people)
   m <- pca(people, scale = TRUE)

   # plot has no labels
   p <- plotResiduals(m$res$cal, show.labels = TRUE)
})

test_that("bug #111 is fixed", {
   data(people)
   X <- people[, -4]
   Y <- people[, 4]

   expect_warning({m <- pls(people[, -4], people[, 4], 7)})
   expect_silent({m <- selectCompNum(m, 7)})
})


######################################
# testing the new parameter cv.scope #
######################################
context("regmodel: tests for cv-scope.")

# simple function for calibration of MLR model
test.cal <- function(x, y, center, scale, ncomp = 1, method = "mlr", cv = FALSE) {

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

# simple function for predictions
test.pred <- function(x, coeffs, ycenter, yscale) {

   yp <- x %*% coeffs

   # unscale predicted y values
   yp <- if (is.numeric(yscale)) sweep(yp, 2, yscale, "*") else yp

   # uncenter predicted y values
   yp <- if (is.numeric(ycenter)) sweep(yp, 2, ycenter, "+") else yp


   dim(yp) <- c(nrow(x), 1, 1)
   dimnames(yp) <- list(
      rownames(x),
      dimnames(coeffs)[[2]],
      dimnames(coeffs)[[3]]
   )

   return (yp)
}

# create data set - we remove sex several other variables to ret rid of collinearity
data(people)
x <- people[, -c(4, 1, 3, 9, 10, 11)]
y <- people[,  4, drop = FALSE]

# manual implementation of cross-validation with global scaling
regcv.global <- function(x, y, center = TRUE, scale = FALSE, cv = 1) {

   n <- nrow(x)

   xcenter <- if (center) apply(x, 2, mean) else rep(0, ncol(x))
   ycenter <- if (center) apply(y, 2, mean) else rep(0, ncol(y))
   xscale <- if (scale) apply(x, 2, sd) else rep(1, ncol(x))
   yscale <- if (scale) apply(y, 2, sd) else rep(1, ncol(y))

   x <- prep.autoscale(x, xcenter, xscale)
   y <- prep.autoscale(y, ycenter, yscale)

   # get matrix with indices for cv segments
   cv_ind <- crossval(cv, nobj = n, resp = y[, 1])
   nseg <- max(cv_ind)
   nrep <- ncol(cv_ind)

   n <- nrow(x)
   y.pred <- array(0, dim = c(n, 1, 1))
   jk.coeffs <- array(0, dim = c(ncol(x), 1, 1, nseg))

   # so far only full cross-validation
   for (i in seq_len(nseg)) {
      ind <- which(cv_ind[, 1] == i)
      if (length(ind) == 0) next

      xc <- x[-ind, , drop = FALSE]
      yc <- y[-ind, , drop = FALSE]
      xt <- x[ ind, , drop = FALSE]

      mi <- test.cal(xc, yc, center = FALSE, scale = FALSE, cv = TRUE)
      jk.coeffs[, 1, 1, i] <- mi$coeffs$values
      y.pred[ind, 1, 1] <- test.pred(xt, mi$coeffs$values, ycenter = ycenter, yscale = yscale)
   }

   return (list(y.pred = y.pred, jk.coeffs = jk.coeffs))
}

# manual implementation of cross-validation with local scaling
regcv.local <- function(x, y, center = TRUE, scale = FALSE, cv = 1) {


   n <- nrow(x)
   y.pred <- array(0, dim = c(n, 1, 1))
   jk.coeffs <- array(0, dim = c(ncol(x), 1, 1, n))

   # get matrix with indices for cv segments
   cv_ind <- crossval(cv, nobj = n, resp = y[, 1])
   nseg <- max(cv_ind)
   nrep <- ncol(cv_ind)

   # so far only full cross-validation
   for (i in seq_len(nseg)) {
      ind <- which(cv_ind[, 1] == i)
      if (length(ind) == 0) next

      xc <- x[-ind, , drop = FALSE]
      yc <- y[-ind, , drop = FALSE]
      xt <- x[ ind, , drop = FALSE]

      xcenter <- if (center) apply(xc, 2, mean) else rep(0, ncol(xc))
      ycenter <- if (center) apply(yc, 2, mean) else rep(0, ncol(yc))
      xscale <- if (scale) apply(xc, 2, sd) else rep(1, ncol(xc))
      yscale <- if (scale) apply(yc, 2, sd) else rep(1, ncol(yc))

      xc <- prep.autoscale(xc, xcenter, xscale)
      yc <- prep.autoscale(yc, ycenter, yscale)
      xt <- prep.autoscale(xt, xcenter, xscale)

      mi <- test.cal(xc, yc, center = FALSE, scale = FALSE, cv = TRUE)
      jk.coeffs[, 1, 1, i] <- mi$coeffs$values
      y.pred[ind, 1, 1] <- test.pred(xt, mi$coeffs$values, ycenter = ycenter, yscale = yscale)
   }

   return (list(y.pred = y.pred, jk.coeffs = jk.coeffs))
}


test_that("global scaling scope works correctly for LOO", {

   center = TRUE
   scale = TRUE
   m <- test.cal(x, y, center = center, scale = scale)
   cvres1 <- regcv.global(x, y, center = center, scale = scale)
   cvres2 <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal, pred.fun = test.pred, cv.scope = 'global')
   expect_equivalent(cvres1$y.pred, cvres2$y.pred)
   expect_equivalent(apply(cvres1$jk.coeffs, c(1, 2, 3), sum), apply(cvres2$jk.coeffs, c(1, 2,3), sum))

   center = TRUE
   scale = FALSE
   m <- test.cal(x, y, center = center, scale = scale)
   cvres1 <- regcv.global(x, y, center = center, scale = scale)
   cvres2 <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal, pred.fun = test.pred, cv.scope = 'global')
   expect_equivalent(cvres1$y.pred, cvres2$y.pred)
   expect_equivalent(apply(cvres1$jk.coeffs, c(1, 2, 3), sum), apply(cvres2$jk.coeffs, c(1, 2,3), sum))

   center = FALSE
   scale = FALSE
   m <- test.cal(x, y, center = center, scale = scale)
   cvres1 <- regcv.global(x, y, center = center, scale = scale)
   cvres2 <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal, pred.fun = test.pred, cv.scope = 'global')
   expect_equivalent(cvres1$y.pred, cvres2$y.pred)
   expect_equivalent(apply(cvres1$jk.coeffs, c(1, 2, 3), sum), apply(cvres2$jk.coeffs, c(1, 2,3), sum))

   center = FALSE
   scale = TRUE
   m <- test.cal(x, y, center = center, scale = scale)
   cvres1 <- regcv.global(x, y, center = center, scale = scale)
   cvres2 <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal, pred.fun = test.pred, cv.scope = 'global')
   expect_equivalent(cvres1$y.pred, cvres2$y.pred)
   expect_equivalent(apply(cvres1$jk.coeffs, c(1, 2, 3), sum), apply(cvres2$jk.coeffs, c(1, 2,3), sum))

})

test_that("local scaling scope works correctly for LOO", {

   center = TRUE
   scale = TRUE
   m <- test.cal(x, y, center = center, scale = scale)
   cvres1 <- regcv.local(x, y, center = center, scale = scale)
   cvres2 <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal, pred.fun = test.pred, cv.scope = 'local')
   expect_equivalent(apply(cvres1$jk.coeffs, c(1, 2, 3), sum), apply(cvres2$jk.coeffs, c(1, 2,3), sum))
   expect_equivalent(cvres1$y.pred, cvres2$y.pred)

   center = TRUE
   scale = FALSE
   m <- test.cal(x, y, center = center, scale = scale)
   cvres1 <- regcv.local(x, y, center = center, scale = scale)
   cvres2 <- crossval.regmodel(m, x, y, cv = 1, cal.fun = test.cal, pred.fun = test.pred, cv.scope = 'local')
   expect_equivalent(apply(cvres1$jk.coeffs, c(1, 2, 3), sum), apply(cvres2$jk.coeffs, c(1, 2,3), sum))
   expect_equivalent(cvres1$y.pred, cvres2$y.pred)
})

test_that("global scaling scope works correctly for venetian blinds", {

   cv <- list("ven", 8)

   center <- TRUE
   scale <- TRUE
   m <- test.cal(x, y, center = center, scale = scale)
   cvres1 <- regcv.global(x, y, center = center, scale = scale, cv = cv)
   cvres2 <- crossval.regmodel(m, x, y, cv = cv, cal.fun = test.cal, pred.fun = test.pred, cv.scope = 'global')
   expect_equivalent(cvres1$y.pred, cvres2$y.pred)
   expect_equivalent(apply(cvres1$jk.coeffs, c(1, 2, 3), sum), apply(cvres2$jk.coeffs, c(1, 2,3), sum))

   center <- TRUE
   scale <- FALSE
   m <- test.cal(x, y, center = center, scale = scale)
   cvres1 <- regcv.global(x, y, center = center, scale = scale, cv = cv)
   cvres2 <- crossval.regmodel(m, x, y, cv = cv, cal.fun = test.cal, pred.fun = test.pred, cv.scope = 'global')
   expect_equivalent(cvres1$y.pred, cvres2$y.pred)
   expect_equivalent(apply(cvres1$jk.coeffs, c(1, 2, 3), sum), apply(cvres2$jk.coeffs, c(1, 2,3), sum))

})

# testing the new parameter cv.scope
context("pls: tests for cv-scope.")

test_that("global scaling scope works correctly for PLS", {

   data(simdata)
   X <- simdata$spectra.c
   Y <- simdata$conc.c[, 3, drop = FALSE]

   cv <- list("ven", 8)
   Xpv <- pcv::pcvpls(X, Y, ncomp = 10, cv = cv)
   m1 <- pls(X, Y, 10, cv = cv, cv.scope = 'global')
   m2 <- pls(X, Y, 10, cv = cv, cv.scope = 'local')
   m3 <- pls(X, Y, 10, x.test = Xpv, y.test = Y)

   #print(rbind(m1$res$cv$rmse, m2$res$cv$rmse, m3$res$test$rmse))
   expect_equivalent(m1$res$cv$rmse, m3$res$test$rmse)
   expect_true(all(abs(m1$res$cv$rmse - m3$res$test$rmse) < 0.000001))
   expect_false(all(abs(m1$res$cv$rmse - m2$res$cv$rmse) < 0.0001))

   cv <- rep(1:4, 25)[sample(100)]
   Xpv <- pcv::pcvpls(X, Y, ncomp = 10, cv = cv)
   m1 <- pls(X, Y, 10, cv = cv, cv.scope = 'global')
   m2 <- pls(X, Y, 10, cv = cv, cv.scope = 'local')
   m3 <- pls(X, Y, 10, x.test = Xpv, y.test = Y)

   #print(rbind(m1$res$cv$rmse, m2$res$cv$rmse, m3$res$test$rmse))
   expect_equivalent(m1$res$cv$rmse, m3$res$test$rmse)
   expect_true(all(abs(m1$res$cv$rmse - m3$res$test$rmse) < 0.0000001))
   expect_false(all(abs(m1$res$cv$rmse - m2$res$cv$rmse) < 0.0001))

})



