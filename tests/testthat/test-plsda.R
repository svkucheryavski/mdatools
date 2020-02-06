# new tests on top

setup({
   pdf(file = tempfile("mdatools-test-plsda-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-plsda-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

## prepare datasets
data(iris)
cal.ind <- c(1:25, 51:75, 101:125)
val.ind <- c(26:50, 76:100, 126:150)

Xc <- iris[cal.ind, 1:4]
Xv <- iris[val.ind, 1:4]

Xc <- mda.exclrows(Xc, c(1, 10, 15))
Xv <- mda.exclrows(Xv, c(2, 8, 16))

cc.all <- iris[cal.ind, 5]
cc.vir <- cc.all == "virginica"
cv.all <- iris[val.ind, 5]
cv.vir <- cv.all == "virginica"

context("plsda: one class model")

test_that("calibration and cross-validation works fine", {
   expect_error(m <- plsda(Xc, cc.vir, 3, cv = 1))
   expect_silent(m <- plsda(Xc, cc.vir, 3, cv = 1, classname = "virginica"))
   expect_silent(plot(m))
   expect_output(summary(m))
   expect_output(print(m))

   cat("\nOuput for PLS-DA one class model\n:")
   summary(m)
   print(m)
})

test_that("predictions work fine", {
   expect_silent(m <- plsda(Xc, cc.vir, 3, cv = 1, classname = "virginica"))

   # prediction works fine for logical vector as reference
   expect_silent(res <- predict(m, Xv, cv.vir))
   expect_silent(plot(res))
   expect_output(summary(res))
   expect_output(print(res))
   cat("\nOutput for predictions with logical vector as referencen")
   summary(res)
   print(res)

   # prediction works fine with factor as reference
   expect_silent(res <- predict(m, Xv, cv.all))
   expect_silent(plot(res))
   expect_silent(plot(res))
   expect_output(summary(res))
   expect_output(print(res))
   cat("\nOutput for predictions with factor as reference\n")
   summary(res)
   print(res)

   # prediction without reference values work fine
   expect_silent(res <- predict(m, Xv))
   expect_silent(plot(res))
   expect_output(summary(res))
   expect_output(print(res))
   cat("\nOutput for predictions with factor as reference\n")
   summary(res)
   print(res)
})

context("plsda: multiple class model")

test_that("calibration and cross-validation works fine", {
   expect_silent(m <- plsda(Xc, cc.all, 3, cv = 1))
   expect_silent(plot(m))
   expect_output(summary(m))
   expect_output(print(m))

   cat("\nOuput for PLS-DA one class model\n:")
   summary(m)
   print(m)
})

test_that("predictions work fine", {
   expect_silent(m <- plsda(Xc, cc.all, 3, cv = 1))

   # logical vector can not be used as reference
   expect_error(res <- predict(m, Xv, cv.vir))

   # prediction works fine with factor as reference
   expect_silent(res <- predict(m, Xv, cv.all))
   expect_silent(plot(res))
   expect_silent(plot(res))
   expect_output(summary(res))
   expect_output(print(res))
   cat("\nOutput for predictions with factor as reference\n")
   summary(res)
   print(res)

   # prediction without reference values work fine
   expect_silent(res <- predict(m, Xv))
   expect_silent(plot(res))
   expect_output(summary(res))
   expect_output(print(res))
   cat("\nOutput for predictions with factor as reference\n")
   summary(res)
   print(res)
})
