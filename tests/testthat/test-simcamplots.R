###############################################
# Tests for simca and simcares class methods  #
###############################################

setup({
   pdf(file = tempfile("mdatools-test-simcamplots-", fileext = ".pdf"))
   sink(tempfile("mdatools-test-simcamplots-", fileext = ".txt"), append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

data(iris)
x1 <- iris[1:25, 1:4]
x2 <- iris[101:125, 1:4]
m1 <- simca(x1, 'setosa', ncomp = 1, scale = T)
m2 <- simca(x2, 'versicolor', ncomp = 1, scale = T)
m <- simcam(list(m1, m2))
plotPredictions(m)

## prepare datasets
data(iris)
ind.test <- seq(2, nrow(iris), by = 2)

x.cal <- iris[-ind.test, 1:4]
x.cal1 <- x.cal[1:25, ]
x.cal1 <- mda.exclrows(x.cal1, c(1, 10))
x.cal2 <- x.cal[26:50, ]
x.cal2 <- mda.exclrows(x.cal2, c(11, 21))
x.cal3 <- x.cal[51:75, ]
x.cal3 <- mda.exclrows(x.cal3, c(6, 16))

classnames <- levels(iris[, 5])
x.test <- iris[ind.test, 1:4]
c.test <- iris[ind.test, 5]

## create models
m1 <- simca(x.cal1, classnames[1], 4, cv = 1, scale = F)
m1 <- selectCompNum(m1, 2)
m2 <- simca(x.cal2, classnames[2], 4, cv = 5, scale = T, lim.type = "chisq", alpha = 0.01)
m2 <- selectCompNum(m2, 3)
m3 <- simca(x.cal3, classnames[3], 4, cv = list("rand", 5, 10), x.test = x.test,
   c.test = c.test, scale = T, lim.type = "ddrobust", alpha = 0.10)
m3 <- selectCompNum(m3, 3)
m3 <- setDistanceLimits(m3, lim.type = "ddrobust", alpha = 0.05)


m <- simcam(list(m1, m2, m3), info = "Test SIMCAM methods")

context("simcam: main plots")

par(mfrow = c(1, 1))
plot.new()
text(0, 0, "SIMCAM - plots for model", pos = 4)


par(mfrow = c(2, 2))
test_that("model distance plot works correctly", {
   expect_silent(plotModelDistance(m))
   expect_silent(plotModelDistance(m, type = "b", show.labels = T, labels = "values"))
   expect_silent(plotModelDistance(m, nc = 2, show.labels = T, labels = "values"))
   expect_silent(plotModelDistance(m, nc = 3, show.labels = T, labels = "values"))
})

par(mfrow = c(2, 2))
test_that("discrimination power plot works correctly", {
   expect_silent(plotDiscriminationPower(m))
   expect_silent(plotDiscriminationPower(m, c(1, 3), type = "b", xticks = 1:4, xticklabels = colnames(x.cal)))
   expect_silent(plotDiscriminationPower(m, c(1, 3), show.labels = T, labels = "values"))
   expect_silent(plotDiscriminationPower(m, c(2, 3), show.labels = T))
})

par(mfrow = c(2, 2))
test_that("Cooman's plot works correctly", {
   expect_silent(plotCooman(m))
   expect_silent(plotCooman(m, c(1, 3)))
   expect_silent(plotCooman(m, c(1, 3), show.labels = T))
   expect_silent(plotCooman(m, c(2, 3), show.labels = T, show.excluded = TRUE))
})

par(mfrow = c(2, 2))
test_that("Predictions plot works correctly", {
   expect_silent(plotPredictions(m))
   expect_silent(plotPredictions(m, nc = 1))
   expect_silent(plotPredictions(m, nc = c(1, 3), show.labels = T))
   expect_silent(plotPredictions(m, nc = c(2, 3), show.labels = T, show.excluded = TRUE))
})

# just output to check in txt file
fprintf("\nSummary and print methods for model\n")
cat("-------------------------------\n")
print(m)
summary(m)

context("simcam: new predictions with reference values")

par(mfrow = c(1, 1))
plot.new()
text(0, 0, "SIMCAM - results for predictions", pos = 4)

res <- predict(m, x.test, c.test)

par(mfrow = c(2, 2))
test_that("Cooman's plot works correctly", {
   expect_silent(plotCooman(res))
   expect_silent(plotCooman(res, c(1, 3)))
   expect_silent(plotCooman(res, c(1, 3), show.labels = T))
   expect_silent(plotCooman(res, c(2, 3), show.labels = T, show.excluded = TRUE))
})

par(mfrow = c(2, 2))
test_that("Predictions plot works correctly", {
   expect_silent(plotPredictions(res))
   expect_silent(plotPredictions(res, nc = 1))
   expect_silent(plotPredictions(res, nc = c(1, 3), show.labels = T))
   expect_silent(plotPredictions(res, nc = c(2, 3), show.labels = T, show.excluded = TRUE))
})


# just output to check in txt file
fprintf("\nSummary and print methods for result with reference values\n")
cat("-------------------------------\n")
print(res)
summary(res)
showPredictions(res)
print(getConfusionMatrix(res))


context("simcam: new predictions without reference values")

par(mfrow = c(1, 1))
plot.new()
text(0, 0, "SIMCAM - results for predictions without reference", pos = 4)

res <- predict(m, x.test)

par(mfrow = c(2, 2))
test_that("Cooman's plot works correctly", {
   expect_silent(plotCooman(res))
   expect_silent(plotCooman(res, c(1, 3)))
   expect_silent(plotCooman(res, c(1, 3), show.labels = T))
   expect_silent(plotCooman(res, c(2, 3), show.labels = T, show.excluded = TRUE))
})

par(mfrow = c(2, 2))
test_that("Predictions plot works correctly", {
   expect_silent(plotPredictions(res))
   expect_silent(plotPredictions(res, nc = 1))
   expect_silent(plotPredictions(res, nc = c(1, 3), show.labels = T))
   expect_silent(plotPredictions(res, nc = c(2, 3), show.labels = T, show.excluded = TRUE))
})


# just output to check in txt file
fprintf("\nSummary and print methods for result with reference values\n")
cat("-------------------------------\n")
print(res)
summary(res)
showPredictions(res)

test_that("confusion matrix returns an error", {
   expect_error(getConfusionMatrix(res))
})
