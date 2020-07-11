library(mdatools)
data(carbs)


x <- t(carbs$S)
y <- x + rbind(
   dnorm(1:ncol(x), 750, 200) * 10000,
   dnorm(1:ncol(x), 750, 100) * 10000,
   dnorm(1:ncol(x), 500, 100) * 10000
)

y.new <- prep.alsbasecorr(y, 5, 0.01)

par(mfrow = c(3, 1))
for (i in 1:3) {
   mdaplotg(list(
      orig = x[i, , drop = FALSE],
      bad = y[i, , drop = FALSE],
      good = y.new[i, , drop = FALSE]
      ), type = "l", lty = c(2, 1, 1), col = c("black", "red", "green"))
}

# take spectra from carbs dataset
spectra = mda.t(carbs$S)
# apply the correction
pspectra = prep.alsbasecorr(spectra, plambda = 3, p = 0.01)
# show the original and the corrected spectra individually
par(mfrow = c(3, 1))
for (i in 1:3) {
   mdaplotg(list(
      original = mda.subset(spectra, i),
      corrected = mda.subset(pspectra, i)
   ), type = "l", col = c("black", "blue"), lwd = c(2, 1), main = rownames(spectra)[i])
}